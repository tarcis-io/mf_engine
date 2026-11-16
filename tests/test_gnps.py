"""Tests for GNPSProcessor — the records-in / records-out annotation builder."""

from __future__ import annotations

import pytest

from mf_engine import Processor
from mf_engine.processors.gnps import (
    _NUMERIC_COLUMNS,
    ANNOTATION_COLUMNS,
    GNPSInput,
    GNPSProcessor,
    GNPSResult,
)


def _sample_input() -> GNPSInput:
    return GNPSInput(
        annotations=[
            {"id": 1, "compound_name": "Skimmianine", "adduct": "[M+H]+",
             "score": 0.95, "mol_formula": "C14H13NO4"},
            {"id": 2, "compound_name": "Glucose", "adduct": "[M+H]+",
             "score": 0.80, "mol_formula": "C6H12O6"},
        ],
        feature_table=[
            {"row ID": 1, "row m/z": 300.0870, "row retention time": 5.5},
            {"row ID": 2, "row m/z": 180.0634, "row retention time": 1.2},
        ],
        library=[
            {"#Scan#": 1, "Compound_Name": "Skimmianine", "MQScore": 0.90,
             "Adduct": "[M+H]+", "molecular_formula": "C14H13NO4",
             "npclassifier_pathway": "Alkaloids",
             "npclassifier_class": "Furoquinoline alkaloids",
             "superclass": "Organoheterocyclic compounds",
             "class": "Quinolines", "subclass": "Furoquinolines"},
        ],
    )


def test_gnps_is_a_processor():
    assert issubclass(GNPSProcessor, Processor)


def test_gnps_inherits_name_from_base():
    # GNPSProcessor defines no __init__, so it inherits the base Processor's
    # optional name (defaulting to the class name). Lock that public behavior in
    # so a future refactor that adds __init__ can't silently drop it.
    assert GNPSProcessor().name == "GNPSProcessor"
    assert GNPSProcessor(name="merge-step").name == "merge-step"


def test_numeric_columns_are_a_subset_of_the_schema():
    # Guards the single-source-of-truth invariant: every column flagged numeric
    # must be a real output column, so a rename in one place can't silently leave
    # _NUMERIC_COLUMNS pointing at a column that no longer exists.
    assert _NUMERIC_COLUMNS <= set(ANNOTATION_COLUMNS)


def test_io_container_defaults_are_empty_and_independent():
    # The dataclasses default to independent empty lists (field(default_factory)),
    # not a single shared mutable default.
    a, b = GNPSInput(), GNPSInput()
    assert a.annotations == a.feature_table == a.library == []
    assert a.annotations is not b.annotations
    assert GNPSResult().annotations == []


def test_duplicate_feature_table_id_keeps_first_row():
    # The feature table is unscored; on a duplicate id the first row wins
    # (consistent with the tie-break used on the scored tables).
    inp = GNPSInput(
        annotations=[{"id": 1, "compound_name": "a", "score": 0.5}],
        feature_table=[
            {"row ID": 1, "row m/z": 100.0, "row retention time": 1.0},
            {"row ID": 1, "row m/z": 999.0, "row retention time": 9.0},  # loses
        ],
        library=[],
    )
    row = {r["Feature ID"]: r for r in GNPSProcessor().process(inp).annotations}[1]
    assert row["m/z"] == pytest.approx(100.0)
    assert row["RT (min)"] == pytest.approx(1.0)


def test_process_returns_records():
    result = GNPSProcessor().process(_sample_input())
    assert isinstance(result, GNPSResult)
    assert isinstance(result.annotations, list)
    assert all(isinstance(r, dict) for r in result.annotations)
    assert len(result.annotations) == 2


def test_output_has_exactly_the_schema_columns():
    result = GNPSProcessor().process(_sample_input())
    for row in result.annotations:
        assert list(row.keys()) == ANNOTATION_COLUMNS


def test_only_passthrough_fields_present():
    # No derived / interpretive / best-of columns should appear.
    result = GNPSProcessor().process(_sample_input())
    row = result.annotations[0]
    for forbidden in (
        "Compound", "Adduct", "Molecular Formula", "Class",   # best-of merges
        "Is Alkaloid", "MSI", "Source",                       # interpretation
        "Genus Distribution", "Plant Part", "Extract",
    ):
        assert forbidden not in row


def test_passthrough_values():
    result = GNPSProcessor().process(_sample_input())
    by_id = {r["Feature ID"]: r for r in result.annotations}

    f1 = by_id[1]
    assert f1["m/z"] == pytest.approx(300.0870, abs=1e-4)      # feature table
    assert f1["RT (min)"] == pytest.approx(5.5, abs=1e-4)
    assert f1["Score MZmine"] == pytest.approx(0.95)           # annotations
    assert f1["Score GNPS"] == pytest.approx(0.90)             # library
    assert f1["Compound MZmine"] == "Skimmianine"
    assert f1["Compound GNPS"] == "Skimmianine"
    assert f1["NPClassifier Pathway"] == "Alkaloids"
    assert f1["Class (GNPS)"] == "Quinolines"

    f2 = by_id[2]
    assert f2["Score MZmine"] == pytest.approx(0.80)
    assert f2["Score GNPS"] is None                            # no GNPS match
    assert f2["Compound GNPS"] == ""


def test_keeps_best_scoring_duplicate_annotation():
    # Two MZmine annotations for the same feature: the higher score is kept.
    inp = _sample_input()
    inp.annotations.append(
        {"id": 1, "compound_name": "Lower match", "adduct": "[M+H]+",
         "score": 0.10, "mol_formula": "CX"}
    )
    result = GNPSProcessor().process(inp)
    f1 = {r["Feature ID"]: r for r in result.annotations}[1]
    assert f1["Compound MZmine"] == "Skimmianine"
    assert f1["Score MZmine"] == pytest.approx(0.95)


def test_tied_scores_keep_first_in_input_order():
    # Two annotations for one feature with equal scores: the stable sort keeps
    # the first in input order, so the result is deterministic.
    inp = GNPSInput(
        annotations=[
            {"id": 7, "compound_name": "First", "score": 0.80},
            {"id": 7, "compound_name": "Second", "score": 0.80},
        ],
        feature_table=[],
        library=[],
    )
    result = GNPSProcessor().process(inp)
    f7 = {r["Feature ID"]: r for r in result.annotations}[7]
    assert f7["Compound MZmine"] == "First"


def test_gnps_only_feature_without_feature_table_row():
    # A feature annotated only by GNPS and absent from the feature table is
    # still included, with blank m/z / RT.
    inp = _sample_input()
    inp.library.append(
        {"#Scan#": 99, "Compound_Name": "Mystery", "MQScore": 0.88,
         "Adduct": "[M+H]+", "molecular_formula": "C10H10",
         "npclassifier_pathway": "", "npclassifier_class": "",
         "superclass": "", "class": "", "subclass": ""}
    )
    result = GNPSProcessor().process(inp)
    by_id = {r["Feature ID"]: r for r in result.annotations}
    assert 99 in by_id
    assert by_id[99]["Compound GNPS"] == "Mystery"
    assert by_id[99]["m/z"] is None


def test_empty_input_yields_no_annotations():
    result = GNPSProcessor().process(
        GNPSInput(
            annotations=[],
            feature_table=[{"row ID": 1, "row m/z": 1.0,
                            "row retention time": 1.0}],
            library=[],
        )
    )
    assert result.annotations == []


def test_fully_empty_input_yields_no_annotations():
    # All three tables empty: no annotated features, so an empty list — not a
    # crash on the empty frames.
    assert GNPSProcessor().process(GNPSInput()).annotations == []


def test_feature_ids_sorted_numerically():
    # Numeric ids should sort as numbers, not lexically ("10" before "2").
    inp = GNPSInput(
        annotations=[
            {"id": 10, "compound_name": "c10", "score": 0.5},
            {"id": 2, "compound_name": "c2", "score": 0.5},
            {"id": 1, "compound_name": "c1", "score": 0.5},
        ],
        feature_table=[],
        library=[],
    )
    ids = [r["Feature ID"] for r in GNPSProcessor().process(inp).annotations]
    assert ids == [1, 2, 10]


def test_mixed_id_types_do_not_raise():
    # A generic caller might supply mixed int/str ids; sorting must not crash.
    inp = GNPSInput(
        annotations=[
            {"id": 2, "compound_name": "num", "score": 0.5},
            {"id": "alpha", "compound_name": "txt", "score": 0.5},
        ],
        feature_table=[],
        library=[],
    )
    result = GNPSProcessor().process(inp)
    ids = [r["Feature ID"] for r in result.annotations]
    assert set(ids) == {2, "alpha"}       # numbers sort before text
    assert ids[0] == 2


def test_ids_merge_across_mismatched_types():
    # The same feature supplied as int in one table and str in another (e.g. a
    # DB source vs. a plain CSV reader) must merge into one row, not split.
    inp = GNPSInput(
        annotations=[
            {"id": 1, "compound_name": "MZ", "score": 0.5},
        ],
        feature_table=[
            {"row ID": "1", "row m/z": 300.0, "row retention time": 5.0},
        ],
        library=[
            {"#Scan#": 1.0, "Compound_Name": "GN", "MQScore": 0.9},
        ],
    )
    result = GNPSProcessor().process(inp)
    assert len(result.annotations) == 1
    row = result.annotations[0]
    assert row["Feature ID"] == 1
    assert row["Compound MZmine"] == "MZ"      # from annotations (int id)
    assert row["Compound GNPS"] == "GN"        # from library (float id)
    assert row["m/z"] == pytest.approx(300.0)  # from feature table (str id)


def _cross_type_same_id_result(rows: list[dict]) -> dict:
    # Run the processor over a single annotations table and return its one row.
    result = GNPSProcessor().process(
        GNPSInput(annotations=rows, feature_table=[], library=[])
    )
    assert len(result.annotations) == 1   # the two rows collapsed into one feature
    return result.annotations[0]


def test_cross_type_same_id_keeps_highest_score():
    # Two rows for the SAME feature within one table, supplied with different id
    # *types* (int 1 and str "1") and different scores. They must collapse to one
    # feature, and the higher-scoring row must win — here the higher score comes
    # second in input order, so a correct result cannot be input-order luck.
    #
    # This pins the ordering of the two dedup stages: scoring sorts the raw rows
    # *before* the ids are canonicalised, so the cross-type pair is still two
    # distinct raw ids when the best-per-id row is chosen, and only unified
    # afterwards. Getting that order wrong would keep the lower-scoring row.
    row = _cross_type_same_id_result(
        [
            {"id": "1", "compound_name": "LOW_str", "score": 0.10},
            {"id": 1, "compound_name": "HIGH_int", "score": 0.90},
        ]
    )
    assert row["Feature ID"] == 1
    assert row["Compound MZmine"] == "HIGH_int"
    assert row["Score MZmine"] == pytest.approx(0.90)


def test_cross_type_same_id_winner_is_independent_of_input_order():
    # Same as above with the input order reversed (higher score now first): the
    # winner must be unchanged, proving the pick is score-driven, not positional.
    row = _cross_type_same_id_result(
        [
            {"id": 1, "compound_name": "HIGH_int", "score": 0.90},
            {"id": "1", "compound_name": "LOW_str", "score": 0.10},
        ]
    )
    assert row["Feature ID"] == 1
    assert row["Compound MZmine"] == "HIGH_int"
    assert row["Score MZmine"] == pytest.approx(0.90)


def test_both_sources_merge_into_single_row():
    # A feature annotated by BOTH MZmine and GNPS yields one row carrying both
    # sources' passthrough values side by side.
    inp = GNPSInput(
        annotations=[
            {"id": 5, "compound_name": "MZname", "adduct": "[M+H]+",
             "score": 0.7, "mol_formula": "C5H5"},
        ],
        feature_table=[
            {"row ID": 5, "row m/z": 123.45, "row retention time": 2.0},
        ],
        library=[
            {"#Scan#": 5, "Compound_Name": "GNname", "MQScore": 0.6,
             "Adduct": "[M+Na]+", "molecular_formula": "C5H5"},
        ],
    )
    result = GNPSProcessor().process(inp)
    assert len(result.annotations) == 1
    row = result.annotations[0]
    assert row["Compound MZmine"] == "MZname"
    assert row["Compound GNPS"] == "GNname"
    assert row["Adduct MZmine"] == "[M+H]+"
    assert row["Adduct GNPS"] == "[M+Na]+"
    assert row["Score MZmine"] == pytest.approx(0.7)
    assert row["Score GNPS"] == pytest.approx(0.6)


def test_all_nan_score_column_does_not_raise():
    # Every annotation for a feature has a missing score: dedup must still pick
    # a single deterministic row rather than crash on the sort.
    inp = GNPSInput(
        annotations=[
            {"id": 1, "compound_name": "first", "score": None},
            {"id": 1, "compound_name": "second", "score": None},
        ],
        feature_table=[],
        library=[],
    )
    result = GNPSProcessor().process(inp)
    assert len(result.annotations) == 1
    assert result.annotations[0]["Compound MZmine"] == "first"


def test_blank_id_is_dropped_not_crashed():
    # A blank id cell (NaN from pandas, or None) is not a usable feature id. It
    # must be dropped rather than crash the merge or collapse every blank-id row
    # into one bogus feature.
    inp = GNPSInput(
        annotations=[
            {"id": float("nan"), "compound_name": "ghost", "score": 0.5},
            {"id": None, "compound_name": "ghost2", "score": 0.5},
            {"id": "nan", "compound_name": "ghost3", "score": 0.5},  # CSV blank
            {"id": 1, "compound_name": "real", "score": 0.5},
        ],
        feature_table=[],
        library=[],
    )
    result = GNPSProcessor().process(inp)
    ids = [r["Feature ID"] for r in result.annotations]
    assert ids == [1]


def test_non_numeric_score_prefers_numeric_winner():
    # A score column mixing real numbers with non-numeric junk must still keep
    # the genuine numeric best, not sort lexically. Non-numeric scores coerce to
    # NaN and lose to any real number.
    inp = GNPSInput(
        annotations=[
            {"id": 1, "compound_name": "junk", "score": "high"},
            {"id": 1, "compound_name": "winner", "score": 0.5},
        ],
        feature_table=[],
        library=[],
    )
    result = GNPSProcessor().process(inp)
    f1 = {r["Feature ID"]: r for r in result.annotations}[1]
    assert f1["Compound MZmine"] == "winner"


def test_all_non_numeric_scores_keep_first_deterministically():
    # When no score parses as a number, every row ties (all NaN) and the stable
    # sort keeps the first in input order — never a lexical "low" > "high" pick.
    inp = GNPSInput(
        annotations=[
            {"id": 1, "compound_name": "first", "score": "high"},
            {"id": 1, "compound_name": "second", "score": "low"},
        ],
        feature_table=[],
        library=[],
    )
    result = GNPSProcessor().process(inp)
    f1 = {r["Feature ID"]: r for r in result.annotations}[1]
    assert f1["Compound MZmine"] == "first"


def test_feature_table_only_feature_is_excluded():
    # The output is the union of features annotated by MZmine OR GNPS. A feature
    # present only in the feature table (annotated by neither source) is not an
    # annotation and must not appear.
    inp = GNPSInput(
        annotations=[{"id": 1, "compound_name": "annotated", "score": 0.5}],
        feature_table=[
            {"row ID": 1, "row m/z": 1.0, "row retention time": 1.0},
            {"row ID": 42, "row m/z": 9.9, "row retention time": 2.0},
        ],
        library=[],
    )
    result = GNPSProcessor().process(inp)
    ids = {r["Feature ID"] for r in result.annotations}
    assert ids == {1}


def test_large_integer_id_keeps_full_precision():
    # Ids beyond 2**53 must not be rounded through float. Supplied as int or as
    # the string a CSV reader would yield, the id round-trips exactly and the
    # two forms merge into one feature.
    big = 9007199254740993  # 2**53 + 1, not representable as float
    inp = GNPSInput(
        annotations=[{"id": big, "compound_name": "MZ", "score": 0.5}],
        feature_table=[],
        library=[{"#Scan#": str(big), "Compound_Name": "GN", "MQScore": 0.9}],
    )
    result = GNPSProcessor().process(inp)
    assert len(result.annotations) == 1
    row = result.annotations[0]
    assert row["Feature ID"] == big
    assert row["Compound MZmine"] == "MZ"
    assert row["Compound GNPS"] == "GN"


def test_mixed_integer_and_float_ids_keep_their_type():
    # When one table carries both integer-valued ids and a genuinely non-integer
    # id (e.g. 1.5), the integer ids must still come out as int — not promoted to
    # float by pandas inferring a single float dtype for the id set. The merge key
    # was always correct; this pins the *emitted* Feature ID type.
    inp = GNPSInput(
        annotations=[
            {"id": 1.5, "compound_name": "half", "score": 0.5},
            {"id": -3, "compound_name": "neg", "score": 0.5},
            {"id": 4, "compound_name": "four", "score": 0.5},
        ],
        feature_table=[],
        library=[],
    )
    by_id = {r["Feature ID"]: r for r in GNPSProcessor().process(inp).annotations}
    assert isinstance(by_id[-3]["Feature ID"], int)
    assert isinstance(by_id[4]["Feature ID"], int)
    assert isinstance(by_id[1.5]["Feature ID"], float)


def test_large_integer_id_keeps_precision_despite_float_sibling():
    # A >2**53 integer id sharing a table with a non-integer id must still keep
    # full precision. Without object-dtype handling pandas would infer float64 for
    # the id set and round 9007199254740993 to ...992.0.
    big = 9007199254740993  # 2**53 + 1
    inp = GNPSInput(
        annotations=[
            {"id": big, "compound_name": "big", "score": 0.5},
            {"id": 1.5, "compound_name": "half", "score": 0.5},
        ],
        feature_table=[],
        library=[],
    )
    ids = {r["Feature ID"] for r in GNPSProcessor().process(inp).annotations}
    assert big in ids  # exact, not 9007199254740992.0


def test_feature_id_zero_is_a_real_feature():
    # id 0 is a valid feature id and must not be confused with a missing/blank id
    # (the missing check is `is None` / NaN, not falsiness).
    inp = GNPSInput(
        annotations=[{"id": 0, "compound_name": "zero", "score": 0.5}],
        feature_table=[],
        library=[],
    )
    ids = [r["Feature ID"] for r in GNPSProcessor().process(inp).annotations]
    assert ids == [0]


def test_process_is_idempotent_and_does_not_mutate_input():
    # Calling process twice yields equal output and leaves the input dataclass
    # untouched — a processor should be a pure function of its input.
    original = [{"id": 1, "compound_name": "a", "score": 0.5}]
    inp = GNPSInput(annotations=[dict(original[0])], feature_table=[], library=[])
    proc = GNPSProcessor()
    first = proc.process(inp).annotations
    second = proc.process(inp).annotations
    assert first == second
    assert inp.annotations == original  # input not mutated
