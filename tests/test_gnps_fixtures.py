"""Integration tests for GNPSProcessor against real MZmine + GNPS exports.

These run the processor end-to-end over actual data files in ``tests/fixtures``
(an MZmine annotations export, a reformatted GNPS feature table, and a GNPS
library-match TSV), complementing the hand-built cases in ``test_gnps.py``.

The fixtures are loaded exactly as a caller would: read each file with pandas
and hand the rows to the processor as ``list[dict]`` via ``to_dict("records")``.
Paths are resolved against this test file's own directory so the suite runs the
same no matter where pytest is invoked from. The fixtures are committed to the
repository, so a missing one is a broken checkout and fails the test rather than
skipping it silently.
"""

from __future__ import annotations

from pathlib import Path
from typing import cast

import pytest

from mf_engine.processors.gnps import (
    ANNOTATION_COLUMNS,
    GNPSInput,
    GNPSProcessor,
)

pd = pytest.importorskip("pandas")

FIXTURES = Path(__file__).parent / "fixtures"

_ANNOTATIONS_CSV = FIXTURES / "mzmine_annotations.csv"
_FEATURE_TABLE_CSV = FIXTURES / "gnps_feature_table_reformated.csv"
_LIBRARY_TSV = FIXTURES / "gnps_merged_results.tsv"


def _records(path: Path, *, sep: str = ",") -> list[dict]:
    """Load a committed fixture as records; a missing one is a broken checkout."""
    assert path.exists(), f"committed fixture missing: {path.name}"
    return cast("list[dict]", pd.read_csv(path, sep=sep).to_dict("records"))


@pytest.fixture(scope="module")
def real_input() -> GNPSInput:
    return GNPSInput(
        annotations=_records(_ANNOTATIONS_CSV),
        feature_table=_records(_FEATURE_TABLE_CSV),
        library=_records(_LIBRARY_TSV, sep="\t"),
    )


@pytest.fixture(scope="module")
def real_rows(real_input: GNPSInput) -> list[dict]:
    return GNPSProcessor().process(real_input).annotations


def test_fixtures_are_consumable_by_the_processor(real_rows: list[dict]) -> None:
    # Smoke test: the real exports round-trip through the records boundary and
    # produce a non-empty annotation table.
    assert real_rows
    assert all(isinstance(r, dict) for r in real_rows)


def test_every_row_has_exactly_the_schema_columns(real_rows: list[dict]) -> None:
    for row in real_rows:
        assert list(row.keys()) == ANNOTATION_COLUMNS


def test_row_count_is_the_union_of_annotated_features(
    real_input: GNPSInput, real_rows: list[dict]
) -> None:
    # The output is the union of features annotated by MZmine OR GNPS. Derive the
    # expected id set straight from the inputs (coercing ids to int the way the
    # processor canonicalises them) and compare — no magic hard-coded count.
    def _ids(records: list[dict], key: str) -> set[int]:
        col = pd.to_numeric(pd.Series([r.get(key) for r in records]), errors="coerce")
        return {int(v) for v in col.dropna()}

    expected = _ids(real_input.annotations, "id") | _ids(real_input.library, "#Scan#")
    got = {r["Feature ID"] for r in real_rows}
    assert got == expected


def test_feature_ids_are_ints_sorted_ascending(real_rows: list[dict]) -> None:
    ids = [r["Feature ID"] for r in real_rows]
    assert all(isinstance(i, int) for i in ids)
    assert ids == sorted(ids)


def test_mz_and_rt_populated_for_features_in_the_feature_table(
    real_input: GNPSInput, real_rows: list[dict]
) -> None:
    # Every annotated feature in this dataset also appears in the feature table,
    # so m/z and RT must be filled (numeric, not None) for all of them.
    ft_ids = {
        int(v)
        for v in pd.to_numeric(
            pd.Series([r.get("row ID") for r in real_input.feature_table]),
            errors="coerce",
        ).dropna()
    }
    for row in real_rows:
        if row["Feature ID"] in ft_ids:
            assert isinstance(row["m/z"], float)
            assert isinstance(row["RT (min)"], float)


def test_both_source_compounds_are_passed_through(real_rows: list[dict]) -> None:
    # The dataset has features annotated by MZmine only, by GNPS only, and by
    # both — so all three populations must be represented.
    only_mz = sum(
        1 for r in real_rows if r["Compound MZmine"] and not r["Compound GNPS"]
    )
    only_gn = sum(
        1 for r in real_rows if r["Compound GNPS"] and not r["Compound MZmine"]
    )
    both = sum(1 for r in real_rows if r["Compound MZmine"] and r["Compound GNPS"])
    assert only_mz > 0
    assert only_gn > 0
    assert both > 0


def test_missing_library_columns_fill_blank_without_raising(
    real_rows: list[dict],
) -> None:
    # The GNPS TSV in this dataset lacks the molecular-formula and ClassyFire
    # taxonomy columns. Those output fields must be present and blank for every
    # row (the "missing column -> blanks, never raise" contract), not absent.
    absent_in_source = [
        "Formula GNPS",
        "Superclass (GNPS)",
        "Class (GNPS)",
        "Subclass (GNPS)",
    ]
    for row in real_rows:
        for col in absent_in_source:
            assert col in row
            assert row[col] == ""


def test_feature_table_only_features_are_excluded_in_real_data(
    real_input: GNPSInput, real_rows: list[dict]
) -> None:
    # At integration scale: features present in the feature table but annotated
    # by neither MZmine nor GNPS must not appear in the output (mirrors the
    # unit-level test_feature_table_only_feature_is_excluded). Only assert the
    # exclusion if the dataset actually contains such features.
    def _ids(records: list[dict], key: str) -> set[int]:
        col = pd.to_numeric(pd.Series([r.get(key) for r in records]), errors="coerce")
        return {int(v) for v in col.dropna()}

    annotated = _ids(real_input.annotations, "id") | _ids(
        real_input.library, "#Scan#"
    )
    feature_only = _ids(real_input.feature_table, "row ID") - annotated
    got = {r["Feature ID"] for r in real_rows}
    if feature_only:
        assert not (feature_only & got)


def test_best_scoring_annotation_kept_per_feature(
    real_input: GNPSInput, real_rows: list[dict]
) -> None:
    # The annotations export carries several rows per feature; the processor must
    # keep the highest MZmine score for each. Check it against the raw maxima.
    best_by_id: dict[int, float] = {}
    for rec in real_input.annotations:
        score = pd.to_numeric(rec.get("score"), errors="coerce")
        if pd.isna(score):
            continue
        fid = pd.to_numeric(rec.get("id"), errors="coerce")
        if pd.isna(fid):
            continue
        fid = int(fid)
        best_by_id[fid] = max(best_by_id.get(fid, float("-inf")), float(score))

    for row in real_rows:
        kept = row["Score MZmine"]
        if kept is None:
            continue
        assert kept == pytest.approx(best_by_id[row["Feature ID"]])
