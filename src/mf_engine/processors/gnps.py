"""GNPSProcessor — build a per-feature annotation table from MZmine + GNPS data.

This processor merges three input tables (MZmine annotations, the GNPS feature
table, and the GNPS library matches) into a single *annotations* table: one
record per feature annotated by either source, carrying the measured m/z and
retention time, the per-source name/adduct/formula/score, and the GNPS taxonomy
fields.

Every output field is a **direct passthrough** of an input column — nothing is
interpreted, scored against a threshold, picked between sources, or otherwise
derived. Any downstream decisions (choosing a best-of name, confidence levels,
compound-class flags, group comparisons, summaries, ...) are left to the caller.

Unlike a file-based pipeline, this processor is **I/O-free**: it accepts and
returns plain Python records (``list[dict]``), one dict per row. That maps
directly onto loaded CSV rows, database result sets, JSON arrays, or
``DataFrame.to_dict("records")``. pandas is used internally for the merge, but
never appears at the boundary.

Feature matching
----------------
Features are joined across the three tables by the **value** of their shared
feature id (``id`` / ``row ID`` / ``#Scan#`` — the first recognised column in
each table). Ids are canonicalised first, so the same feature matches even when
sources disagree on type: ``1`` (int), ``1.0`` (float from a CSV with blanks)
and ``"1"`` (str from a plain CSV reader) are treated as one feature, and
integer-valued ids keep full precision rather than rounding through ``float``.
A feature annotated by only one source is still emitted, with the other
source's columns left blank — the output row set is the **union of features
annotated by either source**. A feature present only in the feature table
(annotated by neither MZmine nor GNPS) is *not* an annotation and is omitted.

Robustness
----------
The processor tolerates the mess real tables carry:

* **Missing id** — a blank id cell (``None``, or the ``NaN`` pandas yields for
  a blank CSV cell) is not a usable key, so that row is dropped rather than
  merged with every other blank-id row.
* **Score column** — the per-source "best match" is the row with the highest
  ``score`` / ``MQScore``. Scores are compared numerically; a non-numeric or
  blank score is treated as missing and loses to any real number. When *no*
  score for a feature parses as a number, the first row in input order is kept,
  so the result is always deterministic.
* **Missing column or table** — a source whose table is empty, or which lacks a
  mapped column, contributes blanks for those fields instead of raising.

Example
-------
    from mf_engine.processors.gnps import GNPSProcessor, GNPSInput

    result = GNPSProcessor().process(
        GNPSInput(
            annotations=annotations_records,
            feature_table=feature_table_records,
            library=library_records,
        )
    )
    result.annotations  # -> list[dict], one per feature
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any

import pandas as pd

from mf_engine.processors.base import Processor

# ---------------------------------------------------------------------------
# Schema
# ---------------------------------------------------------------------------

# Candidate names for the shared feature id in each input table (the first one
# present is used). ``id`` == ``row ID`` == ``#Scan#`` across the three tables.
_ANNOTATION_ID_CANDIDATES = ["id", "row ID", "featureID", "#Scan#"]
_FEATURE_TABLE_ID_CANDIDATES = ["row ID", "id", "featureID"]
_LIBRARY_ID_CANDIDATES = ["#Scan#", "scan", "Scan", "feature_id"]

# The output schema is a single source of truth: output column -> source input
# column, per table. Every output column is a verbatim copy of one input column
# (no derived values). The column order and ANNOTATION_COLUMNS are both derived
# from these mappings, so adding a field in one place can never drift.
_FEATURE_TABLE_FIELDS = {
    "m/z": "row m/z",
    "RT (min)": "row retention time",
}
_ANNOTATION_FIELDS = {
    "Score MZmine": "score",
    "Compound MZmine": "compound_name",
    "Adduct MZmine": "adduct",
    "Formula MZmine": "mol_formula",
}
_LIBRARY_FIELDS = {
    "Score GNPS": "MQScore",
    "Compound GNPS": "Compound_Name",
    "Adduct GNPS": "Adduct",
    "Formula GNPS": "molecular_formula",
    "NPClassifier Pathway": "npclassifier_pathway",
    "NPClassifier Class": "npclassifier_class",
    "Superclass (GNPS)": "superclass",
    "Class (GNPS)": "class",
    "Subclass (GNPS)": "subclass",
}

# Output columns, in order: the feature id, then every mapped field.
ANNOTATION_COLUMNS = [
    "Feature ID",
    *_FEATURE_TABLE_FIELDS,
    *_ANNOTATION_FIELDS,
    *_LIBRARY_FIELDS,
]

# Numeric output columns (coerced to float / None); everything else is text.
_NUMERIC_COLUMNS = {"m/z", "RT (min)", "Score MZmine", "Score GNPS"}

# A record is a single table row as a mapping of column name -> value.
Record = dict[str, Any]
Records = list[Record]


# ---------------------------------------------------------------------------
# I/O containers
# ---------------------------------------------------------------------------


@dataclass
class GNPSInput:
    """The three input tables, each a list of records (one dict per row).

    Column names follow the MZmine / GNPS conventions:
      * annotations:    id, compound_name, adduct, score, mol_formula, ...
      * feature_table:  row ID, row m/z, row retention time, ...
      * library:        #Scan#, Compound_Name, MQScore, Adduct,
                        molecular_formula, npclassifier_*, superclass, class,
                        subclass, ...
    """

    annotations: Records = field(default_factory=list)
    feature_table: Records = field(default_factory=list)
    library: Records = field(default_factory=list)


@dataclass
class GNPSResult:
    """Output of GNPSProcessor: the per-feature annotation records."""

    annotations: Records = field(default_factory=list)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _frame(records: Records) -> pd.DataFrame:
    """Build a DataFrame from records, keeping columns at ``object`` dtype.

    Why object dtype: ``pd.DataFrame`` infers a column's dtype from its values, so
    an id column mixing an integer with any non-integer (e.g. ``1`` and ``1.5``)
    is inferred as ``float64`` *at construction* — before canonicalisation runs —
    silently rounding integers through float (``9007199254740993`` ->
    ``...992.0``). Forcing object dtype keeps every cell as its original Python
    value, so :func:`_canonical_id` sees exact ids. Numeric handling is
    unaffected: scoring and the numeric output columns call
    ``pd.to_numeric(..., errors="coerce")`` themselves, which coerces object input
    fine. An empty record list yields an empty frame (no columns), as before.
    """
    records = list(records)
    return pd.DataFrame(records, dtype=object) if records else pd.DataFrame()


def _first_present(columns: Any, candidates: list[str]) -> str | None:
    """Return the first column name from ``candidates`` present in ``columns``."""
    cols = set(columns)
    for col in candidates:
        if col in cols:
            return col
    return None


# Sentinel returned by :func:`_canonical_id` for a value that is not a usable
# feature id (a missing/blank id, i.e. ``None`` or ``NaN``). Rows carrying it are
# dropped during indexing rather than merged together under a bogus shared key.
_MISSING_ID = object()


def _is_bool(value: Any) -> bool:
    """True for a ``bool``, which — being an ``int`` subclass — must be treated as
    text rather than as the integer ids ``0``/``1``.

    Shared by :func:`_canonical_id` and :func:`_sort_key` so the rule lives in one
    place and the merge key and the sort order can't drift apart.
    """
    return isinstance(value, bool)


def _canonical_id(value: Any) -> Any:
    """Normalise a feature id so the same feature matches across tables.

    Feature ids arrive in whatever type each source happens to use: a DB or
    JSON source yields ``int``, pandas reading a CSV with blanks yields
    ``float`` (``1.0``), and a plain CSV reader yields ``str`` (``"1"``). All of
    these denote the same feature and must collapse to one key, or the merge
    silently splits a feature into separate half-populated rows.

    The rule: an integer-valued id becomes ``int`` (so ``1``, ``1.0`` and
    ``"1"`` all become ``1``, and large ids keep full integer precision rather
    than being rounded through ``float``); a non-integer number becomes
    ``float``; anything else becomes ``str``. A missing id — ``None`` or a
    floating-point ``NaN`` (what pandas yields for a blank id cell) — returns the
    :data:`_MISSING_ID` sentinel so the row can be dropped instead of crashing.

    Because a *string* id is parsed numerically, some text that merely looks
    numeric is canonicalised rather than kept verbatim — intentional, since real
    feature ids are numbers that different readers stringify differently:

    * ``"1e3"`` (and other float-parseable text) becomes the number ``1000`` —
      a scientific-notation id merges with the integer ``1000``.
    * ``"NaN"`` / ``"nan"`` parses to a float ``NaN``, so it is treated as a
      *missing* id and its row is dropped — a genuine text id literally spelled
      "NaN" would therefore disappear. ``"inf"`` / ``"-inf"`` parses to an
      infinity, which is not a usable id and is kept as text (``fallback``) so it
      can't masquerade as a feature. These collisions don't arise with real
      MZmine/GNPS numeric ids; only a non-numeric text id survives untouched.
    """
    if value is None:
        return _MISSING_ID
    if _is_bool(value):  # bool is an int subclass; treat as text
        return str(value)
    # Integer ids stay exact: never route through float, so ids beyond 2**53
    # don't lose precision. This also catches numpy integer scalars (int64,
    # etc.) that pandas yields, which are not Python ``int`` instances.
    if isinstance(value, int):
        return value
    if hasattr(value, "__index__"):  # numpy.int64 and other integral scalars
        return int(value)
    if isinstance(value, str):
        # Try an exact integer parse first (keeps full precision for "12345...").
        try:
            return int(value)
        except ValueError:
            pass
        # Then a float parse for things like "1.0" / "1.5" / "1e3".
        try:
            as_float = float(value)
        except ValueError:
            return value  # genuine text id, kept as-is
        return _from_float(as_float, fallback=value)
    # Any other numeric-ish value (e.g. a numpy/pandas float).
    try:
        as_float = float(value)
    except (TypeError, ValueError):
        return str(value)
    return _from_float(as_float, fallback=str(value))


def _from_float(as_float: float, fallback: Any) -> Any:
    """Canonicalise a parsed float: int when integer-valued, else float.

    NaN (a blank id cell) and ±inf are not usable ids: NaN -> missing sentinel,
    inf -> ``fallback`` kept as text so it can't masquerade as a feature.
    """
    if math.isnan(as_float):
        return _MISSING_ID
    if math.isinf(as_float):
        return str(fallback)
    as_int = int(as_float)
    return as_int if as_float == as_int else as_float


def _canonical_index(index: pd.Index) -> pd.Index:
    """Return ``index`` with every label canonicalised by :func:`_canonical_id`.

    Built at ``object`` dtype so a mix of int and non-integer ids keeps each
    label's exact Python type (and full integer precision) instead of pandas
    coercing the whole index to ``float64``.
    """
    return pd.Index([_canonical_id(v) for v in index], name=index.name, dtype=object)


def _drop_missing_ids(indexed: pd.DataFrame) -> pd.DataFrame:
    """Drop rows whose canonical id is :data:`_MISSING_ID` (missing/blank id)."""
    keep = [label is not _MISSING_ID for label in indexed.index]
    return indexed[keep]


def _unique_by_canonical_id(indexed: pd.DataFrame) -> pd.DataFrame:
    """Canonicalise a just-set index, drop missing ids, keep the first per id.

    Shared by every indexing path (annotations, library, feature table) so they
    stay consistent: canonicalise the labels, drop rows with no usable id, and
    collapse any ids that canonicalise to the same value (e.g. ``1`` and
    ``"1"``) by keeping the first — which upstream has already ordered so the
    best-scoring row wins.
    """
    indexed = indexed.copy()
    indexed.index = _canonical_index(indexed.index)
    indexed = _drop_missing_ids(indexed)
    return indexed[~indexed.index.duplicated(keep="first")]


def _best_per_feature(
    df: pd.DataFrame, id_candidates: list[str], score_col: str | None = None
) -> pd.DataFrame | None:
    """Index ``df`` by feature id, keeping the highest-``score_col`` row per id.

    The sort is stable, so when several rows for a feature share the same score
    the first one in input order is kept — making the result deterministic. With
    ``score_col=None`` (or a column the table lacks) no scoring is applied and
    the first row per id wins — the behaviour wanted for the unscored feature
    table. Returns ``None`` when the table has no recognisable id column (e.g.
    empty).
    """
    id_col = _first_present(df.columns, id_candidates)
    if id_col is None:
        return None
    if score_col is not None and score_col in df.columns:
        # Sort by the score *numerically*: a non-numeric or blank score column
        # would otherwise sort lexically (``"low" > "high"``) and silently keep
        # the wrong row. Non-numeric/blank scores coerce to NaN and sort last,
        # so a real numeric winner is always preferred over them.
        score = pd.to_numeric(df[score_col], errors="coerce")
        ordered = df.iloc[
            score.reset_index(drop=True)
            .sort_values(ascending=False, na_position="last", kind="stable")
            .index
        ]
    else:
        ordered = df
    # Keep the first row per raw id (already highest-scoring), then let the
    # shared helper canonicalise, drop missing ids, and collapse cross-type
    # duplicates (e.g. 1 and "1") — keeping the first, i.e. the best-scoring.
    indexed = ordered.drop_duplicates(subset=[id_col], keep="first").set_index(id_col)
    return _unique_by_canonical_id(indexed)


def _column(
    source: pd.DataFrame | None, feature_ids: pd.Index, src_col: str
) -> pd.Series:
    """Align ``source[src_col]`` to ``feature_ids`` (missing rows/table -> NA)."""
    if source is not None and src_col in source.columns:
        return source[src_col].reindex(feature_ids)
    return pd.Series(pd.NA, index=feature_ids)


def _sort_key(value: Any) -> tuple[int, Any]:
    """Order ids so mixed/incomparable types never raise.

    Numbers sort before non-numbers; within each group the natural order
    applies. Ids are already canonicalised by :func:`_canonical_id` before they
    reach here — so numeric values are ``int``/``float`` and the rest ``str`` —
    which is exactly why this stays consistent with the merge key.
    """
    if _is_bool(value):  # bool is an int subclass; treat as text
        return (1, str(value))
    if isinstance(value, (int, float)):
        return (0, float(value))
    return (1, str(value))


# ---------------------------------------------------------------------------
# The processor
# ---------------------------------------------------------------------------


class GNPSProcessor(Processor):
    """Merge MZmine + GNPS tables into a per-feature annotation table.

    Input and output are plain records (``list[dict]``) via :class:`GNPSInput`
    and :class:`GNPSResult`, so the processor is independent of any storage
    format (CSV, database, JSON, ...).
    """

    def process(self, data: GNPSInput) -> GNPSResult:
        """Build the annotations table from the three input tables."""
        annotations = _frame(data.annotations)
        feature_table = _frame(data.feature_table)
        library = _frame(data.library)

        records = self._build_annotations(annotations, feature_table, library)
        return GNPSResult(annotations=records)

    # -- internal -----------------------------------------------------------

    def _build_annotations(
        self,
        annotations: pd.DataFrame,
        feature_table: pd.DataFrame,
        library: pd.DataFrame,
    ) -> Records:
        # Index each table by feature id. The annotation and library tables keep
        # the highest-scoring row per feature; the feature table is unscored, so
        # the first row per id wins on the off chance of a duplicate id (same
        # path, score_col=None).
        ann = _best_per_feature(annotations, _ANNOTATION_ID_CANDIDATES, "score")
        lib = _best_per_feature(library, _LIBRARY_ID_CANDIDATES, "MQScore")
        feat = _best_per_feature(feature_table, _FEATURE_TABLE_ID_CANDIDATES)

        # The feature set is the union of features annotated by either source.
        ids: set[Any] = set()
        if ann is not None:
            ids |= set(ann.index)
        if lib is not None:
            ids |= set(lib.index)
        if not ids:
            return []
        # Keep the ordered ids as a plain Python list: it is the source of truth
        # for the Feature ID *values and their types*. The pandas index built from
        # it is object-dtype so alignment never coerces a mixed int/float id set
        # to float; but the emitted Feature ID is taken from this list, not from
        # the index round-trip, so each id keeps its exact Python type.
        ordered_ids = sorted(ids, key=_sort_key)
        feature_ids = pd.Index(ordered_ids, name="Feature ID", dtype=object)

        # Assemble the output column-by-column (vectorised); each column is a
        # straight copy of a source column aligned to the feature id index.
        # _column tolerates a missing table/column by filling with NA.
        out = pd.DataFrame(index=feature_ids)
        for source, mapping in (
            (feat, _FEATURE_TABLE_FIELDS),
            (ann, _ANNOTATION_FIELDS),
            (lib, _LIBRARY_FIELDS),
        ):
            for out_col, src_col in mapping.items():
                out[out_col] = _column(source, feature_ids, src_col)

        return self._to_records(ordered_ids, out)

    @staticmethod
    def _to_records(feature_ids: list[Any], out: pd.DataFrame) -> Records:
        """Convert the assembled frame to records with tidy value types.

        ``feature_ids`` is the ordered list of canonical ids (same order as
        ``out``'s rows); it supplies the ``Feature ID`` column directly so each id
        keeps its exact Python type rather than being re-coerced through a pandas
        index. Numeric columns become ``float`` or ``None``; text columns become
        ``str`` or ``""``. Coercion is done per-column (vectorised).
        """
        # Row-position index so each tidied column lines up with feature_ids by
        # position, independent of the (object-dtype) feature id index.
        tidy = pd.DataFrame(index=range(len(feature_ids)))
        for col in out.columns:
            series = out[col].reset_index(drop=True)
            if col in _NUMERIC_COLUMNS:
                # NaN -> None (object dtype so None survives to_dict).
                nums = pd.to_numeric(series, errors="coerce")
                tidy[col] = nums.astype(object).where(nums.notna(), None)
            else:
                # Text columns are assumed string-valued at the source (compound
                # name, adduct, formula, taxonomy). A numeric value slipping into
                # one would stringify verbatim (float 5.0 -> "5.0"); that never
                # happens with real MZmine/GNPS exports, so no numeric guard here.
                tidy[col] = series.where(series.notna(), "").astype(str)

        # to_dict("records") is typed as list[dict[Hashable, Any]]; prepend the
        # Feature ID from the id list and rebuild with str keys so the return
        # matches the declared Records type exactly.
        raw = tidy.to_dict("records")
        return [
            {"Feature ID": fid, **{str(k): v for k, v in row.items()}}
            for fid, row in zip(feature_ids, raw)
        ]
