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

See the README for a usage example and the loading-from-CSV pattern.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any, cast

import numpy as np
import pandas as pd

from mf_engine.processors.base import Processor

# ---------------------------------------------------------------------------
# Schema
# ---------------------------------------------------------------------------

# Candidate names for the shared feature id in each input table, in preference
# order (first present wins). ``id`` == ``row ID`` == ``#Scan#`` across tables.
# The spellings differ per list by design — each is the set of header variants
# that tool actually emits (e.g. MZmine's ``featureID`` vs GNPS's ``feature_id``).
_ANNOTATION_ID_CANDIDATES = ["id", "row ID", "featureID", "#Scan#"]
_FEATURE_TABLE_ID_CANDIDATES = ["row ID", "id", "featureID"]
_LIBRARY_ID_CANDIDATES = ["#Scan#", "scan", "Scan", "feature_id"]

# Output schema, single source of truth: output column -> source input column,
# per table (every output column is a verbatim copy of one input column).
# ANNOTATION_COLUMNS and the column order both derive from these, so they can't
# drift.
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

# Public output schema: the exact keys, in order, of every GNPSResult record
# (feature id, then each mapped field). Re-exported from mf_engine.processors.
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


def _records_to_frame(records: Records) -> pd.DataFrame:
    """Build a DataFrame from records, forcing ``object`` dtype.

    Object dtype keeps every cell as its original Python value, so
    :func:`_canonical_id` sees exact ids instead of a column pandas has already
    coerced to ``float64`` (which would round large integer ids through float at
    construction time). The numeric output paths coerce with ``pd.to_numeric``
    themselves and are unaffected.
    """
    records = list(records)
    frame = pd.DataFrame(records, dtype=object) if records else pd.DataFrame()
    return cast(pd.DataFrame, frame)


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
    # bool is an int subclass; shared by _canonical_id and _sort_key so True/False
    # ids are kept as text, not the integer ids 0/1, in both the key and the sort.
    return isinstance(value, bool)


def _canonical_id(value: Any) -> Any:
    """Normalise a feature id so the same feature matches across tables.

    The rule: an integer-valued id becomes ``int`` (so ``1``, ``1.0`` and ``"1"``
    all collapse to ``1``, at full precision for large ids); a non-integer number
    becomes ``float``; anything else stays ``str``. A missing id (``None`` or the
    ``NaN`` pandas yields for a blank cell) returns :data:`_MISSING_ID` so the row
    is dropped rather than crashing. See the module docstring ("Feature
    matching") for why, and for the edge cases of parsing string ids numerically
    (``"1e3"`` merges with ``1000``; ``"inf"`` stays text; ``"nan"`` is missing).
    """
    if value is None:
        return _MISSING_ID
    if _is_bool(value):
        return str(value)
    # Integer ids stay exact (never routed through float), so ids beyond 2**53
    # keep full precision. __index__ also catches numpy integer scalars.
    if isinstance(value, int):
        return value
    if hasattr(value, "__index__"):
        return int(value)
    if isinstance(value, str):
        # Exact integer parse first (full precision), then float ("1.0"/"1e3").
        try:
            return int(value)
        except ValueError:
            pass
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

    NaN -> missing sentinel; ±inf -> ``fallback`` kept as text (see
    :func:`_canonical_id` for why neither is a usable id).
    """
    if math.isnan(as_float):
        return _MISSING_ID
    if math.isinf(as_float):
        return str(fallback)
    as_int = int(as_float)
    return as_int if as_float == as_int else as_float


def _canonical_index(index: pd.Index) -> pd.Index:
    """Apply :func:`_canonical_id` to every label, as an object-dtype Index."""
    labels = [_canonical_id(v) for v in index]
    return cast(pd.Index, pd.Index(labels, name=index.name, dtype=object))


def _unique_by_canonical_id(indexed: pd.DataFrame) -> pd.DataFrame:
    """Canonicalise a just-set index, drop missing ids, keep the first per id.

    Shared by every indexing path so they stay consistent. Rows whose canonical
    id is :data:`_MISSING_ID` (a missing/blank id) are dropped. Collapsing
    cross-type duplicates (e.g. ``1`` and ``"1"``) keeps the first, which
    upstream has already ordered so the best-scoring row wins.
    """
    indexed = indexed.copy()
    indexed.index = _canonical_index(indexed.index)
    indexed = indexed[[label is not _MISSING_ID for label in indexed.index]]
    return cast(pd.DataFrame, indexed[~indexed.index.duplicated(keep="first")])


def _index_by_feature(
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
        # Rank numerically, not lexically ("low" > "high"). Non-numeric, blank,
        # and non-finite (±inf) scores all coerce to NaN and sort last, so a real
        # finite winner always beats them.
        score = pd.to_numeric(df[score_col], errors="coerce")
        score = score.where(np.isfinite(score))
        ordered = df.iloc[
            score.reset_index(drop=True)
            .sort_values(ascending=False, na_position="last", kind="stable")
            .index
        ]
    else:
        ordered = df
    # Deduplicate on the *canonical* id (not the raw one), keeping the first per
    # id — which the score sort above has made the best-scoring. Canonicalising
    # first is what lets 1 and "1" collapse into one feature.
    indexed = ordered.set_index(id_col)
    return _unique_by_canonical_id(indexed)


def _aligned_column(
    source: pd.DataFrame | None, feature_ids: pd.Index, src_col: str
) -> pd.Series:
    """Align ``source[src_col]`` to ``feature_ids`` (missing rows/table -> NA)."""
    if source is not None and src_col in source.columns:
        return source[src_col].reindex(feature_ids)
    return cast(pd.Series, pd.Series(pd.NA, index=feature_ids))


def _to_text(value: Any) -> str:
    """Render a text-column cell as a string: blank -> "".

    Real exports supply these columns as strings already. The one subtlety is a
    stray numeric: a bare ``str(1.0)`` would leak ``"1.0"``, so an integer-valued
    float is rendered without the trailing ``.0`` (``1.0`` -> ``"1"``), matching
    how the same value reads as a feature id. Non-finite floats render as "".
    """
    if value is None or (isinstance(value, float) and not math.isfinite(value)):
        return ""
    try:
        if pd.isna(value):
            return ""
    except (TypeError, ValueError):
        pass
    if isinstance(value, float) and value.is_integer():
        return str(int(value))
    return str(value)


def _sort_key(value: Any) -> tuple[int, str, Any]:
    """Order ids so mixed/incomparable types never raise.

    Numbers sort (numerically) before non-numbers, each in natural order. The
    numeric value is kept as-is rather than cast to ``float``, so two integer ids
    that differ only beyond 2**53 still order correctly. The tuple shape is
    uniform across groups, but number-vs-string comparisons are gated by the
    first element, so ``int``/``float`` and ``str`` are never compared directly.
    """
    if _is_bool(value):  # bool is an int subclass; order as text
        return (1, str(value), None)
    if isinstance(value, (int, float)):
        return (0, "", value)
    return (1, str(value), None)


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
        annotations = _records_to_frame(data.annotations)
        feature_table = _records_to_frame(data.feature_table)
        library = _records_to_frame(data.library)

        records = self._build_annotations(annotations, feature_table, library)
        return GNPSResult(annotations=records)

    # -- internal -----------------------------------------------------------

    def _build_annotations(
        self,
        annotations: pd.DataFrame,
        feature_table: pd.DataFrame,
        library: pd.DataFrame,
    ) -> Records:
        # Scored tables keep the highest-scoring row per feature; the feature
        # table is unscored (score_col=None -> first row per id wins).
        ann = _index_by_feature(annotations, _ANNOTATION_ID_CANDIDATES, "score")
        lib = _index_by_feature(library, _LIBRARY_ID_CANDIDATES, "MQScore")
        feat = _index_by_feature(feature_table, _FEATURE_TABLE_ID_CANDIDATES)

        # Output row set = union of features annotated by either source (see the
        # module docstring); a feature-table-only feature contributes no id here.
        ids: set[Any] = set()
        if ann is not None:
            ids |= set(ann.index)
        if lib is not None:
            ids |= set(lib.index)
        if not ids:
            return []
        # ordered_ids is the source of truth for the emitted Feature ID values
        # and their exact types; the object-dtype index below is only for
        # aligning the source columns and never re-coerces them.
        ordered_ids = sorted(ids, key=_sort_key)
        feature_ids = pd.Index(ordered_ids, name="Feature ID", dtype=object)

        # Each column is a source column aligned to the feature ids; a missing
        # table/column fills with NA.
        out = pd.DataFrame(index=feature_ids)
        for source, mapping in (
            (feat, _FEATURE_TABLE_FIELDS),
            (ann, _ANNOTATION_FIELDS),
            (lib, _LIBRARY_FIELDS),
        ):
            for out_col, src_col in mapping.items():
                out[out_col] = _aligned_column(source, feature_ids, src_col)

        return self._to_records(ordered_ids, out)

    @staticmethod
    def _to_records(feature_ids: list[Any], out: pd.DataFrame) -> Records:
        """Convert the assembled frame to records with tidy value types.

        ``feature_ids`` (same order as ``out``'s rows) supplies the ``Feature ID``
        column directly, so each id keeps its exact Python type. Numeric columns
        become ``float``/``None``, text columns ``str``/``""``, per column.
        """
        # Row-position index so columns line up with feature_ids by position,
        # independent of the object-dtype feature id index.
        tidy = pd.DataFrame(index=range(len(feature_ids)))
        for col in out.columns:
            series = out[col].reset_index(drop=True)
            if col in _NUMERIC_COLUMNS:
                # Blank/non-numeric/±inf all -> None (object dtype so None
                # survives to_dict); only finite floats are emitted.
                nums = pd.to_numeric(series, errors="coerce")
                tidy[col] = nums.astype(object).where(np.isfinite(nums), None)
            else:
                tidy[col] = [_to_text(v) for v in series]

        # Prepend the Feature ID (from the id list, exact type) and rebuild rows
        # with str keys so the result matches the declared Records type.
        raw = tidy.to_dict("records")
        return [
            {"Feature ID": fid, **{str(k): v for k, v in row.items()}}
            for fid, row in zip(feature_ids, raw, strict=True)
        ]
