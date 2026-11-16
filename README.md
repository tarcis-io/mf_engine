# mf_engine

The **MassFlow engine** — a reusable, importable Python package of *processors*
for MZmine + GNPS metabolomics batches.

## Install

From the package directory (editable install for development):

```bash
pip install -e .
```

With the optional data-processing dependencies (pandas / openpyxl), used by
concrete data processors:

```bash
pip install -e ".[data]"
```

## Usage

```python
from mf_engine import Processor


class UppercaseProcessor(Processor):
    def process(self, data: str) -> str:
        return data.upper()


proc = UppercaseProcessor()
proc.process("hello")  # -> "HELLO"
```

`Processor` is abstract: it cannot be instantiated directly, and every subclass
must implement `process`.

### GNPSProcessor

`GNPSProcessor` merges MZmine + GNPS tables into a per-feature **annotation
table**. It is I/O-free: inputs and outputs are plain records (`list[dict]`),
so it works with loaded CSV rows, database result sets, or JSON — anything that
round-trips through `DataFrame.to_dict("records")`. It returns only the
annotations; any summaries or interpretation are left to the caller.

```python
from mf_engine.processors.gnps import GNPSProcessor, GNPSInput

result = GNPSProcessor().process(
    GNPSInput(
        annotations=annotations_records,      # MZmine annotations rows
        feature_table=feature_table_records,  # GNPS feature table rows
        library=library_records,              # GNPS library match rows
    )
)
result.annotations  # list[dict], one per annotated feature
```

The records boundary means loading real exports is just reading each file and
calling `to_dict("records")`:

```python
import pandas as pd
from mf_engine.processors.gnps import GNPSProcessor, GNPSInput

result = GNPSProcessor().process(
    GNPSInput(
        annotations=pd.read_csv("mzmine_annotations.csv").to_dict("records"),
        feature_table=pd.read_csv("gnps_feature_table.csv").to_dict("records"),
        library=pd.read_csv("gnps_merged_results.tsv", sep="\t").to_dict("records"),
    )
)
```

The output schema — the exact keys, in order, of every result record — is
available as `mf_engine.processors.ANNOTATION_COLUMNS`, handy for building a
header row or validating output.

Key guarantees (the `GNPSProcessor` module docstring has the full detail):

- Every output field is a **direct passthrough** of an input column — no derived
  or interpretive fields (best-of picks, confidence, compound-class flags). Those
  belong to the caller.
- The output is the **union of features annotated by either source**; a feature
  present only in the feature table is omitted. One-source features are returned
  with the other source's columns blank.
- Features join by feature id matched *by value regardless of type* (`1`, `1.0`,
  `"1"` are one feature), with integer ids kept at full precision.
- Per source, the highest-scoring row per feature wins; non-numeric, blank, and
  non-finite scores count as missing. Ties and all-missing fall back to input
  order, so output is deterministic.
- Messy inputs are tolerated, not fatal: a blank id drops just that row; a
  missing column or empty table yields blanks rather than raising.

`GNPSProcessor` needs pandas — install the `data` extra. The base `Processor`
and `import mf_engine` stay dependency-free; pandas is only loaded when a GNPS
symbol is used.

## Version

`mf_engine.__version__` is read from the installed package metadata (the single
source of truth is `pyproject.toml`). When the package is imported from a source
tree that has **not** been installed, the metadata is unavailable and the value
falls back to `"0.0.0"` — so `__version__` matches `pyproject.toml` only once
the package is installed (`pip install -e .`). In practice this means a bare
`src` checkout that was never installed reports `0.0.0`; install the package to
get the real version.

## Development

Install the dev tooling and run the checks:

```bash
pip install -e ".[data,dev]"
pytest            # tests
ruff check .      # lint
mypy              # type-check (config in pyproject.toml drives what's checked)
```

Or run all three at once with the bundled `Makefile`:

```bash
make check        # lint + type-check + tests
```

## Layout

```
mf_engine/                    # repository root
├── pyproject.toml
├── README.md
├── LICENSE
├── Makefile
├── src/
│   └── mf_engine/            # the importable package
│       ├── __init__.py       # exposes `Processor`, `__version__`
│       ├── py.typed          # PEP 561 typing marker
│       └── processors/
│           ├── __init__.py   # exposes `Processor`; lazily exposes GNPS* + ANNOTATION_COLUMNS
│           ├── base.py       # the Processor ABC (stdlib only)
│           └── gnps.py       # GNPSProcessor (needs pandas)
└── tests/
    ├── test_processor.py     # the Processor ABC + package surface
    ├── test_gnps.py          # GNPSProcessor unit tests (hand-built cases)
    ├── test_gnps_fixtures.py # GNPSProcessor integration tests (real exports)
    └── fixtures/             # sample MZmine + GNPS export files
```

This uses the **src layout** so tests always run against the installed package,
not a stray copy on the path.
