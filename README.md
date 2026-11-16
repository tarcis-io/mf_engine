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

In short: every output field is a direct passthrough (no derived or best-of
fields — those are the caller's job); the rows are the union of features
annotated by either source; features join by id matched regardless of type; the
highest-scoring row per feature wins, deterministically; and messy inputs yield
blanks rather than raising. The `GNPSProcessor` module docstring is the
authoritative reference for all of these.

`GNPSProcessor` needs pandas — install the `data` extra. The base `Processor`
and `import mf_engine` stay dependency-free; pandas is only loaded when a GNPS
symbol is used.

## Version

`mf_engine.__version__` is read from the installed package metadata (the single
source of truth is `pyproject.toml`). A source tree that was never installed has
no metadata, so it falls back to `"0.0.0"`; install the package
(`pip install -e .`) to get the real version.

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
    └── fixtures/             # sample MZmine annotations + GNPS feature/library tables
```

This uses the **src layout** so tests always run against the installed package,
not a stray copy on the path.
