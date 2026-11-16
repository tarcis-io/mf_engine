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

Every output field is a **direct passthrough** of an input column (feature id,
m/z, RT, per-source score/compound/adduct/formula, and the GNPS taxonomy
fields) — no derived or interpretive fields such as best-of picks, confidence
levels, compound-class flags, or group comparisons. Those belong to the caller.
Features are joined across the three tables by the value of their shared
feature id. Ids are matched by value regardless of type, so the same feature
lines up whether a source supplies it as `1`, `1.0`, or `"1"` (as a plain CSV
reader would), and integer ids keep full precision rather than rounding through
`float`. A feature annotated by only one source is still returned, with the
other source's columns left blank. The output is the **union of features
annotated by either source**; a feature present only in the feature table
(annotated by neither) is omitted.

When a feature has several rows in a source, the one with the highest
`score` / `MQScore` is kept. Scores are compared numerically, so a non-numeric
or blank score is treated as missing and loses to any real number; when no
score for a feature is numeric, the first row in input order is kept (the
result is always deterministic). Messy inputs are tolerated rather than fatal:
a blank/`NaN` feature id drops just that row, and a missing column or an empty
table contributes blanks for the affected fields instead of raising.

`GNPSProcessor` needs pandas — install the `data` extra. The base `Processor`
and `import mf_engine` stay dependency-free; pandas is only loaded when a GNPS
symbol is used.

## Version

`mf_engine.__version__` is read from the installed package metadata (the single
source of truth is `pyproject.toml`). When the package is imported from a source
tree that has **not** been installed, the metadata is unavailable and the value
falls back to `"0.0.0"` — so `__version__` matches `pyproject.toml` only once
the package is installed (`pip install -e .`).

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
mf_engine/
├── pyproject.toml
├── README.md
├── LICENSE
├── src/
│   └── mf_engine/            # the importable package
│       ├── __init__.py       # exposes `Processor`, `__version__`
│       ├── py.typed          # PEP 561 typing marker
│       └── processors/
│           ├── __init__.py   # exposes `Processor`; lazily exposes GNPS*
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
