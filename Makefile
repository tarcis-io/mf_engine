# Developer convenience targets. Install the dev tooling first:
#     pip install -e ".[data,dev]"
#
# `make check` runs lint, type-check, and tests together.

.PHONY: help install test cov lint typecheck check

help:  ## Show this help
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) \
		| awk 'BEGIN {FS = ":.*?## "}; {printf "  %-12s %s\n", $$1, $$2}'

install:  ## Editable install with data + dev extras
	pip install -e ".[data,dev]"

test:  ## Run the test suite
	pytest

cov:  ## Run the test suite with a coverage report
	pytest --cov=mf_engine

lint:  ## Lint with ruff
	ruff check .

typecheck:  ## Type-check with mypy
	mypy

check: lint typecheck test  ## Run lint, type-check, and tests
