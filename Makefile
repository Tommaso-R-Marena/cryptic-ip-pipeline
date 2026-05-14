.PHONY: install check test test-cov validation download-validation \
        download-yeast index-yeast screen-yeast-fast screen-yeast-full \
        report-yeast clean-small docker lint

PYTHON ?= python
PIP    ?= pip

install:
	$(PIP) install -e .[dev]

check:
	crypticip check-env

test:
	$(PYTHON) -m pytest -q

test-cov:
	$(PYTHON) -m pytest --cov=crypticip --cov-report=term-missing

lint:
	ruff check crypticip tests || true

download-validation:
	crypticip download-validation

validation:
	crypticip validate --config config/validation.yaml

download-yeast:
	crypticip download-proteome --organism yeast --resume

index-yeast:
	crypticip index-proteome --organism yeast

screen-yeast-fast:
	crypticip screen --organism yeast --mode fast --workers 4 --resume

screen-yeast-full:
	crypticip screen --organism yeast --mode full --workers 4 --resume

report-yeast:
	crypticip report --organism yeast

clean-small:
	rm -rf results/screening/_tmp results/reports/_tmp __pycache__ */__pycache__ */*/__pycache__ .pytest_cache

docker:
	docker build -t crypticip:latest .
