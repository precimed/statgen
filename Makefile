.PHONY: install fixtures test test-python test-octave test-matlab

install:
	pip install -e python/

fixtures:
	python tests/fixtures/generate.py

test:
	pytest tests/

test-python:
	pytest tests/ -m "not octave"

test-octave:
	pytest tests/ -m octave

test-matlab:
	STATGEN_MATLAB=1 pytest tests/ -m octave
