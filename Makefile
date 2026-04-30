.PHONY: install fixtures test

install:
	pip install -e python/

fixtures:
	python tests/fixtures/generate.py

test:
	pytest tests/
