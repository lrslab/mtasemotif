PYTHON ?= python

.PHONY: lint test sdist-check package-check

lint:
	$(PYTHON) -m ruff check .

test:
	$(PYTHON) -m pytest -q

sdist-check:
	rm -rf .sdist-check
	$(PYTHON) setup.py sdist --dist-dir .sdist-check
	$(PYTHON) -m pytest -q tests/test_no_rebase_bundling.py tests/test_packaging_release.py

package-check:
	$(MAKE) sdist-check
	rm -rf build dist .sdist-check
	$(PYTHON) -m build
	$(PYTHON) -m twine check dist/*
	$(PYTHON) -m pytest -q tests/test_no_rebase_bundling.py tests/test_packaging_release.py
