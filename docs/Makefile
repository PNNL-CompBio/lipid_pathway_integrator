.PHONY: docs docs-clean docs-serve

docs:
	uv run sphinx-build -W -b html -d _build/doctrees . _build/html

docs-clean:
	rm -rf _build

docs-serve: docs
	python -m http.server --directory _build/html 8000