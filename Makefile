.PHONY : format lint precommit test tox build docs

init:
	uv sync --locked
	uv run pre-commit install

format:
	uvx ruff format

lint:
	uvx ruff check

precommit:
	uv run pre-commit run --all-files

test:
	uv run --no-dev pytest ./tests

tox:
	uv run tox

build:
	uv build

docs:
	uv run --group=docs --directory=docs bash -c "make clean && make html"
