## Building the documentation

1. Install the sphinx requirements:

```console
uv sync --group docs
```

2. Build the documentation as an html site:

```console
make docs
```
or
```console
uv run --group=docs --directory=docs bash -c "make clean && make html"
```

3. Open the built docs (`build/html/index.html`) in a browser

```console
open build/html/index.html
```

4. Run step 2 again and refresh the browser window to rebuild the docs
after making edits.

## Editing the docs

`source/index.rst` will define the "home" page of the documentation site.
Any other pages should be added in the `source/` directory and referenced
using the toctree and cross-references (see the sphinx documentation) to
make them findable by users.
