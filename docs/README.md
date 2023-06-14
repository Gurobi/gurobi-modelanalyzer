## Building the documentation

1. Install the sphinx requirements:

```console
python -m pip install -r requirements.txt
```

2. Build the documentation:

```console
make html
```

3. Open the main docs page (`build/html/index.html`) in a browser

```
open build/html/index.html
```

## Editing the docs

`source/index.rst` will define the "home" page of the documentation site.
Any other pages should be added in the `source/` directory and referenced
using the toctree and cross-references (see the sphinx documentation) to
make them findable by users.
