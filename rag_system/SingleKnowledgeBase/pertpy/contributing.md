# Contributing guide

Scanpy provides extensive [developer documentation][scanpy developer guide], most of which applies to this repo, too.
This document will not reproduce the entire content from there. Instead, it aims at summarizing the most important
information to get you started on contributing.

We assume that you are already familiar with git and with making pull requests on GitHub. If not, please refer
to the [scanpy developer guide][].

## Installing dev dependencies

In addition to the packages needed to _use_ this package, you need additional python packages to _run tests_ and _build
the documentation_. It's easy to install them using `pip`:

```bash
git clone https://github.com/scverse/pertpy.git
cd pertpy
pip install -e ".[dev,test,doc]"
```

```{note}
If you're working on an Apple Silicon machine, the installation is slightly more complex. In that case, follow
the steps described in the Apple Silicon section of the installation guide and replace the last two steps described there with the code above.
```

## Code-style

This project uses [pre-commit][] to enforce consistent code-styles. On every commit, pre-commit checks will either
automatically fix issues with the code, or raise an error message.

To enable pre-commit locally, simply run

```bash
pre-commit install
```

in the root of the repository. Pre-commit will automatically download all dependencies when it is run for the first time.

Alternatively, you can rely on the [pre-commit.ci][] service enabled on GitHub. If you didn't run `pre-commit` before
pushing changes to GitHub it will automatically commit fixes to your pull request, or show an error message.

If pre-commit.ci added a commit on a branch you still have been working on locally, simply use

```bash
git pull --rebase
```

to integrate the changes into yours.
While the [pre-commit.ci][] is useful, we strongly encourage installing and running pre-commit locally first to understand its usage.

Finally, most editors have an _autoformat on save_ feature. Consider enabling this option for [black][black-editors]
and [prettier][prettier-editors].

[black-editors]: https://black.readthedocs.io/en/stable/integrations/editors.html
[prettier-editors]: https://prettier.io/docs/en/editors.html

## Writing tests

```{note}
Remember to first install the package with `pip install '-e[dev,test]'`
```

This package uses the [pytest][] for automated testing. Please [write tests][scanpy-test-docs] for every function added
to the package.

Most IDEs integrate with pytest and provide a GUI to run tests. Alternatively, you can run all tests from the
command line by executing

```bash
pytest
```

in the root of the repository. Continuous integration will automatically run the tests on all pull requests.

[scanpy-test-docs]: https://scanpy.readthedocs.io/en/latest/dev/testing.html#writing-tests

## Publishing a release

### Updating the version number

Before making a release, you need to update the version number. Please adhere to [Semantic Versioning][semver], in brief

> Given a version number MAJOR.MINOR.PATCH, increment the:
>
> 1.  MAJOR version when you make incompatible API changes,
> 2.  MINOR version when you add functionality in a backwards compatible manner, and
> 3.  PATCH version when you make backwards compatible bug fixes.
>
> Additional labels for pre-release and build metadata are available as extensions to the MAJOR.MINOR.PATCH format.

Once you are done, run

```
git push --tags
```

to publish the created tag on GitHub.

### Building and publishing the package on PyPI

Python packages are not distributed as source code, but as _distributions_. The most common distribution format is the so-called _wheel_. To build a _wheel_, run

```bash
python -m build
```

This command creates a _source archive_ and a _wheel_, which are required for publishing your package to [PyPI][]. These files are created directly in the root of the repository.

Before uploading them to [PyPI][] you can check that your _distribution_ is valid by running:

```bash
twine check dist/*
```

and finally publishing it with:

```bash
twine upload dist/*
```

Provide your username and password when requested and then go check out your package on [PyPI][]!

For more information, follow the [Python packaging tutorial][].

It is possible to automate this with GitHub actions, see also [this feature request][pypi-feature-request]
in the cookiecutter-scverse template.

[python packaging tutorial]: https://packaging.python.org/en/latest/tutorials/packaging-projects/#generating-distribution-archives
[pypi-feature-request]: https://github.com/scverse/cookiecutter-scverse/issues/88

## Writing documentation

Please write documentation for new or changed features and use-cases. This project uses [sphinx][] with the following features:

-   the [myst][] extension allows to write documentation in markdown/Markedly Structured Text
-   Google-style docstrings
-   Jupyter notebooks as tutorials through [myst-nb][] (See [Tutorials with myst-nb](#tutorials-with-myst-nb-and-jupyter-notebooks))
-   [Sphinx autodoc typehints][], to automatically reference annotated input and output types
-   Citations (like {cite:p}`Virshup_2023`) can be included with [sphinxcontrib-bibtex](https://sphinxcontrib-bibtex.readthedocs.io/)

See the [scanpy developer docs](https://scanpy.readthedocs.io/en/latest/dev/documentation.html) for more information
on how to write documentation.

### Tutorials with myst-nb and jupyter notebooks

The documentation is set-up to render jupyter notebooks stored in the `docs/tutorials` directory using [myst-nb][].
Currently, only notebooks in `.ipynb` format are supported that will be included with both their input and output cells.

These notebooks come from [pert-tutorials](https://github.com/scverse/pertpy-tutorials) which is a git submodule of pertpy.

#### Hints

-   If you refer to objects from other packages, please add an entry to `intersphinx_mapping` in `docs/conf.py`. Only
    if you do so can sphinx automatically create a link to the external documentation.
-   If building the documentation fails because of a missing link that is outside your control, you can add an entry to
    the `nitpick_ignore` list in `docs/conf.py`

#### Building the docs locally

```bash
cd docs
make html
open _build/html/index.html
```

<!-- Links -->

[scanpy developer guide]: https://scanpy.readthedocs.io/en/latest/dev/index.html
[cookiecutter-scverse-instance]: https://cookiecutter-scverse-instance.readthedocs.io/en/latest/template_usage.html
[github quickstart guide]: https://docs.github.com/en/get-started/quickstart/create-a-repo?tool=webui
[codecov]: https://about.codecov.io/sign-up/
[codecov docs]: https://docs.codecov.com/docs
[codecov bot]: https://docs.codecov.com/docs/team-bot
[codecov app]: https://github.com/apps/codecov
[pre-commit.ci]: https://pre-commit.ci/
[readthedocs.org]: https://readthedocs.org/
[myst-nb]: https://myst-nb.readthedocs.io/en/latest/
[jupytext]: https://jupytext.readthedocs.io/en/latest/
[pre-commit]: https://pre-commit.com/
[anndata]: https://github.com/scverse/anndata
[mudata]: https://github.com/scverse/mudata
[pytest]: https://docs.pytest.org/
[semver]: https://semver.org/
[sphinx]: https://www.sphinx-doc.org/en/master/
[myst]: https://myst-parser.readthedocs.io/en/latest/intro.html
[numpydoc-napoleon]: https://www.sphinx-doc.org/en/master/usage/extensions/napoleon.html
[numpydoc]: https://numpydoc.readthedocs.io/en/latest/format.html
[sphinx autodoc typehints]: https://github.com/tox-dev/sphinx-autodoc-typehints
[pypi]: https://pypi.org/
