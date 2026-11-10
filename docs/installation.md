# Installation

orthoSynAssign can be installed via **Conda/Bioconda** (recommended), via **Pixi**, or from **source** (regular install or editable developer install).

---

## 1. Conda Installation (Recommended)

The simplest way to install orthoSynAssign and its dependencies is through [Bioconda](https://anaconda.org/bioconda/orthosynassign). Create a dedicated conda environment first, then install into it:

```sh
conda create -n orthosynassign
conda activate orthosynassign
conda install bioconda::orthosynassign
```

---

## 2. Pixi Installation

orthoSynAssign defines a [Pixi](https://pixi.sh) workspace in `pyproject.toml` (`[tool.pixi.*]`) for `linux-64`, `osx-64`, and `osx-arm64`.

### Global CLI Installation

```sh
pixi global install -c bioconda orthosynassign
```

### Adding to a Pixi Project

```sh
pixi init my_orthosynassign_project
cd my_orthosynassign_project
pixi add -c bioconda orthosynassign

# Run orthosynassign commands inside the pixi environment
pixi run orthosynassign --help
```

---

## 3. Install from Source

### Clone the repository

Via SSH:

```sh
git clone git@github.com:stajichlab/orthoSynAssign.git
```

Via HTTPS:

```sh
git clone https://github.com/stajichlab/orthoSynAssign.git
```

### Install

Navigate to the project directory and install the package:

```sh
cd orthoSynAssign
pip install .
```

The package bundles a Rust extension module (built via [maturin](https://www.maturin.rs/)) that implements the core synteny-comparison engine, so a Rust toolchain (`>= 1.75`) is required if a prebuilt wheel is not available for your platform.

---

## 4. Install from Latest Stable Release

You can download the source code from the [latest stable release](https://github.com/stajichlab/orthoSynAssign/releases/latest) on GitHub, decompress it, and install via pip:

```sh
tar -zxvf orthoSynAssign-X.X.X.tar.gz
cd orthoSynAssign-X.X.X
pip install .
```

> **Note:** Building from a release tarball still requires a Rust toolchain (`>= 1.75`) if a prebuilt wheel is not available for your platform.

---

## 5. Developer Setup

Developers should clone the project directly and install with the `dev` extra, then set up pre-commit hooks before making commits.

### Option A: Using Pixi (Recommended)

The `pyproject.toml` also defines a Pixi `dev` environment/feature with an editable install and [ruff](https://docs.astral.sh/ruff/) for linting:

```sh
git clone https://github.com/stajichlab/orthoSynAssign.git
cd orthoSynAssign

# Install the dev environment
pixi install --environment dev

# Activate the development environment shell
pixi shell -e dev

# Set up pre-commit hooks
pre-commit install
```

#### Key Pixi Tasks (defined in `pyproject.toml`)

The tasks belong to the `dev` feature, so run them from a `pixi shell -e dev` or with `pixi run -e dev <task>`:

- **Build the Rust extension in place**: `pixi run -e dev build` (runs `maturin develop`)
- **Run tests**: `pixi run -e dev test` (runs `pytest`, which reports coverage for the `orthosynassign` package per `[tool.coverage.run]`)
- **Build the documentation**: `pixi run -e dev docs` (writes `docs/_build/html`)
- **Preview the documentation with live reload**: `pixi run -e dev docs-serve` (runs `sphinx-autobuild`)
- **Lint**: `pixi run -e dev ruff check .` (uses the `ruff` dev dependency; line length 130, rule sets `E`, `F`, `I`)

### Option B: Using pip

```sh
git clone https://github.com/stajichlab/orthoSynAssign.git
cd orthoSynAssign
pip install -e ".[dev]"
pre-commit install
```

The `dev` extra installs `maturin`, `pre-commit`, `pytest`, and `pytest-cov`. The Pixi `dev` environment additionally
installs the `docs` extra (`sphinx`, `myst-parser`, `furo`, `sphinx-copybutton`, `sphinx-autobuild`); with pip, install it
explicitly with `pip install -e ".[dev,docs]"` if you want to build the documentation.

> **Tip:** Even with pip, it is recommended to create an isolated virtual environment (e.g. via `conda create -n orthosynassign-dev python=3.13` or `python -m venv`) before installing, rather than installing into your system/base Python.

---

## System Requirements

- **[Python]**: `>= 3.9, < 3.15`
- **[pyGenomeViz]**: `>= 1.6.0` (used by `orthosynassign-vis`)
- **[Rust]**: `>= 1.75`, used to build the core synteny engine (only needed when installing from source without a prebuilt wheel)

[Python]: https://www.python.org/
[pyGenomeViz]: https://github.com/moshi4/pyGenomeViz
[Rust]: https://rust-lang.org/
