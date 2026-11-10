"""Sphinx configuration for the orthoSynAssign documentation."""

from __future__ import annotations

from importlib.metadata import PackageNotFoundError, version

project = "orthoSynAssign"
author = "Cheng-Hung Tsai, Carolina G. Piña Páez, Jason E. Stajich"
copyright = "2026, Cheng-Hung Tsai, Carolina G. Piña Páez and Jason E. Stajich"

try:
    release = version("orthosynassign")
except PackageNotFoundError:
    release = "unknown"
version = ".".join(release.split(".")[:2])

extensions = [
    "myst_parser",
    "sphinx_copybutton",
]

source_suffix = {".md": "markdown", ".rst": "restructuredtext"}
root_doc = "index"
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

myst_enable_extensions = ["colon_fence", "deflist", "tasklist", "attrs_inline"]
myst_heading_anchors = 3

html_theme = "furo"
html_title = f"orthoSynAssign {release}"
html_static_path = ["_static"]
html_theme_options = {
    "source_repository": "https://github.com/stajichlab/orthoSynAssign/",
    "source_branch": "main",
    "source_directory": "docs/",
}
