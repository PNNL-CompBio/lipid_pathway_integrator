# docs/conf.py
import os
import sys
from pathlib import Path
import tomllib 


# Add source code to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

import lipi


# Read pyproject.toml
pyproject_path = Path(__file__).parent.parent / "pyproject.toml"
with open(pyproject_path, "rb") as f:
    pyproject = tomllib.load(f)

# Extract project metadata
project_data = pyproject["project"]
project = "lipi"
author = ", ".join(project_data.get("authors", [{}])[0].get("name", ""))
release = lipi.__version__
description = project_data.get("description", "")

# Extract URLs
urls = project_data.get("urls", {})
repository_url = urls.get("Repository", "")

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.intersphinx",
    "sphinx.ext.napoleon",
    "numpydoc",
    "sphinx_autodoc_typehints",
]


# Numpydoc configuration
numpydoc_show_class_members = True
numpydoc_use_plots = True


# Autosummary configuration
autosummary_generate = True
autosummary_generate_overwrite = True  # Regenerate even if files exist

# Autodoc configuration
autodoc_default_options = {
    "members": True,
    "member-order": "bysource",
    "undoc-members": False,
    "show-inheritance": True,
}

autodoc_typehints = "description"
autodoc_member_order = "bysource"

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

html_theme = "sphinx_rtd_theme"
html_theme_options = {
    "logo_only": False,
    #"display_version": True,
    "prev_next_buttons_location": "bottom",
    "style_external_links": False,
    #"vcs_url": repository_url,
}
html_static_path = ["_static"]

# Intersphinx mapping for external documentation
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
}