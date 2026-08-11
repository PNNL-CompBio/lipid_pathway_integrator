"""
lipi/__init__.py
Dylan Ross (dylan.ross@pnnl.gov)

    Lipid Pathway Integrator
"""


__version__ = "1.0.2"

# Make submodules discoverable (for docs, IDE autocomplete, etc.)
from . import identifiers
from . import mapping
from . import network
from . import refpath
from . import translation


# __all__ = [
#     "identifiers",
#     "mapping",
#     "network", 
#     "refpath",
#     "translation",
# ]