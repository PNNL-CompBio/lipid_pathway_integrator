"""
lipi/identifiers/ingestion/lipid_maps/__init__.py
Dylan Ross (dylan.ross@pnnl.gov)

    Sub-package for mapping common lipid identifiers to the LIPID MAPS ontology
"""


import os
import pickle
from typing import Tuple, Mapping, Union, Dict, Any
from collections.abc import Iterable
from logging import getLogger

from lipi.refpath.manager import RefPathManager
from lipi.identifiers.ingestion.lipid_maps._remote import fetch_record


#===============================================================================
# Constants

# module-level logger
_LOGGER = getLogger(__name__)

# set some parameters for enumerating the LMSD 
# key is category
# value is how many classes are in the category (indexed from 0)
# for subclasses, we will just pick a suitably high number (18) for all and enumerate up to that
# but only enumerate for subclass indices >0
_CATEGORIES_CLASSES = {
    "FA": 14,
    "GL": 8,
    "GP": 26,
    "SP": 10,
    "ST": 6,
    "PR": 5,
    "SL": 6,
    "PK": 16
}


#===============================================================================
# Helper Function

def _parse_identifiers_from_result(
    result: Dict[str, Any]
) -> Iterable[Tuple[str, Mapping[str, Union[str, int, float, bool]]]] :
    name = n if (n := result.get("name")) is not None else result["sys_name"]
    identifiers = [
        # {"name": name, "tags": tags}
        {"name": name, "tags": {"display": True}},
        {
            "name": (lm_id := result["lm_id"]), 
            "tags": {
                "source": "LIPID MAPS",
                "classification": {
                    "category": result["core"],
                    "class": c if (c := result["main_class"]) is not None else f"? [{lm_id[2:6]}]",
                    "subclass": sc if (sc := result["sub_class"]) is not None else f"? [{lm_id[2:8]}]"
                },
                "is_lipid_group": True
            }
        },
    ]
    # check for ones that are not always there
    if (sn := result.get("sys_name")) is not None and sn != name:
        identifiers.append({"name": sn})
    if (syns := result.get("synonyms")) is not None:
        for syn in syns.split("; "):
            identifiers.append({"name": syn})
    if (kegg := result.get("kegg_id")) is not None:
        identifiers.append({"name": "cpd:" + kegg, "tags": {"source": "KEGG"}})
    if (chebi := result.get("chebi_id")) is not None:
        identifiers.append({"name": int(chebi), "tags": {"source": "ChEBI"}})
    if (pubchem := result.get("pubchem_cid")) is not None:
        identifiers.append({"name": int(pubchem), "tags": {"source": "PubChem"}})
    # the above was already implemented using list of dicts, I would rather not re-write it 
    # so before returning convert the list to a list of tuples
    return [(d["name"], d.get("tags", {})) for d in identifiers]


#===============================================================================
# Main Ingestion Function

def ingest_lipid_maps_classes(
    lmsd_cache_file: str,
    manager: RefPathManager
) -> None :
    """
    """
    _LOGGER.info(
        "%(module)s.%(func)s: ingesting lipid class identifiers from Lipid Maps (LMSD)",
        {
            "module": __name__,
            "func": "fetch_record",
        }
    )

    # make sure RefPathManager is connected
    if not manager.connected:
        raise RuntimeError("reference pathway databse manager is not connected")

    # load the lmsd cache (pickle)
    if not os.path.isfile(lmsd_cache_file):
        # make an empty cache if cache file not found
        lmsd_cache = {}
        _LOGGER.debug(
            "%(module)s.%(func)s: LMSD cache file %(lmsd_cache_file)s not found, creating empty cache",
            {
                "module": __name__,
                "func": "fetch_record",
                "lmsd_cache_file": lmsd_cache_file
            }
        )
    else:
        # otherwise load the existing ones
        with open(lmsd_cache_file, "rb") as pf: 
            lmsd_cache = pickle.load(pf)
            _LOGGER.debug(
                "%(module)s.%(func)s: loaded LMSD cache file %(lmsd_cache_file)s",
                {
                    "module": __name__,
                    "func": "fetch_record",
                    "lmsd_cache_file": lmsd_cache_file
                }
            )

    # enumerate some of the LMSD to get lipid class info
    for cat, n_classes in _CATEGORIES_CLASSES.items():
        for cls in range(n_classes + 1):
            for subcls in range(19):
                if (result := fetch_record(f"LM{cat}{cls:02d}{subcls:02d}0000", lmsd_cache)) != {}:
                    # parse identifiers from each result and add them as an identifier group
                    manager.entities.insert_or_update_from_identifiers(
                        _parse_identifiers_from_result(result)
                    )
                    
    # save the cache in case anything new was added 
    # TODO: Make the cache a more formal data structure with a flag indicating whether it has been updated. 
    #       This way we could check the flag and not bother re-saving if there have been no updates
    with open(lmsd_cache_file, "wb") as pf:
        pickle.dump(lmsd_cache, pf)
        _LOGGER.debug(
            "%(module)s.%(func)s: saved LMSD cache file %(lmsd_cache_file)s",
            {
                "module": __name__,
                "func": "fetch_record",
                "lmsd_cache_file": lmsd_cache_file
            }
        )