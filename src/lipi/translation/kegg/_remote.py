"""
lipi/translation/kegg/_remote.py
Dylan Ross (dylan.ross@pnnl.gov)

    Module for interacting with the KEGG REST API
"""


from typing import Dict, Literal
import time
from logging import getLogger
import uuid
import sys

import requests


#===============================================================================
# Constants

# module-level logger
_LOGGER = getLogger(__name__)

# set a delay (in seconds) before sending HTTP requests to limit server load
_HTTP_REQUEST_DELAY = 0.5


#===============================================================================
# Type annotation for convenience

type KeggCache = Dict[Literal["pathway", "gene", "compound"], Dict[str, str]]


#===============================================================================
# convenience function for empty kegg cache

def get_empty_kegg_cache() -> KeggCache :
    return {
        "pathway": {},
        "gene": {},
        "compound": {}
    }


#===============================================================================
# Main function for fetching info 

def fetch_record(
    identifier: str, 
    record_type: Literal["gene", "compound", "pathway"], 
    cache: Dict[Literal["gene", "compound", "pathway"], Dict[str, str]],
) -> str :
    """
    Fetch a compound or gene record from KEGG using the REST API
    
    Checks the supplied cache for the requested identifier before sending HTTP 
    request to the REST API, and updates the cache after each HTTP request is made

    Parameters
    ----------
    identifier
        gene or compound identifier
    record type
        "compound" or "gene", defines the type of the identifier 
    cache
        cache dict with two sub-dicts (keys "compound" and "gene") that map identifiers to records 

    Returns
    -------
    record
        requested record, either from cache or API
    """
    # random ID for the request
    request_id = str(uuid.uuid4())
    kegg_rest_base = "https://rest.kegg.jp/get/"
    if (record := cache[record_type].get(identifier)) is None:
        _LOGGER.debug(
            "%(module)s.%(func)s: fetching %(record_type)s record with identifier=%(identifier)s [request_id=%(request_id)s]",
            {
                "module": __name__,
                "func": sys._getframe().f_code.co_name,
                "record_type": record_type,
                "identifier": identifier,
                "request_id": request_id
            }
        )
        # delay before sending the request to aviod overloading the server
        time.sleep(_HTTP_REQUEST_DELAY)
        # add /kgml to pathway identifiers to request the XML
        suffix = "/kgml" if record_type == "pathway" else ""
        record = requests.get(kegg_rest_base + identifier + suffix).text
        cache[record_type][identifier] = record
        _LOGGER.debug(
            "%(module)s.%(func)s: [request_id=%(request_id)s] done",
            {
                "module": __name__,
                "func": sys._getframe().f_code.co_name,
                "request_id": request_id
            }
        )
    return record