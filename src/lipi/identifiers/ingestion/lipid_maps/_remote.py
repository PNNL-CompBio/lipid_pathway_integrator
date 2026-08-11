"""
lipi/identifiers/ingestion/lipid_maps/_remote.py
Dylan Ross (dylan.ross@pnnl.gov)

    Module for interacting with the LIPID MAPS REST API
"""


from typing import Dict, Any
import time
import uuid
from logging import getLogger

import requests


#===============================================================================
# Constants

# module-level logger
_LOGGER = getLogger(__name__)

# set a delay (in seconds) before sending HTTP requests to limit server load
_HTTP_REQUEST_DELAY = 0.5


#===============================================================================
# Type annotation for convenience

# map the request url to the parsed JSON response or None if the request was otherwise 
# successful but there were no results
type LipidMapsCache = Dict[str, Dict[str, Any]]


#===============================================================================
# basic function for fetching info from the API

def fetch_record(
    lm_id: str, 
    cache: LipidMapsCache,
) -> Dict[str, Any] :
    """
    Fetch a compound record from LMSD using the REST API
    
    Checks the supplied cache before sending HTTP request to the REST API, 
    and updates the cache after each new HTTP request is made

    Parameters
    ----------
    lm_id
        lipid identifier
    cache
        cache dict that maps request URL to response JSON (or None for empty responses)

    Returns
    -------
    record
        requested record, dict parsed from JSON response or None
    """
    # random ID for the request
    request_id = str(uuid.uuid4())
    url = f"https://www.lipidmaps.org/rest/compound/lm_id/{lm_id}/all"
    if (record := cache.get(url)) is None:
        # not in cache, new request
        _LOGGER.debug(
            "%(module)s.%(func)s: fetching record with lm_id=%(lm_id)s [request_id=%(request_id)s]",
            {
                "module": __name__,
                "func": "fetch_record",
                "lm_id": lm_id
            }
        )
        # delay before sending the request to aviod overloading the server
        time.sleep(_HTTP_REQUEST_DELAY)
        record = requests.get(url).json()
        # convert an empty result to None
        record = {} if record == [] else record
        # cache the result
        cache[url] = record
        _LOGGER.debug(
            "%(module)s.%(func)s: [request_id=%(request_id)s] done",
            {
                "module": __name__,
                "func": "fetch_record",
                "request_id": request_id
            }
        )
    return record
