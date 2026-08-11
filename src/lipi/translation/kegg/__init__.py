"""
lipi/translation/kegg/__init__.py
Dylan Ross (dylan.ross@pnnl.gov)

    Sub-package for ingesting reference pathways from KEGG
"""


from xml.etree import ElementTree
import re
from typing import Dict, Tuple, Set, Union, Optional
from collections.abc import Iterable
import pickle
import os
import sys
from logging import getLogger

from lipi.refpath.manager import RefPathManager
from lipi.translation.kegg._remote import KeggCache, get_empty_kegg_cache, fetch_record


#===============================================================================
# Constants

# module-level logger
_LOGGER = getLogger(__name__)


#===============================================================================
# Helper functions 

def _parse_record(record: str) -> Dict[str, str] :
    """ Parse a KEGG record into key:value pairs """
    parsed = {}
    k = None
    v = []
    for line in record.splitlines():
        if (match := re.match(r"^([A-Z_]+)\s+(.*)", line)):
            # encountered a line defining a new key
            if k is not None:
                # store parsed info for the previous key
                parsed[k] = "\n".join(v)
            # begin accumulating information for the new key
            # parse out the key and value line from the regex groups
            k = match.group(1)
            v = [match.group(2)]
        else:
            # continue to accumulate values under the new key
            v.append(line)
        # store parsed info for the final key
        if k is not None:
            parsed[k] = "\n".join(v)
    return parsed


def _unpack_dblinks(dblinks: Optional[str]) -> Iterable[Tuple[str, str]] :
    """ """
    unpacked = []
    if dblinks is None:
        return unpacked
    for a, *bs in map(lambda s: re.split(r":*\s+", s), re.split(r"\n\s+", dblinks)):
        for b in bs:
            unpacked.append((a, b))
    return unpacked


def _deposit_node_data(
    root: ElementTree.Element,
    kegg_cache: KeggCache,
    man: RefPathManager, 
    pathway_id: int
) -> Dict[int, Iterable[int]] :
    """ """
    kegg_id_to_node_id: Dict[int, Iterable[int]] = {}
    for element in root:
        if element.tag == "entry":
            match element.attrib["type"]:
                case "compound":
                    for cid in re.findall(r"cpd:C[0-9]+", element.attrib["name"]):
                        record = _parse_record(fetch_record(cid, "compound", kegg_cache))
                        identifiers = [
                            (name, {})
                            for name in re.split(r";\s+", record["NAME"])
                        ] + [
                            (ext_id, {"source": source})
                            for source, ext_id in _unpack_dblinks(record.get("DBLINKS"))
                        ] + [
                            (cid, {"source": "KEGG"}),
                        ]
                        # use the first name entry as the display name for the identifier group
                        # by adding the "display"=True tag
                        identifiers[0][1]["display"] = True
                        # create or update the entity
                        entity_id = man.entities.insert_or_update_from_identifiers(identifiers)
                        kegg_id_to_node_id[int(element.attrib["id"])] = [man.nodes.insert(
                            pathway_id,
                            entity_id,
                            # capture the layout from the original reference pathway
                            # assuming here that all "entry" elements have a child "graphics" 
                            # with x and y position attributes
                            x=float(element.find("graphics").attrib["x"]),  # type: ignore
                            y=float(element.find("graphics").attrib["y"])   # type: ignore
                        )]
                case "gene":
                    for gid in re.findall(r"[a-z]+:[0-9a-z_A-Z]+", element.attrib["name"]):
                        record = _parse_record(fetch_record(gid, "gene", kegg_cache))
                        if (symbol := record.get("SYMBOL")) is not None:
                            identifiers = [
                                # simple identifiers get tagged to indicate they are genes/gene products
                                (name, {"is_protein": True})
                                for name in re.split(r",\s+", symbol)
                            ] + [
                                (ext_id, {"source": source})
                                for source, ext_id in _unpack_dblinks(record["DBLINKS"])  
                            ] + [
                                (gid, {"source": "KEGG"})
                            ]
                            # use the first name entry as the display name for the identifier group
                            # by adding the "display"=True tag
                            identifiers[0][1]["display"] = True  # type: ignore
                            entity_id = man.entities.insert_or_update_from_identifiers(identifiers)
                            kegg_id_to_node_id[int(element.attrib["id"])] = [man.nodes.insert(
                                pathway_id,
                                entity_id,
                                # capture the layout from the original reference pathway
                                # assuming here that all "entry" elements have a child "graphics" 
                                # with x and y position attributes
                                x=float(element.find("graphics").attrib["x"]),  # type: ignore
                                y=float(element.find("graphics").attrib["y"])   # type: ignore
                            )]
                case "group":
                    # map this group's pathway ID to a list of node IDs for all of its components
                    # assume:
                    #   - all of the components are already present in path_id_to_node_id
                    #   - the components map straight to single node_ids not groups of node_ids (int not Iterable[int])
                    kegg_id_to_node_id[int(element.attrib["id"])] = [   # type: ignore
                        kegg_id_to_node_id[int(e.attrib["id"])] for e in element.findall("component")
                    ]
    return kegg_id_to_node_id


def _deposit_edge_data(
    root: ElementTree.Element,
    kegg_id_to_node_id: Dict[int, Iterable[int]],
    man: RefPathManager
) : 
    """ """
    for element in root:
        match element.tag:
            case "relation":
                match element.attrib["type"]:
                    case "ECrel":
                        # add two sets of edges:
                        #   1) connect entry1 to the id of the subtype 
                        #   2) connect id of the subtype to entry2
                        subtype_id = int(element.find("subtype").attrib["value"])  # type: ignore
                        entry1_id = int(element.attrib["entry1"])
                        entry2_id = int(element.attrib["entry2"])
                        for subtype_nid in kegg_id_to_node_id[subtype_id]:
                            for entry1_nid in kegg_id_to_node_id[entry1_id]:
                                man.edges.insert_or_ignore(
                                    entry1_nid,
                                    subtype_nid,
                                    "REACTION"
                                )
                            for entry2_nid in kegg_id_to_node_id[entry2_id]:
                                man.edges.insert_or_ignore(
                                    subtype_nid,
                                    entry2_nid,
                                    "REACTION"
                                )
                    case "PPrel" | "GErel":
                        if (subt := element.find("subtype")).attrib["name"] == "compound":  # type: ignore
                            # add two sets of edges:
                            #   1) connect entry1 to the id of the compound
                            #   2) connect id of the compound to entry2
                            icid = int(subt.attrib["value"])  # type: ignore
                            entry1_id = int(element.attrib["entry1"])
                            entry2_id = int(element.attrib["entry2"])
                            for icnid in kegg_id_to_node_id[icid]:
                                for entry1_nid in kegg_id_to_node_id[entry1_id]:
                                    man.edges.insert_or_ignore(
                                        entry1_nid,
                                        subtype_nid,
                                        "REACTION"
                                    )
                                for entry2_nid in kegg_id_to_node_id[entry2_id]:
                                    man.edges.insert_or_ignore(
                                        icnid,
                                        entry2_nid,
                                        "REACTION"
                                    )
                        else:
                            # add single edge between entry1 and entry2
                            for entry1_nid in kegg_id_to_node_id[entry1_id]:
                                for entry2_nid in kegg_id_to_node_id[entry2_id]:
                                    man.edges.insert_or_ignore(
                                        entry1_nid,
                                        entry2_nid,
                                        "REACTION"
                                    )
            case "reaction":
                # add two sets of edges: 
                #   1) connect the substrate(s) to the entity that is associated with this reaction
                #   2) then connect that entity to the product(s)
                for r_nid in kegg_id_to_node_id[int(element.attrib["id"])]:
                    for sube in element.findall("substrate"):
                        for sube_nid in kegg_id_to_node_id[int(sube.attrib["id"])]:
                            man.edges.insert_or_ignore(
                                sube_nid,
                                r_nid,
                                "REACTION"
                            )
                    for prde in element.findall("product"):
                        for prde_nid in kegg_id_to_node_id[int(prde.attrib["id"])]:
                            man.edges.insert_or_ignore(
                                r_nid,
                                prde_nid,
                                "REACTION"
                            )


#===============================================================================
# Main ingestion function

def ingest_kegg_pathway(
    kegg_pathway_identifier: str, 
    kegg_cache_file: str,
    manager: RefPathManager
) -> None :
    """
    """
    _LOGGER.info(
        "%(module)s.%(func)s: ingesting KEGG pathway: %(kegg_pathway_identifier)s",
        {
            "module": __name__,
            "func": sys._getframe().f_code.co_name,
            "kegg_pathway_identifier": kegg_pathway_identifier
        }
    )

    # make sure RefPathManager is connected
    if not manager.connected:
        raise RuntimeError("reference pathway databse manager is not connected")

    # load the KEGG cache (pickle)
    if not os.path.isfile(kegg_cache_file):
        # make an empty cache if cache file not found
        kegg_cache = get_empty_kegg_cache()
        _LOGGER.debug(
            "%(module)s.%(func)s: KEGG cache file %(kegg_cache_file)s not found, creating empty cache",
            {
                "module": __name__,
                "func": sys._getframe().f_code.co_name,
                "kegg_cache_file": kegg_cache_file
            }
        )
    else:
        # otherwise load the existing ones
        with open(kegg_cache_file, "rb") as pf: 
            kegg_cache = pickle.load(pf)
            _LOGGER.debug(
                "%(module)s.%(func)s: loaded KEGG cache file %(kegg_cache_file)s",
                {
                    "module": __name__,
                    "func": sys._getframe().f_code.co_name,
                    "kegg_cache_file": kegg_cache_file
                }
            )

    # fetch pathway record (KMGL XML) from KEGG, either from a cache or HTTP request to REST API
    try: 
        pathway_record = fetch_record(kegg_pathway_identifier, "pathway", kegg_cache)
    except Exception as e:
        # failed to fetch the record, bubble up the exception
        raise RuntimeError(
            f"failed to fetch KGML XML for KEGG pathway record: {kegg_pathway_identifier=}"
        ) from e
    
    # parse the xml 
    # this is the "pathway" element
    # elements under this are "entry" (node) or "reaction"/"relation" (edge) types that define the graph
    root = ElementTree.fromstring(pathway_record)
    
    # insert a pathway entry into the databse, get the database pathway ID
    # NOTE: assume that none of the .get() calls return None
    pathway_id = manager.pathways.insert_or_ignore(
        root.get("title"),        # type: ignore
        source="KEGG", 
        kegg_id=root.get("name"), # type: ignore
        url=root.get("link")      # type: ignore
    )

    # deposit all of the node data into the database
    # create a mapping between the pathway IDs and the database node IDs
    # then deposit all of the edge data into the database
    # track all of the database edge IDs that are added 
    _deposit_edge_data(
        root, 
        _deposit_node_data(
            root, 
            kegg_cache, 
            manager,
            pathway_id
        ), 
        manager
    )

    # save the cache in case anything new was added 
    # TODO: Make the cache a more formal data structure with a flag indicating whether it has been updated. 
    #       This way we could check the flag and not bother re-saving if there have been no updates
    with open(kegg_cache_file, "wb") as pf:
        pickle.dump(kegg_cache, pf)
        _LOGGER.debug(
            "%(module)s.%(func)s: saved KEGG cache file %(kegg_cache_file)s",
            {
                "module": __name__,
                "func": sys._getframe().f_code.co_name,
                "kegg_cache_file": kegg_cache_file
            }
        )