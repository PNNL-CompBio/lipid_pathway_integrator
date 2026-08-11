"""
lipi/translation/yaml.py
Dylan Ross (dylan.ross@pnnl.gov)

    Module for translating pathway data from/to YAML files
"""


import os
import sys
from logging import getLogger
from typing import Dict, Any

import yaml

from lipi.refpath.manager import RefPathManager


#===============================================================================
# Constants

# module-level logger
_LOGGER = getLogger(__name__)


#===============================================================================
# helper functions

def _validate_yaml_data(
    yaml_data: Dict[Any, Any]
) -> None :
    """
    ensure the YAML has the expected structure, raise a RuntimeError if not
    """
    # TODO: implement me


def _add_entities(
    yaml_data: Dict[Any, Any],
    manager: RefPathManager
) -> Dict[int, int] :
    """
    adds all of the entities from the entities section of the YAML
    returns a dictionary mapping internal entity IDs in the YAML to 
    the corresponding reference pathway database entity ID
    """
    internal_to_db_entity_ids = {}
    for internal_entity_id, identifiers in yaml_data["entities"].items():
        # create an entity in the database, get its ID
        db_entity_id = manager.entities.insert_or_update_from_identifiers(
            [
                (identifier["name"], identifier["tags"])
                for identifier in identifiers
            ]
        )
        # map internal to YAML and database entity IDs
        internal_to_db_entity_ids[internal_entity_id] = db_entity_id
    return internal_to_db_entity_ids


def _add_nodes(
    yaml_data: Dict[Any, Any],
    internal_to_db_entity_ids: Dict[int, int],
    pathway_id: int,
    manager: RefPathManager
) -> Dict[int, int] :
    """
    adds all of the nodes from the nodes section of the YAML
    returns a dictionary mapping internal node IDs in the YAML to 
    the corresponding reference pathway database node ID
    """
    internal_to_db_node_ids = {}
    for internal_node_id, node_data in yaml_data["nodes"].items():
        # insert a node into the database, get its ID
        db_node_id = manager.nodes.insert(
            pathway_id,
            # convert internal YAML entity ID to database entity ID
            internal_to_db_entity_ids[node_data["entity"]],
            float(node_data["x"]),
            float(node_data["y"])
        )
        # map internal to YAML and database entity IDs
        internal_to_db_node_ids[internal_node_id] = db_node_id
    return internal_to_db_node_ids


def _add_edges(
    yaml_data: Dict[Any, Any],
    internal_to_db_node_ids: Dict[int, int],
    manager: RefPathManager
) -> None :
    """
    adds all of the edges from the edges section of the YAML
    """
    for edge_data in yaml_data["edges"]:
        # insert an edge into the database
        # no need to store the edge ID 
        _ = manager.edges.insert_or_ignore(
            # convert internal YAML node IDs to database entity IDs
            internal_to_db_node_ids[edge_data["source_node"]],
            internal_to_db_node_ids[edge_data["target_node"]],
            edge_data["type"],
        )


#===============================================================================
# main ingestion function

def ingest_yaml_pathway(
    yaml_file: str,
    name: str,
    manager: RefPathManager
) -> int :
    """
    ingest a pathway defined in a YAML file

    Parameters
    ----------
    yaml_file
        path to input YAML file
    name
        assign a name for the pathway
    
    Returns
    -------
    pathway_id
        assigned pathway identifier
    """
    _LOGGER.info(
        "%(module)s.%(func)s: ingesting pathway from YAML file: %(yaml_file)s",
        {
            "module": __name__,
            "func": sys._getframe().f_code.co_name,
            "yaml_file": yaml_file
        }
    )
    # load data from the input YAML file and validate it
    if not os.path.isfile(yaml_file):
        raise FileNotFoundError(yaml_file)
    with open(yaml_file, "r") as yf:
        yaml_data = yaml.safe_load(yf)
    _validate_yaml_data(yaml_data)
    # insert a pathway entry into the database, get the database pathway ID
    pathway_id = manager.pathways.insert_or_ignore(
        name, 
        source="YAML", 
    )
    # insert the entities
    internal_to_db_entity_ids = _add_entities(
        yaml_data, 
        manager
    )
    # insert the nodes
    internal_to_db_node_ids = _add_nodes(
        yaml_data,
        internal_to_db_entity_ids,
        pathway_id,
        manager
    )
    # insert the edges
    _add_edges(
        yaml_data,
        internal_to_db_node_ids,
        manager
    )
    # return the pathway ID from the database
    return pathway_id
    


# TODO: export a pathway to YAML
