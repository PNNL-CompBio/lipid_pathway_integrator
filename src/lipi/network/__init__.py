"""
lipid_pathway_integrator/network/__init__.py
Dylan Ross (dylan.ross@pnnl.gov)

    Sub-package for defining the network-based pathway data model
"""


import json
from typing import List, Iterable, Union, Set, Tuple

from networkx import DiGraph

from lipi.refpath.manager import RefPathManager
from lipi.identifiers import Identifier, IdentifierGroup
from lipi.network.nodes import Protein, Metabolite, Lipid, LipidGroup


#===============================================================================
# Helper Functions

def _create_igroup(
    entity_id: int,
    names: Iterable[str],
    tags: Iterable[str]
) -> IdentifierGroup :
    """ create an IdentifierGroup instance from lists of names and tags """
    identifiers: List[Identifier] = []
    display_name = None
    for name, tags in zip(names, tags):
        # update display name for the group if any of the identifiers have the 
        # "display": True tag
        if json.loads(tags).get("display"):
            display_name = name
        identifiers.append(Identifier(name, tags))
        # if none of the identifiers had the "display": True tag, then just use the name of the 
        # first one in the group as the display name for the group
        if display_name is None:
            display_name = identifiers[0].name
    assert type(display_name) is str
    return IdentifierGroup(display_name, identifiers, entity_id=entity_id)


def _parse_fa_restrictions(
    original_fas: List[str]
) -> Set[Tuple[int, int]] :
    """
    parse the original FA restriction specification from a list of strings like 
    ["C:U", ...] into a set of tuples with ints {(C, U), ...}
    """
    return set([
        tuple(map(int, fa_str.split(":")))
        for fa_str in original_fas
    ])  # type: ignore (I know this will only ever be a 2-tuple of ints, not arbirary len)


def _create_node_data(
    node_type: str,
    igroup: IdentifierGroup,
) -> Union[Protein, Metabolite, Lipid, LipidGroup] :
    """  """
    match node_type:
        case "PROTEIN":
            return Protein(igroup)
        case "METABOLITE":
            return Metabolite(igroup)
        case "LIPID":
            # TODO: populate extra fields like FA sum composition and 
            #       individual FAs
            return Lipid(igroup)
        case "LIPID_GROUP":
            # find and unpack the LIPID MAPS classification info
            # search for the tags string that contains "is_lipid_group", this should also have the
            # LIPID MAPS classification info
            # also load/parse any FA restrictions that have been defined 
            # (require_fas, exclude_fas, exact_fas)
            classification_tags = None
            require_fas = set()
            exclude_fas = set()
            exact_fas = set()
            for identifier in igroup.components:
                if "is_lipid_group" in identifier.tags:
                    parsed_tags = json.loads(identifier.tags)
                    classification_tags = parsed_tags["classification"]
                    # parse any FA restrictions that were specified
                    if (req_fas := parsed_tags.get("require_fas")) is not None:
                        require_fas = _parse_fa_restrictions(req_fas)
                    if (exc_fas := parsed_tags.get("exclude_fas")) is not None:
                        exclude_fas = _parse_fa_restrictions(exc_fas)
                    if (exa_fas := parsed_tags.get("exact_fas")) is not None:
                        exact_fas = _parse_fa_restrictions(exa_fas)
            if classification_tags is None:
                # not sure how this would ever happen, but handle the possibility
                raise ValueError(
                    "_create_node_data: node_type is LIPID_GROUP but no identifiers with "
                    "'is_lipid_group' were found in the identifier group"
                )
            return LipidGroup(
                igroup,
                classification_tags["category"],
                classification_tags["class"],
                classification_tags["subclass"],
                require_fas=require_fas,
                exclude_fas=exclude_fas,
                exact_fas=exact_fas
            )
        case _:
            raise ValueError(f"_create_node_data: unrecognized node type: {node_type}")


#===============================================================================
# Main Pathway Loading Function

def load_pathway(
    manager: RefPathManager, 
    pathway_id: int
) -> DiGraph :
    """ load a specified pathway (by pathway ID) into a `networkx.DiGraph` """
    # TODO: This creates a full copy of each node data instance. It would be better to have
    #       like a dictionary or something that can map potentially multiple references to a 
    #       single instance of the node data
    # initialize graph (directed)
    g = DiGraph()
    # fetch node/edge ids
    node_ids, edge_ids = manager.pathways.fetch_node_and_edge_ids(pathway_id)
    # add nodes
    for node_id, node_type, entity_id, ig_names, ig_tags, x, y in manager.nodes.fetch_data(node_ids):
        # create identifier group
        igroup = _create_igroup(entity_id, ig_names, ig_tags)
        # create node data object
        data = _create_node_data(node_type, igroup)
        # add node to the graph
        # even though node type and name are accessible through the node data object,
        # add them as additional direct node attributes for convenience
        g.add_node(node_id, data=data, type=node_type, label=igroup.name, x=x, y=y)
    # add edges
    for src_node_id, tgt_node_id, edge_type in manager.edges.fetch_data(edge_ids):
        g.add_edge(src_node_id, tgt_node_id, type=edge_type)
    # return the graph
    return g
