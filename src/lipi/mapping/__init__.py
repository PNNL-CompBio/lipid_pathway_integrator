"""
lipid_pathway_integrator/mapping/__init__.py
Dylan Ross (dylan.ross@pnnl.gov)

    Sub-package with utilities for mapping omics data into pathway data model
"""


from typing import Literal, Optional, Union, Set, Tuple
from collections.abc import Iterable
from logging import getLogger

from networkx import DiGraph
from lipidimea._lipidlib.parser import parse_lipid_name
from lipidimea._lipidlib.lipids import (
    Lipid as LipidIMEA_Lipid, 
    LipidWithChains as LipidIMEA_LipidWithChains
)

from lipi.identifiers import IdentifierGroup
from lipi.stats import Statistic
from lipi.network.nodes import Lipid, LipidGroup


#===============================================================================
# CONSTANTS

# module-level logger
_LOGGER = getLogger(__name__)


#===============================================================================
# Helper Functions

def _try_parse_lipid(
    igroup: IdentifierGroup
) -> Optional[Union[LipidIMEA_Lipid, LipidIMEA_LipidWithChains]] :
    """
    Go through each of the identifiers in the group, try to parse into a 
    `LipidIMEA_Lipid` or `LipidWithChains` object. Return the first successfully 
    parsed object, or None if none were able to be parsed.

    Note
    ----
    The `LipidIMEA_Lipid` and `LipidIMEA_LipidWithChains` objects come from 
    `lipidimea` and are renamed on import with the prepended LipidIMEA_ to 
    avoid confusion with the `Lipid` and `LipidGroup` node-associated data 
    objects from this codebase.
    """
    for identifier in igroup.components:
        # We could add some more conditions to omit identifiers that have 
        # associated tags, but for now I think we can just rely on the fact
        # that those are very unlikely to be successfully parsed and will 
        # be skipped anyway.
        if (lipid := parse_lipid_name(identifier.name)) is not None: 
            return lipid
    # if we make it through the whole identifier group without successfully 
    # parsing a lipid name, then return None
    return None


def _get_fas_from_lipidimea_lipid(
    lipid: Union[LipidIMEA_Lipid, LipidIMEA_LipidWithChains]
) -> Set[Tuple[int, int]] :
    # NOTE: Check for LipidWithChains first since it is a subclass of Lipid
    #       so isinstance(lipid, LipidIMEA_Lipid) will return True which is 
    #       not what we want.
    if isinstance(lipid, LipidIMEA_LipidWithChains):
        return set([
            (c, u) 
            for c, u in zip(lipid.fa_carbon_chains, lipid.fa_unsat_chains)
        ])
    if isinstance(lipid, LipidIMEA_Lipid):
        return set()
    # fallthrough, should never get here
    # TODO: change this to an exception, asserts can be ignored at runtime under
    #       certain conditions and are only useful with debugging enabled
    assert False


def _fatty_acid_restrictions_satisfied(
    lipid: Union[LipidIMEA_Lipid, LipidIMEA_LipidWithChains],
    lipid_group: LipidGroup
) -> bool :
    """
    return a boolean indicating wether a lipid (LipidIMEA Lipid or LipidWithChains) 
    satisfies whatever fatty acid restrictions have been defined for a LipidGroup
    """
    incoming_fas = _get_fas_from_lipidimea_lipid(lipid)
    # check excluded FAs
    if len(incoming_fas & lipid_group.exclude_fas) > 0:
        # contains an excluded FA, not satisfied
        return False
    # check required FAs
    if len(incoming_fas & lipid_group.require_fas) < len(lipid_group.require_fas):
        # does not contain all required FAs, not satisfied
        return False
    # check exact FAs
    if len(lipid_group.exact_fas) > 0 and incoming_fas != lipid_group.exact_fas:
        # exact FAs were specified but incoming FAs do not match, not satisfied
        return False
    # if none of the above checks failed, the restrictions are satisfied
    return True


#===============================================================================
# Main Mapping Functions

# TODO: This function might be a little big. Could be good to factor some 
#       parts down into separate helper functions. The transcriptomics/
#       proteomics/metabolomics branch and the lipidomics branch would be 
#       good candidates to split out into two separate helper functions

def map_omics_data_to_pathway(
    pathway: DiGraph,
    omics_type: Literal["transcriptomics", "proteomics", "metabolomics", "lipidomics"],
    igroups: Iterable[IdentifierGroup],
    stats: Iterable[Statistic],
) -> None :
    """
    map a set of omics features (as iterables of identifier groups 
    and corresponding statistics) to a pathway

    Parameters
    ----------
    pathway
        a pathway instance (`networkx.DiGraph`)
    omics_type
        transcriptomics, proteomics, metabolomics, or lipidomics
    igroups
        iterable of identifier groups for the omics features
    stats
        iterable of statistics for the omics features

    Note
    ----
    Automatically adds the "omics_type" tag to all of the input `Statistic`s
    """
    match omics_type:
        case "transcriptomics" | "proteomics" | "metabolomics" :
            # map the omics type to the corresponding node type 
            # that should recieve the statistic
            node_type = {
                "transcriptomics": "PROTEIN",
                "proteomics": "PROTEIN",
                "metabolomics": "METABOLITE",
            }[omics_type]
            for igroup, stat in zip(igroups, stats):
                for node_id, node in pathway.nodes(data=True):
                    if node["type"] == node_type and node["data"].igroup == igroup:
                        _LOGGER.debug(
                            "%(module)s.%(func)s: %(igroup)s matched to node: %(node)s (node_id=%(node_id)d)",
                            {
                                "module": __name__,
                                "func": "map_omics_data_to_pathway",
                                "igroup": igroup,
                                "node": node["data"],
                                "node_id": node_id
                            }
                        )
                        # add omics_type tag to the Statistic before adding
                        node["data"].stats.add(stat.add_tags(omics_type=omics_type))
        case "lipidomics":
            # lipids can map to either individual lipid nodes or lipid groups
            for igroup, stat in zip(igroups, stats):
                # add omics_type tag to each Statistic
                stat = stat.add_tags(omics_type=omics_type)
                # first, ensure we can parse out a lipid species from at least
                # one of the identifiers in the identifier group
                if (lipid := _try_parse_lipid(igroup)) is not None:
                    for node_id, node in pathway.nodes(data=True):
                        match node["type"]:
                            case "LIPID":
                                # Re-parse the node into a Lipid/LipidWithChains to 
                                # compare it to the omics feature. Since it is already
                                # associated with a LIPID node type, we assume that 
                                # parsing will be successful (i.e. not None)
                                node_lipid = _try_parse_lipid(node["data"].igroup)
                                assert node_lipid is not None
                                # make the comparison
                                if (
                                    lipid.lmaps_category == node_lipid.lmaps_category
                                    and lipid.lmaps_class == node_lipid.lmaps_class
                                    and lipid.lmaps_subclass == node_lipid.lmaps_subclass
                                    and lipid.fa_carbon == node_lipid.fa_carbon
                                    and lipid.fa_unsat == node_lipid.fa_unsat
                                    and (
                                        _get_fas_from_lipidimea_lipid(lipid)
                                            == _get_fas_from_lipidimea_lipid(node_lipid)
                                    )
                                ): 
                                    _LOGGER.debug(
                                        "%(module)s.%(func)s: %(lipid)s matched to node: %(node)s (node_id=%(node_id)d)",
                                        {
                                            "module": __name__,
                                            "func": "map_omics_data_to_pathway",
                                            "lipid": lipid,
                                            "node": node["data"],
                                            "node_id": node_id
                                        }
                                    )
                                    node["data"].stats.add(stat)
                            case "LIPID_GROUP":
                                # see if the incoming lipid feature fits in the LIPID_GROUP
                                # based on classification info and fatty acid restrictions (if specified)
                                lipid_group: LipidGroup = node["data"]
                                fa_restrictions_flag = _fatty_acid_restrictions_satisfied(lipid, lipid_group)
                                if (
                                    lipid.lmaps_category == lipid_group.lm_category
                                    and lipid.lmaps_class == lipid_group.lm_class
                                    and lipid.lmaps_subclass == lipid_group.lm_subclass
                                    and fa_restrictions_flag
                                ): 
                                    _LOGGER.debug(
                                        "%(module)s.%(func)s: %(lipid)s matched to node: %(node)s (node_id=%(node_id)d)",
                                        {
                                            "module": __name__,
                                            "func": "map_omics_data_to_pathway",
                                            "lipid": lipid,
                                            "node": lipid_group,
                                            "node_id": node_id
                                        }
                                    )
                                    # NOTE: In the following block we are getting a reference
                                    #       to a node (in matched_single_lipid_node) which is
                                    #       contained within the set in lipid_group.lipids. We
                                    #       are then modifying that node using the reference
                                    #       matched_single_lipid_node, which means we are indirectly
                                    #       modifying the lipid_group.
                                    # Create an instance of an individual Lipid node and 
                                    # add that into the LIPID_GROUP node's container for 
                                    # lipids (LipidGroup.lipids attribute). Check if one 
                                    # already exists first, and add one if not. 
                                    # Map the statistic onto that individual Lipid node.
                                    matched_single_lipid_node = None
                                    for node_lipid in lipid_group.lipids:
                                        # we already know classification info is fine from
                                        # the LIPID_GROUP level, just need to check that the
                                        # sum composition matches and individual FAs if known
                                        if (
                                            lipid.fa_carbon == node_lipid.fa_carbon
                                            and lipid.fa_unsat == node_lipid.fa_unsat
                                            and (
                                                _get_fas_from_lipidimea_lipid(lipid) == node_lipid.fas
                                            )
                                        ): 
                                            # found a matching node, grab a reference
                                            matched_single_lipid_node = node_lipid
                                            break
                                    # if we did not find a matching contained lipid node,
                                    # then create and add a new one
                                    if matched_single_lipid_node is None:
                                        matched_single_lipid_node = Lipid(
                                            IdentifierGroup.from_name(str(lipid)),
                                            fa_carbon=lipid.fa_carbon,
                                            fa_unsat=lipid.fa_unsat,
                                            fas=_get_fas_from_lipidimea_lipid(lipid),
                                        )
                                        lipid_group.lipids.append(matched_single_lipid_node)
                                    # finally: map the statistic to the matched lipid node
                                    # inside of the lipid_group via the reference we have
                                    matched_single_lipid_node.stats.add(stat)
                            # case _: ignore all other node types