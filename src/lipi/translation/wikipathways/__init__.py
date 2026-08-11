"""
lipid_pathway_integrator/translation/wikipathways/__init__.py
Dylan Ross (dylan.ross@pnnl.gov)

    Sub-package for ingesting reference pathways from WikiPathways
"""


import re
from xml.etree import ElementTree
from typing import Dict, Any

import networkx as nx

from lipi.utility import Source
from lipi.identifiers import Identifier
from lipi.identifiers.factory import name_to_igroup
from lipi.network import Pathway
from lipi.network.nodes import Protein, Metabolite, Lipid


#===============================================================================
# Constants

_XREF_SOURCE_TO_SOURCE_ENUM = {
    "Uniprot-TrEMBL": Source.UNIPROT,
    "LIPID MAPS": Source.LIPID_MAPS,
    "Enzyme Nomenclature": Source.EC,
    "ChEBI": Source.CHEBI,
    "HMDB": Source.HMDB,
    "PubChem-compound": Source.PUBCHEM,
    "Ensembl": Source.ENSEMBL,
    "CAS": Source.CAS,
    "Wikidata": Source.WIKIDATA,
    "EMBL": Source.EMBL,
    "Entrez Gene": Source.ENTREZ_GENE,
}


#===============================================================================
# Helper functions 

def _load_inital_wp_graph(gpml: str) -> nx.DiGraph :
    """ ? """
    gpml = re.sub(r'xmlns="[A-Za-z0-9:/.]+"', "", gpml)
    root = ElementTree.fromstring(gpml)
    wpg = nx.DiGraph()
    data_node_eid_to_gid = {}
    group_gid_to_eid = {}
    for element in root:
        print(element)
        match element.tag:
            case "DataNode":
                if (eid := element.attrib.get("GraphId")) is not None:
                    name = element.attrib["TextLabel"]
                    ntype = element.attrib.get("Type")
                    xrefs = {e.attrib["Database"]: e.attrib["ID"] for e in element.findall("Xref")}
                    wpg.add_node(eid, name=name, ntype=ntype, gpml="DataNode", xrefs=xrefs)
                    if (gid := element.attrib.get("GroupRef")) is not None:
                        data_node_eid_to_gid[eid] = gid
            case "Group":
                eid = element.attrib["GraphId"]
                gid = element.attrib["GroupId"]
                wpg.add_node(eid, gpml="Group")
                group_gid_to_eid[gid] = eid
            case "Interaction":
                eid = element.attrib["GraphId"]
                points = [e.attrib.get("GraphRef") for e in element.find("Graphics").findall("Point")]  # type: ignore
                if len(points) == 2:
                    ifrom, ito = points
                    wpg.add_node(eid, gpml="Interaction")
                    if ifrom is not None:
                        wpg.add_edge(ifrom, eid)
                    if ito is not None:
                        wpg.add_edge(eid, ito)
                    for anchor in element.find("Graphics").findall("Anchor"):  # type: ignore
                        wpg.add_edge(anchor.attrib["GraphId"], eid)
    # connect data nodes to group nodes, requires mapping between group ids and element ids
    # add an edge attribute to distinguish connection of group members to group node from other relations
    for eid, gid in data_node_eid_to_gid.items():
        wpg.add_edge(eid, group_gid_to_eid[gid], group_member=True)
    return wpg


def _collapse_interaction_nodes(wpg: nx.DiGraph) -> None :
    """ ? """
    for node in wpg.copy():
        # must check that the node is still in the graph first, because it could have been removed as a part
        # of a double interaction node
        if node in wpg:
            upstream_int = [
                n 
                for n in wpg.predecessors(node) 
                if wpg.nodes[n].get("gpml") in ["Interaction", None]
            ]
            upstream_nonint = [
                n 
                for n in wpg.predecessors(node) 
                if wpg.nodes[n].get("gpml") in ["DataNode", "Group"]
            ]
            downstream_int = [
                n 
                for n in wpg.successors(node) 
                if wpg.nodes[n].get("gpml") in ["Interaction", None]
            ]
            downstream_nonint = [
                n 
                for n in wpg.successors(node) 
                if wpg.nodes[n].get("gpml") in ["DataNode", "Group"]
            ]
            match (gpml := wpg.nodes[node].get("gpml")):
                case "Interaction":
                    # if the interaction connects non-Interaction nodes directly upstream and downstream
                    # directly connect all up- and downstream nodes
                    if upstream_int == [] and downstream_int == []:
                        for us in upstream_nonint:
                            for ds in downstream_nonint:
                                wpg.add_edge(us, ds)
                        wpg.remove_node(node)
                case None:  # this is an Anchor
                    if len(upstream_int) == 1 and len(downstream_int) == 1:
                        # a set of nodes consisting of an anchor node that separates two interaction nodes up- 
                        # and downstream can be fully collapsed, connecting the down-upstream nodes to the 
                        # down-downstream nodes via the up-upstream nodes, for example:
                        # DataNodeB->Int1->Anchor->Int2->DataNodeC   becomes   DataNodeA->DataNodeB->DataNodeC
                        #                 DataNodeA-^
                        # collect references to up-upstream, down-upstream, and down-downstream non-interaction nodes
                        upups = [
                            n 
                            for n in wpg.predecessors(upstream_int[0]) 
                            if wpg.nodes[n].get("gpml") in ["DataNode", "Group"]
                        ]
                        downups = [
                            n 
                            for n in wpg.predecessors(downstream_int[0]) 
                            if wpg.nodes[n].get("gpml") in ["DataNode", "Group"]
                        ]
                        downdowns = [
                            n 
                            for n in wpg.successors(downstream_int[0]) 
                            if wpg.nodes[n].get("gpml") in ["DataNode", "Group"]
                        ]
                        # add new edges to connect up-upstream to down-downstream nodes 
                        # and down-upstream to down-downstream nodes
                        for uu in upups:
                            for du in downups:
                                wpg.add_edge(du, uu)
                            for dd in downdowns:
                                wpg.add_edge(uu, dd)
                        # delete the the Interaction nodes from this set
                        wpg.remove_nodes_from([upstream_int[0], node, downstream_int[0]])


def _ungroup_nodes(wpg: nx.DiGraph) -> None :
    """ ? """
    for node in wpg.copy():
        if wpg.nodes[node].get("gpml") == "Group":
            group_members = [p for p in wpg.predecessors(node) if wpg.get_edge_data(p, node).get("group_member")]
            upstream = [p for p in wpg.predecessors(node) if not wpg.get_edge_data(p, node).get("group_member")]
            downstream = list(wpg.successors(node))
            for group_member in group_members:
                for us in upstream:
                    wpg.add_edge(us, group_member)
                for ds in downstream:
                    wpg.add_edge(group_member, ds)
            wpg.remove_node(node)


def _ingest_gpml(
        gpml: str
    ) -> nx.DiGraph :
    """ ? """
    # in-place patch the GPML to remove the annoying namespace that is used to prefix all of the tags
    wpg = _load_inital_wp_graph(gpml)
    # first remove any floating anchor nodes (no nodes upstream)
    for node in wpg.copy():
        if wpg.nodes[node].get("gpml") is None and len(list(wpg.predecessors(node))) == 0:
            wpg.remove_node(node)
    # collapse interaction nodes
    _collapse_interaction_nodes(wpg)
    # remove any unconnected nodes
    wpg.remove_nodes_from(nx.isolates(wpg.copy()))
    # ungroup grouped nodes
    _ungroup_nodes(wpg)
    return wpg


def _translate_wp_graph(wpg: nx.DiGraph) -> nx.DiGraph :
    """ ? """
    # translate nodes first, indexed by identifier from wikipathways graph
    wpg_nodes_to_igroups = {}
    for node in wpg:
        if (name := wpg.nodes[node].get("name")) is not None:
            ig = name_to_igroup(name)
            for xref_source, xref_id in wpg.nodes[node]["xrefs"].items():
                if xref_source != "":
                    ig += Identifier(xref_id, _XREF_SOURCE_TO_SOURCE_ENUM[xref_source])
            wpg_nodes_to_igroups[node] = ig
            # check the translated nodes for potenially redundant identifier groups and combine if necessary
            added = False
            for n, i in wpg_nodes_to_igroups.copy().items():
                if i == ig:
                    # combine identifier groups if matching one already present
                    wpg_nodes_to_igroups[n] += ig
                    # copy the entry
                    wpg_nodes_to_igroups[node] = wpg_nodes_to_igroups[n]
                    added = True
            if not added:
                # add a new one if there was not an existing igroup that matched
                wpg_nodes_to_igroups[node] = ig
    # now build up the new graph based on edges from the wikipathways graph
    ng = nx.DiGraph()
    for wpnode_start, wpnode_end in wpg.edges: 
        nnode_start = None
        ig_start = wpg_nodes_to_igroups[wpnode_start]
        match wpg.nodes[wpnode_start]["ntype"]:
            case "Metabolite":
                nnode_start = Metabolite(ig_start)
            case "GeneProduct":
                nnode_start = Protein(ig_start)
        nnode_end = None
        ig_end = wpg_nodes_to_igroups[wpnode_end]
        match wpg.nodes[wpnode_end]["ntype"]:
            case "Metabolite":
                nnode_end = Metabolite(ig_end)
            case "GeneProduct":
                nnode_end = Protein(ig_end)
        if nnode_start is not None and nnode_end is not None:
            ng.add_edge(nnode_start, nnode_end)
    # return the translated graph
    return ng


#===============================================================================
# Main ingestion function

def ingest_wikipathways_pathway(
        gpml: str
    ) -> Pathway :
    """
    Parameters
    ----------
    gpml
        pathway in GPML as a string
    """
    # initial ingestion of graph data in WikiPathways format
    wpg = _ingest_gpml(gpml)

    # translate wikipathways graph into internal data model
    ng = _translate_wp_graph(wpg)

    # initialize and return Pathway instance
    return Pathway("pathway identifier", ng, Source.WIKIPATHWAYS)



