"""
lipid_pathway_integrator/visualization.py
Dylan Ross (dylan.ross@pnnl.gov)

    Define functions for visualizing pathway data
"""


from typing import Optional, Tuple
from collections.abc import Collection, Mapping, Iterable

import networkx as nx
import numpy as np
from matplotlib import pyplot as plt, colors, cm
from matplotlib.axes import Axes
from matplotlib.patches import ArrowStyle, Wedge, Circle

from lipi.stats import Statistic


#===============================================================================
# CONSTANTS

# define colors for different plot elements
_COLORS = {
    "edge": "#666666",
    "protein": "#DA6666",
    "metabolite": "#9C98E3",
    "lipid_group": "#78DF9A",
    "lipid": "#DCDF78"
}

#===============================================================================
# PLOTTING FUNCTIONS


def plot_pathway_layout(
    pathway: nx.DiGraph,
    node_size: float = 256,
    node_pos: Optional[Mapping[int, Collection[float]]] = None, 
    label_font_size: int = 6,
    label_font_family: str = "sans-serif",
    add_labels: bool = True,
    ax: Optional[Axes] = None,
    figsize: Optional[Tuple[float, float]] = None,
    tight_layout: bool = True,
    show: bool = True,
    close: bool = True,
    figname: Optional[str] = None
) -> Axes :
    """
    Parameters
    ----------
    pathway
    node_size
    node_pos
    label_font_size
    ax
    figsize
    tight_layout
    show
    close
    figname
    """
    if ax is None: 
        # if an Axes instance was not passed in, create one
        if figsize is None:
            # if no figure size is specified, make a guess based on size of pathway
            sz = min(max(len(pathway) / 4 + 0.5, 1.25), 12)
            figsize = (sz, sz)
        fig, ax = plt.subplots(figsize=figsize)
    if node_pos is None:
        # if no node position mapping was provided, use kamada kawai layout
        node_pos = nx.kamada_kawai_layout(pathway)
    # draw metabolite nodes
    nx.draw_networkx_nodes(
        pathway, 
        node_pos,
        nodelist=[n for n in pathway.nodes if pathway.nodes[n]["type"] == "PROTEIN"],
        node_shape="h",
        node_size=node_size,
        node_color=_COLORS["protein"], 
        ax=ax, 
        edgecolors=_COLORS["edge"]
    )
    nx.draw_networkx_nodes(
        pathway, 
        node_pos,
        nodelist=[n for n in pathway.nodes if pathway.nodes[n]["type"] == "METABOLITE"],
        node_shape="s",
        node_size=0.8 * node_size,
        node_color=_COLORS["metabolite"], 
        ax=ax, 
        edgecolors=_COLORS["edge"]
    )
    nx.draw_networkx_nodes(
        pathway, 
        node_pos,
        nodelist=[n for n in pathway.nodes if pathway.nodes[n]["type"] == "LIPID_GROUP"],
        node_shape="o",
        node_size=node_size,
        node_color=_COLORS["lipid_group"], 
        ax=ax, 
        edgecolors=_COLORS["edge"]
    )
    nx.draw_networkx_nodes(
        pathway, 
        node_pos,
        nodelist=[n for n in pathway.nodes if pathway.nodes[n]["type"] == "LIPID"],
        node_shape="^",
        node_size=1.1 * node_size,
        node_color=_COLORS["lipid"], 
        ax=ax, 
        edgecolors=_COLORS["edge"]
    )
    # label nodes
    if add_labels:
        nx.draw_networkx_labels(
            pathway, 
            node_pos, 
            labels={n: f"{n}\n{pathway.nodes[n]["label"]}" for n in pathway}, 
            horizontalalignment="left",
            ax=ax, 
            font_size=label_font_size,
            clip_on=False, 
            font_family=label_font_family,
        )
    # draw edges
    nx.draw_networkx_edges(
        pathway, 
        node_pos, 
        ax=ax, 
        arrowstyle="-|>", 
        arrowsize=8,
        node_size=node_size, 
        edge_color=_COLORS["edge"]
    )
    # finish plot
    ax.axis("off")
    ax.set_aspect("equal")
    # final optional tidying and showing and/or saving to file
    if tight_layout:
        plt.tight_layout()
    if figname is not None:
        plt.savefig(figname, bbox_inches="tight", dpi=400, transparent=True)
    if show:
        plt.show()
    if close: 
        plt.close()
    return ax


def _v2c(
    stat: None | Statistic | Collection[Statistic], 
    vmin: float, 
    vmax: float
) -> str : 
    """
    map a stat value to a color using the "coolwarm" colormap
    """
    if stat is None:
        return "#FFFFFF"
    elif type(stat) is Statistic:
        value = stat.value
    else:
        # collection of statistics
        values = [_.value for _ in stat]  # type: ignore
        if len(values) == 0:
            return "#FFFFFF"
        value = np.mean(values)
    return colors.to_hex(
        cm.ScalarMappable(
            norm=colors.Normalize(vmin=vmin, vmax=vmax), 
            cmap="coolwarm"
        ).to_rgba(value)  # type: ignore
    )


def _scale_radius(
    stat: Statistic, 
    r: float, 
    vmin: float, 
    vmax: float
) -> float :
    """ 
    return a fraction of the radius that is proportional to the difference between 
    the stat value and vmin or vmax (depends what side of the midpoint between vmin
    and vmax that the value falls on)
    """
    half_diff = (vmax - vmin) / 2
    midpoint = vmin + half_diff
    return min(1, abs(stat.value - midpoint) / half_diff) * r


def _draw_lipid_pie(
    ax: Axes, 
    stats, 
    x: float, 
    y: float, 
    r: float,
    vmin: float,
    vmax: float,
    scale_radius: bool
) -> None :
    """
    Helper function to draw a little pie chart for the stats from the lipid nodes
    """
    n = len(stats)
    if n > 0:
        slice_degrees = 360. / n
        start_degrees = 90.
        for stat in sorted(stats, key=lambda s: s.value):
            scaled_r = _scale_radius(stat, r, vmin, vmax) if scale_radius else r
            ax.add_patch(
                Wedge(
                    (x, y), 
                    scaled_r, 
                    start_degrees, 
                    start_degrees + slice_degrees, 
                    fc=_v2c(stat, vmin, vmax),
                    #ec=_COLORS["edge"],
                    lw=0.5
                )
            )
            start_degrees += slice_degrees
            ax.plot([x, x], [y, y + r * 1.15], ls="-", c=_COLORS["edge"], lw=0.5)
    ax.add_patch(
        Circle(
            (x, y),
            r,
            fc="none",
            ec=_COLORS["edge"],
            lw=0.5
        )
    )


def plot_pathway_stats(
    pathway, 
    vmin: float = -1,
    vmax: float = 1,
    scale_radius: bool = False,
    node_pos: Optional[Mapping[int, Collection[float]]] = None, 
    label_font_size: int = 6,
    label_font_family: str = "sans-serif",
    label_rotation: Optional[float] = None,
    lipid_pie_r: float = 10,
    ax: Optional[Axes] = None,
    figsize: Optional[Tuple[float, float]] = None,
    tight_layout: bool = True,
    show: bool = True,
    close: bool = True,
    figname: Optional[str] = None,
    contains_fas: Optional[Iterable[Tuple[int, int]]] = None
) -> Axes :
    """
    """
    # ============ SETUP ================
    if ax is None: 
        # if an Axes instance was not passed in, create one
        if figsize is None:
            # if no figure size is specified, make a guess based on size of pathway
            sz = min(max(len(pathway) / 4. + 0.5, 1.25), 12.)
            figsize = (sz, sz)
        else: 
            sz = max(figsize)
        fig, ax = plt.subplots(figsize=figsize)
    if node_pos is None:
        # if no node position mapping was provided, use kamada kawai layout
        node_pos = nx.kamada_kawai_layout(pathway)
    # ============ EDGES ================
    nx.draw_networkx_edges(
        pathway, 
        node_pos, 
        ax=ax, 
        arrowstyle=ArrowStyle.CurveB(),  # type: ignore
        edge_color=_COLORS["edge"],
        node_size=555,
        arrowsize=6
    )  # type: ignore
    # ============ RNA/PROTEIN NODES ================
    nx.draw_networkx_nodes(
        pathway, 
        node_pos, 
        nodelist=[n for n in pathway.nodes if pathway.nodes[n]["type"] == "PROTEIN"],
        node_color=[
            # color is average of all stats entries with omics_type="proteomics" tag
            _v2c(
                node["data"].find_stats_with_tags(omics_type="proteomics"),
                vmin,
                vmax
            ) 
            for node_id, node in pathway.nodes(data=True) 
            if node["type"] == "PROTEIN"
        ],
        edgecolors=_COLORS["edge"],
        linewidths=0.5,
        node_size=256,
        node_shape="h"
    )
    nx.draw_networkx_nodes(
        pathway, 
        node_pos, 
        nodelist=[n for n in pathway.nodes if pathway.nodes[n]["type"] == "PROTEIN"],
        node_color=[
            # color is average of all stats entries with omics_type="transcriptomics" tag
            _v2c(
                node["data"].find_stats_with_tags(omics_type="transcriptomics"),
                vmin, 
                vmax
            ) 
            for node_id, node in pathway.nodes(data=True) 
            if node["type"] == "PROTEIN"
        ],
        edgecolors=_COLORS["edge"],
        linewidths=0.5,
        node_size=64,
        node_shape="h"
    )
    # ============ METABOLITE NODES ================
    nx.draw_networkx_nodes(
        pathway, 
        node_pos, 
        nodelist=[n for n in pathway.nodes if pathway.nodes[n]["type"] == "METABOLITE"],
        node_color=[
            # color is average of all stats entries with omics_type="metabolomics" tag
            _v2c(
                node["data"].find_stats_with_tags(omics_type="metabolomics"), 
                vmin, 
                vmax
            )
            for node_id, node in pathway.nodes(data=True) 
            if node["type"] == "METABOLITE"
        ],
        edgecolors=_COLORS["edge"],
        linewidths=0.5,
        node_size=120,
        node_shape="s"
    )
    # ============ LIPID NODES ================
    # lipid groups
    for node_id, node in pathway.nodes(data=True):
        if node["type"] == "LIPID_GROUP":
            # collect all of the relevant Statistics from the contained lipids
            stats = list(node["data"].find_filtered_lipid_stats_with_tags(
                omics_type="lipidomics", 
                contains_fas=contains_fas
            ))
            x, y = node_pos[node_id]
            _draw_lipid_pie(ax, stats, x, y, lipid_pie_r, vmin, vmax, scale_radius)
    # single lipid nodes
    nx.draw_networkx_nodes(
        pathway, 
        node_pos, 
        nodelist=[n for n in pathway.nodes if pathway.nodes[n]["type"] == "LIPID"],
        node_color=[
            # color is average of all stats entries with omics_type="lipidomics" tag
            _v2c(
                node["data"].find_stats_with_tags(omics_type="lipidomics"),
                vmin, 
                vmax
            ) 
            for node_id, node in pathway.nodes(data=True) 
            if node["type"] == "LIPID"
        ],
        edgecolors=_COLORS["edge"],
        linewidths=0.5,
        node_size=120,
        node_shape="^"
    )
    # ============ NODE LABELS ================
    lbls = nx.draw_networkx_labels(
        pathway, 
        node_pos,
        {
            #node_id: f"{node_id}\n{node["label"]}"
            node_id: node["label"]
            for node_id, node in pathway.nodes(data=True)
        },
        font_size=label_font_size,
        font_family=label_font_family,
        horizontalalignment="left",
        verticalalignment="bottom", 
        clip_on=False
    )
    # optionally rotate the text labels
    if label_rotation is not None:
        for _, text in lbls.items():
            text.set_verticalalignment("center")
            text.set_rotation_mode("anchor")
            text.set_rotation(label_rotation)
    # ============ FINISH PLOT ================
    ax.axis("off")
    ax.set_aspect("equal")
    # final optional tidying and showing and/or saving to file
    if tight_layout:
        plt.tight_layout()
    if figname is not None:
        plt.savefig(figname, bbox_inches="tight", dpi=400, transparent=True)
    if show:
        plt.show()
    if close: 
        plt.close()
    return ax

