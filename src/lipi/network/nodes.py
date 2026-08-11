"""
lipid_pathway_integrator/network/nodes.py
Dylan Ross (dylan.ross@pnnl.gov)

    Node class definitions
"""


from dataclasses import dataclass, field
from collections.abc import Iterable, Iterator
from typing import Optional, Set, List, Union, Tuple
import json
import itertools

from lipi.identifiers import IdentifierGroup
from lipi.stats import Statistic


#==============================================================================
# Simple nodes (_Node and subclasses) 

@dataclass
class _Node:
    """ 
    base class for node-associated data
    
    Attributes
    ----------
    igroup : `IdentifierGroup`
        the group of identifiers that are associated with this node
    stats : `Set[Statistic]`
        Set of primary statistics (`lipid_pathway_integrator.stats.Statistc`) associated with this node
    """
    igroup: IdentifierGroup
    stats: Set[Statistic] = field(default_factory=set) 
    
    def __repr__(self):
        return f"{self.__class__.__name__}({self.igroup.name} | {len(self.stats)} mapped statistics)"
    
    def find_stats_with_tags(self,
        **tags : Union[str, int, float, bool] 
    ) -> Iterator[Statistic] :
        """
        Generator that yields any `Statistic`s in `self.stats` that contain the 
        specified tags (specified tags only need to be a subset of all tags contained
        in the `Statistic`, it does not need to be an exact match)
        """
        for stat in self.stats:
            # The <= operator on dictionary views treats them as sets and checks 
            # if the left side is a subset of the right side.
            if tags.items() <= json.loads(stat.tags).items():
                yield stat


@dataclass
class Metabolite(_Node):
    """
    Node data type for a metabolite
    
    Attributes
    ----------
   igroup : `IdentifierGroup`
        the group of identifiers that are associated with this node
    stats : `Iterable[Statistic]`
        collection of primary statistics (`lipid_pathway_integrator.stats.Statistc`) associated with this node
    """


@dataclass
class Protein(_Node):
    """ 
    Node data type for a protein

    Attributes
    ----------
    igroup : `IdentifierGroup`
        the group of identifiers that are associated with this node
    stats : `Iterable[Statistic]`
        collection of primary statistics (`lipid_pathway_integrator.stats.Statistc`) associated with this node
    """

@dataclass
class Lipid(_Node):
    """ 
    Node type for a single lipid species
    
    Attributes
    ----------
    igroup : `IdentifierGroup`
        the group of identifiers that are associated with this node
    stats : `Iterable[Statistic]`
        collection of primary statistics (`lipid_pathway_integrator.stats.Statistc`) associated with this node
    fa_carbon, fa_unsat : `int`
        sum composition of the lipid species
    fas : `set(tuple(int, int))`
        individual FAs (if known), empty set if individual FAs are not known
    """
    fa_carbon: Optional[int] = None
    fa_unsat: Optional[int] = None
    fas: Set[Tuple[int, int]] = field(default_factory=set)


#==============================================================================
# Group nodes

class LipidGroup:
    """
    Node encapsulating a group of multiple lipids (`Lipid`), all under a 
    specified category, class, and subclass (using Lipid Maps ontology)
    
    Attributes
    ----------
    name : `str`
        lipid name
    lm_category : `str`
    lm_class : `str`
    lm_subclass : `str` or `None`
        Lipid MAPS classification, category and classification are required 
        but subclass can also be provided when relevant
    lipids : ``list(Lipid)``
        individual lipid species that belong to this group
    require_fas : ``set(tuple(int, int))``
    exclude_fas : ``set(tuple(int, int))``
    exact_fas : ``set(tuple(int, int))``
    """

    def __init__(self,
        igroup: IdentifierGroup,
        lm_category: str,
        lm_class: str,
        lm_subclass: Optional[str],
        lipids: Iterable[Lipid] = [], 
        require_fas: Set[Tuple[int, int]] = set(),
        exclude_fas: Set[Tuple[int, int]] = set(),
        exact_fas: Set[Tuple[int, int]] = set()
    ):
        self.__igroup: IdentifierGroup = igroup
        self.__lm_category: str = lm_category
        self.__lm_class: str = lm_class
        self.__lm_subclass: Optional[str] = lm_subclass
        self.__lipids: List[Lipid] = list(lipids) 
        self.__require_fas: Set[Tuple[int, int]] = require_fas
        self.__exclude_fas: Set[Tuple[int, int]] = exclude_fas
        self.__exact_fas: Set[Tuple[int, int]] = exact_fas

    def __repr__(self):
        return (
            f"LipidGroup({self.igroup} "
            f"[{self.lm_category} > {self.lm_class} > {self.lm_subclass}] "
            f"{self._gen_fa_restrictions_str()})"
        )

    # --- attribute getters ---

    @property
    def igroup(self) -> IdentifierGroup : return self.__igroup

    @property
    def lm_category(self) -> str : return self.__lm_category

    @property
    def lm_class(self) -> str : return self.__lm_class

    @property
    def lm_subclass(self) -> Optional[str] : return self.__lm_subclass

    @property
    def lipids(self) -> List[Lipid] : return self.__lipids

    @property
    def require_fas(self) -> Set[Tuple[int, int]] : return self.__require_fas

    @property
    def exclude_fas(self) -> Set[Tuple[int, int]] : return self.__exclude_fas

    @property
    def exact_fas(self) -> Set[Tuple[int, int]] : return self.__exact_fas

    # @property
    # def stats(self) -> Set[Statistic] : return set().union(*[lipid.stats for lipid in self.lipids])

    # --- methods ---

    def _gen_fa_restrictions_str(self) -> str :
        """
        helper function to generate a concise string representation of whatever 
        FA restrictions are defined for this LipidGroup. Only include them if there
        are any FAs actually defined otherwise don't bother
        """
        fa_str = ""
        if len(self.require_fas) > 0:
            fa_str += " require_fas=" + "|".join([f"{c}:{u}" for c, u in self.require_fas])
        if len(self.exclude_fas) > 0:
            fa_str += " exclude_fas=" + "|".join([f"{c}:{u}" for c, u in self.exclude_fas])
        if len(self.exact_fas) > 0:
            fa_str += " exact_fas=" + "|".join([f"{c}:{u}" for c, u in self.exact_fas])
        return fa_str

    def find_filtered_lipid_stats_with_tags(self,
        contains_fas: Optional[Iterable[Tuple[int, int]]] = None,
        **stat_tags: Union[str, float, int, bool]
    ) -> Iterator[Statistic] :
        """
        Generator that yields any `Statistic`s from `Lipid` nodes contained in `self.lipids` 
        that contain the specified tags (specified tags only need to be a subset of all tags 
        present in the `Statistic`, it does not need to be an exact match). The lipids that 
        the `Statistic`s are gathered from can be optionally be filtered by FA content.

        Parameters
        ----------
        [contains_fas]
            iterable of tuples (carbons, unsaturations) for specific fatty acids that must be
            present in the gathered `Statistic`s
        """
        # collect all of the relevant Statistics from the contained lipids
        yield from itertools.chain.from_iterable([
            lipid.find_stats_with_tags(**stat_tags)
            for lipid in self.lipids
            if (
                contains_fas is None
                or len(set(contains_fas) & lipid.fas) > 0
            )
        ])