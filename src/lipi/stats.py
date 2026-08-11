"""
lipid_pathway_integrator/stats.py
Dylan Ross (dylan.ross@pnnl.gov)

    Define `Statistic` dataclass
"""


from dataclasses import dataclass
import json
from typing import Union, Tuple, List

import polars as pl

from lipi.identifiers import IdentifierGroup


@dataclass(frozen=True)
class Statistic:
    """ 
    dataclass for an individual statistic
    
    The statistic has a value (`float`) and tags and dictionary with optional metadata 
    (a JSON string) to assign things like sources This allows relevant statistics to be 
    selected from collections when multiple are present.
    """
    value: float
    tags: str

    @staticmethod
    def from_kwargs(
        value: float, 
        **kwargs: Union[str, float, int, bool]
    ) -> Statistic : 
        """
        initialize a Statistic using kwargs to set tags

        Parameters
        ----------
        value
            value for the statistic
        **kwargs
            set tags for the statistic using kwargs,
            they will be packed into a JSON string
        """
        return Statistic(value, json.dumps(kwargs))
    
    def add_tags(self,
        **kwargs: Union[str, float, int, bool]
    ) -> Statistic :
        """
        Add additional tags to a `Statistic`, as kwargs, returns
        the new `Statistic` instance. Any newly added tags that 
        were already present will have their values overridden.
        """
        return Statistic.from_kwargs(
            self.value,            
            **(kwargs | json.loads(self.tags))
        )


def load_stats_from_csv(
    csv_f: str,
    **tags
) -> Tuple[List[IdentifierGroup], List[Statistic]] :
    """
    Load some statistics from a `.csv` file returning a collection of 
    `IdentifierGroups` and `Statistics` that can be used with the 
    `lipid_pathway_integrator.mapping.map_omics_data_to_pathway` function. 
    
    Expect two columns: 
    - metabolite/lipid/protein identifier
    - statistic value

    The statistic can be any continuous value (_e.g._ fold-change, z-score, p-value, ...)

    Keyword args are added as tags to the resulting `Statistic`s. The `stat_type` tag
    is recommended to distinguish, for instance, a z-score from a p-value assigned to 
    the same node.

    Parameters
    ----------
    csv_f
        filename of input CSV
    tags
        kwargs specifying tags to assign to the `Statistic`s

    Returns
    -------
    identifier_groups
        list of identifier groups (from the names in column 1 of the CSV)
    statistics
        list of statistics (from the values in column 2 of the CSV)
    """
    igroups, stats = [], []
    for name, value in pl.read_csv(csv_f).iter_rows():
        igroups.append(IdentifierGroup.from_name(name))
        stats.append(Statistic.from_kwargs(
            value, 
            **tags
        ))
    return igroups, stats
