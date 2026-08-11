"""
lipiidentifiers/__init__.py
Dylan Ross (dylan.ross@pnnl.gov)

    Utility for dealing with various eqivalent identifiers.

    The main component is a class, `IdentifierGroup`, which represents groups 
    of equivalent identifiers (class `Identifier`)
"""


from dataclasses import dataclass
from typing import Set, Optional, Union
from collections.abc import Iterable
import json


@dataclass(frozen=True)
class Identifier:
    """ 
    dataclass for an individual identifier
    
    Identifiers have a name and dictionary with optional tags (JSON string) 
    to assign things like sources
    """
    name: str
    tags: str

    @staticmethod
    def from_kwargs(
        name: str, 
        **kwargs: Union[str, float, int, bool]
    ) -> Identifier : 
        """
        initialize a Statistic using kwargs to set tags

        Parameters
        ----------
        name
            name for the identifier group
        **kwargs
            set tags for the identifier using kwargs,
            they will be packed into a JSON string
        """
        return Identifier(name, json.dumps(kwargs))
    

class IdentifierGroup:
    """
    A collection of `Identifier`s that are all equivalent to one another
    
    Attributes
    ----------
    name : `str` 
        A single name that can serve as a default name for this group.
    components : `set(Identifier)`
        The collection of equivalent individual identifiers
    entity_id: `int or None`
        optional entity ID associated with the identifier
    """

    def __init__(self,
        display_name: str, 
        identifiers: Iterable[Identifier],
        entity_id: Optional[int] = None
    ) :
        """  
        Initialize an `IdentifierGroup` from an iterable of equivalent 
        individual `Identifier`s

        Parameters
        ----------
        display_name
            a single display name for the group
        identifiers
            equivalent identifiers, must not be empty
        [entity_id]
            associate this identifier with an entity in reference pathway
        """
        self.__name: str = display_name
        self.__entity_id: Optional[int] = entity_id
        if len(identifiers) < 1:  # type: ignore
            raise ValueError("identifiers must have at least 1 element")
        self.__components: Set[Identifier] = set(identifiers)

    def __eq__(self, other):
        match other:
            case IdentifierGroup():
                # run self.contains on all of the components from the other group
                return any(map(self.contains, other.components))
            case Identifier():
                return self.contains(other)
            case str():
                # for simple string assume no tags and create an identifier
                return self.contains(Identifier(other, "{}"))
        return NotImplemented
    
    def __repr__(self):
        return (
            f"IdentifierGroup(igroup_id={self.entity_id}, "
            f"{", ".join(map(str, self.components))})"
        )
    
    def __str__(self): return self.name

    # --- attribute getters ---

    @property
    def components(self) -> Set[Identifier] : return self.__components

    @property
    def name(self) -> str : return self.__name

    @property
    def entity_id(self) -> Optional[int] : return self.__entity_id
        
    @property
    def size(self) -> int : return len(self.components)

    # --- static methods ---

    @staticmethod
    def from_name(
        name: str
    ) -> IdentifierGroup :
        """
        Convenience method to create an identifer group from a single name without tags
        """
        return IdentifierGroup(name, [Identifier(name, "{}")])

    # --- methods ---        

    def contains(self, 
            identifier: Identifier
        ) -> bool :
        """ Test if an identifier is contained in this group """
        tags = json.loads(identifier.tags)
        # ignore "display" tag 
        _ = tags.pop("display", None)
        # ignore "is_protein" tag
        _ = tags.pop("is_protein", None)
        for component in self.components:
            if component.name == identifier.name:
                # matched name, check tags
                ctags = json.loads(component.tags)
                # ignore "display" tag
                _ = ctags.pop("display", None)
                # ignore "is_protein" tag
                _ = ctags.pop("is_protein", None)
                # ignore
                if tags == ctags:
                    return True
        # fallthrough, no match found
        return False
   

