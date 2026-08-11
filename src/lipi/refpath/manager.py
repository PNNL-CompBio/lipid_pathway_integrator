"""
lipid_pathway_integrator/refpath/manager.py
Dylan Ross (dylan.ross@pnnl.gov)

    Module with interface for interacting with reference pathway database
"""

import os
import sqlite3
import json
from importlib.resources import read_text
import warnings
from typing import Optional, Union, Literal, Tuple, Set, Mapping
from collections.abc import Iterable

from lipidimea._lipidlib.parser import parse_lipid_name

from lipi import __version__ as LIPID_PATHWAY_INTEGRATOR_VERSION


#===============================================================================
# Constants

# load the database schema, built into this package as a resource
_SCHEMA = read_text("lipid_pathway_integrator._resources", "refpath_schema.sql")


#===============================================================================
# Helper Functions

def _insert_qry(
    table: str, 
    n_binds: int
) -> str :
    """ Generate an insert query with the specified number of bindings """
    return f"""--beginsql
        INSERT INTO {table} VALUES ({("?," * n_binds)[:-1]});
    --endsql"""


def _store_package_version(
    cur: sqlite3.Cursor
) -> None :
    """
    Store the version of the `lipid_pathway_integrator` codebase that was used to initialize the database
    """
    cur.execute(_insert_qry("Versions", 2), ("lipid_pathway_integrator", LIPID_PATHWAY_INTEGRATOR_VERSION))


def _name_is_parsable_as_lipid(name: str) -> bool :
    """
    Return a boolean flag indicating whether a name was successfully parseable using
    parsing function from lipidimea, meaning that it should be associated with the 
    "LIPID" node type
    """
    # discard the resulting lipid object, 
    # just return a flag indicating successful parsing or failure
    return parse_lipid_name(name) is not None


#===============================================================================
# Main Database Interface Class

class RefPathManager:
    """
    Interface for interacting with the reference pathway database

    Attributes
    ----------
    db_path
        path to database file
    connected
        a boolean flag indicating whether the connection to the 
        database is still active
    """

    def __init__(self,
        refpath_db_file: Optional[str] = None
    ) : 
        """ 
        Connect to an existing reference pathway database, or create one in memory

        Parameters
        ----------
        refpath_db_file
            path to the reference pathway database file or `None` to create one in-memory
        """
        # store database path
        self.__db_path: str = refpath_db_file if refpath_db_file is not None else ":memory:"
        # ensure the database file exists already
        if self.db_path != ":memory:" and not os.path.isfile(self.db_path):
            raise FileNotFoundError(self.db_path)
        # connect to the database
        self.__con: sqlite3.Connection = sqlite3.connect(self.db_path)
        cur: sqlite3.Cursor = self.__con.cursor()
        self.__connected: bool = True
        # enforce foreign key constraints
        cur.execute("PRAGMA foreign_keys=1")
        # if in-memory database, then initialize with the schema
        if self.db_path == ":memory:":
            cur.executescript(_SCHEMA)
            _store_package_version(cur)
        # Check the codebase version, this is important because the database schema can change.
        self._check_codebase_version()
        # connect the helper classes with grouped methods
        self._connect_helpers()

    # --- static methods ---

    @staticmethod
    def create(
        refpath_file: str,
        overwrite: bool = False
    ) -> None :
        """
        creates a sqlite database for molecular signature library using built in schema

        raises a RuntimeError if the database already exists (unless overwrite flag set)

        Parameters
        ----------
        siglib_file 
            filename/path of the signature library database
        overwrite
            if the database file already exists and this flag is True, then overwrite existing database 
            and do not raise the RuntimeError
        """
        # see if the file exists
        if os.path.exists(refpath_file):
            if overwrite:
                os.remove(refpath_file)
            else:
                msg = f"create: database file ({refpath_file}) already exists"
                raise FileExistsError(msg)
        # ensure the directory for the database file exists, create it if not
        if (dir := os.path.dirname(refpath_file)) != "":
            # No need to mkdir if the dirname is just "" (i.e. when connecting to a file in the
            # current workding directory)
            os.makedirs(dir, exist_ok=True)
        # initial connection creates the DB
        con = sqlite3.connect(refpath_file)  
        cur = con.cursor()
        # Read the schema file from package resources
        cur.executescript(_SCHEMA)
        # keep track of the codebase version used to initialize the database
        _store_package_version(cur)
        # save and close the database
        con.commit()
        con.close()

    # --- attribute getters ---

    @property
    def db_path(self) -> str : return self.__db_path

    @property
    def connected(self) -> bool : return self.__connected

    # --- methods: miscellaneous --- 

    def _connect_helpers(self):
        """ connect all of the helper classes with grouped methods """
        # Reinitialize helpers with new connection
        self.identifiers = self._IdentifierMethods(self.__con)
        self.entities = self._EntityMethods(self.__con, self.identifiers)
        self.nodes = self._NodeMethods(self.__con)
        self.edges = self._EdgeMethods(self.__con)
        self.pathways = self._PathwayMethods(self.__con)

    def get_read_only_connection(self) -> sqlite3.Connection :
        """ expose a read-only database connection """
        # open a separate connection in read-only mode
        return sqlite3.connect(f"file:{self.__db_path}?mode=ro", uri=True)

    # pass through database connection commit and close
    def commit(self): self.__con.commit()
    
    def close(self): 
        self.__con.close()
        self.__connected = False

    def reconnect(self):
        """ enable reconnecting to a previously close database """
        if not self.connected:
            self.__con = sqlite3.connect(self.db_path)
            self.__connected = True
            # reinitialize helpers with new database connection
            self._connect_helpers()

    def _check_codebase_version(self):
        """ ensure the version is matched between the codebase and database, warn if not """
        cur = self.__con.cursor()
        qry = """--beginsql
            SELECT version FROM Versions WHERE component='lipid_pathway_integrator';
        --endsql"""
        if (db_ver := cur.execute(qry).fetchone()[0]) != LIPID_PATHWAY_INTEGRATOR_VERSION:
            msg = f"RefPathManager: version mismatch between codebase ({LIPID_PATHWAY_INTEGRATOR_VERSION}) and database ({db_ver})"
            warnings.warn(msg)

    # --- methods: identifiers  ---

    class _IdentifierMethods:
        """ helper class to group methods related to Identifiers table """

        def __init__(self, con: sqlite3.Connection):
            # capture the private database connection in the outer RefPathManager object
            self.__con = con
        
        def get_entity_id(self,
            name: str,
            **tags: Union[str, int, float, bool]
        ) -> Optional[int] :
            """ try to lookup identifier group ID for an identifier entry, if no matching entry exists returns None """
            cur = self.__con.cursor()
            jtags = json.dumps(tags)
            qry = """--beginsql
                SELECT entity_id FROM Identifiers WHERE name=? AND tags=?
            --endsql"""
            if (entity_id := cur.execute(qry, (name, jtags)).fetchone()) is not None:
                return entity_id[0]
            return None
        
        def insert_or_ignore(self,
            name: str,
            entity_id: int,
            **tags: Union[str, int, float, bool]
        ) -> None : 
            """ tries to insert an entry into Identifiers table, silently ignores non-unique entries  """
            cur = self.__con.cursor()
            jtags = json.dumps(tags)
            try: 
                cur.execute(_insert_qry("Identifiers", 3), (name, jtags, entity_id))
            except sqlite3.IntegrityError as ie:
                if "UNIQUE" not in str(ie):
                    # re-raise any errors besides a uniqueness constraint violation
                    # which just gets silently ignored
                    raise ie

    # --- methods: entities ---

    class _EntityMethods:
        """ helper class to group methods related to Entities table """

        def __init__(self, 
            con: sqlite3.Connection,
            identifiers: RefPathManager._IdentifierMethods
        ):
            # capture the private database connection from the outer RefPathManager object
            self.__con = con
            # capture the identifier methods object from the outer RefPathManager object
            self.__identifiers = identifiers

        def insert(self,
            entity_type: Literal["PROTEIN", "METABOLITE", "LIPID", "LIPID_GROUP"]
        ) -> int : 
            """ inserts an entry into Entities, returns entity_id """
            cur = self.__con.cursor()
            cur.execute(_insert_qry("Entities", 2), (None, entity_type))
            assert cur.lastrowid is not None
            return cur.lastrowid
        
        def _infer_entity_type_from_identifiers(self,
            identifiers: Iterable[
                Tuple[
                    str, 
                    Mapping[str, Union[str, int, float, bool]]
                ]
            ], 
        ) -> Literal["PROTEIN", "METABOLITE", "LIPID", "LIPID_GROUP"] :
            """ 
            infer the node type associated with a group of identifiers based on presence 
            of various identifier tags 
            """
            # flags for different tags that indicate node type they should associate with
            is_protein = False
            is_lipid_group = False
            is_lipid = False
            for name, tags in identifiers:
                # take these flag values directly from the tags dict if present
                is_protein = is_protein or v if (v := tags.get("is_protein")) is not None else is_protein 
                is_lipid_group = is_lipid_group or v if (v := tags.get("is_lipid_group")) is not None else is_lipid_group  
                # is_lipid flag is inferred based on whether a name is parsable as a lipid (using pygoslin)
                # only try to parse normal names, not external identifiers, i.e. ignore any 
                # that have the "source" tag
                is_lipid = is_lipid or (tags.get("source") is None and _name_is_parsable_as_lipid(name))
            # determine the node type based on the flags
            match (is_protein, is_lipid_group, is_lipid):
                case (True, False, False) | (True, False, True):
                    # the is_protein flag takes precedence over the is_lipid flag
                    return "PROTEIN"
                case (False, True, False) | (False, True, True):
                    # the is_lipid_group flag takes precedence over the is_lipid flag
                    return "LIPID_GROUP"
                case (False, False, True):
                    return "LIPID"
                case (False, False, False):
                    # if none of the flags are set, assume metabolite
                    return "METABOLITE"
                case _:
                    # any other combination that is not explicitly accounted for is an error
                    # this may signal an improper grouping of unrelated identifiers
                    raise RuntimeError(
                        "RefPathManager._Entities._infer_entity_type_from_identifiers: "
                        "bad combination of flags: "
                        f"{is_protein=}, {is_lipid_group=}, {is_lipid=}"
                    )

        def insert_or_update_from_identifiers(self,
            identifiers: Iterable[
                Tuple[
                    str, 
                    Mapping[str, Union[str, int, float, bool]]
                ]
            ], 
        ) -> int :
            """
            Take a group of identifiers, first check if there are any existing identifier
            entries and grab the entity_id from them, if not, make a new entity and 
            add all of the identifiers using that. Returns the entity id
            
            Returns
            -------
            entity_id
            """
            # consume once in case it's a generator
            # because we need to iterate over it twice
            identifiers = list(identifiers)
            entity_id = None
            for name, tags in identifiers:
                if (existing_entity_id := self.__identifiers.get_entity_id(name, **tags)) is not None:
                    entity_id = existing_entity_id
            if entity_id is None:
                entity_id = self.insert(
                    self._infer_entity_type_from_identifiers(identifiers)
                )
            # now we have a valid entity_id, make sure to add all of the identifiers
            # (already existing entries silently skipped)
            for name, tags in identifiers:
                self.__identifiers.insert_or_ignore(name, entity_id, **tags)
            return entity_id    
        
    # === nodes ===

    class _NodeMethods:
        """ helper class to group methods related to Nodes table """

        def __init__(self, con: sqlite3.Connection):
            # capture the private database connection from the outer RefPathManager object
            self.__con = con

        def insert(self,
            pathway_id: int,
            entity_id: int,
            x: Optional[float] = None,
            y: Optional[float] = None
        ) -> int :
            """
            insert a node 

            Parameters
            ----------
            pathway_id
                pathway to add node to
            entity_id
                the underlying entity
            x, y
                optional position
            
            Returns
            -------
            node_id
            """
            cur = self.__con.cursor()
            cur.execute(
                _insert_qry("Nodes", 5),
                (None, entity_id, pathway_id, x, y)
            )
            assert cur.lastrowid is not None
            return cur.lastrowid
        
        def fetch_data(self,
            node_ids: Iterable[int]
        ) -> Iterable[Tuple[int, str, int, Iterable[str], Iterable[str], Optional[float], Optional[float]]] :
            """
            fetch node data from specified node IDs

            Yields
            ------
            node_id
            node_type
            entity_id
            entity_names
            entity_tags
            x
            y
            """
            # the node_ids parameter is an arbitrary iterable of ints, this could be a set in some cases
            # but execute does not like a set for query parameters, so convert to a list just in case
            # this also makes sure len(node_ids) works as needed 
            node_ids = list(node_ids)
            cur = self.__con.cursor()
            qry = """--beginsql
                SELECT 
                    node_id,
                    type,
                    entity_id,
                    GROUP_CONCAT(name, '|'),
                    GROUP_CONCAT(tags, '|'),
                    x,
                    y
                FROM 
                    Nodes
                    JOIN Entities USING(entity_id)
                    JOIN Identifiers USING(entity_id)
                WHERE 
                    node_id IN ({})
                GROUP BY 
                    node_id
            --endsql""".format(','.join('?' * len(node_ids)))
            for node_id, node_type, entity_id, entity_names, entity_tags, x, y in cur.execute(qry, node_ids):
                yield (
                    node_id, 
                    node_type, 
                    entity_id,
                    entity_names.split("|"),
                    entity_tags.split("|"),
                    x, 
                    y
                )

    # === edges ===

    class _EdgeMethods:
        """ helper class to group methods related to Edges table """

        def __init__(self, con: sqlite3.Connection):
            # capture the private database connection from the outer RefPathManager object
            self.__con = con
        
        def insert_or_ignore(self,
            source_node_id: int,
            target_node_id: int,
            edge_type: Literal["REACTION", "ACTIVATION", "INHIBITION"]
        ) -> int :
            """ try to insert an entry into the Edges table, returns edge_id """
            cur = self.__con.cursor()
            try: 
                cur.execute(_insert_qry("Edges", 4), (None, source_node_id, target_node_id, edge_type))
                assert cur.lastrowid is not None
                return cur.lastrowid
            except sqlite3.IntegrityError as ie:
                if "UNIQUE" in str(ie):
                    # uniqueness constraint violation
                    # lookup the node ID based on name and tags
                    qry = """--beginsql
                        SELECT edge_id FROM Edges WHERE source_node=? AND target_node=? AND type=?;
                    --endsql"""
                    # assume that the lookup is successful because we got a uniqueness constraint violation
                    return cur.execute(qry, (source_node_id, target_node_id, edge_type)).fetchone()[0]
                else:
                    # re-raise any other errors
                    raise ie
    
        def fetch_data(self,
            edge_ids: Iterable[int]
        ) -> Iterable[Tuple[int, int, str]] :
            """
            fetch edge data from specified edge IDs

            Yields
            ------
            source_node_id
            target_node_id
            edge_type
            """
            # the edge_ids parameter is an arbitrary iterable of ints, this could be a set in some cases
            # but execute does not like a set for query parameters, so convert to a list just in case
            # this also makes sure len(edge_ids) works as needed 
            edge_ids = list(edge_ids)
            cur = self.__con.cursor()
            qry = """--beginsql
                SELECT 
                    source_node,
                    target_node,
                    type
                FROM 
                    Edges
                WHERE 
                    edge_id IN ({})
            --endsql""".format(','.join('?' * len(edge_ids)))
            yield from cur.execute(qry, edge_ids)

    # === pathways ===

    class _PathwayMethods:
        """ helper class to group methods related to Pathways table """

        def __init__(self, con: sqlite3.Connection):
            # capture the private database connection from the outer RefPathManager object
            self.__con = con
    
        def insert_or_ignore(self,
            name: str,
            **tags: Union[str, int, float, bool]
        ) -> int : 
            """ 
            tries to insert an entry into Pathways table, returns pathway_id 
            does not insert and returns existing pathway_id if matching entry
            was already present
            """
            cur = self.__con.cursor()
            jtags = json.dumps(tags)
            try: 
                cur.execute(_insert_qry("Pathways", 3), (None, name, jtags))
                assert cur.lastrowid is not None
                return cur.lastrowid
            except sqlite3.IntegrityError as ie:
                if "UNIQUE" in str(ie):
                    # uniqueness constraint violation
                    # lookup the pathway ID based on name and tags
                    qry = """--beginsql
                        SELECT pathway_id FROM Pathways WHERE name=? AND tags=?;
                    --endsql"""
                    # assume that the lookup is successful because we got a uniqueness constraint violation
                    return cur.execute(qry, (name, jtags)).fetchone()[0]
                # re-raise any other errors
                raise ie
            
        def fetch_node_and_edge_ids(self,
            pathway_id: int
        ) -> Tuple[Set[int], Set[int]] :
            """
            fetch all edge and node IDs assiciated with a specified pathway

            Returns
            -------
            edge_ids
            node_ids
            """
            cur = self.__con.cursor()
            qry_n = """--beginsql
                SELECT 
                    node_id 
                FROM 
                    Nodes
                WHERE 
                    pathway_id=?;
            --endsql"""
            node_ids = set([_[0] for _ in cur.execute(qry_n, (pathway_id,)).fetchall()])
            q_marks = ','.join('?' * len(node_ids))
            qry_e = """--beginsql
                SELECT
                    edge_id
                FROM 
                    Edges
                WHERE 
                    source_node IN ({})
                    AND target_node IN ({})
            --endsql""".format(q_marks, q_marks)
            # this looks a little icky, there's probably a better way to do this...
            edge_ids = set([_[0] for _ in cur.execute(qry_e, list(node_ids) + list(node_ids)).fetchall()])
            # combine source and target node IDs 
            return node_ids, edge_ids
    