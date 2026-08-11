-- 
CREATE TABLE Versions (
    -- component, required field
    component TEXT NOT NULL,  
    -- version, required field
    version TEXT NOT NULL,  
    -- ensures UNIQUE (component, version)
    PRIMARY KEY (component, version) 
);


-- 
CREATE TABLE Entities (
    -- unique identifier for the entity
    entity_id INTEGER PRIMARY KEY,
    -- type of node (restricted to certain values)
    type TEXT NOT NULL,
    -- restricts type to specific values
    CHECK (type IN ('PROTEIN', 'METABOLITE', 'LIPID', 'LIPID_GROUP'))
);


--
CREATE TABLE Identifiers (
    -- name, required field 
    -- (part of composite primary key)
    name TEXT NOT NULL,  
    -- tags stored as JSON; default is empty object "{}"
    -- (part of composite primary key)
    tags TEXT NOT NULL DEFAULT '{}',
    -- link to an entity, all identifier entries must map to 
    -- exactly one entity
    entity_id INT NOT NULL,
    -- use name and tags as composite primary key
    -- ensures uniqueness for the combination of name and tags
    PRIMARY KEY (name, tags),
    -- enforces integrity with Entities
    FOREIGN KEY (entity_id) REFERENCES Entities (entity_id)
        ON UPDATE CASCADE
        ON DELETE CASCADE   
);


--
CREATE TABLE Pathways (
    -- unique identifier as the primary key
    pathway_id INTEGER PRIMARY KEY,  
    -- name, required field
    name TEXT NOT NULL,  
    -- tags stored as JSON; default is empty object "{}"
    tags TEXT NOT NULL DEFAULT '{}',
    -- ensures uniqueness for the combination of name and tags
    UNIQUE (name, tags)  
);


CREATE TABLE Nodes (
    node_id INTEGER PRIMARY KEY,
    -- foreign key referencing the entity_id column in Entities
    entity_id INTEGER NOT NULL,
    -- foreign key referencing the pathway_id column in Pathways
    pathway_id INTEGER NOT NULL,
    x REAL DEFAULT NULL,
    y REAL DEFAULT NULL,
    -- enforces integrity with Entities table
    FOREIGN KEY (entity_id) REFERENCES Entities (entity_id)
        ON UPDATE CASCADE 
        ON DELETE CASCADE,
    -- enforces integrity with Pathways table
    FOREIGN KEY (pathway_id) REFERENCES Pathways (pathway_id)
        ON UPDATE CASCADE 
        ON DELETE CASCADE
);


-- 
CREATE TABLE Edges (
    -- Unique identifier for the edge
    edge_id INTEGER PRIMARY KEY,
    -- foreign key referencing the node_id column in Nodes for the source node
    source_node INTEGER NOT NULL,
    -- foreign key referencing the node_id column in Nodes for the target node
    target_node INTEGER NOT NULL,
    -- type of edge (restricted to certain values)
    type TEXT NOT NULL,
    -- ensures uniqueness for the combination of source_node, target_node, and type
    UNIQUE (source_node, target_node, type),
    -- enforces integrity with the Nodes table for source_node
    FOREIGN KEY (source_node) REFERENCES Nodes (node_id)
        ON UPDATE CASCADE
        ON DELETE CASCADE,
    -- enforces integrity with the Nodes table for target_node
    FOREIGN KEY (target_node) REFERENCES Nodes (node_id)
        ON UPDATE CASCADE
        ON DELETE CASCADE,
    -- restricts type to specific values
    CHECK (type IN ('REACTION', 'ACTIVATION', 'INHIBITION'))
);