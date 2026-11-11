LiPI
=================

Lipid Pathway Integrator

- Unifies lipid ontology-aware representation and pathway-based multi-omics integration in a single framework
- Ingests existing reference pathways from KEGG and supports defining custom pathways via YAML
- Pathways are represented as editable graph-structured data specifically designed to accommodate individual molecular lipid species
- Feature-level statistics from any -omics data type can be mapped to pathway nodes and visualized to support interpretation

Installation
------------

.. note::

   Requires Python version 3.14+

Examples are provided for different installation methods using either ``pip`` or `uv <https://docs.astral.sh/uv/>`_.

from PyPI
~~~~~~~~~

.. code-block:: sh

   # install (pip)
   pip install lipid-pathway-integrator

   # install (uv)
   uv add lipid-pathway-integrator

from source
~~~~~~~~~~~

.. code-block:: sh

   # clone repository
   git clone https://github.com/PNNL-CompBio/lipid_pathway_integrator.git

   # install (pip)
   pip install ./lipid_pathway_integrator/

   # install (uv)
   uv add ./lipid_pathway_integrator/

setting up for local development
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: sh

   # clone repository
   git clone https://github.com/PNNL-CompBio/lipid_pathway_integrator.git
   cd ./lipid_pathway_integrator/

   # initialize the development environment (uv)
   uv sync


.. toctree::
   :maxdepth: 6
   :caption: Contents:

   api/index

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`