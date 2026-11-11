# LiPI
<b>Li</b>pidomics   <b>P</b>athway   <b>I</b>ntegrator

- Unifies lipid ontology-aware representation and pathway-based multi-omics integration in a single framework
- Ingests existing reference pathways from KEGG and supports defining custom pathways via YAML
- Pathways are represented as editable graph-structured data specifically designed to accommodate individual molecular lipid species
- Feature-level statistics from any -omics data type can be mapped to pathway nodes and visualized to support interpretation


## Installation

> [!NOTE]
> Requires Python version 3.14+

Examples are provided for different installation methods using either `pip` or [`uv`](https://docs.astral.sh/uv/)

### from PyPI

```sh
# install (pip)
pip install lipid-pathway-integrator

# install (uv)
uv add lipid-pathway-integrator
```

### from source

```sh
# clone repository
git clone https://github.com/PNNL-CompBio/lipid_pathway_integrator.git

# install (pip)
pip install ./lipid_pathway_integrator/

# install (uv)
uv add ./lipid_pathway_integrator/
```

### setting up for local development

```sh
# clone repository
git clone https://github.com/PNNL-CompBio/lipid_pathway_integrator.git
cd ./lipid_pathway_integrator/

# initialize the development environment (uv)
uv sync
```


## Documentation
- [https://lipi.readthedocs.io/](https://lipi.readthedocs.io/)

## Demos
- [GBM (CPTAC3)](demos/cptac_3_gbm/demo.ipynb)
- [AML (CPTAC4)](demos/cptac_4_aml/demo.ipynb)

## License & Disclaimer
- [license (BSD)](license.txt)
- [disclaimer](disclaimer.txt)

## Contributors:
- Dylan Ross (dylan.ross@pnnl.gov)

