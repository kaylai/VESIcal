```{image} img/header_transparent.png
```

# What is VESIcal?

VESIcal is a framework for thermodynamic modeling of magmatic volatiles written in Python. As such, VESIcal provides a standard way to interact with multiple (currrently seven) published volatile solubility models. This allows a user to:

- Run automatic caluclations on large datasets
- Easily compare models using their own data
- Interrogate choices made by model authors
- Transform geochemical data (silicate liquids and H-O-C fluids): convert between units, normalize compositions
- Make plots like isobar diagrams and degassing paths

## Some Common Use Cases

- Calculate pressures from melt inclusions
- Calculate and plot magma degassing paths
- Calculate equilibrium state of liquid-vapor system
- Calculate the density and viscosity of silicate liquids

## Installation

See the {doc}`install` section for detailed instructions and dependencies.

### Ways to use VESIcal

In general, VESIcal can be accessed and used in a variety of ways. From most flexible (advanced) to least flexible (novice), these are:

- local installation of the VESIcal library (get the full VESIcal experience)
- through the ENKI server (<http://enki-portal.org/> recommended for most users)
- through a jupyter notebook version of the VESIcal manuscript (<https://mybinder.org/v2/gh/kaylai/vesical-binder/HEAD?filepath=Manuscript.ipynb>)
- via the web-app (<https://vesical.anvil.app/>)

### Local installation

VESIcal can be installed with pip:

```bash
pip install VESIcal
```

Always use the most up-to-date version of the code:

```bash
pip install VESIcal --upgrade
```

## GitHub

Download the VESIcal source code, create issues, and more.

<https://github.com/kaylai/VESIcal>

## Peer-reviewed Manuscripts

Instructions on how to cite VESIcal are {doc}`here </about>`

### Interactive versions of manuscripts

> - [Manuscript on Binder](https://mybinder.org/v2/gh/kaylai/vesical-binder/HEAD?filepath=Manuscript.ipynb) - Direct link to interactive VESIcal Part I manuscript
> - [Jupyter Hub](https://mybinder.org/v2/gh/kaylai/vesical-binder/HEAD) - Jupyter Notebook hub with VESIcal:

### PDF versions of manuscripts

> - [VESIcal Part I: An Open-Source Thermodynamic Model Engine for Mixed Volatile (H2O-CO2) Solubility in Silicate Melts](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2020EA001584)
> - [VESIcal Part II: A critical approach to volatile solubility modelling using an open-source Python3 engine](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2021EA001932)

- {ref}`genindex`

% * :ref:`modindex`

% * :ref:`search`

```{toctree}
:caption: Getting Started
:hidden: true
:maxdepth: 2

install
models
```

```{toctree}
:caption: Tutorials
:hidden: true
:maxdepth: 2

quick_tutorials
tutorials
advanced_tutorials
youtube
integration
```

```{toctree}
:caption: Reference
:hidden: true
:maxdepth: 2

codedoc
changelog
news
about
license
```

```{toctree}
:caption: Support
:hidden: true
:maxdepth: 2

support
faq
```

```{toctree}
:caption: Community Resources
:hidden: true
:maxdepth: 2

workshops
GitHub Repo <https://github.com/kaylai/VESIcal>
Link Tree <https://linktr.ee/VESIcal>
```
