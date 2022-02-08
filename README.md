# VESIcal
A generalized python library for calculating and plotting various things related to mixed volatile (H2O-CO2) solubility in silicate melts.

[![Documentation Status](https://readthedocs.org/projects/vesical/badge/?version=latest)](https://vesical.readthedocs.io/en/latest/?badge=latest)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5095382.svg)](https://doi.org/10.5281/zenodo.5095382)

## Documentation
Check here first for all your VESIcal questions! And be sure to read the manuscripts.

   - Read all of our documentation, inlcuding quickstart guides here: https://vesical.readthedocs.io/en/latest/
   - Check our our YouTube channel for videos on how to use VESIcal here: https://www.youtube.com/channel/UCpvCCs5KMXzOxXWm0seF8Qw

### Interactive versions of manuscripts

   - Direct link to interactive VESIcal Part I manuscript: [![Manuscript on Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/kaylai/vesical-binder/HEAD?filepath=Manuscript.ipynb)
   - Jupyter Notebook hub with VESIcal: [![Manuscript on Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/kaylai/vesical-binder/HEAD)

### PDF versions of manuscripts

   - [VESIcal Part I: An Open-Source Thermodynamic Model Engine for Mixed Volatile (H2O-CO2) Solubility in Silicate Melts](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2020EA001584)
   - [VESIcal Part II: A critical approach to volatile solubility modelling using an open-source Python3 engine](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2021EA001932)


## Installation and online use

In general, VESIcal can be accessed and used in a variety of ways. From most flexible (advanced) to least flexible (novice), these are:

   - local installation of the VESIcal library
   - through the ENKI server (recommended for most users)
   - through a [jupyter notebook version of the VESIcal manuscript](https://mybinder.org/v2/gh/kaylai/vesical-binder/HEAD?filepath=Manuscript.ipynb).
   - via the web-app [https://vesical.anvil.app/](https://vesical.anvil.app/)

### VESIcal on the ENKI server

VESIcal requires installation of not only the VESIcal library but also some other python libraries, one of which is a bit tricky to install (ENKI/thermoengine aka the engine behind MELTS). But, we have a solution! All dependencies and the latest version of VESIcal are all installed on the ENKI server, within a Jupyter Notebook Hub. Steps to use VESIcal on the ENKI server are:

   1. Create a (free) GitLab account, which you'll use to sign into ENKI here: https://gitlab.com/users/sign_up
   2. Email ENKI PI Mark Ghiorso at ghiorso@ofm-research.org with your GitLab username and requet access to the ENKI server.
   3. Access the ENKI Production Server by going to http://enki-portal.org/ and clicking "SERVERS" > "PRODUCTION SERVER"
   4. Sign in with your GitLab credentials: You are now in your own jupyter notebook workspace! You can upload and create files here. They won't be accessible to anyone else. 
   5. Click the green "CLOSE THIS SCREEN" button
   6. Create a new notebook by clicking the blue plus button and then selecting Python3 under Notebook. Or select from the menu File > New > Notebook
   7. Be sure to import VESIcal as v at the top of your file, and now you are ready to get to work!

See video tutorials on our ReadTheDocs page for more.

### Installing locally

**Important! Thermoengine must be installed!**
 Please see below for details on how to install thermoengine, the python implementation of MELTS/MagmaSat  

First, obtain Python3.x if you do not already have it installed. If you are new to python, we recommend installing it via [anaconda3](https://www.anaconda.com/products/individual). VESIcal can be installed with one line. Open a terminal and type the following:

```
pip install VESIcal
```

Check that the installation worked by entering the following lines into a terminal:

```
python
import VESIcal as v
```

If no output is returned, VESIcal has installed properly! You will very likely, however, see a warning telling you that no module named 'thermoengine' could be found. The installation you performed via pip attempts to install all dependencies (other libraries that VESIcal requires), but thermoengine is not available via pip and so must be manually installed.

Dependencies that should automatically be installed for you are:

   - pandas
   - numpy
   - matplotlib
   - cycler
   - abc
   - scipy
   - sys
   - sympy
   - copy

If any warnings related to these libraries appear, try installing them as you did VESIcal: with 'pip install [package]'.

### Installing thermoengine

Thermoengine is the ENKI implementation of MELTS (MagmaSat), which is the default solubility model implemented in VESIcal. VESIcal cannot be run without thermoengine at this time, however a VESIcal-lite that does not include MagmaSat is planned. To install thermoengine, please refer to the ENKI documentation here: https://gitlab.com/ENKI-portal/ThermoEngine.

In almost all cases you will need to install thermoengine using docker. The thermoengine devs have kindly put together a docker image for you. We suggest you follow those instructions here: https://gitlab.com/ENKI-portal/ThermoEngine/-/tree/master/#running-a-container-image-locally.

## Updating

To upgrade to the most recent version of VESIcal, type the following into terminal:

```
pip install VESIcal --upgrade
```

## Contributing
Issues are tracked on [GitHub](https://github.com/kaylai/VESIcal/issues).

Patches may be submitted via a [Github pull request](https://github.com/kaylai/VESIcal/pulls). All changes should include tests (VESIcal uses python's unittest library) and pass [flake8](https://pypi.org/project/flake8/).
