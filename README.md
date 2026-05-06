# VESIcal
A generalized python library for calculating and plotting various things related to mixed volatile (H2O-CO2) solubility in silicate melts.

[![Documentation Status](https://readthedocs.org/projects/vesical/badge/?version=latest)](https://vesical.readthedocs.io/en/latest/?badge=latest)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5095382.svg)](https://doi.org/10.5281/zenodo.5095382)

## Documentation
Check here first for all your VESIcal questions! And be sure to read the manuscripts.

   - Read all of our documentation, inlcuding quickstart guides here: https://vesical.readthedocs.io/en/latest/
   - Check our our YouTube channel for videos on how to use VESIcal here: https://www.youtube.com/channel/UCpvCCs5KMXzOxXWm0seF8Qw

### Interactive versions of manuscripts

   - Curvenote/AGU Notebooks Now! implementation of VESIcal Part I manuscript: [Manuscript on Curvenote](https://agu.curve.space/articles/NN0001)
   - Direct link to interactive jupyter notebook version of VESIcal Part I manuscript: [Manuscript on Binder](https://agu-binder.curvenote.dev/user/2be900e9-fb5d-4-9778d16a48c.zip-4wrcztow/lab/tree/Manuscript.ipynb?token=EzBUfh6US4qFq4UW0MSkYA)
   - Jupyter Notebook hub with VESIcal: [![Manuscript on Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/kaylai/vesical-binder/HEAD)

### PDF versions of manuscripts

   - [VESIcal Part I: An Open-Source Thermodynamic Model Engine for Mixed Volatile (H2O-CO2) Solubility in Silicate Melts](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2020EA001584)
   - [VESIcal Part II: A critical approach to volatile solubility modelling using an open-source Python3 engine](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2021EA001932)


## Installation and online use

In general, VESIcal can be accessed and used in a variety of ways. From most flexible (advanced) to least flexible (novice), these are:

   - local installation of the VESIcal library (we have made the thermoengine install much simpler for macOS, see below!)
   - through the ENKI server (recommended for non-macOS users)
   - through a [jupyter notebook version of the VESIcal manuscript](https://agu-binder.curvenote.dev/user/2be900e9-fb5d-4-9778d16a48c.zip-4wrcztow/lab/tree/Manuscript.ipynb?token=EzBUfh6US4qFq4UW0MSkYA).
   - via the web-app [https://vesical.anvil.app/](https://vesical.anvil.app/)

### Installing locally

**Important! Thermoengine must be installed for the default MagmaSat model to function.**
 Please see below for details on how to install thermoengine, the python implementation of MELTS/MagmaSat.

#### VESIcal without thermoengine
Thermoengine is the ENKI implementation of MELTS (MagmaSat), which is the default solubility model implemented in VESIcal. You can install VESIcal with a simple `pip install VESIcal` on your local machine and run everything in VESIcal *except for the MagmaSat model*. 

```
pip install VESIcal
```

Check that the installation worked by entering the following lines into a terminal:

```
python
import VESIcal as v
```

If no output is returned, VESIcal has installed properly! Optionally, type `v.__version__`, and the VESIcal version number will print to the terminal. You will very likely, however, see a warning telling you that no module named 'thermoengine' could be found. The installation you performed via pip attempts to install all dependencies (other libraries that VESIcal requires), but thermoengine is not available via pip and so must be manually installed. Dependencies that should automatically be installed for you are listed in the requirements.txt file in the repo root. If any warnings related to those libraries appear, try installing them as you did VESIcal: with `pip install [package]`.

If thermoengine is not installed, you will see a warning after running `import VESIcal as v` notifying you as such. The only other noteable caveat is that you must pass the name of the model you want to use for every calculation, since MagmaSat is the default. Do so like this (with VESIcal imported as `v`):
```python
v.calculate_saturation_pressure(sample=mysample, temperature=mytemperature, model="IaconoMarziano")
```

### Installing thermoengine
After installing both themoengine and VESIcal, you can use all VESIcal sub-models including MagmaSat to your heart's content.

#### Simplified instructions for installing thermoengine on macOS
We have forked the original thermoengine v1 GitLab repository and made the few necessary edits for a working install on macOS. Find those instructions in the forked repo's README here: https://gitlab.com/kaylai/thermo-engine-for-mac

#### For Windows
To install thermoengine, please refer to the ENKI documentation here: https://gitlab.com/ENKI-portal/ThermoEngine.

In almost all cases you will need to install thermoengine using docker. The thermoengine devs have kindly put together a docker image for you. We suggest you follow those instructions here: https://gitlab.com/ENKI-portal/ThermoEngine/-/tree/master/#running-a-container-image-locally.

#### For Linux
The Windows installation instructions are tested to work for Ubuntu, but you might also try the simplified macOS instructions as well. That's untested, but there's a good chance it just works. If you try it and it does, let us know!

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

#### How to install Python for newcomers
First, obtain Python 3.10 or newer if you do not already have it installed. If you are new to Python, we recommend installing it via [Miniforge](https://github.com/conda-forge/miniforge), a free, community-maintained Python distribution that includes the `conda` package manager and uses the [conda-forge](https://conda-forge.org/) channel by default.

## Updating

To upgrade to the most recent version of VESIcal, type the following into terminal:

```
pip install VESIcal --upgrade
```

## Contributing
Issues are tracked on [GitHub](https://github.com/kaylai/VESIcal/issues).

Patches may be submitted via a [Github pull request](https://github.com/kaylai/VESIcal/pulls). All changes should include tests (VESIcal uses python's unittest library) and pass [flake8](https://pypi.org/project/flake8/).
