============
Installation
============

.. contents::
   :depth: 2
   :local:

VESIcal can be accessed and used in a variety of ways. From most flexible (local installation; advanced), quite flexible but slightly cumbersome (via the ENKI jupyter protal), quite flexible but potentially overly cumbersome (manuscript binder) to least flexible (the web app; best for absolute beginners)

----------------------

1. Installing locally
#####################

- Local installation of the VESIcal library
		- The best and most complete version of VESIcal
		- Recommended for users comfortable installing ENKI thermoengine (MagmaSat) via the instructions below
		- Recommended for all users that don't need to use MagmaSat 

.. important:: 
	**Thermoengine must be installed to use the MagmaSat model**
 
 	Please see below for details on how to install thermoengine, the python implementation of MELTS/MagmaSat. VESIcal will run fine without this, but you will not be able to use the default model MagmaSat. In this case, simply specify which model you would like to use (any model other than MagmaSat, see :doc:`models` for a complete list) in any calculation you perform with the argument `model="<your-desired-model>"`.

 	Nota bene: We made MagmaSat the default model for a reason. In many cases, it is the best model for the job. However, there are plenty of reasons one would wish to use a different model, which is why we have made them all available to you. We strongly recommend reading the `VESIcal Part II manuscript <https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2021EA001932>`_, which details the pros, cons, follies, and pitfalls of all models implemented in VESIcal. We also recommend using VESIcal's tools (like `calib_plot()`) to determine which model will best suit your needs.

1.1 Installing thermoengine
---------------------------

1.1.1 On Mac and probably Linux, too
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
**Grab the code and get detailed installation instructions here:** 
`https://gitlab.com/kaylai/thermo-engine-for-mac <https://gitlab.com/kaylai/thermo-engine-for-mac>`_

Thermoengine is the ENKI implementation of MELTS (MagmaSat), which is the default solubility model implemented in VESIcal. To streamline the installation process without failed builds, we have forked the original thermoengine GitLab repo into one tested on MacOS (Intel and Apple Silicon).

1.1.2 On Windows
^^^^^^^^^^^^^^^^
On Windows you will need to install thermoengine using docker. The thermoengine devs have kindly put together a docker image for you. We suggest you follow those instructions here: `https://gitlab.com/ENKI-portal/ThermoEngine/-/tree/master/#running-a-container-image-locally <https://gitlab.com/ENKI-portal/ThermoEngine/-/tree/master/#running-a-container-image-locally>`_.

Installing on Windows x64? Check out Liam Peterson's instructions on `Installing and Running thermoengine on Windows x64 <https://github.com/kaylai/VESIcal/raw/master/docs/thermoengine_local_install_Windowsx64.docx>`_

1.2 Installing VESIcal via pip
------------------------------
On all systems, simply:

.. code-block:: python

	pip install VESIcal

Check that the installation worked by entering the following lines into a terminal:

.. code-block:: python

	python
	import VESIcal as v

If no output is returned, VESIcal has installed properly! The installation you performed via pip attempts to install all dependencies (other libraries that VESIcal requires). To use VESIcal's default model MagmaSat, thermoengine must be installed separately but is not available via pip and so must be manually installed. VESIcal can be used without thermoengine, but a model must be specified for each calculation performed.

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

1.2.1 Updating
^^^^^^^^^^^^^^

To upgrade to the most recent version of VESIcal, type the following into terminal:

.. code-block:: python

	pip install VESIcal --upgrade

----------------------

2. VESIcal on the ENKI server
#############################

- Recommended for users who are less comfortable in a command line and prefer and "plug and play" experience.
- Requires a GitLab account and joining the ENKI server run by Mark Ghiorso

VESIcal requires installation of not only the VESIcal library but also some other python libraries, one of which is a bit tricky to install (ENKI/thermoengine aka the engine behind MELTS). But, we have a solution! All dependencies and the latest version of VESIcal are all installed on the ENKI server, within a Jupyter Notebook Hub. Steps to use VESIcal on the ENKI server are:

	1. Create a (free) GitLab account, which you'll use to sign into ENKI here: (`https://gitlab.com/users/sign_up <https://gitlab.com/users/sign_up>`_)
	2. Email ENKI PI Mark Ghiorso at ghiorso@ofm-research.org with your GitLab username and requet access to the ENKI server.
	3. Access the ENKI Production Server by going to `http://enki-portal.org/ <http://enki-portal.org/>`_ and clicking "SERVER"
	4. Sign in with your GitLab credentials: You are now in your own jupyter notebook workspace! You can upload and create files here. They won't be accessible to anyone else. 
	5. Click the green "CLOSE THIS SCREEN" button
	6. Create a new notebook by clicking the blue plus button and then selecting Python3 under Notebook. Or select from the menu File > New > Notebook
	7. Be sure to import VESIcal as v at the top of your file, and now you are ready to get to work!

See this video tutorial on accessing the ENKI server for more:

.. raw:: html

	<iframe width="560" height="315" src="https://www.youtube.com/embed/jUshguhFpjk" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

|

----------------------


3. Interactive Manuscripts
##########################

**Full details on the** :doc:`manuscripts` **page.**

You can access interactive versions of the VESIcal Part I manuscript via a Binder (a mechanism to bundle a jupyter notebook so that it can be accessed in the web browser without installing anything). 

It works... some times.

----------------------

4. VESIcal Web App
##################

**Go straight to the web app here:**
`https://vesical.anvil.app/ <https://vesical.anvil.app/>`_


- Best for absolute beginners to coding, particuarly in Python. 
- Upload a file, click a button, get your results!
- Does not currently offer the full functionality of VESIcal
- Does not support customizing functions and calculations

.. tip::
	If you encounter any issues, :doc:`get in touch </support>`.

----------------------
