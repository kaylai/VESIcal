============
Installation
============

Many Ways to Use VESIcal
########################
VESIcal requires installation of not only the VESIcal library but also some other python libraries, one of which is a bit tricky to install (ENKI/thermoengine aka the engine behind MELTS). But, we have a solution! All dependencies and the latest version of VESIcal are all installed on the ENKI server, within a Jupyter Notebook Hub. Steps to use VESIcal on the ENKI server are:

	1. Create a (free) GitLab account, which you'll use to sign into ENKI here: (`https://gitlab.com/users/sign_up <https://gitlab.com/users/sign_up>`_)
	2. Email ENKI PI Mark Ghiorso at ghiorso@ofm-research.org with your GitLab username and requet access to the ENKI server.
	3. Access the ENKI Production Server by going to `http://enki-portal.org/ <http://enki-portal.org/>`_ and clicking "SERVERS" > "PRODUCTION SERVER"
	4. Sign in with your GitLab credentials: You are now in your own jupyter notebook workspace! You can upload and create files here. They won't be accessible to anyone else. 
	5. Click the green "CLOSE THIS SCREEN" button
	6. Create a new notebook by clicking the blue plus button and then selecting Python3 under Notebook. Or select from the menu File > New > Notebook
	7. Be sure to import VESIcal as v at the top of your file, and now you are ready to get to work!

In general, VESIcal can be accessed and used in a variety of ways. From most flexible (advanced) to least flexible (novice), these are:

	- local installation of the VESIcal library
	- through a jupyter notebook hosted online (the VESIcal manuscript will be in the form of a jupyter notebook once published)
	- via the web-app (`https://vesical.anvil.app/ <https://vesical.anvil.app/>`_)

This documentation mainly serves to help users interact with VESIcal as a python library. That is, the user will instal VESIcal onto their machine (or run VESIcal in a jupyter notebook) and type python code to execute VESIcal.

Installing locally
##################

**Important! Thermoengine must be installed!**
 Please see below for details on how to install thermoengine, the python implementation of MELTS/MagmaSat  

First, obtain Python3.x if you do not already have it installed. If you are new to python, we recommend installing it via `anaconda3 <https://www.anaconda.com/products/individual>`_. VESIcal can be installed with one line. Open a terminal and type the following:

.. code-block:: python

   pip install VESIcal

Check that the installation worked by entering the following lines into a terminal:

.. code-block:: python

   python
   import VESIcal as v

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

Installing thermoengine
#######################

Thermoengine is the ENKI implementation of MELTS (MagmaSat), which is the default solubility model implemented in VESIcal. VESIcal cannot be run without thermoengine at this time, however a VESIcal-lite that does not include MagmaSat is planned. To install thermoengine, please refer to the ENKI documentation here: `https://gitlab.com/ENKI-portal/ThermoEngine <https://gitlab.com/ENKI-portal/ThermoEngine>`_.

In almost all cases you will need to install thermoengine using docker. The thermoengine devs have kindly put together a docker image for you. We suggest you follow those instructions here: `https://gitlab.com/ENKI-portal/ThermoEngine/-/tree/master/#running-a-container-image-locally <https://gitlab.com/ENKI-portal/ThermoEngine/-/tree/master/#running-a-container-image-locally>`_.

Updating
########

To upgrade to the most recent version of VESIcal, type the following into terminal:

.. code-block:: python

   pip install VESIcal --upgrade
