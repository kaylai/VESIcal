============
Installation
============

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

thermoengine is the ENKI implementation of MELTS (MagmaSat), which is the backbone of the entire VESIcal library. VESIcal cannot be run without thermoengine at this time, however a VESIcal-lite that does not include MagmaSat is planned. To install thermoengine, please refer to the ENKI documentation here: `https://gitlab.com/ENKI-portal/ThermoEngine <https://gitlab.com/ENKI-portal/ThermoEngine>`_.

Updating
========

To upgrade to the most recent version of VESIcal, type the following into terminal:

.. code-block:: python

   pip install VESIcal --upgrade