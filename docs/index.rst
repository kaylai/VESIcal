##################
VESIcal Quickstart
##################

.. image:: img/header_transparent.png

Installation
------------
See the :doc:`install` section for detailed instructions and dependencies.

Ways to use VESIcal
^^^^^^^^^^^^^^^^^^^
In general, VESIcal can be accessed and used in a variety of ways. From most flexible (advanced) to least flexible (novice), these are:

- local installation of the VESIcal library (get the full VESIcal experience)
- through the ENKI server (http://enki-portal.org/ recommended for most users)
- through a jupyter notebook version of the VESIcal manuscript (https://mybinder.org/v2/gh/kaylai/vesical-binder/HEAD?filepath=Manuscript.ipynb)
- via the web app (https://vesical.anvil.app/). The web app sometimes crashes. If so, it should restart within 1 minute. If it does not, please email kayla.iacovino at nasa dot gov!

Local installation
^^^^^^^^^^^^^^^^^^
VESIcal can be installed with pip:

.. code-block:: bash

    pip install VESIcal

Always use the most up-to-date version of the code:

.. code-block:: bash

   pip install VESIcal --upgrade

GitHub
------
Download the VESIcal source code, create issues, and more.

https://github.com/kaylai/VESIcal

Peer-reviewed Manuscripts
-------------------------

Instructions on how to cite VESIcal are :doc:`here </about>`

Interactive versions of manuscripts
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
   - Direct link to interactive VESIcal Part I manuscript: `Manuscript on Binder <https://mybinder.org/v2/gh/kaylai/vesical-binder/HEAD?filepath=Manuscript.ipynb>`_
   - Jupyter Notebook hub with VESIcal: `Manuscript on Binder <https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/kaylai/vesical-binder/HEAD>`_

PDF versions of manuscripts
^^^^^^^^^^^^^^^^^^^^^^^^^^^
   - `VESIcal Part I: An Open-Source Thermodynamic Model Engine for Mixed Volatile (H2O-CO2) Solubility in Silicate Melts <https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2020EA001584>`_
   - `VESIcal Part II: A critical approach to volatile solubility modelling using an open-source Python3 engine <https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2021EA001932>`_

News
----
VESIcal can now be used without needing to install thermoengine!
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Without thermoengine installed, MagmaSat (the default model) cannot be used, however all other models, all plotting capability, and the thermo package can be used. Simply make sure you are using version 1.2.0 or higher following the simple installation instructions below. Then be sure to explicitly tell VESIcal which model you want to use (any model other than MagmaSat. See :doc:`models` for a complete list) by passing `model="some-model-name"` when performing a calculation.

For example:

.. code-block:: python

   v.calculate_saturation_pressure(sample=<your-sample-here>, model="<some-model-name-here>")

If thermoengine is not installed and you import VESIcal, you will be warned that you won't be able to use MagmaSat, but everything else will work as expected. Remember that if you do not pass a model name, VESIcal will default to MagmaSat, and an error will be generated telling you that you need to pass a model name.

Indices and tables
^^^^^^^^^^^^^^^^^^

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. toctree::
    :maxdepth: 2
    :caption: Contents:

    install
    changelog
    about
    models
    quick_tutorials
    tutorials
    advanced_tutorials
    youtube
    integration
    workshops
    codedoc
    faq
    support
    license
