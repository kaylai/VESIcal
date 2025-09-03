.. image:: img/header_transparent.png

###############################################
VESIcal: Open-source volatile solubility engine
###############################################

VESIcal is a framework for thermodynamic modeling of magmatic volatiles written in Python. As such, VESIcal provides a standard way to interact with multiple (currrently seven) published volatile solubility models. This allows a user to:
- Run automatic caluclations on large datasets
- Easily compare models using their own data
- Interrogate choices made by model authors
- Transform geochemical data (silicate liquids and H-O-C fluids): convert between units, normalize compositions
- Make plots like isobar diagrams and degassing paths

Some Common Use Cases
---------------------
- Calculate pressures from melt inclusions
- Calculate and plot magma degassing paths
- Calculate equilibrium state of liquid-vapor system
- Calculate the density and viscosity of silicate liquids

Installation
------------
See the :doc:`install` section for detailed instructions and dependencies.

More Important Stuff
--------------------
- `VESIcal github <https://github.com/kaylai/VESIcal>`_
- :doc:`Peer-reviewed manuscripts </manuscripts>`


.. Indices and tables
.. ^^^^^^^^^^^^^^^^^^

.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`

.. toctree::
   :hidden:
   :maxdepth: 2
   :caption: Getting Started

   install
   models

.. toctree::
   :hidden:
   :maxdepth: 2
   :caption: Tutorials

   quick_tutorials
   tutorials
   advanced_tutorials
   youtube
   integration

.. toctree::
   :hidden:
   :maxdepth: 2
   :caption: Reference

   manuscripts
   codedoc
   faq
   changelog
   support
   license

.. toctree::
   :hidden:
   :maxdepth: 2
   :caption: Community Resources

   workshops
   GitHub Repo <https://github.com/kaylai/VESIcal>
   Link Tree <https://linktr.ee/VESIcal>
