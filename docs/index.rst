##################
VESIcal Quickstart
##################

.. image:: img/header_transparent.png

News
----
VESIcal can now be used without needing to install thermoengine!
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Without thermoengine installed, MagmaSat (the default model) cannot be used, however all other models, all plotting capability, and the thermo package can be used. Simply make sure you are using version 1.2.0 or higher following the simple installation instructions below. Then be sure to explicitly tell VESIcal which model you want to use (any model other than MagmaSat. See :doc:`models` for a complete list) by passing `model="some-model-name"` when performing a calculation.

For example:

.. code-block:: python

   v.calculate_saturation_pressure(sample=<your-sample-here>, model="<some-model-name-here>")

If thermoengine is not installed and you import VESIcal, you will be warned that you won't be able to use MagmaSat, but everything else will work as expected. Remember that if you do not pass a model name, VESIcal will default to MagmaSat, and an error will be generated telling you that you need to pass a model name.

Installation
------------
VESIcal can be installed with pip:

.. code-block:: bash

    pip install VESIcal

Always use the most up-to-date version of the code:

.. code-block:: bash

   pip install VESIcal --upgrade

See the :doc:`install` section for detailed instructions and dependencies.

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


Indices and tables
^^^^^^^^^^^^^^^^^^

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
