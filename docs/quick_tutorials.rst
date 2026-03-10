###############
Quick Reference
###############

Calculations
============

.. list-table::
   :header-rows: 1
   :widths: 30 35 35

   * - Calculation
     - Function
     - Required args
   * - :ref:`Saturation pressure <calc-satp>`
     - ``calculate_saturation_pressure()``
     - ``sample``, ``temperature``
   * - :ref:`Dissolved volatiles <calc-dissolved>`
     - ``calculate_dissolved_volatiles()``
     - ``sample``, ``temperature``, ``pressure``, ``X_fluid``
   * - :ref:`Equilibrium fluid <calc-eqfluid>`
     - ``calculate_equilibrium_fluid_comp()``
     - ``sample``, ``temperature``, ``pressure``
   * - :ref:`Isobars & isopleths <calc-isobars>`
     - ``calculate_isobars_and_isopleths()``
     - ``sample``, ``temperature``, ``pressure_list``, ``isopleth_list``
   * - :ref:`Degassing path <calc-degassing>`
     - ``calculate_degassing_path()``
     - ``sample``, ``temperature``
   * - :ref:`Liquid density <calc-density>`
     - ``calculate_liquid_density()``
     - ``sample``, ``temperature``, ``pressure``
   * - :ref:`Liquid viscosity <calc-viscosity>`
     - ``calculate_liquid_viscosity()``
     - ``sample``, ``temperature``

.. tip::

   **Single sample:** ``v.calculate_*(...).result`` — append ``.result`` to get values.

   **Batch:** ``myfile.calculate_*(...)`` — called on a ``BatchFile`` object, no ``.result`` needed.

Models
======
Use ``v.get_model_names()`` to print a list of all models.

.. list-table::
    :header-rows: 0

    * - Mixed
      - ``MagmaSat`` (default), ``ShishkinaIdealMixing``, ``Dixon``, ``IaconoMarziano``, ``Liu``
    * - H2O
      - ``ShishkinaWater``,  ``DixonWater``,  ``IaconoMarzianoWater``,  ``MooreWater``, ``LiuWater``, 
    * - CO2
      - ``ShishkinaCarbon``, ``DixonCarbon``, ``IaconoMarzianoCarbon``, ``AllisonCarbon``, ``AllisonCarbon_sunset``, ``AllisonCarbon_sfvf``, ``AllisonCarbon_erebus``, ``AllisonCarbon_vesuvius``, ``AllisonCarbon_etna``, ``AllisonCarbon_stromboli``, ``LiuCarbon``
    
    

Arguments
=========

.. list-table::
   :header-rows: 1
   :widths: 20 15 65

   * - Argument
     - Type
     - Description
   * - ``sample``
     - Sample
     - A ``Sample`` object. Required for single-sample calculations; omit for batch.
   * - ``temperature``
     - float or str
     - Temperature in degrees C. Pass a column name (str) for batch calculations.
   * - ``pressure``
     - float or str
     - Pressure in bars. Pass a column name (str) for batch calculations.
   * - ``X_fluid``
     - float or str
     - Mole fraction of CO2 in the fluid (0 to 1). Pass a column name (str) for batch.
   * - ``model``
     - str
     - Model name. Default: ``'MagmaSat'``. See ``v.get_model_names()``.
   * - ``verbose``
     - bool
     - Return extra values (single-sample only). Default: ``False``.
   * - ``print_status``
     - bool
     - Print progress (batch only). Default: ``False``.

Sample
======

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Task
     - Syntax
   * - Create manually
     - ``v.Sample({'SiO2': 77.3, 'Al2O3': 12.6, ...})``
   * - With normalization
     - ``v.Sample({...}, default_normalization='standard')``
   * - With input units
     - ``v.Sample({...}, units='mol_oxides')``
   * - With output units
     - ``v.Sample({...}, default_units='mol_oxides')``
   * - Extract from file
     - ``myfile.get_sample_composition('SampleName', asSampleClass=True)``
   * - Get composition
     - ``my_sample.get_composition()``
   * - Get normalized
     - ``my_sample.get_composition(normalization='standard')``
   * - Update composition
     - ``my_sample.change_composition(new_comp)``

BatchFile
=========

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Task
     - Syntax
   * - Import CSV or Excel file
     - ``v.BatchFile('file.extension')``
   * - Specific sheet (name)
     - ``v.BatchFile('file.xlsx', sheet_name="Sheet2")``
   * - Specific sheet (index)
     - ``v.BatchFile('file.xlsx', sheet_name=0)``
   * - With normalization
     - ``v.BatchFile('file.xlsx', default_normalization='standard')``
   * - Specify input units (default is wt% oxides)
     - ``v.BatchFile('file.xlsx', units='mol_oxides')``
   * - Specify output units (default is wt% oxides)
     - ``v.BatchFile('file.xlsx', default_units='mol_oxides')``
   * - From DataFrame
     - ``v.BatchFile(filename=None, dataframe=my_df)``
   * - Get all data
     - ``myfile.get_data()``
   * - Get normalized data
     - ``myfile.get_data(normalization='standard')``
   * - Extract one sample (defaults to return as dict, pass ``asSampleClass=True`` to return as Sample() object.)
     - ``myfile.get_sample_composition('Label', asSampleClass=True)``


Normalization
=============

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Value
     - Description
   * - ``"none"`` (default)
     - No normalization is applied.
   * - ``"standard"``
     - Normalizes all oxides (including volatiles) to 100 wt%.
   * - ``"fixedvolatiles"``
     - Normalizes to 100 wt%, but H2O and CO2 stay fixed while other oxides are reduced proportionally.
   * - ``"additionalvolatiles"``
     - Normalizes non-volatile oxides to 100 wt%. Original H2O and CO2 are added on top, so the total may exceed 100%.

:ref:`See normalization examples below <norm-examples>`

Units
=====

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Value
     - Description
   * - ``"wtpt_oxides"`` (default)
     - Weight percent oxides.
   * - ``"mol_oxides"``
     - Mol fraction oxides.
   * - ``"mol_cations"``
     - Mol fraction cations.
   * - ``"mol_singleO"``
     - Mol fraction on a single-oxygen basis.

.. admonition:: **units** vs **default_units**

   ``units`` tells VESIcal what your **input** data are in.

   ``default_units`` tells VESIcal what units to **return** data in.

:ref:`See units examples below <units-examples>`

Saving
======

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Method
     - Description
   * - :ref:`save_excel() <save-excel>`
     - Save to a ``.xlsx`` file with one sheet per calculation.
   * - :ref:`save_csv() <save-csv>`
     - Save to one CSV file per calculation.

----------

Examples
========

.. _setup:

Setup
-----

.. code-block:: python

   import VESIcal as v

.. _create-sample:

Create a Sample
^^^^^^^^^^^^^^^

.. code-block:: python

   my_sample = v.Sample({
       'SiO2': 77.3, 'TiO2': 0.08, 'Al2O3': 12.6,
       'Fe2O3': 0.207, 'Cr2O3': 0.0, 'FeO': 0.473,
       'MnO': 0.0, 'MgO': 0.03, 'NiO': 0.0, 'CoO': 0.0,
       'CaO': 0.43, 'Na2O': 3.98, 'K2O': 4.88,
       'P2O5': 0.0, 'H2O': 6.5, 'CO2': 0.05
   })

.. autoclass:: VESIcal.sample_class.Sample

Extract a sample from an imported file:

.. code-block:: python

   extracted_sample = myfile.get_sample_composition('SampleOne', asSampleClass=True)

.. _import-file:

Import an Excel or CSV File
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   myfile = v.BatchFile('path/to/your/file.xlsx')

   # specific sheet by name
   myfile = v.BatchFile('path/to/your/file.xlsx', sheet_name="NameOfYourSheet")

   # specific sheet by index (0-based)
   myfile = v.BatchFile('path/to/your/file.xlsx', sheet_name=0)

.. autoclass:: VESIcal.batchfile.BatchFile

----------

.. _norm-examples:

Normalization
-------------

On import:

.. code-block:: python

   my_sample = v.Sample({...}, default_normalization='standard')
   myfile = v.BatchFile('file.xlsx', default_normalization='standard')

On existing data:

.. code-block:: python

   normed = mysample.get_composition(normalization="standard")
   mysample.change_composition(normed)

Normalize an entire BatchFile:

.. code-block:: python

   my_normed_data = myfile.get_data(normalization="standard")
   myNewData = v.BatchFile(filename=None, dataframe=my_normed_data)

----------

.. _units-examples:

Units
-----

.. code-block:: python

   # input is mol fraction oxides, output as mol fraction oxides
   my_sample = v.Sample({...}, units='mol_oxides', default_units='mol_oxides')

   # input is wt% (default), convert output to mol fraction oxides
   myfile = v.BatchFile('file.xlsx', default_units='mol_oxides')

----------

.. _tips:

Tips and Tricks
---------------

Pull arguments from a file:

.. code-block:: python

   myfile.calculate_dissolved_volatiles(
       temperature="MyTemps",    # column name in your file
       pressure="SomePs",        # column name in your file
       X_fluid=0.35              # single value applied to all samples
   ).result

Choose a model:

.. code-block:: python

   v.get_model_names()  # list all available models

   v.calculate_saturation_pressure(
       sample=my_sample, temperature=900,
       model='ShishkinaIdealMixing'
   ).result

----------

.. _calc-dissolved:

Dissolved Volatile Concentrations
----------------------------------
``v.calculate_dissolved_volatiles(sample, temperature, pressure, X_fluid)``

.. code-block:: python

   # single sample
   v.calculate_dissolved_volatiles(
       sample=my_sample, temperature=1000, pressure=2000, X_fluid=0.5
   ).result

   # batch
   myfile.calculate_dissolved_volatiles(temperature=1000, pressure=2000, X_fluid=0.5)

----------

.. _calc-eqfluid:

Equilibrium Fluid Compositions
-------------------------------
``v.calculate_equilibrium_fluid_comp(sample, temperature, pressure)``

.. code-block:: python

   # single sample
   v.calculate_equilibrium_fluid_comp(
       sample=my_sample, temperature=1000, pressure=2000
   ).result

   # batch
   myfile.calculate_equilibrium_fluid_comp(temperature=1000, pressure=2000)

----------

.. _calc-satp:

Saturation Pressures
---------------------
``v.calculate_saturation_pressure(sample, temperature)``

.. code-block:: python

   # single sample
   v.calculate_saturation_pressure(sample=my_sample, temperature=1000).result

   # batch
   myfile.calculate_saturation_pressure(temperature=1000)

----------

.. _calc-isobars:

Isobars and Isopleths
-----------------------
``v.calculate_isobars_and_isopleths(sample, temperature, pressure_list, isopleth_list)``

Single-sample only.

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Extra argument
     - Description
   * - ``pressure_list``
     - List of pressures (bars) at which to calculate isobars.
   * - ``isopleth_list``
     - List of X_fluid values at which to calculate isopleths.

.. code-block:: python

   isobars, isopleths = v.calculate_isobars_and_isopleths(
       sample=my_sample,
       temperature=1000,
       pressure_list=[500, 1000, 2000],
       isopleth_list=[0.25, 0.5, 0.75]
   ).result

   fig, ax = v.plot(isobars=isobars, isopleths=isopleths)
   v.show()

For custom plotting, get smoothed data as a DataFrame:

.. code-block:: python

   smoothed = v.vplot.smooth_isobars_and_isopleths(isobars=isobars, isopleths=isopleths)

----------

.. _calc-degassing:

Degassing Paths
----------------
``v.calculate_degassing_path(sample, temperature)``

Single-sample only.

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Extra argument
     - Description
   * - ``init_vapor``
     - Initial percent of exsolved fluid (e.g., ``2.0`` for 2%). Default: ``0.0``.
   * - ``fractionate_vapor``
     - Fraction of fluid removed at each step. ``0.0`` = closed (default), ``1.0`` = fully open, ``0.2`` = 20% removed per step.

.. code-block:: python

   # closed system (default)
   closed = v.calculate_degassing_path(sample=my_sample, temperature=1000).result

   # closed system with 2% initial fluid
   init = v.calculate_degassing_path(sample=my_sample, temperature=1000, init_vapor=2.0).result

   # fully open system
   opened = v.calculate_degassing_path(sample=my_sample, temperature=1000, fractionate_vapor=1.0).result

   # partially open (20% removed per step)
   partial = v.calculate_degassing_path(sample=my_sample, temperature=1000, fractionate_vapor=0.2).result

.. code-block:: python

   fig, ax = v.plot(
       degassing_paths=[closed, init, opened, partial],
       degassing_path_labels=["Closed", "2% Initial Fluid", "Open", "Partly Open"]
   )
   v.show()

----------

.. _calc-density:

Liquid Density
---------------
``v.calculate_liquid_density(sample, temperature, pressure)``

.. code-block:: python

   # single sample
   v.calculate_liquid_density(sample=my_sample, temperature=1000, pressure=2000).result

   # batch
   myfile.calculate_liquid_density(temperature=1000, pressure=2000)

----------

.. _calc-viscosity:

Liquid Viscosity
-----------------
``v.calculate_liquid_viscosity(sample, temperature)``

.. code-block:: python

   # single sample
   v.calculate_liquid_viscosity(sample=my_sample, temperature=1000).result

   # batch
   myfile.calculate_liquid_viscosity(temperature=1000)

----------

.. _save-excel:

Save to Excel
--------------

.. code-block:: python

   dissolved = myfile.calculate_dissolved_volatiles(temperature=900, pressure=1000, X_fluid=0.5)
   SatP = myfile.calculate_saturation_pressure(temperature=900)

   myfile.save_excel("myoutput.xlsx", calculations=[dissolved, SatP])

   # with custom sheet names
   myfile.save_excel("myoutput.xlsx",
       calculations=[dissolved, SatP],
       sheet_names=["Dissolved", "Saturation Pressures"]
   )

.. _save-csv:

Save to CSV
------------

.. code-block:: python

   myfile.save_csv(
       filename=["my_dissolved_output.csv", "my_SatP_output.csv"],
       calculations=[dissolved, SatP]
   )
