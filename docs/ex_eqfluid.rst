##########################################
Calculating equilibrium fluid compositions
##########################################
.. contents::

The :py:meth:`VESIcal.calculate_equilibrium_fluid_comp()` function calculates the composition of a fluid phase in equilibrium with a given silicate melt with known pressure, temperature, and dissolved H2O and CO2 concentrations. The calculation is performed simply by calculating the equilibrium state of the given sample at the given conditions and determining if that melt is fluid saturated. If the melt is saturated, fluid composition and mass are reported back. If the calculation finds that the melt is not saturated at the given pressure and temperature, values of 0.0 will be returned for the H2O and CO2 concentrations in the fluid.

**Method structure:**

Single sample:

.. code-block:: python

	def calculate_equilibrium_fluid_comp(self, sample, temperature, pressure, verbose=False).result

ExcelFile batch process:

.. code-block:: python

	def calculate_equilibrium_fluid_comp(self, temperature, pressure, print_status=False, model='MagmaSat')

**Required inputs:**

:py:meth:`sample`: *Only for single-sample calculations*. The composition of a sample. A single sample may be passed as a dictionary of values, with compositions of oxides in wt%.

:py:meth:`temperature`, :py:meth:`pressure`: the temperature in degrees C and the pressure in bars. Temperature and pressure of the sample or samples must be passed unless an ExcelFile object with a column for temperature and/or pressure is passed to sample. If, alternatively, the user wishes to use temperature or pressure information in their ExcelFile object, the title of the column containing temperature or pressure data should be passed in quotes (as a string) to temperature or pressure respectively. Note for batch calculations that if temperature or pressure information exists in the ExcelFile but a single numerical value is defined for one or both of these variables, both the original information plus the values used for the calculations will be returned.

**Optional inputs:**

:py:meth:`verbose`: *Only for single-sample calculations.* Default value is False. If set to True, additional parameters are returned in a dictionary: H2O and CO2 concentrations in the fluid, mass of the fluid in grams, and proportion of the fluid in the system in wt%.

:py:meth:`print_status`: *Only for ExcelFile batch calcualtions*. The default value is False. If True is passed, the progress of the calculation will be printed to the terminal. 

**Calculated outputs:**
If a single sample is passed to sample, a dictionary with keys ‘H2O’ and ‘CO2’ is returned (plus additional variables ‘FluidMass_grams’ and ‘FluidProportion_wtper’ if verbose is set to True).

If mutliple samples are passed as an ExcelFile object, a pandas DataFrame is returned with sample information plus calculated equilibrium fluid compositions, mass of the fluid in grams, and proportion of the fluid in the system in wt%. Pressure (in bars) and Temperature (in degrees C) columns are always returned.

For an entire dataset
=====================
Import an Excel file
--------------------

.. code-block:: python

	myfile = v.ExcelFile('example_data.xlsx')
	myfile.data

.. csv-table:: Output
   :file: tables/example_data.csv
   :header-rows: 1

Do the calculation
------------------

.. code-block:: python

	eqfluid = myfile.calculate_equilibrium_fluid_comp(temperature=900.0, pressure=200.0)
	eqfluid

.. csv-table:: Output
   :file: tables/eqfluid.csv
   :header-rows: 1

For a single sample
===================

Extract a single sample from your dataset
-----------------------------------------

.. code-block:: python

	SampleName = 'BT-ex'
	extracted_bulk_comp = myfile.get_sample_oxide_comp(SampleName)

Do the calculation
------------------

.. code-block:: python

	v.calculate_equilibrium_fluid_comp(sample=extracted_bulk_comp, temperature=900.0, pressure=200.0).result

.. code-block:: python

	{'H2O': 0.994834052532463, 'CO2': 0.00516594746753703}






