###############################
Calculating saturation presures
###############################
.. contents::

The :py:meth:`VESIcal.calculate_saturation_pressure()` function calculates the minimum pressure at which a given silicate melt with known temperature and H2O and CO2 concentrations would be saturated with fluid. This is calcualted by finding the pressure at which the smallest amount of vapor is present. This function also calculates the composition of the vapor in equilibrium with the melt at those conditions.

The function works by calculating the equilibrium state of the given melt at very high pressure (2,000 MPa) and then decreasing the pressure in steps of 100 MPa until the mass of vapor is >0 grams. At this point, the pressure space is narrowed and searched in steps of 10 MPa and then in steps of 1 MPa until the saturation pressure is found.

**Method structure:**

Single sample:

.. code-block:: python

	def calculate_saturation_pressure(self, sample, temperature, verbose=False).result

ExcelFile batch process:

.. code-block:: python

	def calculate_saturation_pressure(self, temperature, print_status=False)

**Required inputs:**

:py:meth:`sample`: *Only for single-sample calculations*. The composition of a sample. A single sample may be passed as a dictionary of values, with compositions of oxides in wt%.

:py:meth:`temperature`: The temperature in degres C. For ExcelFile batch calculations, if temperature information is present in the ExcelFile (e.g., as a column with unique temperature values for each sample), this can be accessed by passing the column name in quotes to the temperature variable.

**Optional inputs:**

:py:meth:`verbose`: *Only for single-sample calculations.* Default value is False. If set to True, additional parameters are returned in a dictionary: saturation pressure in bars, H2O and CO2 concentrations in the fluid, mass of the fluid in grams, and proportion of the fluid in the system in wt%.

:py:meth:`print_status`: *Only for ExcelFile batch calcualtions*. The default value is False. If True is passed, the progress of the calculation will be printed to the terminal. 

**Calculated outputs:**
If a single sample is passed to sample, the saturation pressure in bars is returned as a numerical value (float) (plus additional variables ‘XH2O_fl’, ‘XCO2_fl’, ‘FluidMass_grams’, and ‘FluidProportion_wtper’ if verbose is set to True).

If mutliple samples are passed as an ExcelFile object, a pandas DataFrame is returned with sample information plus calculated saturation pressures, equilibrium fluid compositions, mass of the fluid in grams, and proportion of the fluid in the system in wt%. Temperature (in degrees C) is always returned.

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

	satPs = myfile.calculate_saturation_pressure(temperature=925.0)
	satPs

.. code-block:: python

	Calculating sample BT-ex
	Calculating sample TVZMa-ex
	Calculating sample TVZOh-ex
	Calculating sample Oh48-FTIR1-MI1-a
	Calculating sample Oh48-FTIR1-MI1-b
	Calculating sample Oh48-FTIR1-MI1-IRc
	Calculating sample Oh50-4.1
	Calculating sample Oh50-4.2
	Calculating sample Oh49-4.1
	Calculating sample Oh49-4.2
	Calculating sample Ma55-5a.1
	Calculating sample Ma57-3b.2
	Calculating sample Ma57-3c.1
	Calculating sample Ma57-3c.2
	Done!

.. csv-table:: Output
   :file: tables/satP.csv
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

	v.calculate_saturation_pressure(sample=extracted_bulk_comp, temperature=925.0).result

.. code-block:: python

	2310.0
