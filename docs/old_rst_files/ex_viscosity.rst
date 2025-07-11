##############################
Calculating liquid viscosities
##############################
.. contents::

The :py:meth:`VESIcal.calculate_liquid_viscosity()` function calculates the viscosity of the silicate liquid given composition and temperature. The function uses the model of Giordano et al (2008). No other models are currently available for this calculation.

**Method Structure**

Single sample:

.. code-block:: python

	def calculate_liquid_viscosity(self, sample, temperature).result

BatchFile process:

.. code-block:: python

	def calculate_liquid_viscosity(self, temperature)

**Required inputs:**

``sample``: *Only for single-sample calculations*. The composition of a sample as Sample class.

``temperature``: The temperature in degres C. For BatchFile calculations, if temperature information is present in the file (e.g., as a column with unique temperature values for each sample), this can be accessed by passing the column name in quotes to the temperature variable.

**Calculated outputs:**
The Log viscosity of the liquid in Pa*s, rounded to 4 dp.

For an entire dataset
=====================
Import a data file
------------------

.. code-block:: python

	myfile = v.BatchFile('example_data.xlsx')
	myfile.get_data()

.. csv-table:: Output
   :file: tables/example_data.csv
   :header-rows: 1

Do the calculation
------------------

.. code-block:: python

	viscosities = myfile.calculate_liquid_viscosity(temperature=900)
	viscosities

.. csv-table:: Output
   :file: tables/viscosities.csv
   :header-rows: 1

For a single sample
===================

Extract a single sample from your dataset
-----------------------------------------

.. code-block:: python

	SampleName = 'BT-ex'
	extracted_bulk_comp = myfile.get_sample_composition(SampleName, asSampleClass=True)

Do the calculation
------------------

.. code-block:: python

	v.calculate_liquid_viscosity(sample=extracted_bulk_comp, temperature=900).result

.. code-block:: python

	3.9209
