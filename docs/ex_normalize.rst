#################################
Normalizing and Transforming Data
#################################
.. contents::

Before performing model calculations on your data, it may be desired to normalize the input composition to a total of 100 wt%. VESIcal has multiple methods for normalizing sample data using various routines. Each of the normalization routines can be accessed by the user at any time to normalize either a signle sample or all samples in a BatchFile object.

The standard normalization type returns the composition normalized to 100%, including any volatiles. The FixedVolatiles type the oxides to 100%, but volatiles remain fixed while other major element oxides are reduced proporitonally so that the total is 100 wt%. The AdditionalVolatiles type normalizes oxides to 100% assuming the sample is volatile-free. If H\ :subscript:`2`\ O or CO\ :subscript:`2` concentrations are passed to the function, their un-normalized values will be retained in addition to the normalized non-volatile oxides, summing to >100%.

.. code-block:: python

	import VESIcal as v

Normalizing an entire dataset
=============================
Import an Excel file
--------------------

.. code-block:: python

	myfile = v.BatchFile('example_data.xlsx')
	myfile.get_data()

.. csv-table:: Output
   :file: tables/example_data.csv
   :header-rows: 1

Standard Normalization
----------------------
Returns the composition normalized to 100%, including any volatiles.

.. code-block:: python

	standard = myfile.get_data(normalization='standard')
	standard

.. csv-table:: Output
   :file: tables/NormStandard.csv
   :header-rows: 1

The variable standards is a pandas DataFrame object, not a VESIcal.BatchFile object. To apply this normalization to your data, you will need to create a new BatchFile object (or overwrite myfile):

.. code-block:: python

	my_normalized_file = v.BatchFile_from_DataFrame(standard)

or

.. code-block:: python

	my_normalized_file = v.BatchFile(filename=None, dataframe=standard)

FixedVolatiles Normalization
----------------------------
Normalizes the oxides to 100%, but volatiles remain fixed while other major element oxides are reduced proporitonally so that the total is 100 wt%.

.. code-block:: python

	fixed_vols = myfile.get_data(normalization='fixedvolatiles')
	fixed_vols

.. csv-table:: Output
   :file: tables/NormFixedVolatiles.csv
   :header-rows: 1

AdditionalVolatiles Normalization
---------------------------------
Normalizes oxides to 100% assuming the sample is volatile-free. If H_2O or CO_2 concentrations are passed to the function, their un-normalized values will be retained in addition to the normalized non-volatile oxides, summing to >100%.

.. code-block:: python

	additional_vols = myfile.get_data(normalization='additionalvolatiles')
	additional_vols

.. csv-table:: Output
   :file: tables/NormAdditionalVolatiles.csv
   :header-rows: 1

Normalize a single sample composition
=====================================
Extract a single sample from your dataset
-----------------------------------------
Here, a composition is extracted from a BatchFile object and returned as a Sample object. Set asSampleClass=False to return as a dictionary.

.. code-block:: python

	SampleName = 'BT-ex'
	extracted_bulk_comp = myfile.get_sample_composition(SampleName, asSampleClass=True)

The normalization type can be passed to get_sample_composition directly:

.. code-block:: python
	extracted_bulk_comp = myfile.get_sample_composition(SampleName, normalization=<normalization-type>, asSampleClass=True)

Or, normalization can be done to any Sample object, as shown below.

Standard Normalization
----------------------
In the following three examples, the normalized composition is returned as a dictionary, not as a Sample object.

.. code-block:: python

	single_standard = extracted_bulk_comp.get_composition(normalization='standard')
	single_standard

.. code-block:: python

	{'SiO2': 73.3693079617533,
	 'TiO2': 0.07573605983148728,
	 'Al2O3': 11.833759348669886,
	 'Fe2O3': 0.1959670548139733,
	 'Cr2O3': 0.0,
	 'FeO': 0.44778945375366846,
	 'MnO': 0.0,
	 'MgO': 0.028401022436807727,
	 'NiO': 0.0,
	 'CoO': 0.0,
	 'CaO': 0.4070813215942441,
	 'Na2O': 3.7678689766164917,
	 'K2O': 4.619899649720724,
	 'P2O5': 0.0,
	 'H2O': 5.2068541134147495,
	 'CO2': 0.04733503739467954}

FixedVolatiles Normalization
----------------------------
.. code-block:: python

	single_fixed = extracted_bulk_comp.get_composition(normalization='fixedvolatiles')
	single_fixed

.. code-block:: python

	{'SiO2': 73.1402378097522,
	 'TiO2': 0.07549960031974419,
	 'Al2O3': 11.79681254996003,
	 'Fe2O3': 0.19535521582733809,
	 'Cr2O3': 0.0,
	 'FeO': 0.4463913868904875,
	 'MnO': 0.0,
	 'MgO': 0.02831235011990407,
	 'NiO': 0.0,
	 'CoO': 0.0,
	 'CaO': 0.405810351718625,
	 'Na2O': 3.756105115907274,
	 'K2O': 4.6054756195043955,
	 'P2O5': 0.0,
	 'CO2': 0.05,
	 'H2O': 5.5}

AdditionalVolatiles Normalization
---------------------------------
.. code-block:: python

	single_additional = extracted_bulk_comp.get_composition(normalization='additionalvolatiles')
	single_additional

.. code-block:: python

	{'SiO2': 77.4380495603517,
	 'TiO2': 0.07993605115907274,
	 'Al2O3': 12.490007993605113,
	 'Fe2O3': 0.20683453237410068,
	 'Cr2O3': 0.0,
	 'FeO': 0.4726219024780175,
	 'MnO': 0.0,
	 'MgO': 0.029976019184652272,
	 'NiO': 0.0,
	 'CoO': 0.0,
	 'CaO': 0.4296562749800159,
	 'Na2O': 3.9768185451638685,
	 'K2O': 4.8760991207034365,
	 'P2O5': 0.0,
	 'H2O': 5.5,
	 'CO2': 0.05}


