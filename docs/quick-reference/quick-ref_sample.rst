=========================
Creating a VESIcal Sample
=========================
VESIcal allows you to perform calculations on a single sample, either defined as a dictionary or pulled from a supplied Excel or CSV file, or on an entire dataset, supplied by an Excel or CSV file. To define a single sample manually, create a dictionary with oxide names as keys and oxide concentrations as values. Note that you are defining a bulk system here, so the H2O and CO2 concentrations refer to those oxides in the entire melt +/- fluid system, not just dissolved in the melt.

Define a single sample manually, like this:

.. code-block:: python

	my_sample = v.Sample({'SiO2':  77.3, 
         'TiO2':   0.08, 
         'Al2O3': 12.6, 
         'Fe2O3':  0.207,
         'Cr2O3':  0.0, 
         'FeO':    0.473, 
         'MnO':    0.0,
         'MgO':    0.03, 
         'NiO':    0.0, 
         'CoO':    0.0,
         'CaO':    0.43, 
         'Na2O':   3.98, 
         'K2O':    4.88, 
         'P2O5':   0.0, 
         'H2O':    6.5,
         'CO2':    0.05})

You can also extract a sample from a provided Excel or CSV file and save it to a variable. In this example, we have imported an Excel file to a variable named `myfile` and wish to extract the sample named 'SampleOne' (also see Import an Excel or CSV File below):

.. code-block:: python

	extracted_sample = myfile.get_sample_composition('SampleOne', asSampleClass=True)