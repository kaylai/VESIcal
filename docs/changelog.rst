*******************************
ChangeLog: What's new in v 1.0?
*******************************

Some minor but important changes to the VESIcal code were made between version 0.9 and 1.0. These are mainly bug fixes and editing of the code itself to be fully flake8 compliant. We have also implemented much more rigorous testing routines, important if you are interested in forking and developing for VESIcal. No functionality changes were implemented during this upgrade, so if your code worked in version 0.9, it should continue to work in 1.0. If you are upgrading code originally written for a version <0.9, please see below for important updates.


Upgrading your code to v 0.9+
#############################
In early 2021 the VESIcal code went through some major structural changes. From the user's perspective, not a lot of the funcationality has changed, except for loads of new features and some key changes to fundamental function calls. Those are detailed here.

In a couple of instances, code written with versions <0.9 will not execute properly in 0.9+. Major changes are few but key:

	- ExcelFile() has changed to BatchFile()
	- ExcelFile.data has changed to BatchFile.get_data()
	- Samples must now be in the form of VESIcal's Sample class before being passed to a calculation
	- To extract a sample from a BatchFile is now done as BatchFile.get_sample_composition() instead of ExcelFile.get_sample_oxide_comp()

Try this first: Quick and dirty update
**************************************
As a first pass, simply try making these changes:

	1. Change all instances of ExcelFile to BatchFile
	2. Change all instances of .data to .get_data()
	3. Change all instances of get_sample_oxide_comp() to get_sample_composition(). Note that if you wish to then pass the extracted sample to a calculation, it will need to be in the form of a Sample class. Do so with get_sample_composition(<your-sample-name>, asSampleClass=True)
	4. Construct a sample as a Sample class from scratch like so:

.. code-block:: python

	mysample = v.Sample({'SiO2': 77.5,
		 'TiO2': 0.08,
		 'Al2O3': 12.5,
		 'Fe2O3': 0.207,
		 'Cr2O3': 0.0,
		 'FeO': 0.473,
		 'MnO': 0.0,
		 'MgO': 0.03,
		 'NiO': 0.0,
		 'CoO': 0.0,
		 'CaO': 0.43,
		 'Na2O': 3.98,
		 'K2O': 4.88,
		 'P2O5': 0.0,
		 'H2O': 5.5,
		 'CO2': 0.05})

If your code continues to throw errors, please refer to the guides in this documentation, which have been updated to reflect changes made for version 0.9. If all else fails, give us a shout: kayla.iacovino@nasa.gov.