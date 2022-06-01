*******************************
ChangeLog: What's new in v 1.2?
*******************************
A small but mighty update comes to VESIcal v1.2.0. VESIcal can now be used without the installation of the thermoengine library (which can only be built from source on a mac or run via a docker image on a PC or linux machine). MagmaSat, VESIcal's default model, requires thermoengine to run. We chose MagmaSat as the default model for a reason (it's usually the best model to use), but there are plenty of cases where access to the other models and not to MagmaSat is desired.

The usage of VESIcal remains unchanged in the new version. Simply install using pip (see :doc:`install`) and get to work. If you do not have thermoengine installed, remember to specify which model to use in each calculation performed.

Previous version history
########################

Version 1.1.1
^^^^^^^^^^^^^
Fixed bug where users could not calculate or plot degassing paths for compositions with only one volatile. Fixed bug where degassing paths in MagmaSat were not calculated down to 1 bar.

Version 1.1.0
^^^^^^^^^^^^^
New in version 1.1.0 is VESIcal's thermo package. This package begins to introduce more calculations that a user might wish to perform. In this version, we have added the ability to calculate the density of a liquid (using DensityX, Iacovino and Till, 2019) and the viscosity of a liquid (Giordano et al., 2008). Check out the tutorials to learn how to use these functions, which use the same syntax as all other core calculations.

Version 1.0.4
^^^^^^^^^^^^^
Diverted Duan H2O driver warnings to the void.

Version 1.0.3
^^^^^^^^^^^^^
Batch calculations for :py:meth:`calculate_equilibrium_fluid_comp()` with the MagmaSat (default) model can now take the ``verbose`` argument. If passed as ``verbose=True``, the calculation will return the mass and weight proportion of fluid in equilibrium with the system. Previously this only worked for single sample calculations.


Version 1.0
###########
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
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
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