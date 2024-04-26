*********
ChangeLog
*********

What's new in v 1.2?
####################
A small but mighty update comes to VESIcal v1.2.0. VESIcal can now be used without the installation of the thermoengine library (which can only be built from source on a mac or run via a docker image on a PC or linux machine). MagmaSat, VESIcal's default model, requires thermoengine to run. We chose MagmaSat as the default model for a reason (it's usually the best model to use), but there are plenty of cases where access to the other models and not to MagmaSat is desired.

The usage of VESIcal remains unchanged in the new version. Simply install using pip (see :doc:`install`) and get to work. If you do not have thermoengine installed, remember to specify which model to use in each calculation performed.

Version history
###############
Version 1.2.10
^^^^^^^^^^^^^^
Our attempts to remove DuanDriver error messages thrown by MagmaSat ended up causing other things to break. This update removes DuanDriver "fixes". Expect to see more DuanDriver print messages when running MagmaSat if you got used to them not being there!

Version 1.2.9
^^^^^^^^^^^^^
This fixes the bug outlined in issue https://github.com/kaylai/VESIcal/issues/188, or rather provides a workaround for the bug in thermoengine that was causing it. Also fixed some flake8 complaints, and made some of the unit tests more robust to machine error and random thermoengine errors (see issue https://github.com/kaylai/VESIcal/issues/197). The testing routine passes on the ENKI server.

Version 1.2.8
^^^^^^^^^^^^^
Fixed bug in calculate_degassing_path() where the default final_pressure() value was set to 100 bars instead of 1 bar. Added code to silence unhelpful "Duan Driver" warning messages from thermoengine. Made significant updates to the unittest testing routines.

Version 1.2.7
^^^^^^^^^^^^^
Fixed a bug where a cation concentration was not returned for a Sample composition given in terms of oxides and vice versa.

Version 1.2.6
^^^^^^^^^^^^^
Made a minor adjustment to how pandas dataframes are handled that makes VESIcal work with Pandas 2.0. VESIcal will continue to work with Pandas 1.

Version 1.2.5
^^^^^^^^^^^^^
Fixed a bug where different models would return different data types. Any float is now returned as np.float rather than python's native float. Also removed some (but not all, unfortunately) unnecessary MagmaSat warnings.

Version 1.2.4
^^^^^^^^^^^^^
Fixed bug where a degassing path could not be calculated for some models using a numpy array as a pressure input. Now, a user can pass, for example np.arange(450, 200, -1) to calculate a degassing path from 450 to 200 bars in steps of 1 bar.

Version 1.2.3
^^^^^^^^^^^^^
Critical bug fix. With the update to matplotlib v3.6.0, styles were renamed. Namely, the "seaborn-colorblind" style used in VESIcal was updated to "seaborn-v0_8-colorblind". This change caused VESIcal to not import at all if matplotlib >=3.6.0 was installed. Try/Except block was added to first try "seaborn-colorblind" and then except to use "seaborn-v0_8-colorblind" in order to prevent VESIcal from breaking if older or newer versions of matplotlib are used.

Version 1.2.2
^^^^^^^^^^^^^
Fixed bug in MagmaSat model. To be internally consistent between MagmaSat calculations and between MagmaSat and other models, we now force "fixedvolatiles" normalization on any sample when performing MagmaSat calculations. This makes VESIcal's behavior consistent with the behavior of the MagmaSat app, which may treat input H2O and CO2 contents as either dissolved in the melt or in the system (melt + fluid). This fixes a bug where a user might calculate the saturation pressure, then when calculating the equilibrium fluid composition at that pressure, VESIcal would say the composition is undersaturated. Results in very small differences to the equilibrium_fluid_comp() results, within error of the model (on the order of 0.01 mole fraction difference for H2O or CO2 in the fluid). The saturation pressure results are unchanged as they already performed a "fixedvolatiles" normalization.

Version 1.2.1
^^^^^^^^^^^^^
Fixed small bug in the Iacono-Marziano model. In most cases this will result in either a very small difference to the calculated results or no difference at all.

Version 1.2.0
^^^^^^^^^^^^^
VESIcal can now be installed without the thermoengine library. Note that the MagmaSat model will not work without thermoengine installed (which we recommend in most cases).

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
