##############################################
Calculating and plotting isobars and isopleths
##############################################
.. contents::

Calculate isobars and isopleths
===============================

In this example, we demonstrate how isobars (lines of constant pressure) and isopleths (lines of constant fluid composition) can be calculated for any one composition with :py:meth:`VESIcal.calculate_isobars_and_isopleths()`. A single melt composition can be extracted from a loaded excel file, or a composition can be enterred by hand an stored within a dictionary. Due to computational intensity, isobars and isopleths can only be computed for one sample composition at a time.

Once a single composition is defined, conditions over which to calculate isobars and isopleths must be specified. The generated plot is isothermal, so only one temperature can be chosen. Isobars can be calculated for any number of pressures, passed as a list.

The calculation is performed by iterating through possible concentrations of H2O and CO2 and calculating the equilibrium state for the system. The iteration begins at a fixed H2O concentration, increasing the CO2 concentration in steps of 0.1 wt% until a fluid phase is stable. The H2O concentration is then increased by 0.5 wt% and CO2 is again increased from 0 until a fluid phase is stable. This process is repeated for H2O values ranging from 0–15 wt%. The H2O and CO2 concentrations from each system for which a fluid phase was found to be stable are saved and written to a pandas DataFrame, which is returned upon completion of the calculation.

**Method structure:**

*Only single-sample calculations:*

.. code-block:: python

	def calculate_isobars_and_isopleths(self, sample, temperature, pressure_list, isopleth_list, print_status=False).result

**Required inputs:**

``sample``: The composition of a sample. A single sample may be passed as a dictionary of values, with compositions of oxides in wt%.

``temperature``: The temperature in degres C. 

``pressure_list``: A list of all pressures in bars at which to calculate isobars.

**Optional inputs:**

``isopleth_list``: The default value is None in which case only isobars will be calcualted. A list of all fluid composition values, in mole fraction H2O (XH2Ofluid), at which to calculate isopleths. Values can range from 0–1.

``print_status``: The default value is False. If True is passed, the progress of the calculation will be printed to the terminal. 

**Calculated outputs:**
The function returns two pandas DataFrames: the first has isobar data, and the second has isopleth data. Columns in the isobar dataframe are ‘Pressure’, ‘H2Omelt’, and ‘CO2melt’, correpsonding to pressure in bars and dissolved H2O and CO2 in the liquid in wt%. Columns in the isopleth dataframe are ‘XH2O_fl’, ‘H2O_liq’, and ‘CO2_liq’, corresponding to XH2Ofluid and dissolved H2O and CO2 in the liquid in wt%.

Import a data file and extract a single sample
----------------------------------------------

.. code-block:: python

	myfile = v.BatchFile('../manuscript/example_data.xlsx')
	SampleName = 'BT-ex'
	extracted_bulk_comp = myfile.get_sample_composition(SampleName, asSampleClass=True)

Do the calculation
------------------

.. code-block:: python

	isobars, isopleths = v.calculate_isobars_and_isopleths(sample=extracted_bulk_comp,
                                            	temperature=1000.0,
                                            	pressure_list=[1000.0, 2000.0, 3000.0],
                                            	isopleth_list=[0.25,0.5,0.75],
                                            	print_status=True).result


Plot isobars and isopleths
==========================

Once isobars and isopleths are calculated, they can be plotted in an H2O versus CO2 plot with one simple function call. The plot will be printed directly in the notebook or, if the code is run as script in a command line, the plot will appear it its own window, at which point it can be saved as an image file. The :py:meth:`VESIcal.plot_isobars_and_isopleths()` function takes in a pandas DataFrame with calculated isobar and isopleth information (e.g., output from calculate_isobars_and_isopleths) and plots data as isobars (lines of constant pressure) and isopleths (lines of constant fluid composition). These lines represent the saturation pressures of the melt composition used to calculate the isobar and isopleth information.

**Method structure:** 

*Only single-sample calculations.* 

.. code-block:: python
	
	plot(isobars=None, isopleths=None, degassing_paths=None, custom_H2O=None, custom_CO2=None,
		 isobar_labels=None, isopleth_labels=None, degassing_path_labels=None, custom_labels=None, **kwargs)

**Required inputs for plotting isobars and isopleths:**

``isobars``: Pandas DataFrame object containing isobar information as calculated by calculate_isobars_and_isopleths.

``isopleths``: Pandas DataFrame object containing isopleth information as calculated by calculate_isobars_and_isopleths.

**Calculated outputs:**
The function returns a matplotlib object with the x-axis as H2O wt% in the melt and y-axis as CO2 wt% in the melt. Isobars, or lines of constant pressure at which the sample magma composition is saturated, and isopleths, or lines of constant fluid composition at which the sample magma composition is saturated, are plotted.

Plot your data
--------------

.. code-block:: python

	fig, ax = v.plot(isobars=isobars, isopleths=isopleths)
	v.show()

.. image:: img/ex_isobarsandisopleths_img1.png
   :width: 600

Custom plotting
---------------
The user may wish to apply custom formatting to the plot rather than relying on our built-in plotting method. The output from `calculate_isobars_and_isopleths()` can look a bit messy, which is why our `plot_isobars_and_isopleths()` function incorporates curve smoothing so that plots look pretty. In order to facilitate easy custom plotting, we have implemented a second method, `smooth_isobars_and_isopleths()` which takes in isobars and isopleths calulated with VESIcal, performs smoothing on them using numpy, and then returns a pandas DataFrame with isobars and isopleths ready for plotting by the user.

Following the example above, we could get smoothed isobars and/or isopleths with the following:

.. code-block:: python

	smoothed_isobars, smoothed_isopleths = v.vplot.smooth_isobars_and_isopleths(isobars, isopleths)

The method can also do isobars only or isopleths only, like:

.. code-block:: python

	smoothed_isobars = v.vplot.smooth_isobars_and_isopleths(isobars)
	smoothed_isopleths = v.vplot.smooth_isobars_and_isopleths(isopleths=isopleths)

Below is essentially the same code that our `plot()` and `smooth_isobars_and_isopleths()` methods use to perform smoothing with numpy. Executing the code below will produce the same output as a call to `plot()`. Feel free to grab this code and edit it to customize to your hearts desire.

.. code-block:: python

	import matplotlib.pyplot as plt
	import numpy as np

	# #-----------------------CODE TO PLOT ISOBARS AND ISOPLETHS-----------------#
	P_vals = isobars.Pressure.unique() #grab all unique pressure values
	XH2O_vals = isopleths.XH2O_fl.unique() #grab all unique XH2O values
	isobars_lists = isobars.values.tolist() #transform isobar dataframe into a list of lists
	isopleths_lists = isopleths.values.tolist() #transform isopleth dataframe into a list of lists

	# add zero values to volatiles list
	isobars_lists.append([0.0, 0.0, 0.0, 0.0])

	# draw the figure
	fig, ax1 = plt.subplots()
	plt.xlabel('H2O wt%')
	plt.ylabel('CO2 wt%')

	# do some data smoothing
	for pressure in P_vals:
	    Pxs = [item[1] for item in isobars_lists if item[0] == pressure]
	    Pys = [item[2] for item in isobars_lists if item[0] == pressure]

	    try:
	        np.seterr(divide='ignore', invalid='ignore') #turn off numpy warning
	        ## calcualte polynomial
	        Pz = np.polyfit(Pxs, Pys, 3)
	        Pf = np.poly1d(Pz)

	        ## calculate new x's and y's
	        Px_new = np.linspace(Pxs[0], Pxs[-1], 50)
	        Py_new = Pf(Px_new)

	        # Plot some stuff
	        ax1.plot(Px_new, Py_new)
	    except:
	        ax1.plot(Pxs, Pys)

	for Xfl in XH2O_vals:
	    Xxs = [item[1] for item in isopleths_lists if item[0] == Xfl]
	    Xys = [item[2] for item in isopleths_lists if item[0] == Xfl]

	    try:
	        ## calcualte polynomial
	        Xz = np.polyfit(Xxs, Xys, 2)
	        Xf = np.poly1d(Xz)

	        ## calculate new x's and y's
	        Xx_new = np.linspace(Xxs[0], Xxs[-1], 50)
	        Xy_new = Xf(Xx_new)

	        # Plot some stuff
	        ax1.plot(Xx_new, Xy_new, ls='dashed', color='k')
	    except:
	        ax1.plot(Xxs, Xys, ls='dashed', color='k')

	labels = P_vals
	ax1.legend(labels)

	np.seterr(divide='warn', invalid='warn') #turn numpy warning back on

	plt.show()

