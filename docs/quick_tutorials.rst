###############
Quick Reference
###############
.. contents::

General Notes, Tips, and Tricks
===============================

Single-sample vs. batch processing
----------------------------------
VESIcal allows you to perform calculations on a single sample, either defined as a dictionary or pulled from a supplied Excel file, or on an entire dataset, supplied by an Excel file. To define a single sample manually, create a dictionary with oxide names as keys and oxide concentrations as values. Note that you are defining a bulk system here, so the H2O and CO2 concentrations refer to those oxides in the entire melt +/- fluid system, not just dissolved in the melt.

Define a single sample manually, like this:

.. code-block:: python

	my_sample = {'SiO2':  77.3, 
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
         'CO2':    0.05}

You can also extract a sample from a provided Excel file and save it to a variable. In this example, we have imported an Excel file to a variable named `myfile` and wish to extract the sample named 'SampleOne' (also see Import an Excel File below):

.. code-block:: python

	extracted_sample = myfile.get_sample_oxide_comp('SampleOne')


Using models other than MagmaSat
--------------------------------
MagmaSat (i.e., MELTS v.1.2.0) is the default model for all function calls. But, one of the great powers of VESIcal is the ability to use any of the supplied models for any function call. You can get a list of all available models by typing:

.. code-block:: python

	get_model_names()

which returns a list of model names, as strings.

You can then pass any one of those model names to any calculation, both for batch and single-sample calculations, where `<your_sample>` is a variable (not a string). For example:

.. code-block:: python

	calculate_saturation_pressure(sample=<your_sample>,
					temperature=<your_temp>,
					model='Shishkina').result

Pull arguments (P, T, X_fluid) from Excel file
----------------------------------------------
For any batch calcultions that take `pressure`, `temperature`, or `X_fluid` arguments, those arguments can either be defined directly in the function call, in which case the one value will be applied to all samples, or the arguments can be passed from the Excel file. For example, let's say we have an Excel file, which we've imported into VESIcal and named `myfile`, which contains compositional data, pressure, and temperature values for all of our samples. Our column with temperature values is named "MyTemps", and our column with pressure values is named "SomePs". We will apply one value for X_fluid to the whole dataset. Note that, even if a column of values for X_fluid exists in our Excel file, the following call will ignore it and instead use the value provided for all samples.

.. code-block:: python

	myfile.calculate_dissolved_volatiles(temperature="MyTemps",
						pressure="SomePs",
						X_fluid=0.35).result

print_status
------------
You can print the progress of any batch calcultion by adding

.. code-block:: python

	print_status=True

as an argument to the function call.

verbose
-------
You can make any single sample calculation return extra computed values by adding

.. code-block:: python

	verbose=True

as an argument to the function call. The values returned depend upon the calculation being performed.

----------

Import an Excel file
====================
You can import an excel file containing compositional data describing your samples using the `ExcelFile` class. Your file should have each sample in a separate row, with data in terms of oxides. You can pass the optional argument `input_type` if oxide concentrations are not in wt% (options are 'wtpercent' (default), 'molpercent', and 'molfrac'). You can pass the optional argument 'label' to define the column title referring to the column containing sample names. The default value is 'Label'.

.. code-block:: python

	ExcelFile('path/to/your/file.xlsx')

You'll want to save this ExcelFile object to a variable. Do that like this:

.. code-block:: python

	myfile = ExcelFile('path/to/your/file.xlsx')

If your excel file has multiple sheets, you can specify which sheet to import. Note that you can only import one sheet at a time.

.. code-block:: python

	myfile = ExcelFile('path/to/your/file.xlsx', sheet_name="SameOfYourSheet")

You can also specify the sheet name by it's number (e.g. the 1st, 2nd, 3rd... sheet in the file) as:

.. code-block:: python

	myfile = ExcelFile('path/to/your/file.xlsx', sheet_name=0) #import the first sheet
	myotherfile = ExcelFile('path/to/your/file.xlsx', sheet_name=4) #import the fifth sheet

----------

Save Your Calculations to an Excel File
=======================================
Once you have performed some calculations and have assigned their outputs to variables, you can write all of your data to an excel file. Let's assume you have imported a file and written it to a variable called `myfile`. You then performed two calculations: `calculate_dissolved_volatiles()` and `calculate_saturation_pressure()`. You've written those outputs to teh variables `dissolved` and `SatP`, respectively. Here's how you would save these data to an excel file. What gets created is a .xlsx file with the first sheet containing your originally input data, the second sheet containing the dissolved data, and the third sheet containing the SatP data.

.. code-block:: python

	myfile.save_excelfile("myoutput.xlsx", calculations=[dissolved, SatP])

Optionally, you can tell VESIcal what to name your new sheets in your new excel file:

.. code-block:: python

	myfile.save_excelfile("myoutput.xlsx", calculations=[dissolved, SatP], sheet_name=["My dissolved data", "My saturation data"])

Normalize and Transform Data
============================

Standard normalization
----------------------
Returns the composition normalized to 100%, including any volatiles.

Normalize an entire dataset:

.. code-block:: python

	normalize(<your_excelfile_object>)

For example:

.. code-block:: python

	normalize(myfile)

Or normalize a single sample,where `<your_sample>` is a variable (not a string):

.. code-block:: python

	normalize(<your_sample>)

FixedVolatiles Normalization
----------------------------
Normalizes the oxides to 100%, but volatiles remain fixed while other major element oxides are reduced proporitonally so that the total is 100 wt%.

.. code-block:: python

	normalize_FixedVolatiles(<your_excelfile_or_sample>)

AdditionalVolatiles Normalization
---------------------------------
Normalizes oxides to 100% assuming the sample is volatile-free. If H2O or CO2  concentrations are passed to the function, their un-normalized values will be retained in addition to the normalized non-volatile oxides, summing to >100%.

.. code-block:: python

	normalize_AdditionalVolatiles(<your_excelfile_or_sample>)

----------

Calculate Dissolved Volatile Concentrations
===========================================
For an entire dataset, where `myfile` is an ExcelFile object:

.. code-block:: python

	myfile.calculate_dissolved_volatiles(temperature=<your_temp>, 
						pressure=<your_pressure>, 
						X_fluid=<your_X_fluid>)

Or for a single sample, where `<your_sample>` is a variable (not a string):

.. code-block:: python

	calculate_dissolved_volatiles(sample=<your_sample>, 
					temperature=<your_temp>, 
					pressure=<your_pressure>, 
					X_fluid=<your_X_fluid>).result

----------

Calculate Equilibrium Fluid Compositions
========================================
For an entire dataset, where `myfile` is an ExcelFile object:

.. code-block:: python

	myfile.calculate_equilibrium_fluid_comp(temperature=<your_temp>, 
						pressure=<your_pressure>)

Or for a single sample, where `<your_sample>` is a variable (not a string):

.. code-block:: python

	calculate_equilibrium_fluid_comp(sample=<your_sample>, 
					temperature=<your_temp>, 
					pressure=<your_pressure>).result

----------

Calculate Saturation Pressures
==============================
For an entire dataset, where `myfile` is an ExcelFile object:

.. code-block:: python

	myfile.calculate_saturation_pressure(temperature=<your_temp>)

Or for a single sample, where `<your_sample>` is a variable (not a string):

.. code-block:: python

	calculate_saturation_pressure(sample=<your_sample>, 
					temperature=<your_temp>).result

----------

Calculate and Plot Isobars and Isopleths
========================================
You can only do this for a single sample. First, calculate the isobars and isopleths like so, where `<your_sample>` is a variable (not a string):

.. code-block:: python

	isobars, isopleths = v.calculate_isobars_and_isopleths(sample=<your_sample>, 
                                            temperature=<your_temp>,
                                            pressure_list=[<pressure1>, <pressure2>, <pressure3>],
                                            isopleth_list=[<isopleth1>, <isopleth2>, <isopleth3>].result

Then, you can very easily plot your newly calculated isobars and isopleths, like so:

.. code-block:: python

	plot(isobars=isobars, isopleths=isopleths)

You may wish to do some custom plotting of your isobar and isopleth data without relying on our built-in plot function. However, the raw isobars and isopleths output by the calculate method are a bit messy. `plot_isobars_and_isopleths()` has curve smoothing built-in. We have also implemented the same smoothing in a separate method, called `smooth_isobars_and_isopleths()` which takes isobars and/or isopleths as inputs and returns a pandas DataFrame with smoothed data ready for plotting. Use that function like so:

.. code-block:: puython

	smooth_isobars_and_isopleths(isobars=isobars, isopleths=isopleths)

----------

Calculate and Plot Degassing Paths
==================================
You can only do this for a single sample. First, calculate the degassing path. 

Closed-system
-------------
This example shows the default degassing path, which is closed system degassing with 0% initial fluid. Here, `<your_sample>` is a variable (not a string)

.. code-block:: python

	degass_closed = calculate_degassing_path(sample=<your_sample>,
					temperature=<your_temp>).result

Closed-system with initial fluid
--------------------------------
You might wish to calculate a degassing path for a closed-system, but where your initial magma already contains some percentage of exsolved fluid. In this case, use the `init_vapor` argument. In this example, we calculate the degassing path with 2% initial fluid, where `<your_sample>` is a variable (not a string):

.. code-block:: python

	degass_init = calculate_degassing_path(sample=<your_sample>,
					temperature=<your_temp>,
					init_vapor=2.0).result

Open-system
-----------
You may with to calculate an open or partially open system degassing path. This is acheived using the `fractionate_vapor` argument. A value of 1.0 is a completely open system, in which 100% of the fluid is removed at each calculation step. A value of 0.2 would represent a partially open system, in which 20% of the fluid is removed at each calculation step. 

A completely open system, where `<your_sample>` is a variable (not a string):

.. code-block:: python

	degass_open = calculate_degassing_path(sample=<your_sample>,
					temperature=<your_temp>,
					fractionate_vapor=1.0).result

A partially open system, where 20% of vapor is fractionated at each calculation step, where `<your_sample>` is a variable (not a string):

.. code-block:: python

	degass_partly_open = calculate_degassing_path(sample=<your_sample>,
					temperature=<your_temp>,
					fractionate_vapor=0.2).result

You can then easily plot your newly calculated degassing paths like so:

.. code-block:: python

	plot(degassing_paths=[degass_closed, degass_init, degass_open, degass_partly_open],
            degassing_path_labels=["Closed System", "2% Initial Fluid", "Open System", "Partly Open System"])



