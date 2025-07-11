===============
Quick Reference
===============
.. #TODO make this actually succinct and move longer discussions elsewhere... need happy medium between super quick and "verbose tutorials" which should be called walkthroughs or something.
Need a quick reminder on how to run a particular calculation? You're in the right place. Need a more thorough walkthrough? Check out our :doc:`tutorials`.

.. toctree::
   :maxdepth: 2
   :caption: Working with a Sample

   /quick-reference/quick-ref_sample
   /quick-reference/quick-ref_normalization

.. toctree::
   :maxdepth: 2
   :caption: Working with a BatchFile

   /quick-reference/quick-ref_batchfile

.. toctree::
   :maxdepth: 2
   :caption: Data Transformations

   /quick-reference/quick-ref_converting-units

.. contents::

Importing VESIcal
=================
Start with this:

.. code-block:: python

	import VESIcal as v


General Notes, Tips, and Tricks
===============================


Choosing a model
----------------
MagmaSat (i.e., MELTS v.1.2.0) is the default model for all function calls. But, one of the great powers of VESIcal is the ability to use any of the supplied models for any function call. You can get a list of all available models by typing:

.. code-block:: python

	v.get_model_names()

which returns a list of model names, as strings.

You can then pass any one of those model names to any calculation, both for batch and single-sample calculations, where `<your_sample>` is a variable (not a string). For example:

.. code-block:: python

	v.calculate_saturation_pressure(sample=<your_sample>,
					temperature=<your_temp>,
					model='ShishkinaIdealMixing').result

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

Running VESIcal's Core Calculations
===================================

Calculate Dissolved Volatile Concentrations
-------------------------------------------
For an entire dataset, where `myfile` is an BatchFile object:

.. code-block:: python

	myfile.calculate_dissolved_volatiles(temperature=<your_temp>, 
						pressure=<your_pressure>, 
						X_fluid=<your_X_fluid>)

Or for a single sample, where `<your_sample>` is a variable (not a string):

.. code-block:: python

	v.calculate_dissolved_volatiles(sample=<your_sample>, 
					temperature=<your_temp>, 
					pressure=<your_pressure>, 
					X_fluid=<your_X_fluid>).result

----------

Calculate Equilibrium Fluid Compositions
----------------------------------------
For an entire dataset, where `myfile` is an BatchFile object:

.. code-block:: python

	myfile.calculate_equilibrium_fluid_comp(temperature=<your_temp>, 
						pressure=<your_pressure>)

Or for a single sample, where `<your_sample>` is a variable (not a string):

.. code-block:: python

	v.calculate_equilibrium_fluid_comp(sample=<your_sample>, 
					temperature=<your_temp>, 
					pressure=<your_pressure>).result

----------

Calculate Saturation Pressures
------------------------------
For an entire dataset, where `myfile` is an BatchFile object:

.. code-block:: python

	myfile.calculate_saturation_pressure(temperature=<your_temp>)

Or for a single sample, where `<your_sample>` is a variable (not a string):

.. code-block:: python

	v.calculate_saturation_pressure(sample=<your_sample>, 
					temperature=<your_temp>).result

----------

Calculate and Plot Isobars and Isopleths
----------------------------------------
You can only do this for a single sample. First, calculate the isobars and isopleths like so, where `<your_sample>` is a variable (not a string):

.. code-block:: python

	isobars, isopleths = v.calculate_isobars_and_isopleths(sample=<your_sample>, 
                                            temperature=<your_temp>,
                                            pressure_list=[<pressure1>, <pressure2>, <pressure3>],
                                            isopleth_list=[<isopleth1>, <isopleth2>, <isopleth3>].result

Then, you can very easily plot your newly calculated isobars and isopleths, like so:

.. code-block:: python

	fig, ax = v.plot(isobars=isobars, isopleths=isopleths)
	show()

You may wish to do some custom plotting of your isobar and isopleth data without relying on our built-in plot function. However, the raw isobars and isopleths output by the calculate method are a bit messy. `plot_isobars_and_isopleths()` has curve smoothing built-in. We have also implemented the same smoothing in a separate method, called `smooth_isobars_and_isopleths()` which takes isobars and/or isopleths as inputs and returns a pandas DataFrame with smoothed data ready for plotting. Use that function like so:

.. code-block:: python

	v.vplot.smooth_isobars_and_isopleths(isobars=isobars, isopleths=isopleths)

----------

Calculate and Plot Degassing Paths
----------------------------------
You can only do this for a single sample. First, calculate the degassing path. 

Closed-system
^^^^^^^^^^^^^
This example shows the default degassing path, which is closed system degassing with 0% initial fluid. Here, `<your_sample>` is a variable (not a string)

.. code-block:: python

	degass_closed = v.calculate_degassing_path(sample=<your_sample>,
					temperature=<your_temp>).result

Closed-system with initial fluid
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
You might wish to calculate a degassing path for a closed-system, but where your initial magma already contains some percentage of exsolved fluid. In this case, use the `init_vapor` argument. In this example, we calculate the degassing path with 2% initial fluid, where `<your_sample>` is a variable (not a string):

.. code-block:: python

	degass_init = v.calculate_degassing_path(sample=<your_sample>,
					temperature=<your_temp>,
					init_vapor=2.0).result

Open-system
^^^^^^^^^^^
You may with to calculate an open or partially open system degassing path. This is acheived using the `fractionate_vapor` argument. A value of 1.0 is a completely open system, in which 100% of the fluid is removed at each calculation step. A value of 0.2 would represent a partially open system, in which 20% of the fluid is removed at each calculation step. 

A completely open system, where `<your_sample>` is a variable (not a string):

.. code-block:: python

	degass_open = v.calculate_degassing_path(sample=<your_sample>,
					temperature=<your_temp>,
					fractionate_vapor=1.0).result

A partially open system, where 20% of vapor is fractionated at each calculation step, where `<your_sample>` is a variable (not a string):

.. code-block:: python

	degass_partly_open = v.calculate_degassing_path(sample=<your_sample>,
					temperature=<your_temp>,
					fractionate_vapor=0.2).result

You can then easily plot your newly calculated degassing paths like so:

.. code-block:: python

	fig, ax = v.plot(degassing_paths=[degass_closed, degass_init, degass_open, degass_partly_open],
            		degassing_path_labels=["Closed System", "2% Initial Fluid", "Open System", "Partly Open System"])
    v.show()

------------

Running thermo Calculations
===================================

Calculate Liquid Density
-------------------------------------------
For an entire dataset, where `myfile` is an BatchFile object:

.. code-block:: python

	myfile.calculate_liquid_density(temperature=<your_temp>, 
						pressure=<your_pressure>)

Or for a single sample, where `<your_sample>` is a variable (not a string):

.. code-block:: python

	v.calculate_liquid_density(sample=<your_sample>, 
					temperature=<your_temp>, 
					pressure=<your_pressure>).result

----------

Calculate Liquid Viscosity
-------------------------------------------
For an entire dataset, where `myfile` is an BatchFile object:

.. code-block:: python

	myfile.calculate_liquid_viscosity(temperature=<your_temp>)

Or for a single sample, where `<your_sample>` is a variable (not a string):

.. code-block:: python

	v.calculate_liquid_viscosity(sample=<your_sample>, 
					temperature=<your_temp>).result

----------

Save Your Calculations to an Excel or CSV File
==============================================
Once you have performed some calculations and have assigned their outputs to variables, you can write all of your data to an excel or CSV file or files. Let's assume you have imported a file and written it to a variable called `myfile`. You then performed two calculations: `calculate_dissolved_volatiles()` and `calculate_saturation_pressure()`. You've written those outputs to teh variables `dissolved` and `SatP`, respectively. Here's how you would save these data to an excel file. What gets created is a .xlsx file with the first sheet containing your originally input data, the second sheet containing the dissolved data, and the third sheet containing the SatP data.

.. code-block:: python

	myfile.save_excel("myoutput.xlsx", calculations=[dissolved, SatP])

Optionally, you can tell VESIcal what to name your new sheets in your new excel file:

.. code-block:: python

	myfile.save_excel("myoutput.xlsx", calculations=[dissolved, SatP], sheet_names=["My dissolved data", "My saturation data"])

If instead you wish to save these calculations to CSV files, you can do so as:

.. code-block:: python

	myfile.save_csv(filenames=[my_dissolved_output.csv", "my_SatP_output.csv"], calculations=[dissolved, SatP])

Your calculations will be saved to two CSV files: one for each calculation.

