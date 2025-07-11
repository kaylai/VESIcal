=============================
Creating and Using BatchFiles
=============================
.. contents::

Import a CSV or Excel file
--------------------------
You can import an Excel or CSV file containing compositional data describing your samples using the `BatchFile` class. Your file should have each sample in a separate row, with data in terms of oxides. You can pass the optional argument `input_type` if oxide concentrations are not in wt% (options are 'wtpercent' (default), 'molpercent', and 'molfrac'). You can pass the optional argument 'label' to define the column title referring to the column containing sample names. The default value is 'Label'.

.. code-block:: python

	v.BatchFile('path/to/your/file.xlsx')

You'll want to save this BatchFile object to a variable. Do that like this:

.. code-block:: python

	myfile = v.BatchFile('path/to/your/file.xlsx')

If your excel file has multiple sheets, you can specify which sheet to import. Note that you can only import one sheet at a time.

.. code-block:: python

	myfile = v.BatchFile('path/to/your/file.xlsx', sheet_name="SameOfYourSheet")

You can also specify the sheet name by it's number (e.g. the 1st, 2nd, 3rd... sheet in the file) as:

.. code-block:: python

	myfile = v.BatchFile('path/to/your/file.xlsx', sheet_name=0) #import the first sheet
	myotherfile = v.BatchFile('path/to/your/file.xlsx', sheet_name=4) #import the fifth sheet

Pull arguments (P, T, X_fluid) from a file
------------------------------------------
If you have pressure and temperature data for each sample in an imported dataset, you can tell VESIcal to use those values in a BatchFile calculation.

..  figure:: /img/temp-press-from-file.png

    Example dataset containing unique pressure and temperature values for each sample in columns titled "Press" and "Temp", respectively.

.. code-block:: python

	myfile.calculate_dissolved_volatiles(temperature="Temp",
						pressure="Press",
						X_fluid=0.35).result

For any batch calcultions that take `pressure`, `temperature`, or `X_fluid` arguments, those arguments can either be defined directly in the function call, in which case the one value will be applied to all samples, or the arguments can be passed from the batch file. For example, let's say we have an Excel file, which we've imported into VESIcal and named `myfile`, which contains compositional data, pressure, and temperature values for all of our samples. Our column with temperature values is named "Temp", and our column with pressure values is named "Press". We will apply one value for X_fluid to the whole dataset. Note that, even if a column of values for X_fluid exists in our Excel file, the following call will ignore it and instead use the value provided for all samples.