import pandas as pd
import numpy as np
import os
import sys
import warnings as w

from core import *

def rename_duplicates(df, suffix='-duplicate-'):
	appendents = (suffix + df.groupby(level=0).cumcount().astype(str).replace('0','')).replace(suffix, '')
	return df.set_index(df.index.astype(str) + appendents)

def status_bar(percent, sample_name, barLen=20):
	"""
	Prints a status bar to the terminal.

	percent: float
		Percent value of progress from 0 to 1

	barLen: int
		Length of bar to print
	""" 
	sys.stdout.write("\r")
	sys.stdout.write("[{:<{}}] {:.0f}%".format("=" * int(barLen * percent), barLen, percent * 100))
	sys.stdout.write("  Working on sample " + str(sample_name))
	if percent == 1.0:
		sys.stdout.write("\n")
	sys.stdout.flush()

# ---------- BATCHFILE CLASS --------- #
class BatchFile(object):
	"""A batch file with sample names and oxide compositions

	Attributes
	----------
		filename: str
			Path to the batch file, e.g., "my_file.xlsx". This always needs to be passed, even if the user is passing a pandas DataFrame
			rather than an batch file. If passing a DataFrame, filename should be set to None. File can be excel file (.xlsx) or .csv.

		sheet_name: str
			OPTIONAL. For Excel files. Default value is 0 which gets the first sheet in the batch spreadsheet file. This implements the pandas.
			read_excel() sheet_name parameter. But functionality to read in more than one sheet at a time (e.g., pandas.read_excel(sheet_name=None))
			is not yet imlpemented in VESIcal. From the pandas 1.0.4 documentation:
				Available cases:
					- Defaults to 0: 1st sheet as a DataFrame
					- 1: 2nd sheet as a DataFrame
					- "Sheet1": Load sheet with name “Sheet1”

		file_type: str
			OPTIONAL. Default is 'excel', which denotes that passed file has extension .xlsx. Other option is 'csv', which denotes that
			the passed file has extension .csv.

		input_type: str
			OPTIONAL. Default is 'wtpercent'. String defining whether the oxide composition is given in wt percent
			("wtpercent", which is the default), mole percent ("molpercent"), or mole fraction ("molfrac").

		label: str
			OPTIONAL. Default is 'Label'. Name of the column within the passed file referring to sample names.

		dataframe: pandas DataFrame
			OPTIONAL. Default is None in which case this argument is ignored. This argument is used when the user wishes to turn
			a pandas DataFrame into an BatchFile object, for example when user data is already in python rather than being imported
			from a file. In this case set `dataframe` equal to the dataframe object being passed in. If using this option, pass
			None to filename.
	"""
	def __init__(self, filename, sheet_name=0, file_type='excel', input_type='wtpercent', label='Label', dataframe=None, **kwargs):
		"""Return a BatchFile object whose parameters are defined here."""
		self.input_type = input_type

		file_name, file_extension = os.path.splitext(filename)
		if file_extension == '.xlsx' or file_extension == '.xls':
			file_type = 'excel'
		if file_extension == '.csv':
			file_type = 'csv'

		if isinstance(sheet_name, str) or isinstance(sheet_name, int):
			pass
		else:
			raise InputError("If sheet_name is passed, it must be of type str or int. Currently, VESIcal cannot import more than one sheet at a time.")

		if dataframe is not None:
			data = dataframe
			data = self.try_set_index(data, label)
		else:
			if file_type == 'excel':
				data = pd.read_excel(filename, sheet_name=sheet_name)
				data = self.try_set_index(data, label)
			elif file_type == 'csv':
				data = pd.read_csv(filename)
				data = self.try_set_index(data, label)
			else:
				raise InputError("file_type must be one of \'excel\' or \'csv\'.")

		data = rename_duplicates(data) #handle any duplicated sample names
		data = data.dropna(how='all') #drop any rows that are all NaNs
		data = data.fillna(0) #fill in any missing data with 0's

		if 'model' in kwargs:
			w.warn("You don't need to pass a model here, so it will be ignored. You can specify a model when performing calculations on your dataset (e.g., calculate_dissolved_volatiles())",RuntimeWarning,stacklevel=2)

		if 'norm' in kwargs:
			w.warn("We noticed you passed a norm argument here. This does nothing. You can normalize your BatchFile and save it to a new variable name after import using normalize(BatchFileObject). See the documentation for more info.",RuntimeWarning,stacklevel=2)

		total_iron_columns = ["FeOt", "FeOT", "FeOtot", "FeOtotal", "FeOstar", "FeO*"]
		for name in total_iron_columns:
			if name in data.columns:
				if 'FeO' in data.columns:
					for row in data.itertuples():
						if data.at[row.Index, "FeO"] == 0 and data.at[row.Index, name] > 0:
							w.warn("Sample " + str(row.Index) + ": " + str(name) + " value of " + str(data.at[row.Index, name]) + " used as FeO. Fe2O3 set to 0.0.",RuntimeWarning,stacklevel=2)
							data.at[row.Index, "Fe2O3"] = 0.0
							data.at[row.Index, "FeO"] = data.at[row.Index, name]
				else:
					w.warn("Total iron column " + str(name) + " detected. This column will be treated as FeO. If Fe2O3 data are not given, Fe2O3 will be 0.0. In future, an option to calcualte FeO/Fe2O3 based on fO2 will be implemented.",RuntimeWarning,stacklevel=2)
					data['FeO'] = data[name]

		for oxide in oxides:
			if oxide in data.columns:
				pass
			else:
				data[oxide] = 0.0

		if input_type == "wtpercent":
			pass
		if input_type == "molpercent":
			data = mol_to_wtpercent(data)
		if input_type == "molfrac":
			data = mol_to_wtpercent(data)

		self.data = data

	def try_set_index(self, dataframe, label):
		"""
		Method to handle setting the index column in an BatchFile object. If no column is passed that matches the default index name,
		then this method will attempt to choose the 'best' column that the user might want to serve as an index column.

		Parameters
		----------
		dataframe: pandas DataFrame

		label: str
			Name of the column within the passed Excel file referring to sample names.
		"""
		_dataframe = dataframe.copy()
		try:
			_dataframe = _dataframe.set_index(label)
		except:
			label_found = False
			for col in _dataframe.columns:
				if col in oxides:
					pass
				else:
					_dataframe = _dataframe.set_index(col)
					label_found = True
					w.warn("No Label column given, so column '" + str(col) + "' was chosen for you. To choose your own, set label='<column-name>'.",RuntimeWarning,stacklevel=2)
					break
			if label_found == False:
				_dataframe.index.name = 'Label'
				w.warn("No Label column given, so one was created for you. To choose your own, set label='<column-name>'.",RuntimeWarning,stacklevel=2)

		return _dataframe

	def preprocess_sample(self,sample):
		"""
		Adds 0.0 values to any oxide data not passed.

		Parameters
		----------
		sample: pandas DataFrame
			self.data composition of samples in wt% oxides

		Returns
		-------
		pandas DataFrame
		"""
		for oxide in oxides:
			if oxide in self.data.columns:
				pass
			else:
				self.data[oxide] = 0.0

		return sample

	def get_sample_oxide_comp(self, samplename, norm='none'):
		"""
		Returns oxide composition of a single sample from a user-imported file as a dictionary

		Parameters
		----------
		samplename: string
			Name of the desired sample

		norm_style: string
			OPTIONAL. Default value is 'standard'. This specifies the style of normalization applied to the sample.

			'standard' normalizes the entire input composition (including any volatiles) to 100%.

			'fixedvolatiles' normalizes oxides to 100%, including volatiles. The volatile
			wt% will remain fixed, whilst the other major element oxides are reduced proportionally
			so that the total is 100 wt%.

			'additionalvolatiles' normalizes oxides to 100%, assuming it is volatile-free. If
			H2O or CO2 are passed to the function, their un-normalized values will be retained
			in addition to the normalized non-volatile oxides, summing to >100%.

			'none' returns the value-for-value un-normalized composition.

		Returns
		-------
		dictionary
			Composition of the sample as oxides
		"""
		if norm == 'none' or norm == 'standard' or norm == 'fixedvolatiles' or norm == 'additionalvolatiles':
			pass
		else:
			raise InputError('norm must be either none, standard, fixedvolatiles, or additionalvolatiles.')

		data = self.data
		my_sample = pd.DataFrame(data.loc[samplename])
		sample_dict = (my_sample.to_dict()[samplename])
		sample_oxides = {}
		for item, value in sample_dict.items():
			if item in oxides:
				sample_oxides.update({item: value})

		if norm == 'standard':
			return normalize(sample_oxides)
		if norm == 'fixedvolatiles':
			return normalize_FixedVolatiles(sample_oxides)
		if norm == 'additionalvolatiles':
			return normalize_AdditionalVolatiles(sample_oxides)
		if norm == 'none':
			return sample_oxides


	def save_excel(self, filename, calculations, sheet_names=None):
		"""
		Saves data calculated by the user in batch processing mode (using the BatchFile class methods) to an organized
		Excel file, with the original user data plus any calculated data.

		Parameters
		----------
		filename: string
			Name of the file. Extension (.xlsx) should be passed along with the name itself, all in quotes (e.g., 'myfile.xlsx').

		calculations: pandas DataFrame or list of pandas DataFrames
			A single DataFrame or list of DataFrames (e.g., calculated outputs from any of the core BatchFile functions: 
			calculate_dissolved_volatiles, calculate_equilibrium_fluid_comp, and calculate_saturation_pressure). If None, 
			only the original user data will be saved.

		sheet_names: None, string, or list
			OPTIONAL. Default value is None. Allows user to set the name of the sheet or sheets written to the Excel file.

		Returns
		-------
			Creates and saves an Excel file with data from each calculation saved to its own sheet.
		"""
		if isinstance(calculations, list):
			if isinstance(sheet_names, list) or sheet_names is None:
				pass
			else:
				raise InputError("If calculations is passed as list, sheet_names must also be list of same length")
		elif calculations == None:
			pass
		else:
			calculations = [calculations]

		with pd.ExcelWriter(filename) as writer:
			self.data.to_excel(writer, 'Original_User_Data')
			if isinstance(calculations, list):
				if sheet_names is None:
					for n, df in enumerate(calculations):
						df.to_excel(writer, 'Calc%s' % n)
				elif isinstance(sheet_names, list):
					pass
				else:
					sheet_names = [sheet_names]
				if isinstance(sheet_names, list):
					if len(sheet_names) == len(calculations):
						pass
					else:
						raise InputError("calculations and sheet_names must have the same length")

					for i in range(len(calculations)):
						if isinstance(sheet_names[i], str):
							calculations[i].to_excel(writer, sheet_names[i])
						else:
							raise InputError("if sheet_names is passed, it must be list of strings")
			elif calculations == None:
				pass
		return print("Saved " + str(filename))

	def save_csv(self, filenames, calculations, **kwargs):
		"""
		Saves data calculated by the user in batch processing mode to a comma-separated values (csv) file. Mirros the pandas.to_csv()
		method. Any argument that can be passed to pandas.csv() can be passed here. One csv file will be saved for each calculation
		passed.

		Parameters
		----------
		filenames: string or list of strings
			Name of the file. Extension (.csv) should be passed along with the name itself, all in quotes (e.g., 'myfile.csv'). The number
			of calculations passed must match the number of filenames passed. If passing more than one, should be passed as a list.

		calculations: pandas DataFrame or list of pandas DataFrames
			A single variable or list of variables containing calculated outputs from any of the core BatchFile functions:
			calculate_dissolved_volatiles, calculate_equilibrium_fluid_comp, and calculate_saturation_pressure.

		Returns
		-------
			Creates and saves a CSV file or files with data from each calculation saved to its own file.
		"""
		if type(filenames) != list:
			filenames = [filenames]
		if type(calculations) != list:
			calculations = [calculations]
		if len(filenames) != len(calculations):
			raise InputError("calculations and filenames must have the same length")

		for i in range(len(filenames)):
			calculations[i].to_csv(filenames[i], **kwargs)
			print ("Saved " + str(filenames[i]))

	