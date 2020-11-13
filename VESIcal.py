"""
VESIcal

A generalized python library for calculating and plotting various things related to mixed volatile (H2O-CO2) solubility in silicate melts.
"""

__version__ = "0.1.4"
__author__ = 'Kayla Iacovino, Simon Matthews, and Penny Wieser'

#--------------TURN OFF WARNINGS-------------#
import warnings as w
w.filterwarnings("ignore", message="rubicon.objc.ctypes_patch has only been tested ")
w.filterwarnings("ignore", message="The handle")

#-----------------IMPORTS-----------------#
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
from cycler import cycler
from abc import ABC, abstractmethod
from scipy.optimize import root_scalar
from scipy.optimize import root
from scipy.optimize import minimize
import sys
import sympy
from copy import copy
# import anvil_server

#--------------MELTS preamble---------------#
from thermoengine import equilibrate
# instantiate thermoengine equilibrate MELTS instance
melts = equilibrate.MELTSmodel('1.2.0')

# Suppress phases not required in the melts simulation
phases = melts.get_phase_names()
for phase in phases:
	melts.set_phase_inclusion_status({phase: False})
melts.set_phase_inclusion_status({'Fluid': True, 'Liquid': True})


#----------DEFINE SOME CONSTANTS-------------#
oxides = ['SiO2', 'TiO2', 'Al2O3', 'Fe2O3', 'Cr2O3', 'FeO', 'MnO', 'MgO', 'NiO', 'CoO', 'CaO', 'Na2O', 'K2O', 'P2O5',
		  'H2O', 'CO2']
anhydrous_oxides = ['SiO2', 'TiO2', 'Al2O3', 'Fe2O3', 'Cr2O3', 'FeO', 'MnO', 'MgO', 'NiO', 'CoO', 'CaO', 'Na2O', 'K2O', 'P2O5']
volatiles = ['H2O', 'CO2']
oxideMass = {'SiO2': 28.085+32, 'MgO': 24.305+16, 'FeO': 55.845+16, 'CaO': 40.078+16, 'Al2O3': 2*26.982+16*3, 'Na2O': 22.99*2+16,
			 'K2O': 39.098*2+16, 'MnO': 54.938+16, 'TiO2': 47.867+32, 'P2O5': 2*30.974+5*16, 'Cr2O3': 51.996*2+3*16,
			 'NiO': 58.693+16, 'CoO': 28.01+16, 'Fe2O3': 55.845*2+16*3,
			 'H2O': 18.02, 'CO2': 44.01}
CationNum = {'SiO2': 1, 'MgO': 1, 'FeO': 1, 'CaO': 1, 'Al2O3': 2, 'Na2O': 2,
			 'K2O': 2, 'MnO': 1, 'TiO2': 1, 'P2O5': 2, 'Cr2O3': 2,
			 'NiO': 1, 'CoO': 1, 'Fe2O3': 2, 'H2O': 2, 'CO2': 1}
OxygenNum = {'SiO2': 2, 'MgO': 1, 'FeO': 1, 'CaO': 1, 'Al2O3': 3, 'Na2O': 1,
			 'K2O': 1, 'MnO': 1, 'TiO2': 2, 'P2O5': 5, 'Cr2O3': 3,
			 'NiO': 1, 'CoO': 1, 'Fe2O3': 3, 'H2O': 1, 'CO2': 2}
CationCharge = {'SiO2': 4, 'MgO': 2, 'FeO': 2, 'CaO': 2, 'Al2O3': 3, 'Na2O': 1,
			 'K2O': 1, 'MnO': 2, 'TiO2': 4, 'P2O5': 5, 'Cr2O3': 3,
			 'NiO': 2, 'CoO': 2, 'Fe2O3': 3, 'H2O': 1, 'CO2': 4}
CationMass = {'SiO2': 28.085, 'MgO': 24.305, 'FeO': 55.845, 'CaO': 40.078, 'Al2O3': 26.982, 'Na2O': 22.990,
			 'K2O': 39.098, 'MnO': 54.938, 'TiO2': 47.867, 'P2O5': 30.974, 'Cr2O3': 51.996,
			 'NiO': 58.693, 'CoO': 28.01, 'Fe2O3': 55.845, 'H2O': 2, 'CO2': 12.01}
oxides_to_cations = {'SiO2': 'Si', 'MgO': 'Mg', 'FeO': 'Fe', 'CaO': 'Ca', 'Al2O3': 'Al', 'Na2O': 'Na',
			 'K2O': 'K', 'MnO': 'Mn', 'TiO2': 'Ti', 'P2O5': 'P', 'Cr2O3': 'Cr',
			 'NiO': 'Ni', 'CoO': 'Co', 'Fe2O3': 'Fe3', 'H2O': 'H', 'CO2': 'C', 'F': 'F'}
cations_to_oxides = {'Si': 'SiO2', 'Mg': 'MgO', 'Fe': 'FeO', 'Ca': 'CaO', 'Al': 'Al2O3', 'Na': 'Na2O',
			 'K': 'K2O', 'Mn': 'MnO', 'Ti': 'TiO2', 'P': 'P2O5', 'Cr': 'Cr2O3',
			 'Ni': 'NiO', 'Co': 'CoO', 'Fe3': 'Fe2O3', 'H': 'H2O', 'C': 'CO2'}


#----------DEFINE SOME EXCEPTIONS--------------#

class Error(Exception):
	"""Base class for exceptions in this module."""
	pass


class InputError(Error):
	"""Exception raised for errors in the input.

	Attributes:
		expression -- input expression in which the error occurred
		message -- explanation of the error
	"""

	def __init__(self, message):
		self.message = message

class SaturationError(Error):
	"""Exception raised for errors thrown when a sample does not reach saturation.

	Attributes:
		expression -- input expression in which the error occurred
		message -- explanation of the error
	"""

	def __init__(self, message):
		self.message = message


#----------DEFINE CUSTOM PLOTTING FORMATTING------------#
style = "seaborn-colorblind"
plt.style.use(style)
plt.rcParams["mathtext.default"] = "regular"
plt.rcParams["mathtext.fontset"] = "dejavusans"
mpl.rcParams['patch.linewidth'] = 1
mpl.rcParams['axes.linewidth'] = 1
plt.rcParams['axes.titlesize'] = 20
plt.rcParams['axes.labelsize'] = 18
plt.rcParams['xtick.labelsize'] = 14
plt.rcParams['ytick.labelsize'] = 14
plt.rcParams['legend.fontsize'] = 14
mpl.rcParams['lines.markersize'] = 10

#Define color cycler based on plot style set here
the_rc = plt.style.library[style] #get style formatting set by plt.style.use()
color_list = the_rc['axes.prop_cycle'].by_key()['color'] * 10 #list of colors by hex code
color_cyler = the_rc['axes.prop_cycle'] #get the cycler

def printTable(myDict):
	""" Pretty print a dictionary (as pandas DataFrame)

	Parameters
	----------
	myDict: dict
		A dictionary

	Returns
	-------
	pandas DataFrame
		The input dictionary converted to a pandas DataFrame
	"""
	_myDict = myDict.copy()
	try:
		oxidesum = sum(_myDict[oxide] for oxide in oxides)
		_myDict.update({"Sum oxides": oxidesum})
	except:
		pass
	table = pd.DataFrame([v for v in _myDict.values()], columns = ['value'],
						 index = [k for k in _myDict.keys()])
	return table

#----------DEFINE SOME UNIVERSAL INFORMATIVE METHODS--------------#
def get_model_names():
	"""
	Returns all available model names as a list of strings.
	"""
	model_names = []
	for key, value in default_models.items():
		model_names.append(key)

	return model_names

#----------DEFINE SOME BASIC DATA TRANSFORMATION METHODS-----------#

def mol_to_wtpercent(sample):
	"""
	Takes in a pandas DataFrame containing multi-sample input or a dictionary containing single-sample input
	and returns a pandas DataFrame object with oxide values converted from mole percent to wt percent.

	Parameters
	----------
	oxides: pandas DataFrame object or dictionary
		Variable name referring to the pandas DataFrame object that contains user-imported data or a dictionary
		for single-sample input.
	"""
	data = sample

	if isinstance(sample, pd.DataFrame):
		for key, value in oxideMass.items():
			data.loc[:, key] *= value

		data["MPOSum"] = sum([data[oxide] for oxide in oxides])

		for oxide in oxides:
			data.loc[:, oxide] /= data['MPOSum']
			data.loc[:, oxide] *= 100
		del data['MPOSum']

	elif isinstance(sample, dict):
		for oxide in oxides:
			if oxide in data.keys():
				pass
			else:
				data[oxide] = 0.0

		data = {oxide:  data[oxide] for oxide in oxides}

		for key, value in oxideMass.items():
			data.update({key: (data[key] * value)})
		MPOSum = sum(data.values())
		for key, value in data.items():
			data.update({key: 100 * value / MPOSum})

	return data

def wtpercentOxides_to_molCations(oxides):
	"""Takes in a pandas Series containing major element oxides in wt%, and converts it
	to molar proportions of cations (normalised to 1).

	Parameters
	----------
	oxides         dict or pandas Series
		Major element oxides in wt%.

	Returns
	-------
	dict or pandas Series
		Molar proportions of cations, normalised to 1.
	"""
	molCations = {}
	_oxides = oxides.copy()
	if type(oxides) == dict:
		oxideslist = list(_oxides.keys())
	elif type(oxides) == pd.core.series.Series:
		oxideslist = list(_oxides.index)
	else:
		raise InputError("The composition input must be a pandas Series or dictionary.")
	for ox in oxideslist:
		cation = oxides_to_cations[ox]
		molCations[cation] = CationNum[ox]*_oxides[ox]/oxideMass[ox]

	if type(oxides) == pd.core.series.Series:
		molCations = pd.Series(molCations)
		molCations = molCations/molCations.sum()
	else:
		total = np.sum(list(molCations.values()))
		for ox in oxideslist:
			cation = oxides_to_cations[ox]
			molCations[cation] = molCations[cation]/total

	return molCations

def wtpercentOxides_to_molOxides(oxides):
	""" Takes in a pandas Series or dict containing major element oxides in wt%, and converts it
	to molar proportions (normalised to 1).

	Parameters
	----------
	oxides         dict or pandas Series
		Major element oxides in wt%

	Returns
	-------
	dict or pandas Series
		Molar proportions of major element oxides, normalised to 1.
	"""
	molOxides = {}
	_oxides = oxides.copy()
	if type(oxides) == dict or type(oxides) == pd.core.series.Series:
		if type(oxides) == dict:
			oxideslist = list(oxides.keys())
		elif type(oxides) == pd.core.series.Series:
			oxideslist = list(oxides.index)

		for ox in oxideslist:
			molOxides[ox] = _oxides[ox]/oxideMass[ox]

		if type(oxides) == pd.core.series.Series:
			molOxides = pd.Series(molOxides)
			molOxides = molOxides/molOxides.sum()
		else:
			total = np.sum(list(molOxides.values()))
			for ox in oxideslist:
				molOxides[ox] = molOxides[ox]/total

		return molOxides

	elif isinstance(sample, pd.DataFrame):
		data = sample
		for key, value in oxideMass.items():
			data.loc[:, key] /= value

		data["MPOSum"] = sum([data[oxide] for oxide in oxides])

		for oxide in oxides:
			data.loc[:, oxide] /= data['MPOSum']
		del data['MPOSum']

		return data

	else:
		raise InputError("The composition input must be a pandas Series or dictionary.")

def wtpercentOxides_to_molSingleO(oxides,exclude_volatiles=False):
	""" Takes in a pandas Series containing major element oxides in wt%, and constructs
	the chemical formula, on a single oxygen basis.

	Parameters
	----------
	oxides         dict or pandas Series
		Major element oxides in wt%

	Returns
	-------
	dict or pandas Series
		The chemical formula of the composition, on a single oxygen basis. Each element is
		a separate entry in the Series.
	"""
	molCations = {}
	_oxides = oxides.copy()
	if type(oxides) == dict:
		oxideslist = list(oxides.keys())
	elif type(oxides) == pd.core.series.Series:
		oxideslist = list(oxides.index)
	else:
		raise InputError("The composition input must be a pandas Series or dictionary.")

	total_O = 0.0
	for ox in oxideslist:
		if exclude_volatiles == False or (ox != 'H2O' and ox != 'CO2'):
			cation = oxides_to_cations[ox]
			molCations[cation] = CationNum[ox]*oxides[ox]/oxideMass[ox]
			total_O += OxygenNum[ox]*oxides[ox]/oxideMass[ox]

	if type(oxides) == pd.core.series.Series:
		molCations = pd.Series(molCations)
		molCations = molCations/total_O
	else:
		# total = np.sum(list(molCations.values()))
		for ox in oxideslist:
			if exclude_volatiles == False or (ox != 'H2O' and ox != 'CO2'):
				cation = oxides_to_cations[ox]
				molCations[cation] = molCations[cation]/total_O

	return molCations

def wtpercentOxides_to_formulaWeight(sample,exclude_volatiles=False):
	""" Converts major element oxides in wt% to the formula weight (on a 1 oxygen basis).
	Parameters
	----------
	sample     dict or pandas Series
		Major element oxides in wt%.
	exclude_volatiles 	bool
		If True H2O and CO2 will be excluded from the formula weight calculation.

	Returns
	-------
	float
		The formula weight of the composition, on a one oxygen basis.
	"""
	if type(sample) == dict:
		_sample = pd.Series(sample.copy())
	elif type(sample) != pd.core.series.Series:
		raise InputError("The composition input must be a pandas Series or dictionary.")
	else:
		_sample = sample.copy()

	cations = wtpercentOxides_to_molSingleO(_sample,exclude_volatiles=exclude_volatiles)
	if type(cations) != dict:
		cations = dict(cations)

	# if exclude_volatiles == True:
	# 	if 'C' in cations:
	# 		cations.pop('C')
	# 	if 'H' in cations:
	# 		cations.pop('H')
	# 	newsum = 0
	# 	for cation in cations:
	# 		newsum += OxygenNum[cations_to_oxides[cation]]
	# 	for cation in cations:
	# 		cations[cation] = cations[cation]/newsum

	FW = 15.999
	for cation in list(cations.keys()):
		FW += cations[cation]*CationMass[cations_to_oxides[cation]]
	return FW

#----------DATA TRANSFORMATION FOR PANDAS DATAFRAMES---------#
def fluid_molfrac_to_wt(data, H2O_colname='XH2O_fl_VESIcal', CO2_colname='XCO2_fl_VESIcal'):
	"""
	Takes in a pandas dataframe object and converts only the fluid composition from mole fraction to wt%, leaving the melt composition
	in tact. The user must specify the names of the XH2O_fl and XCO2_fl columns.

	Parameters
	----------
	data: pandas DataFrame
		Sample composition(s) containing columns for H2O and CO2 concentrations in the fluid.

	H2O_colname: str
		OPTIONAL. The default value is 'XH2O_fl', which is what is returned by ExcelFile() core calculations.
		String containing the name of the column corresponding to the H2O concentration in the fluid, in mol fraction.

	CO2_colname: str
		OPTIONAL. The default value is 'XCO2_fl', which is what is returned by ExcelFile() core calculations.
		String containing the name of the column corresponding to the CO2 concentration in the fluid, in mol fraction.

	Returns
	-------
	pandas DataFrame
		Original data passed plus newly calculated values are returned.
	"""
	convData = data.copy()

	MPO_H2O_list = []
	MPO_CO2_list = []
	for index, row in convData.iterrows():
		MPO_H2O_list.append(row[H2O_colname] * oxideMass["H2O"])
		MPO_CO2_list.append(row[CO2_colname] * oxideMass["CO2"])

	convData["MPO_H2O"] = MPO_H2O_list
	convData["MPO_CO2"] = MPO_CO2_list
	convData["H2O_fl_wt"] = 100 * convData["MPO_H2O"] / (convData["MPO_H2O"] + convData["MPO_CO2"])
	convData["CO2_fl_wt"] = 100 * convData["MPO_CO2"] / (convData["MPO_H2O"] + convData["MPO_CO2"])

	del convData["MPO_H2O"]
	del convData["MPO_CO2"]

	return convData

def fluid_wt_to_molfrac(data, H2O_colname='H2O_fl_wt', CO2_colname='CO2_fl_wt'):
	"""
	Takes in a pandas dataframe object and converts only the fluid composition from wt% to mole fraction, leaving the melt composition
	in tact. The user must specify the names of the H2O_fl_wt and CO2_fl_wt columns.

	Parameters
	----------
	data: pandas DataFrame
		DataFrame containing columns for H2O and CO2 concentrations in the fluid.

	H2O_colname: str
		OPTIONAL. The default value is 'H2O_fl_wt', which is what is returned by ExcelFile() core calculations.
		String containing the name of the column corresponding to the H2O concentration in the fluid, in wt%.

	CO2_colname: str
		OPTIONAL. The default value is 'CO2_fl_wt', which is what is returned by ExcelFile() core calculations.
		String containing the name of the column corresponding to the CO2 concentration in the fluid, in wt%.

	Returns
	-------
	pandas DataFrame
		Original data passed plus newly calculated values are returned.
	"""
	convData = data.copy()

	MPO_H2O_list = []
	MPO_CO2_list = []
	for index, row in convData.iterrows():
		MPO_H2O_list.append(row[H2O_colname] / oxideMass["H2O"])
		MPO_CO2_list.append(row[CO2_colname] / oxideMass["CO2"])

	convData["MPO_H2O"] = MPO_H2O_list
	convData["MPO_CO2"] = MPO_CO2_list
	convData["XH2O_fl"] = convData["MPO_H2O"] / (convData["MPO_H2O"] + convData["MPO_CO2"])
	convData["XCO2_fl"] = convData["MPO_CO2"] / (convData["MPO_H2O"] + convData["MPO_CO2"])

	del convData["MPO_H2O"]
	del convData["MPO_CO2"]

	return convData

def get_oxides(sample):
	"""
	Returns a sample composition with only compositional oxide data, removing any extranneous data.
	Useful when passing a self-defined sample (e.g. dict or pandas Series) to a some VESIcal function.

	Parameters
	----------
	sample: pandas Series, dictionary
		A sample composition plus other sample information

	Returns
	-------
	Sample passed as > Returned as
			pandas Series > pandas Series
			dictionary > dictionary

				Sample composition with extranneous information removed.
	"""

	clean = {oxide:  sample[oxide] for oxide in oxides}

	if isinstance(sample, dict):
		return clean
	if isinstance(sample, pd.core.series.Series):
		return pd.Series(clean)

def rename_duplicates(df, suffix='-duplicate-'):
	appendents = (suffix + df.groupby(level=0).cumcount().astype(str).replace('0','')).replace(suffix, '')
	return df.set_index(df.index.astype(str) + appendents)

#----------DEFINE SOME NORMALIZATION METHODS-----------#

def normalize(sample):
	"""Normalizes an input composition to 100%. This is the 'standard' normalization routine.

	Parameters
	----------
	sample:    pandas Series, dictionary, pandas DataFrame, or ExcelFile object
		A single composition can be passed as a dictionary. Multiple compositions can be passed either as
		a pandas DataFrame or an ExcelFile object. Compositional information as oxides must be present.

	Returns
	-------
	Sample passed as > Returned as
		pandas Series > pandas Series
		dictionary > dictionary
		pandas DataFrame > pandas DataFrame
		ExcelFile object > pandas DataFrame

			Normalized major element oxides.
	"""
	def single_normalize(sample):
		single_sample = sample
		return {k: 100.0 * v / sum(single_sample.values()) for k, v in single_sample.items()}

	def multi_normalize(sample):
		multi_sample = sample.copy()
		multi_sample["Sum"] = sum([multi_sample[oxide] for oxide in oxides])
		for column in multi_sample:
			if column in oxides:
				multi_sample[column] = 100.0*multi_sample[column]/multi_sample["Sum"]

		del multi_sample["Sum"]
		return multi_sample

	if isinstance(sample, dict):
		_sample = sample.copy()
		return single_normalize(_sample)
	elif isinstance(sample, pd.core.series.Series):
		_sample = pd.Series(sample.copy())
		sample_dict = sample.to_dict()
		return pd.Series(single_normalize(sample_dict))
	elif isinstance(sample, ExcelFile):
		_sample = sample
		data = _sample.data
		return multi_normalize(data)
	elif isinstance(sample, pd.DataFrame):
		return multi_normalize(sample)


def normalize_FixedVolatiles(sample):
	""" Normalizes major element oxides to 100 wt%, including volatiles. The volatile
	wt% will remain fixed, whilst the other major element oxides are reduced proportionally
	so that the total is 100 wt%.

	Parameters
	----------
	sample:    pandas Series, dictionary, pandas DataFrame, or ExcelFile object
		Major element oxides in wt%

	Returns
	-------
	Sample passed as > Returned as
		pandas Series > pandas Series
		dictionary > dictionary
		pandas DataFrame > pandas DataFrame
		ExcelFile object > pandas DataFrame

			Normalized major element oxides.
	"""
	def single_FixedVolatiles(sample):
		normalized = pd.Series({},dtype=float)
		volatiles = 0
		if 'CO2' in list(_sample.index):
			volatiles += _sample['CO2']
		if 'H2O' in list(_sample.index):
			volatiles += _sample['H2O']

		for ox in list(_sample.index):
			if ox != 'H2O' and ox != 'CO2':
				normalized[ox] = _sample[ox]

		normalized = normalized/np.sum(normalized)*(100-volatiles)

		if 'CO2' in list(_sample.index):
			normalized['CO2'] = _sample['CO2']
		if 'H2O' in list(_sample.index):
			normalized['H2O'] = _sample['H2O']

		return normalized

	def multi_FixedVolatiles(sample):
		multi_sample = sample.copy()
		multi_sample["Sum_anhy"] = sum([multi_sample[oxide] for oxide in anhydrous_oxides])
		multi_sample["Sum_vols"] = sum([multi_sample[vol] for vol in volatiles])
		for column in multi_sample:
			if column in anhydrous_oxides:
				multi_sample[column] = 100.0*multi_sample[column]/multi_sample["Sum_anhy"]
				multi_sample[column] = multi_sample[column] / (100.0/(100.0-multi_sample["Sum_vols"]))
		del multi_sample["Sum_anhy"]
		del multi_sample["Sum_vols"]

		return multi_sample

	if isinstance(sample, dict):
		_sample = pd.Series(sample.copy())
		return single_FixedVolatiles(_sample).to_dict()
	elif isinstance(sample, pd.core.series.Series):
		_sample = pd.Series(sample.copy())
		return single_FixedVolatiles(_sample)
	elif isinstance(sample, ExcelFile):
		_sample = sample
		data = _sample.data
		return multi_FixedVolatiles(data)
	elif isinstance(sample, pd.DataFrame):
		return multi_FixedVolatiles(sample)
	else:
		raise InputError("The composition input must be a pandas Series or dictionary for single sample \
							or a pandas DataFrame or ExcelFile object for multi-sample.")


def normalize_AdditionalVolatiles(sample):
	"""Normalises major element oxide wt% to 100%, assuming it is volatile-free. If
	H2O or CO2 are passed to the function, their un-normalized values will be retained
	in addition to the normalized non-volatile oxides, summing to >100%.

	Parameters
	----------
	sample:    pandas Series, dictionary, pandas DataFrame, or ExcelFile object
		Major element oxides in wt%

	Returns
	-------
	Sample passed as > Returned as
		pandas Series > pandas Series
		dictionary > dictionary
		pandas DataFrame > pandas DataFrame
		ExcelFile object > pandas DataFrame

			Normalized major element oxides.
	"""

	def single_AdditionalVolatiles(sample):
		normalized = pd.Series({}, dtype=float)
		for ox in list(_sample.index):
			if ox != 'H2O' and ox != 'CO2':
				normalized[ox] = _sample[ox]

		normalized = normalized/np.sum(normalized)*100
		if 'H2O' in _sample.index:
			normalized['H2O'] = _sample['H2O']
		if 'CO2' in _sample.index:
			normalized['CO2'] = _sample['CO2']

		return normalized

	def multi_AdditionalVolatiles(sample):
		multi_sample = sample.copy()
		multi_sample["Sum"] = sum([multi_sample[oxide] for oxide in anhydrous_oxides])
		for column in multi_sample:
			if column in anhydrous_oxides:
				multi_sample[column] = 100.0*multi_sample[column]/multi_sample["Sum"]

		del multi_sample["Sum"]
		return multi_sample

	if isinstance(sample, dict):
		_sample = pd.Series(sample.copy())
		return single_AdditionalVolatiles(_sample).to_dict()
	elif isinstance(sample, pd.core.series.Series):
		_sample = pd.Series(sample.copy())
		return single_AdditionalVolatiles(sample)
	elif isinstance(sample, ExcelFile):
		_sample = sample
		data = _sample.data
		return multi_AdditionalVolatiles(data)
	elif isinstance(sample, pd.DataFrame):
		return multi_AdditionalVolatiles(sample)
	else:
		raise InputError("The composition input must be a pandas Series or dictionary for single sample \
							or a pandas DataFrame or ExcelFile object for multi-sample.")

#------------DEFINE MAJOR CLASSES-------------------#
class ExcelFile(object):
	"""An excel file with sample names and oxide compositions

	Attributes
	----------
		filename: str
			Path to the excel file, e.g., "my_file.xlsx". This always needs to be passed, even if the user is passing a pandas DataFrame
			rather than an Excel file. If passing a DataFrame, filename should be set to None.

		sheet_name: str
			OPTIONAL. Default value is 0 which gets the first sheet in the excel spreadsheet file. This implements the pandas.
			read_excel() sheet_name parameter. But functionality to read in more than one sheet at a time (e.g., pandas.read_excel(sheet_name=None))
			is not yet imlpemented in VESIcal. From the pandas 1.0.4 documentation:
				Available cases:
					- Defaults to 0: 1st sheet as a DataFrame
					- 1: 2nd sheet as a DataFrame
					- "Sheet1": Load sheet with name “Sheet1”

		input_type: str
			OPTIONAL. Default is 'wtpercent'. String defining whether the oxide composition is given in wt percent
			("wtpercent", which is the default), mole percent ("molpercent"), or mole fraction ("molfrac").

		label: str
			OPTIONAL. Default is 'Label'. Name of the column within the passed Excel file referring to sample names.

	kwargs
	------
		dataframe: pandas DataFrame
				OPTIONAL. Default is None in which case this argument is ignored. This argument is used when the user wishes to turn
				a pandas DataFrame into an ExcelFile object, for example when user data is already in python rather than being imported
				from an Excel file. If using this option, pass None to filename.
	"""

	def __init__(self, filename, sheet_name=0, input_type='wtpercent', label='Label', **kwargs):
		"""Return an ExcelFile object whoes parameters are defined here."""

		if isinstance(sheet_name, str) or isinstance(sheet_name, int):
			pass
		else:
			raise InputError("If sheet_name is passed, it must be of type str or int. Currently, VESIcal cannot import more than one sheet at a time.")

		self.input_type = input_type

		if 'dataframe' in kwargs:
			data = kwargs['dataframe']
			data = data.rename_axis('Label')
		else:
			data = pd.read_excel(filename, sheet_name=sheet_name)

			try:
				data = data.set_index(label)
			except:
				label_found = False
				for col in data.columns:
					if col in oxides:
						pass
					else:
						data = data.set_index(col)
						label_found = True
						w.warn("No Label column given, so column '" + str(col) + "' was chosen for you. To choose your own, set label='<column-name>'.",RuntimeWarning,stacklevel=2)
						break
				if label_found == False:
					data.index.name = 'Label'
					w.warn("No Label column given, so one was created for you. To choose your own, set label='<column-name>'.",RuntimeWarning,stacklevel=2)

		data = rename_duplicates(data) #handle any duplicated sample names
		data = data.fillna(0) #fill in any missing data with 0's

		if 'model' in kwargs:
			w.warn("You don't need to pass a model here, so it will be ignored. You can specify a model when performing calculations on your dataset (e.g., calculate_dissolved_volatiles())",RuntimeWarning,stacklevel=2)

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
		Returns oxide composition of a single sample from a user-imported excel file as a dictionary

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

	def get_XH2O_fluid(self, sample, temperature, pressure, H2O, CO2):
		"""An internally used function to calculate fluid composition.

		Parameters
		----------
		sample: dictionary
			Sample composition in wt% oxides

		temperature: float
			Temperature in degrees C.

		pressure: float
			Pressure in bars

		H2O: float
			wt% H2O in the system

		CO2: float
			wt% CO2 in the system

		Returns
		-------
		float
			Mole fraction of H2O in the H2O-CO2 fluid

		"""
		pressureMPa = pressure / 10.0

		bulk_comp = {oxide:  sample[oxide] for oxide in oxides}
		bulk_comp["H2O"] = H2O
		bulk_comp["CO2"] = CO2
		feasible = melts.set_bulk_composition(bulk_comp)

		output = melts.equilibrate_tp(temperature, pressureMPa, initialize=True)
		(status, temperature, pressureMPa, xmlout) = output[0]
		fluid_comp = melts.get_composition_of_phase(xmlout, phase_name='Fluid', mode='component')
			#NOTE mode='component' returns endmember component keys with values in mol fraction.

		if "Water" in fluid_comp:
			H2O_fl = fluid_comp["Water"]
		else:
			H2O_fl = 0.0

		# if H2O_fl == 0:
		#   raise SaturationError("Composition not fluid saturated.")

		return H2O_fl

	def save_excelfile(self, filename, calculations, sheet_name=None):
		"""
		Saves data calculated by the user in batch processing mode (using the ExcelFile class methods) to an organized
		excel file, with the original user data plus any calculated data.

		Parameters
		----------
		filename: string
			Name of the file. Extension (.xlsx) should be passed along with the name itself, all in quotes (e.g., 'myfile.xlsx').

		calculations: pandas DataFrame or list of pandas DataFrames
			A single variable or list of variables containing calculated outputs from any of the core ExcelFile functions: 
			calculate_dissolved_volatiles, calculate_equilibrium_fluid_comp, and calculate_saturation_pressure.

		sheet_name: None, string, or list
			OPTIONAL. Default value is None. Allows user to set the name of the sheet or sheets written to the Excel file.

		Returns
		-------
		Excel File
			Creates and saves an Excel file with data from each calculation saved to its own sheet.
		"""
		if isinstance(calculations, list):
			if isinstance(sheet_name, list) or sheet_name is None:
				pass
		else:
			calculations = [calculations]
		with pd.ExcelWriter(filename) as writer:
			self.data.to_excel(writer, 'Original_User_Data')
			if sheet_name is None:
				for n, df in enumerate(calculations):
					df.to_excel(writer, 'Calc%s' % n)
			elif isinstance(sheet_name, list):
				pass
			else:
				sheet_name = [sheet_name]
			if isinstance(sheet_name, list):
				if len(sheet_name) == len(calculations):
					pass
				else:
					raise InputError("calculations and sheet_name must have the same length")

				for i in range(len(calculations)):
					if isinstance(sheet_name[i], str):
						calculations[i].to_excel(writer, sheet_name[i])
					else:
						raise InputError("if sheet_name is passed, it must be list of strings")

		return print("Saved " + str(filename))

	def calculate_dissolved_volatiles(self, temperature, pressure, X_fluid=1, print_status=True, model='MagmaSat', record_errors=False, **kwargs):
		"""
		Calculates the amount of H2O and CO2 dissolved in a magma at the given P/T conditions and fluid composition. Fluid composition
		will be matched to within 0.0001 mole fraction.

		Parameters
		----------
		temperature: float, int, or str
			Temperature, in degrees C. Can be passed as float, in which case the
			passed value is used as the temperature for all samples. Alternatively, temperature information for each individual
			sample may already be present in the ExcelFile object. If so, pass the str value corresponding to the column
			title in the ExcelFile object.

		presure: float, int, or str
			Pressure, in bars. Can be passed as float or int, in which case the
			passed value is used as the pressure for all samples. Alternatively, pressure information for each individual
			sample may already be present in the ExcelFile object. If so, pass the str value corresponding to the column
			title in the ExcelFile object.

		X_fluid: float, int, or str
			OPTIONAL: Default value is 1. The mole fraction of H2O in the H2O-CO2 fluid. X_fluid=1 is a pure H2O fluid. X_fluid=0 is a
			pure CO2 fluid. Can be passed as a float or int, in which case the passed value is used as the X_fluid for all samples.
			Alternatively, X_fluid information for each individual sample may already be present in the ExcelFile object. If so, pass
			the str value corresponding to the column title in the ExcelFile object.

		print_status: bool
			OPTIONAL: The default value is True, in which case the progress of the calculation will be printed to the terminal.
			If set to False, nothing will be printed. MagmaSat calculations tend to be slow, and so a value of True is recommended
			for most use cases.

		model: string
			OPTIONAL: Default is 'MagmaSat'. Any other model name can be passed here.

		record_errors: bool
			OPTIONAL: If True, any errors arising during the calculation will be recorded as a column.

		Returns
		-------
		pandas DataFrame
			Original data passed plus newly calculated values are returned.
		"""
		data = self.preprocess_sample(self.data)
		dissolved_data = data.copy()

		if isinstance(temperature, str):
			file_has_temp = True
			temp_name = temperature
		elif isinstance(temperature, float) or isinstance(temperature, int):
			file_has_temp = False
		else:
			raise InputError("temp must be type str or float or int")

		if isinstance(pressure, str):
			file_has_press = True
			press_name = pressure
		elif isinstance(pressure, float) or isinstance(pressure, int):
			file_has_press = False
		else:
			raise InputError("pressure must be type str or float or int")

		if isinstance(X_fluid, str):
			file_has_X = True
			X_name = X_fluid
		elif isinstance(X_fluid, float) or isinstance(X_fluid, int):
			file_has_X = False
			if X_fluid != 0 and X_fluid !=1:
				if X_fluid < 0.001 or X_fluid > 0.999:
					raise InputError("X_fluid is calculated to a precision of 0.0001 mole fraction. \
									 Value for X_fluid must be between 0.0001 and 0.9999.")
		else:
			raise InputError("X_fluid must be type str or float or int")

		H2Ovals = []
		CO2vals = []
		warnings = []
		errors = []
		if model in get_models(models='mixed'):
			for index, row in dissolved_data.iterrows():
				try:
					if file_has_temp == True:
						temperature = row[temp_name]

					if file_has_press == True:
						pressure = row[press_name]

					if file_has_X == True:
						X_fluid = row[X_name]

					bulk_comp = {oxide:  row[oxide] for oxide in oxides}
					calc = calculate_dissolved_volatiles(sample=bulk_comp, pressure=pressure, temperature=temperature,
																	X_fluid=(X_fluid, 1-X_fluid), model=model,
																	silence_warnings=True, **kwargs)
					H2Ovals.append(calc.result['H2O_liq'])
					CO2vals.append(calc.result['CO2_liq'])
					warnings.append(calc.calib_check)
					errors.append('')
				except Exception as inst:
					H2Ovals.append(np.nan)
					CO2vals.append(np.nan)
					warnings.append('Calculation Failed.')
					errors.append(sys.exc_info()[0])
			dissolved_data["H2O_liq_VESIcal"] = H2Ovals
			dissolved_data["CO2_liq_VESIcal"] = CO2vals

			if file_has_temp == False:
				dissolved_data["Temperature_C_VESIcal"] = temperature
			if file_has_press == False:
				dissolved_data["Pressure_bars_VESIcal"] = pressure
			if file_has_X == False:
				dissolved_data["X_fluid_input_VESIcal"] = X_fluid
			dissolved_data["Model"] = model
			dissolved_data["Warnings"] = warnings
			if record_errors == True:
				dissolved_data["Errors"] = errors

			return dissolved_data

		elif model == 'MagmaSat':
			XH2Ovals = []
			XCO2vals = []
			FluidProportionvals = []
			for index, row in dissolved_data.iterrows():
				if print_status == True:
					print("Calculating sample " + str(index))

				if file_has_temp == True:
					temperature = row[temp_name]
				if temperature <= 0:
					H2Ovals.append(np.nan)
					CO2vals.append(np.nan)
					XH2Ovals.append(np.nan)
					XCO2vals.append(np.nan)
					FluidProportionvals.append(np.nan)
					warnings.append("Sample skipped. Bad temperature.")
					errors.append(sys.exc_info()[0])
					w.warn("Temperature for sample " + str(index) + " is <=0. Skipping sample.", stacklevel=2)

				if file_has_press == True:
					pressure = row[press_name]
				if temperature >0 and pressure <= 0:
					H2Ovals.append(np.nan)
					CO2vals.append(np.nan)
					XH2Ovals.append(np.nan)
					XCO2vals.append(np.nan)
					FluidProportionvals.append(np.nan)
					warnings.append("Sample skipped. Bad pressure.")
					errors.append(sys.exc_info()[0])
					w.warn("Pressure for sample " + str(index) + " is <=0. Skipping sample.", stacklevel=2)

				if file_has_X == True:
					X_fluid = row[X_name]
				if temperature >0 and pressure >0 and X_fluid <0:
					H2Ovals.append(np.nan)
					CO2vals.append(np.nan)
					XH2Ovals.append(np.nan)
					XCO2vals.append(np.nan)
					FluidProportionvals.append(np.nan)
					warnings.append("Sample skipped. Bad X_fluid.")
					errors.append(sys.exc_info()[0])
					w.warn("X_fluid for sample " + str(index) + " is <0. Skipping sample.", stacklevel=2)

				if temperature >0 and pressure >0 and X_fluid >1:
					H2Ovals.append(np.nan)
					CO2vals.append(np.nan)
					XH2Ovals.append(np.nan)
					XCO2vals.append(np.nan)
					FluidProportionvals.append(np.nan)
					warnings.append("Sample skipped. Bad X_fluid.")
					errors.append(sys.exc_info()[0])
					w.warn("X_fluid for sample " + str(index) + " is >1. Skipping sample.", stacklevel=2)

				if temperature > 0 and pressure > 0 and X_fluid >=0 and X_fluid <=1:
					try:
						bulk_comp = {oxide:  row[oxide] for oxide in oxides}
						calc = calculate_dissolved_volatiles(sample=bulk_comp, pressure=pressure, temperature=temperature,
																		X_fluid=X_fluid, model=model, silence_warnings=True,
																		verbose=True)
						H2Ovals.append(calc.result['H2O_liq'])
						CO2vals.append(calc.result['CO2_liq'])
						XH2Ovals.append(calc.result['XH2O_fl'])
						XCO2vals.append(calc.result['XCO2_fl'])
						FluidProportionvals.append(calc.result['FluidProportion_wt'])
						warnings.append(calc.calib_check)
						errors.append('')

					except Exception as inst:
						H2Ovals.append(np.nan)
						CO2vals.append(np.nan)
						XH2Ovals.append(np.nan)
						XCO2vals.append(np.nan)
						FluidProportionvals.append(np.nan)
						warnings.append('Calculation Failed.')
						errors.append(sys.exc_info()[0])

			dissolved_data["H2O_liq_VESIcal"] = H2Ovals
			dissolved_data["CO2_liq_VESIcal"] = CO2vals

			if file_has_temp == False:
				dissolved_data["Temperature_C_VESIcal"] = temperature
			if file_has_press == False:
				dissolved_data["Pressure_bars_VESIcal"] = pressure
			if file_has_X == False:
				dissolved_data["X_fluid_input_VESIcal"] = X_fluid
			dissolved_data["Model"] = model
			dissolved_data["Warnings"] = warnings
			if record_errors == True:
				dissolved_data["Errors"] = errors

			return dissolved_data
		else:
			XH2Ovals = []
			XCO2vals = []
			FluidProportionvals = []
			for index, row in dissolved_data.iterrows():
				if file_has_temp == True:
					temperature = row[temp_name]
				if file_has_press == True:
					pressure = row[press_name]
				if file_has_X == True:
					X_fluid = row[X_name]
				bulk_comp = {oxide:  row[oxide] for oxide in oxides}
				if 'Water' in model:
					try:
						calc = calculate_dissolved_volatiles(sample=bulk_comp, pressure=pressure, temperature=temperature,
															 X_fluid=X_fluid, model=model, silence_warnings=True)
						H2Ovals.append(calc.result)
						warnings.append(calc.calib_check)
					except:
						H2Ovals.append(0)
						warnings.append('Calculation Failed #001')
				if 'Carbon' in model:
					try:
						calc = calculate_dissolved_volatiles(sample=bulk_comp, pressure=pressure, temperature=temperature,
															 X_fluid=X_fluid, model=model, silence_warnings=True)
						CO2vals.append(calc.result)
						warnings.append(calc.calib_check)
					except:
						CO2vals.append(0)
						warnings.append('Calculation Failed #002')

			if 'Water' in model:
				dissolved_data["H2O_liq_VESIcal"] = H2Ovals
			if 'Carbon' in model:
				dissolved_data["CO2_liq_VESIcal"] = CO2vals
			if file_has_temp == False:
				dissolved_data["Temperature_C_VESIcal"] = temperature
			if file_has_press == False:
				dissolved_data["Pressure_bars_VESIcal"] = pressure
			if file_has_X == False:
				dissolved_data["X_fluid_input_VESIcal"] = X_fluid
			dissolved_data["Model"] = model
			dissolved_data["Warnings"] = warnings

			return dissolved_data

	def calculate_equilibrium_fluid_comp(self, temperature, pressure=None, print_status=False, model='MagmaSat', **kwargs):
		"""
		Returns H2O and CO2 concentrations in wt% or mole fraction in a fluid in equilibrium with the given sample(s) at the given P/T condition.

		Parameters
		----------
		sample: ExcelFile object
			Compositional information on samples in oxides.

		temperature: float, int, or str
			Temperature, in degrees C. Can be passed as float, in which case the
			passed value is used as the temperature for all samples. Alternatively, temperature information for each individual
			sample may already be present in the ExcelFile object. If so, pass the str value corresponding to the column
			title in the  ExcelFile object.

		presure: float, int, or str
			Pressure, in bars. Can be passed as float or int, in which case the
			passed value is used as the pressure for all samples. Alternatively, pressure information for each individual
			sample may already be present in the ExcelFile object. If so, pass the str value corresponding to the column
			title in the ExcelFile object.

		model: string
			OPTIONAL: Default is 'MagmaSat'. Any other model name can be passed here.

		Returns
		-------
		pandas DataFrame
			Original data passed plus newly calculated values are returned.
		"""
		data = self.preprocess_sample(self.data)
		fluid_data = data.copy()

		if isinstance(temperature, str):
			file_has_temp = True
			temp_name = temperature
		elif isinstance(temperature, float) or isinstance(temperature, int):
			file_has_temp = False
		else:
			raise InputError("temp must be type str or float or int")

		if isinstance(pressure, str):
			file_has_press = True
			press_name = pressure
		elif isinstance(pressure, float) or isinstance(pressure, int) or type(pressure) == type(None):
			file_has_press = False
		else:
			raise InputError("pressure must be type str or float or int")

		H2Ovals = []
		CO2vals = []
		warnings = []
		if model in get_models(models='mixed') or model == "MooreWater":
			for index, row in fluid_data.iterrows():
				try:
					if file_has_temp == True:
						temperature = row[temp_name]
					if file_has_press == True:
						pressure = row[press_name]
					bulk_comp = {oxide:  row[oxide] for oxide in oxides}
					calc = calculate_equilibrium_fluid_comp(sample=bulk_comp, pressure=pressure, temperature=temperature,
															model=model, silence_warnings=True, **kwargs)

					H2Ovals.append(calc.result['H2O'])
					CO2vals.append(calc.result['CO2'])
					if calc.result['H2O'] == 0 and calc.result['CO2'] == 0:
						warnings.append(calc.calib_check + "Sample not saturated at these conditions")
					else:
						warnings.append(calc.calib_check)
				except:
					H2Ovals.append(np.nan)
					CO2vals.append(np.nan)
					warnings.append("Calculation Failed.")
			fluid_data["XH2O_fl_VESIcal"] = H2Ovals
			fluid_data["XCO2_fl_VESIcal"] = CO2vals
			if file_has_temp == False:
				fluid_data["Temperature_C_VESIcal"] = temperature
			if file_has_press == False:
				fluid_data["Pressure_bars_VESIcal"] = pressure
			fluid_data["Model"] = model
			fluid_data["Warnings"] = warnings

			return fluid_data
		elif model == 'MagmaSat':
			for index, row in fluid_data.iterrows():
				if print_status == True:
					print("Calculating sample " + str(index))

				if file_has_temp == True:
					temperature = row[temp_name]
				if temperature <= 0:
					H2Ovals.append(np.nan)
					CO2vals.append(np.nan)
					warnings.append("Calculation skipped. Bad temperature.")
					w.warn("Temperature for sample " + str(index) + " is <=0. Skipping sample.", stacklevel=2)

				if file_has_press == True:
					pressure = row[press_name]
				if temperature >0 and pressure <= 0:
					H2Ovals.append(np.nan)
					CO2vals.append(np.nan)
					warnings.append("Calculation skipped. Bad pressure.")
					w.warn("Pressure for sample " + str(index) + " is <=0. Skipping sample.", stacklevel=2)

				if temperature > 0 and pressure > 0:
					try:
						bulk_comp = {oxide:  row[oxide] for oxide in oxides}
						calc = calculate_equilibrium_fluid_comp(sample=bulk_comp, pressure=pressure, temperature=temperature, model=model, silence_warnings=True)

						H2Ovals.append(calc.result['H2O'])
						CO2vals.append(calc.result['CO2'])
						if calc.result['H2O'] == 0 and calc.result['CO2'] == 0:
							warnings.append(calc.calib_check + "Sample not saturated at these conditions")
						else:
							warnings.append(calc.calib_check)
					except:
						H2Ovals.append(np.nan)
						CO2vals.append(np.nan)
						warnings.append("Calculation Failed.")

			fluid_data["XH2O_fl_VESIcal"] = H2Ovals
			fluid_data["XCO2_fl_VESIcal"] = CO2vals
			if file_has_temp == False:
				fluid_data["Temperature_C_VESIcal"] = temperature
			if file_has_press == False:
				fluid_data["Pressure_bars_VESIcal"] = pressure
			fluid_data["Model"] = model
			fluid_data["Warnings"] = warnings

			return fluid_data

		else:
			saturated = []
			for index, row in fluid_data.iterrows():
				try:
					if file_has_temp == True:
						temperature = row[temp_name]
					if file_has_press == True:
						pressure = row[press_name]
					bulk_comp = {oxide:  row[oxide] for oxide in oxides}
					calc = calculate_equilibrium_fluid_comp(sample=bulk_comp, pressure=pressure, temperature=temperature, model=model, silence_warnings=True)
					saturated.append(calc.result)
					warnings.append(calc.calib_check)
				except:
					saturated.append(np.nan)
					warnings.append("Calculation Failed.")
			fluid_data["Saturated_VESIcal"] = saturated
			if file_has_temp == False:
				fluid_data["Temperature_C_VESIcal"] = temperature
			if file_has_press == False:
				fluid_data["Pressure_bars_VESIcal"] = pressure
			fluid_data["Model"] = model
			fluid_data["Warnings"] = warnings

			return fluid_data

	def calculate_saturation_pressure(self, temperature, print_status=True, model='MagmaSat', **kwargs):
		"""
		Calculates the saturation pressure of multiple sample compositions in the ExcelFile.

		Parameters
		----------
		temperature: float, int, or str
			Temperature at which to calculate saturation pressures, in degrees C. Can be passed as float or int, in which case the
			passed value is used as the temperature for all samples. Alternatively, temperature information for each individual
			sample may already be present in the passed ExcelFile object. If so, pass the str value corresponding to the column
			title in the passed ExcelFile object.

		print_status: bool
			OPTIONAL: The default value is True, in which case the progress of the calculation will be printed to the terminal.
			If set to False, nothing will be printed. MagmaSat calculations tend to be slow, and so a value of True is recommended
			more most use cases.

		model: string
			OPTIONAL: Default is 'MagmaSat'. Any other model name can be passed here.

		Returns
		-------
		pandas DataFrame object
			Values returned are saturation pressure in bars, the mass of fluid present, and the composition of the
			fluid present.
		"""
		data = self.preprocess_sample(self.data)
		satp_data = data.copy()

		if isinstance(temperature, str):
			file_has_temp = True
			temp_name = temperature
		elif isinstance(temperature, float) or isinstance(temperature, int):
			file_has_temp = False
		else:
			raise InputError("temperature must be type str or float or int")

		if model != 'MagmaSat':
			satP = []
			warnings = []
			piStar = []
			for index, row in satp_data.iterrows():
				try:
					if file_has_temp == True:
						temperature = row[temp_name]
					bulk_comp = {oxide:  row[oxide] for oxide in oxides}
					calc = calculate_saturation_pressure(sample=bulk_comp, temperature=temperature,
														 model=model, silence_warnings=True, **kwargs)
					satP.append(calc.result)
					warnings.append(calc.calib_check)
					if model == 'Shishkina':
						piStar.append(default_models['Shishkina'].models[1].PiStar(bulk_comp))
				except:
					satP.append(np.nan)
					warnings.append("Calculation Failed")
			satp_data["SaturationP_bars_VESIcal"] = satP
			if file_has_temp == False:
				satp_data["Temperature_C_VESIcal"] = temperature
			satp_data["Model"] = model
			satp_data["Warnings"] = warnings
			if len(piStar)>0:
				satp_data['PiStar_VESIcal'] = piStar

			return satp_data

		else:
			satP = []
			flmass = []
			flH2O = []
			flCO2 = []
			flsystem_wtper = []
			warnings = []
			for index, row in satp_data.iterrows():
				if print_status == True:
					print("Calculating sample " + str(index))

				if file_has_temp == True:
					temperature = row[temp_name]
				if temperature <= 0:
					satP.append(np.nan)
					flmass.append(np.nan)
					flsystem_wtper.append(np.nan)
					flH2O.append(np.nan)
					flCO2.append(np.nan)
					warnings.append("Calculation skipped. Bad temperature.")
					w.warn("Temperature for sample " + str(index) + " is <=0. Skipping sample.", stacklevel=2)

				if temperature > 0:
					try:
						bulk_comp = {oxide:  row[oxide] for oxide in oxides}
						calc = calculate_saturation_pressure(sample=bulk_comp, temperature=temperature, model=model, verbose=True, silence_warnings=True)
						satP.append(calc.result["SaturationP_bars"])
						flmass.append(calc.result["FluidMass_grams"])
						flsystem_wtper.append(calc.result["FluidProportion_wt"])
						flH2O.append(calc.result["XH2O_fl"])
						flCO2.append(calc.result["XCO2_fl"])
						warnings.append(calc.calib_check)
					except:
						satP.append(np.nan)
						flmass.append(np.nan)
						flsystem_wtper.append(np.nan)
						flH2O.append(np.nan)
						flCO2.append(np.nan)
						warnings.append("Calculation Failed")

			satp_data["SaturationP_bars_VESIcal"] = satP
			if file_has_temp == False:
				satp_data["Temperature_C_VESIcal"] = temperature
			satp_data["XH2O_fl_VESIcal"] = flH2O
			satp_data["XCO2_fl_VESIcal"] = flCO2
			satp_data["FluidMass_grams_VESIcal"] = flmass
			satp_data["FluidSystem_wt_VESIcal"] = flsystem_wtper
			satp_data["Model"] = model
			satp_data["Warnings"] = warnings

			if print_status == True:
				print("Done!")

			return satp_data

class CalibrationRange(object):
	""" The CalibrationRange object allows the range of allowable parameters to be specified and
	used in checking and reporting of the results.
	"""
	def __init__(self, parameter_name, value, checkfunction=None, units='', model_name='',
				fail_msg='',fail_dict={}, pass_msg='', pass_dict={}, description_msg='', description_dict={}):
		self.parameter_name = parameter_name
		self.value = value
		self.checkfunction = checkfunction
		self.units = units
		self.model_name = model_name
		self.fail_msg = (copy(fail_msg), copy(fail_dict))
		self.pass_msg = (copy(pass_msg), copy(pass_dict))
		self.description_msg = (copy(description_msg), copy(description_dict))

	def check(self,parameters):
		"""Method for checking whether parameters satisfy the calibration range."""
		if self.parameter_name in parameters:
			return self.checkfunction(self.value,parameters[self.parameter_name])
		else:
			return None

	def string(self,parameters,report_nonexistance=True):
		"""Returns a string statement of the calibration check"""
		if type(parameters) == type(None):
			msgdict = self.description_msg[1]
			if type(self.value) == float or type(self.value) == int:
				msgdict['calib_val'] = self.value
			elif type(self.value) == list or type(self.value) == tuple or type(self.value) == np.ndarray:
				for i in range(len(self.value)):
					msgdict['calib_val'+str(i)] = self.value[i]
			if 'param_name' not in msgdict:
				msgdict['param_name'] = self.parameter_name
			if 'units' not in msgdict:
				msgdict['units'] = self.units
			if 'model_name' not in msgdict:
				msgdict['model_name'] = self.model_name
			return self.description_msg[0].format(**msgdict)
		else:
			check = self.check(parameters)
			if check == True:
				msgdict = self.pass_msg[1]
				msgdict['param_val'] = parameters[self.parameter_name]
				if type(self.value) == float or type(self.value) == int:
					msgdict['calib_val'] = self.value
				elif type(self.value) == list or type(self.value) == tuple or type(self.value) == np.ndarray:
					for i in range(len(self.value)):
						msgdict['calib_val'+str(i)] = self.value[i]
				if 'param_name' not in msgdict:
					msgdict['param_name'] = self.parameter_name
				if 'units' not in msgdict:
					msgdict['units'] = self.units
				if 'model_name' not in msgdict:
					msgdict['model_name'] = self.model_name
				return self.pass_msg[0].format(**msgdict)
			elif check == False:
				msgdict = self.fail_msg[1]
				msgdict['param_val'] = parameters[self.parameter_name]
				if type(self.value) == float or type(self.value) == int:
					msgdict['calib_val'] = self.value
				elif type(self.value) == list or type(self.value) == tuple or type(self.value) == np.ndarray:
					for i in range(len(self.value)):
						msgdict['calib_val'+str(i)] = self.value[i]
				if 'param_name' not in msgdict:
					msgdict['param_name'] = self.parameter_name
				if 'units' not in msgdict:
					msgdict['units'] = self.units
				if 'model_name' not in msgdict:
					msgdict['model_name'] = self.model_name
				return self.fail_msg[0].format(**msgdict)
			else:
				if report_nonexistance == True:
					return "A value for {} was not provided.".format(self.parameter_name)
				else:
					return ''

# class old_CalibrationRange(object):
# 	""" The CalibrationRange object allows the range of allowable parameters to be specified and
# 	used in checking and reporting of the results.
# 	"""
# 	def __init__(self,parameter_name,value,unit='',modelname='',explanation_string=None,
# 				 parameter_string=None,value_fmt="{:.1f}"):
# 		self.parameter_name = parameter_name
# 		self.value = value
# 		self.value_fmt = value_fmt
# 		self.model_name = modelname
# 		self.unit = unit
# 		self.explanation_string = explanation_string
# 		if parameter_string is not None:
# 			self.parameter_string = parameter_string
# 		else:
# 			self.parameter_string = parameter_name
#
# 	@abstractmethod
# 	def check(self,parameters):
# 		"""Method for checking whether parameters satisfy the calibration range."""
# 		return True
#
# 	@abstractmethod
# 	def string(self,parameters):
# 		"""Returns a string statement of the calibration check"""
# 		return 'No string return defined. '

class Model(object):
	"""The model object implements a volatile solubility model. It is composed
	of the methods needed to evaluate :func:`VESIcal.calculate_dissolved_volatiles`,
	:func:`VESIcal.calculate_equilibrium_fluid_comp`, and :func:`calculate_saturation_pressure`. The
	fugacity and activity models for the volatiles species must be specified,
	defaulting to ideal.
	"""

	def __init__(self):
		self.set_volatile_species(None)
		self.set_fugacity_model(fugacity_idealgas())
		self.set_activity_model(activity_idealsolution())
		self.set_calibration_ranges([])
		self.set_solubility_dependence(False)

	def set_volatile_species(self,volatile_species):
		if type(volatile_species) == str:
			volatile_species = [volatile_species]
		elif type(volatile_species) != list:
			raise InputError("volatile_species must be a str or list.")
		self.volatile_species = volatile_species

	def set_fugacity_model(self,fugacity_model):
		self.fugacity_model = fugacity_model

	def set_activity_model(self,activity_model):
		self.activity_model = activity_model

	def set_calibration_ranges(self,calibration_ranges):
		self.calibration_ranges = calibration_ranges

	def set_solubility_dependence(self,solubility_dependence):
		self.solubility_dependence = solubility_dependence

	@abstractmethod
	def calculate_dissolved_volatiles(self,**kwargs):
		pass

	@abstractmethod
	def calculate_equilibrium_fluid_comp(self,**kwargs):
		pass

	@abstractmethod
	def calculate_saturation_pressure(self,**kwargs):
		pass

	@abstractmethod
	def preprocess_sample(self,**kwargs):
		pass

	# @abstractmethod
	def check_calibration_range(self,parameters,report_nonexistance=True):
		""" Checks whether the given parameters are within the ranges defined by the
		CalibrationRange objects for the model and its fugacity and activity models. An empty
		string will be returned if all parameters are within the calibration range. If a
		parameter is not within the calibration range, a description of the problem will be
		returned in the string.

		Parameters
		----------
		parameters 	dict
			Dictionary keys are the names of the parameters to be checked, e.g., pressure
			temperature, SiO2, etc. Values are the values of each parameter. A complete set
			need not be given.

		Returns
		-------
		str
			String description of any parameters falling outside of the calibration range.
		"""
		s = ''
		for cr in self.calibration_ranges:
			if cr.check(parameters) == False:
				s += cr.string(parameters,report_nonexistance)
		for cr in self.fugacity_model.calibration_ranges:
			if cr.check(parameters) == False:
				s += cr.string(parameters,report_nonexistance)
		for cr in self.activity_model.calibration_ranges:
			if cr.check(parameters) == False:
				s += cr.string(parameters,report_nonexistance)
		return s

	def get_calibration_range(self):
		""" Returns a string describing the calibration ranges defined by the CalibrationRange
		objects for each model, and its associated fugacity and activity models.

		Returns
		-------
		str
			String description of the calibration range objects."""
		s = ''
		for cr in self.calibration_ranges:
			s += cr.string(None)
		for cr in self.fugacity_model.calibration_ranges:
			s += cr.string(None)
		for cr in self.activity_model.calibration_ranges:
			s += cr.string(None)
		return s


class FugacityModel(object):
	""" The fugacity model object is for implementations of fugacity models
	for individual volatile species, though it may depend on the mole
	fraction of other volatile species. It contains all the methods required
	to calculate the fugacity at a given pressure and mole fraction.
	"""

	def __init__(self):
		self.set_calibration_ranges([])

	def set_calibration_ranges(self,calibration_ranges):
		self.calibration_ranges = calibration_ranges

	@abstractmethod
	def fugacity(self,pressure,**kwargs):
		"""
		"""

	# @abstractmethod
	def check_calibration_range(self,parameters,report_nonexistance=True):
		s = ''
		for cr in self.calibration_ranges:
			if cr.check(parameters) == False:
				s += cr.string(parameters,report_nonexistance)
		return s



class activity_model(object):
	""" The activity model object is for implementing activity models
	for volatile species in melts. It contains all the methods required to
	evaluate the activity.
	"""
	def __init__(self):
		self.set_calibration_ranges([])

	def set_calibration_ranges(self,calibration_ranges):
		self.calibration_ranges = calibration_ranges

	@abstractmethod
	def activity(self,X,**kwargs):
		"""
		"""

	# @abstractmethod
	def check_calibration_range(self,parameters,report_nonexistance=True):
		s = ''
		for cr in self.calibration_ranges:
			if cr.check(parameters) == False:
				s += cr.string(parameters,report_nonexistance)
		return s



class Calculate(object):
	""" The Calculate object is a template for implementing user-friendly methods for
	running calculations using the volatile solubility models. All Calculate methods
	have a common workflow- sample is read in, preprocessed, the calculation is performed,
	the calibration range is checked, and the results stored.
	"""
	def __init__(self,sample,model='MagmaSat',silence_warnings=False,preprocess_sample=False,**kwargs):
		if model == 'MagmaSat':
			self.model = MagmaSat()
		elif type(model) == str:
			self.model = default_models[model]
		else:
			self.model = model

		self.sample = sample.copy()
		if preprocess_sample == True:
			self.sample = self.model.preprocess_sample(self.sample)

		self.result = self.calculate(sample=self.sample,**kwargs)
		self.calib_check = self.check_calibration_range(sample=self.sample,**kwargs)

		if self.calib_check is not None and silence_warnings == False:
			if self.calib_check != '':
				w.warn(self.calib_check,RuntimeWarning)

	@abstractmethod
	def calculate(self):
		""" """

	@abstractmethod
	def check_calibration_range(self):
		""" """

#-------------DEFAULT CALIBRATIONRANGE OBJECTS---------------#FIND CALIBRATION RANGES

def crf_EqualTo(calibval,paramval):
	return calibval == paramval
crmsg_EqualTo_pass = "The {param_name} ({param_val:.1f} {units}) is equal to {calib_val:.1f} {units} as required by the calibration range of the {model_name} model. "
crmsg_EqualTo_fail ="{param_name} ({param_val:.1f} {units}) is outside the calibration range of the {model_name} model ({calib_val:.1f} {units}). "
crmsg_EqualTo_description = "The {model_name} model is calibrated for {param_name} equal to {calib_val:.1f} {units}. "
crmsg_EqualTo_fail_AllisonTemp="All calculations for {model_name} are performed at 1200 C (inputted Temp={param_val:.1f} {units}). Allison et al. (2019) suggest the results are likely applicable between 1000-1400°C). "

def crf_GreaterThan(calibval,paramval):
	return paramval >= calibval
crmsg_GreaterThan_pass = "The {param_name} ({param_val:.1f} {units}) is greater than {calib_val:.1f} {units} as required by the calibration range of the {model_name} model. "
crmsg_GreaterThan_fail = "The {param_name} is outside the calibration range of {model_name} ({param_val:.1f}<{calib_val:.1f} {units}. "
crmsg_GreaterThan_description = "The {model_name} model is calibrated for {param_name} greater than {calib_val:.1f} {units}. "
# Warnings for specific models
crmsg_GreaterThan_fail_ShWaterSi = "{param_name}{units} exceeds the lower calibration limit of {model_name} ({calib_val:.1f} {units}, based on the minimum concentration minus 5% in the calibration dataset). " # This check warns users if SiO2<40 in the water model
Dix_40_fail="{param_name} ({param_val:.1f} {units})<40 wt%, which is the lower calibration limit of the Dixon (1997, Pi-SiO2 simpl.) Model. VESIcal has performed calculations assuming SiO2=40wt%. "
crmsg_GreaterThan_fail_Allison="{param_name} ({param_val:.1f} {units}) is less than the lowest P experiment in the calibration dataset ({calib_val:.1f} {units}). "


def crf_LessThan(calibval,paramval):
	return paramval <= calibval
crmsg_LessThan_pass = "The {param_name} ({param_val:.1f} {units}) is less than {calib_val:.1f} {units} as required by the calibration range of the {model_name} model. "
crmsg_LessThan_fail =  "The {param_name} is outside the calibration range of {model_name} ({param_val:.1f}>{calib_val:.1f} {units}. "
crmsg_LessThan_description = "The {model_name} model is calibrated for {param_name} less than {calib_val:.1f} {units}. "
# Warnings for specific model
crmsg_LessThan_fail_ShWaterSi = "{param_name} ({param_val:.1f} {units}) exceeds the upper limit of {calib_val:.1f} {units} suggested by Shishkina et al. for their H2O model. " # This check warns users if SiO2>65 wt% , the upper limit suggested by Shishkina et al. (2014)
Dix_1000bar_fail="{param_name} exceeds 1000 bar, which Iacono-Marziano et al. (2012) suggest as an upper calibration limit of the Dixon (1997, Pi-SiO2 simpl.) Model, "
Dix_2000bar_fail="as well as the upper calibration limit of 2000 bar suggested by Lesne et al. (2011), "
Dix_5000bar_fail="and the upper calibration limit of 5000 bar suggested by Newman and Lowenstern, (2002). "
Dix_49_fail="{param_name} ({param_val:.1f} {units})>49 wt%, which is the calibration limit of the Dixon (1997, Pi-SiO2 simpl.) Model. VESIcal has performed calculations assuming SiO2=49wt% for this sample. "
crmsg_LessThan_fail_Allison="{param_name} ({param_val:.1f} {units}) exceeds the upper limit of {calib_val:.1f} {units} suggested by Allison et al. (2019) - Their spreadsheet would return 7000 bar for this input. "
crmsg_LessThan_fail_AllisonH2O="{param_name} ({param_val:.1f} {units}) is > {param_val:.1f} {units}: this model does not account for the effect of H$_2$O on volatile solubility. VESIcal allows you to combine Allison Carbon with a variety of H$_2$O models.  "

def crf_Between(calibval,paramval):
	return paramval >= calibval[0] and paramval <= calibval[1]
crmsg_Between_pass = "The {param_name} ({param_val:.1f} {units}) is between {calib_val0:.1f} and {calib_val1:.1f} {units} as required by the calibration range of the {model_name} model. "
crmsg_Between_fail = "{param_name} ({param_val:.1f} {units}) is outside the calibration range of the {model_name} model ({calib_val0:.1f}-{calib_val1:.1f} {units}). "
crmsg_Between_description = "The {model_name} model is calibrated for {param_name} between {calib_val0:.1f} and {calib_val1:.1f} {units}. "
# Different wording for compositional checks (loosing "the")
crmsg_BC_pass = "{param_name} ({param_val:.1f} {units}) is between {calib_val0:.1f} and {calib_val1:.1f} {units} as required by the calibration range of the {model_name} model. "
crmsg_BC_fail = "{param_name} ({param_val:.1f} {units}) is outside the calibration range ({calib_val0:.1f}-{calib_val1:.1f} {units}). "
# Warnings for specific model
crmsg_BC_fail_ShT1 = "{param_name} ({param_val:.1f} {units}) is outside the calibration range of {model_name} ({calib_val0:.1f}-{calib_val1:.1f} {units}). Note, the authors recomend that this model is optimally calibrated between 1150-1250C. " # Warning for Shishkina temperature
crmsg_BC_fail_ShC1="{param_name} ({param_val:.1f} {units}) is outside the calibration range of {model_name} (calculated from the max and minimum concentrations in the calibration dataset +-5%; {calib_val0:.1f}-{calib_val1:.1f} {units})." # Warning for Shishkina SiO2
crmsg_BC_DixonT="{param_name} ({param_val:.1f} {units}) lies more than 200°C away from the temperature the Dixon (1997, Pi-SiO2 simpl.) model was calibrated for (the range suggested by Newman and Lowenstern, 2002; VolatileCalc). " # Warning for Dixon temperature
crmsg_BC_IMT="{param_name} ({param_val:.1f} {units}) is outside the broad range suggested by Iacono-Marziano ({calib_val0:.1f}-{calib_val1:.1f} {units}, although they note that this model is best calibrated at 1200-1300C). " # Warning for Iacono-Marziano temperature
crmsg_BC_IMP="{param_name} ({param_val:.1f} {units}) is outside the broad range suggested by Iacono-Marziano ({calib_val0:.1f}-{calib_val1:.1f} {units}, although most  calibration experiments were conducted at <5000 bars).  " # Warning for Iacono-Marziano pressure
crmsg_Between_AllisonTemp="{param_name} ({param_val:.1f} {units} is outside the recomended temperature range for {model_name} (1000-1400°C). "

# Defining compositional ranges for Liu - based on the Max value of the calibration dataset +-5% of that value
LiuWater_SiO2Max=82
LiuWater_SiO2Min=71
LiuWater_TiO2Max=0.21
LiuWater_TiO2Min=0
LiuWater_Al2O3Max=14.2
LiuWater_Al2O3Min=11.5
LiuWater_FeOMax=1.5
LiuWater_FeOMin=0
LiuWater_MgOMax=0.18
LiuWater_MgOMin=0
LiuWater_CaOMax=1.2
LiuWater_CaOMin=0
LiuWater_Na2OMax=4.9
LiuWater_Na2OMin=3.2
LiuWater_K2OMax=6.0
LiuWater_K2OMin=3.4
def crf_LiuWaterComp(calibval=None,sample={}):
	SiTest = sample['SiO2'] >= LiuWater_SiO2Max and sample['SiO2'] <= LiuWater_SiO2Min
	TiTest = sample['TiO2'] >= LiuWater_TiO2Max and sample['TiO2'] <= LiuWater_TiO2Min
	AlTest = sample['Al2O3'] >= LiuWater_Al2O3Min and sample['Al2O3'] <= LiuWater_Al2O3Max
	FeTest = sample['FeO'] >= LiuWater_FeOMin and sample['FeO'] <= LiuWater_FeOMax
	MgTest = sample['MgO'] >= LiuWater_MgOMin and sample['MgO'] <= LiuWater_MgOMax
	CaTest = sample['CaO'] >= LiuWater_CaOMin and sample['CaO'] <= LiuWater_CaOMax
	NaTest = sample['Na2O'] >= LiuWater_Na2OMin and sample['Na2O'] <= LiuWater_Na2OMax
	KTest = sample['K2O'] >= LiuWater_K2OMin and sample['K2O'] <= LiuWater_K2OMax
	return all([SiTest, TiTest, AlTest, FeTest, MgTest, CaTest, NaTest, KTest])

LiuCarbon_SiO2Max=82
LiuCarbon_SiO2Min=73
LiuCarbon_TiO2Max=0.12
LiuCarbon_TiO2Min=0.07
LiuCarbon_Al2O3Max=13.7
LiuCarbon_Al2O3Min=11.9
LiuCarbon_FeOMax=1.1
LiuCarbon_FeOMin=0.36
LiuCarbon_MgOMax=0.08
LiuCarbon_MgOMin=0.05
LiuCarbon_CaOMax=0.6
LiuCarbon_CaOMin=0.23
LiuCarbon_Na2OMax=4.4
LiuCarbon_Na2OMin=3.9
LiuCarbon_K2OMax=5.0
LiuCarbon_K2OMin=4.0

def crf_LiuCarbonComp(calibval=None,sample={}):
	SiTest = sample['SiO2'] >= LiuCarbon_SiO2Min and sample['SiO2'] <= LiuCarbon_SiO2Max
	TiTest = sample['TiO2'] >= LiuCarbon_TiO2Min and sample['TiO2'] <= LiuCarbon_TiO2Max
	AlTest = sample['Al2O3'] >= LiuCarbon_Al2O3Min and sample['Al2O3'] <= LiuCarbon_Al2O3Max
	FeTest = sample['FeO'] >= LiuCarbon_FeOMin and sample['FeO'] <= LiuCarbon_FeOMax
	MgTest = sample['MgO'] >= LiuCarbon_MgOMin and sample['MgO'] <= LiuCarbon_MgOMax
	CaTest = sample['CaO'] >= LiuCarbon_CaOMin and sample['CaO'] <= LiuCarbon_CaOMax
	NaTest = sample['Na2O'] >= LiuCarbon_Na2OMin and sample['Na2O'] <= LiuCarbon_Na2OMax
	KTest = sample['K2O'] >= LiuCarbon_K2OMin and sample['K2O'] <= LiuCarbon_K2OMax
	return all([SiTest, TiTest, AlTest, FeTest, MgTest, CaTest, NaTest, KTest])

#Max and min of combined dataset, as neither dataset is entirely complete due to publications which are hard to track down.
Liu_SiO2Min=min([LiuWater_SiO2Min, LiuCarbon_SiO2Min])
Liu_SiO2Max=max([LiuWater_SiO2Max, LiuCarbon_SiO2Max])
Liu_TiO2Min=min([LiuWater_TiO2Min, LiuCarbon_TiO2Min])
Liu_TiO2Max=max([LiuWater_TiO2Max, LiuCarbon_TiO2Max])
Liu_Al2O3Min=min([LiuWater_Al2O3Min, LiuCarbon_Al2O3Min])
Liu_Al2O3Max=max([LiuWater_Al2O3Max, LiuCarbon_Al2O3Max])
Liu_FeOMin=min([LiuWater_FeOMin, LiuCarbon_FeOMin])
Liu_FeOMax=max([LiuWater_FeOMax, LiuCarbon_FeOMax])
Liu_MgOMin=min([LiuWater_MgOMin, LiuCarbon_MgOMin])
Liu_MgOMax=max([LiuWater_MgOMax, LiuCarbon_MgOMax])
Liu_CaOMin=min([LiuWater_CaOMin, LiuCarbon_CaOMin])
Liu_CaOMax=max([LiuWater_CaOMax, LiuCarbon_CaOMax])
Liu_Na2OMin=min([LiuWater_Na2OMin, LiuCarbon_Na2OMin])
Liu_Na2OMax=max([LiuWater_Na2OMax, LiuCarbon_Na2OMax])
Liu_K2OMin=min([LiuWater_K2OMin, LiuCarbon_K2OMin])
Liu_K2OMax=max([LiuWater_K2OMax, LiuCarbon_K2OMax])


def crf_LiuComp(calibval=None,sample={}):
	SiTest = sample['SiO2'] >= Liu_SiO2Min and sample['SiO2'] <= Liu_SiO2Max
	TiTest = sample['TiO2'] >= Liu_TiO2Min and sample['TiO2'] <= Liu_TiO2Max
	AlTest = sample['Al2O3'] >= Liu_Al2O3Min and sample['Al2O3'] <= Liu_Al2O3Max
	FeTest = sample['FeO'] >= Liu_FeOMin and sample['FeO'] <= Liu_FeOMax
	MgTest = sample['MgO'] >= Liu_MgOMin and sample['MgO'] <= Liu_MgOMax
	CaTest = sample['CaO'] >= Liu_CaOMin and sample['CaO'] <= Liu_CaOMax
	NaTest = sample['Na2O'] >= Liu_Na2OMin and sample['Na2O'] <= Liu_Na2OMax
	KTest = sample['K2O'] >= Liu_K2OMin and sample['K2O'] <= Liu_K2OMax
	return all([SiTest, TiTest, AlTest, FeTest, MgTest, CaTest, NaTest, KTest])

crmsg_LiuComp_pass = "The sample appears to be similar in composition to the rhyolites and haplogranites used to calibrate the Liu et al. model."
crmsg_LiuWaterComp_fail = " These calibration limits were selected based on the minimum and maximum values of these oxides (+-5%) in the Water calibration dataset. As the Liu et al. model incorperates no term for compositional dependence, users must take extreme care when extrapolating this model to compositions which differ significantly from the haplogranites and rhyolites in the calibration dataset. These warnings are simply a guide; we suggest that users carefully compare their major element data to the calibration dataset to check for suitability "
crmsg_LiuCarbonComp_fail = " These calibration limits were selected based on the minimum and maximum values of these oxides (+-5%) in the Carbon calibration dataset. As the Liu et al. model incorperates no term for compositional dependence, users must take extreme care when extrapolating this model to compositions which differ significantly from the haplogranites and rhyolites in the calibration dataset. These warnings are simply a guide; we suggest that users carefully compare their major element data to the calibration dataset to check for suitability "
crmsg_LiuComp_fail = " These calibration limits were selected based on the minimum and maximum values of these oxides (+-5%) in the combined Water and Carbon calibration dataset. As the Liu et al. model incorperates no term for compositional dependence, users must take extreme care when extrapolating this model to compositions which differ significantly from the haplogranites and rhyolites in the calibration dataset. These warnings are simply a guide; we suggest that users carefully compare their major element data to the calibration dataset to check for suitability "
crmsg_LiuComp_description = "The Liu et al. model is suitable for haplogranites and rhyolites."

# ALLISON COMPOSITIONAL LIMITS -DEFINED AS MIN AND MAX OF CALIBRATION DATASET (-5% AND +5% RESPECTIVELY)
Allison_SiO2_sfvf=[50.0, 55.97]
Allison_TiO2_sfvf=[1.08, 1.25]
Allison_Al2O3_sfvf=[15.76, 18.15]
Allison_FeO_sfvf=[7.07, 8.06]
Allison_MgO_sfvf=[5.89, 7.09]
Allison_CaO_sfvf=[8.80, 9.91]
Allison_Na2O_sfvf=[3.05, 3.53]
Allison_K2O_sfvf=[1.25, 1.49]

def crf_AllisonComp_sfvf(calibval=None,sample={}):
	SiTest = sample['SiO2'] >= Allison_SiO2_sfvf[0] and sample['SiO2'] <= Allison_SiO2_sfvf[1]
	TiTest = sample['TiO2'] >= Allison_TiO2_sfvf[0] and sample['TiO2'] <= Allison_TiO2_sfvf[1]
	AlTest = sample['Al2O3'] >= Allison_Al2O3_sfvf[0] and sample['Al2O3'] <= Allison_Al2O3_sfvf[1]
	FeTest = sample['FeO'] >= Allison_FeO_sfvf[0] and sample['FeO'] <= Allison_FeO_sfvf[1]
	MgTest = sample['MgO'] >= Allison_MgO_sfvf[0] and sample['MgO'] <= Allison_MgO_sfvf[1]
	CaTest = sample['CaO'] >= Allison_CaO_sfvf[0] and sample['CaO'] <= Allison_CaO_sfvf[1]
	NaTest = sample['Na2O'] >= Allison_Na2O_sfvf[0] and sample['Na2O'] <= Allison_Na2O_sfvf[1]
	KTest = sample['K2O'] >= Allison_K2O_sfvf[0] and sample['K2O'] <= Allison_K2O_sfvf[1]
	return all([SiTest, TiTest, AlTest, FeTest, MgTest, CaTest, NaTest, KTest])


Allison_SiO2_sunset=[45.72, 50.62]
Allison_TiO2_sunset=[1.75, 1.96]
Allison_Al2O3_sunset=[15.62, 17.45]
Allison_FeO_sunset=[9.13, 10.42]
Allison_MgO_sunset=[8.13, 9.19]
Allison_CaO_sunset=[9.56, 10.59]
Allison_Na2O_sunset=[3.29, 3.64]
Allison_K2O_sunset=[0.77, 0.86]

def crf_AllisonComp_sunset(calibval=None,sample={}):
	SiTest = sample['SiO2'] >= Allison_SiO2_sunset[0] and sample['SiO2'] <= Allison_SiO2_sunset[1]
	TiTest = sample['TiO2'] >= Allison_TiO2_sunset[0] and sample['TiO2'] <= Allison_TiO2_sunset[1]
	AlTest = sample['Al2O3'] >= Allison_Al2O3_sunset[0] and sample['Al2O3'] <= Allison_Al2O3_sunset[1]
	FeTest = sample['FeO'] >= Allison_FeO_sunset[0] and sample['FeO'] <= Allison_FeO_sunset[1]
	MgTest = sample['MgO'] >= Allison_MgO_sunset[0] and sample['MgO'] <= Allison_MgO_sunset[1]
	CaTest = sample['CaO'] >= Allison_CaO_sunset[0] and sample['CaO'] <= Allison_CaO_sunset[1]
	NaTest = sample['Na2O'] >= Allison_Na2O_sunset[0] and sample['Na2O'] <= Allison_Na2O_sunset[1]
	KTest = sample['K2O'] >= Allison_K2O_sunset[0] and sample['K2O'] <= Allison_K2O_sunset[1]
	return all([SiTest, TiTest, AlTest, FeTest, MgTest, CaTest, NaTest, KTest])

Allison_SiO2_erebus=[45.96, 51.04]
Allison_TiO2_erebus=[2.67, 3.00]
Allison_Al2O3_erebus=[18.31, 20.63]
Allison_FeO_erebus=[7.46, 9.37]
Allison_MgO_erebus=[3.03, 3.43]
Allison_CaO_erebus=[6.58, 7.42]
Allison_Na2O_erebus=[5.8, 6.49]
Allison_K2O_erebus=[2.75, 3.13]

def crf_AllisonComp_erebus(calibval=None,sample={}):
	SiTest = sample['SiO2'] >= Allison_SiO2_erebus[0] and sample['SiO2'] <= Allison_SiO2_erebus[1]
	TiTest = sample['TiO2'] >= Allison_TiO2_erebus[0] and sample['TiO2'] <= Allison_TiO2_erebus[1]
	AlTest = sample['Al2O3'] >= Allison_Al2O3_erebus[0] and sample['Al2O3'] <= Allison_Al2O3_erebus[1]
	FeTest = sample['FeO'] >= Allison_FeO_erebus[0] and sample['FeO'] <= Allison_FeO_erebus[1]
	MgTest = sample['MgO'] >= Allison_MgO_erebus[0] and sample['MgO'] <= Allison_MgO_erebus[1]
	CaTest = sample['CaO'] >= Allison_CaO_erebus[0] and sample['CaO'] <= Allison_CaO_erebus[1]
	NaTest = sample['Na2O'] >= Allison_Na2O_erebus[0] and sample['Na2O'] <= Allison_Na2O_erebus[1]
	KTest = sample['K2O'] >= Allison_K2O_erebus[0] and sample['K2O'] <= Allison_K2O_erebus[1]
	return all([SiTest, TiTest, AlTest, FeTest, MgTest, CaTest, NaTest, KTest])

Allison_SiO2_vesuvius=[46.65, 53.29]
Allison_TiO2_vesuvius=[0.93, 1.13]
Allison_Al2O3_vesuvius=[13.75, 16.28]
Allison_FeO_vesuvius=[4.98, 7.48]
Allison_MgO_vesuvius=[6.41, 7.76]
Allison_CaO_vesuvius=[11.12, 14.16]
Allison_Na2O_vesuvius=[1.74, 2.07]
Allison_K2O_vesuvius=[5.48, 6.35]

def crf_AllisonComp_vesuvius(calibval=None,sample={}):
	SiTest = sample['SiO2'] >= Allison_SiO2_vesuvius[0] and sample['SiO2'] <= Allison_SiO2_vesuvius[1]
	TiTest = sample['TiO2'] >= Allison_TiO2_vesuvius[0] and sample['TiO2'] <= Allison_TiO2_vesuvius[1]
	AlTest = sample['Al2O3'] >= Allison_Al2O3_vesuvius[0] and sample['Al2O3'] <= Allison_Al2O3_vesuvius[1]
	FeTest = sample['FeO'] >= Allison_FeO_vesuvius[0] and sample['FeO'] <= Allison_FeO_vesuvius[1]
	MgTest = sample['MgO'] >= Allison_MgO_vesuvius[0] and sample['MgO'] <= Allison_MgO_vesuvius[1]
	CaTest = sample['CaO'] >= Allison_CaO_vesuvius[0] and sample['CaO'] <= Allison_CaO_vesuvius[1]
	NaTest = sample['Na2O'] >= Allison_Na2O_vesuvius[0] and sample['Na2O'] <= Allison_Na2O_vesuvius[1]
	KTest = sample['K2O'] >= Allison_K2O_vesuvius[0] and sample['K2O'] <= Allison_K2O_vesuvius[1]
	return all([SiTest, TiTest, AlTest, FeTest, MgTest, CaTest, NaTest, KTest])


Allison_SiO2_etna=[46.03, 52.97]
Allison_TiO2_etna=[1.61, 1.89]
Allison_Al2O3_etna=[15.87, 18.24]
Allison_FeO_etna=[6.75, 10.21]
Allison_MgO_etna=[5.9, 7]
Allison_CaO_etna=[9.38, 11.99]
Allison_Na2O_etna=[3.4, 3.94]
Allison_K2O_etna=[1.69, 2.25]

def crf_AllisonComp_etna(calibval=None,sample={}):
	SiTest = sample['SiO2'] >= Allison_SiO2_etna[0] and sample['SiO2'] <= Allison_SiO2_etna[1]
	TiTest = sample['TiO2'] >= Allison_TiO2_etna[0] and sample['TiO2'] <= Allison_TiO2_etna[1]
	AlTest = sample['Al2O3'] >= Allison_Al2O3_etna[0] and sample['Al2O3'] <= Allison_Al2O3_etna[1]
	FeTest = sample['FeO'] >= Allison_FeO_etna[0] and sample['FeO'] <= Allison_FeO_etna[1]
	MgTest = sample['MgO'] >= Allison_MgO_etna[0] and sample['MgO'] <= Allison_MgO_etna[1]
	CaTest = sample['CaO'] >= Allison_CaO_etna[0] and sample['CaO'] <= Allison_CaO_etna[1]
	NaTest = sample['Na2O'] >= Allison_Na2O_etna[0] and sample['Na2O'] <= Allison_Na2O_etna[1]
	KTest = sample['K2O'] >= Allison_K2O_etna[0] and sample['K2O'] <= Allison_K2O_etna[1]
	return all([SiTest, TiTest, AlTest, FeTest, MgTest, CaTest, NaTest, KTest])

Allison_SiO2_stromboli=[47.23, 55.15]
Allison_TiO2_stromboli=[0.74, 0.94]
Allison_Al2O3_stromboli=[14.87, 17.73]
Allison_FeO_stromboli=[5.08, 7.51]
Allison_MgO_stromboli=[7.52, 9.26]
Allison_CaO_stromboli=[11.99, 13.46]
Allison_Na2O_stromboli=[2.29, 2.67]
Allison_K2O_stromboli=[1.79, 2.17]

def crf_AllisonComp_stromboli(calibval=None,sample={}):
	SiTest = sample['SiO2'] >= Allison_SiO2_stromboli[0] and sample['SiO2'] <= Allison_SiO2_stromboli[1]
	TiTest = sample['TiO2'] >= Allison_TiO2_stromboli[0] and sample['TiO2'] <= Allison_TiO2_stromboli[1]
	AlTest = sample['Al2O3'] >= Allison_Al2O3_stromboli[0] and sample['Al2O3'] <= Allison_Al2O3_stromboli[1]
	FeTest = sample['FeO'] >= Allison_FeO_stromboli[0] and sample['FeO'] <= Allison_FeO_stromboli[1]
	MgTest = sample['MgO'] >= Allison_MgO_stromboli[0] and sample['MgO'] <= Allison_MgO_stromboli[1]
	CaTest = sample['CaO'] >= Allison_CaO_stromboli[0] and sample['CaO'] <= Allison_CaO_stromboli[1]
	NaTest = sample['Na2O'] >= Allison_Na2O_stromboli[0] and sample['Na2O'] <= Allison_Na2O_stromboli[1]
	KTest = sample['K2O'] >= Allison_K2O_stromboli[0] and sample['K2O'] <= Allison_K2O_stromboli[1]
	return all([SiTest, TiTest, AlTest, FeTest, MgTest, CaTest, NaTest, KTest])




crmsg_AllisonComp_pass = "The sample appears to be similar in composition to the compositional dataset for the selected Carbon model of Allison et al. (2019)."
crmsg_AllisonComp_fail = " These calibration limits were selected based on the minimum and maximum values of these oxides (+-5%) in the calibration dataset. As the Allison et al. model incorperates no term for compositional dependence, users must take extreme care when extrapolating this model to compositions which differ significantly from the calibration dataset. These warnings are simply a guide; we suggest that users carefully compare their major element data to the calibration dataset to check for suitability "
crmsg_AllisonComp_description = "The Allison et al. (2019) Carbon model is defined for 6 different alkali compositions."


#-------------FUGACITY MODELS--------------------------------#

class fugacity_idealgas(FugacityModel):
	""" An instance of FugacityModel for an ideal gas.
	"""

	def fugacity(self,pressure,X_fluid=1.0,**kwargs):
		""" Returns the fugacity of an ideal gas, i.e., the partial pressure.

		Parameters
		----------
		pressure    float
			Total pressure of the system, in bars.
		X_fluid     float
			The mole fraction of the species in the vapour phase.

		Returns
		-------
		float
			Fugacity (partial pressure) in bars
		"""
		return pressure*X_fluid



class fugacity_KJ81_co2(FugacityModel):
	""" Implementation of the Kerrick and Jacobs (1981) EOS for mixed fluids. This class
	will return the properties of the CO2 component of the mixed fluid.
	"""
	def __init__(self):
		self.set_calibration_ranges([CalibrationRange('pressure',20000.0,crf_LessThan,'bar','Kerrick and Jacobs (1981) EOS',
													  fail_msg=crmsg_LessThan_fail, pass_msg=crmsg_LessThan_pass, description_msg=crmsg_LessThan_description),
									 CalibrationRange('temperature',1050,crf_LessThan,'oC','Kerrick and Jacobs (1981) EOS',
									 				  fail_msg=crmsg_LessThan_fail, pass_msg=crmsg_LessThan_pass, description_msg=crmsg_LessThan_description)])

	def fugacity(self,pressure,temperature,X_fluid,**kwargs):
		""" Calculates the fugacity of CO2 in a mixed CO2-H2O fluid. Above 1050C,
		it assumes H2O and CO2 do not interact, as the equations are not defined
		beyond this point.

		Parameters
		----------
		pressure    float
			Total pressure of the system in bars.
		temperature     float
			Temperature in degC
		X_fluid     float
			Mole fraction of CO2 in the fluid.

		Returns
		-------
		float
			fugacity of CO2 in bars
		"""
		if X_fluid == 0:
			return 0
		elif temperature >= 1050.0:
			return pressure*np.exp(self.lnPhi_mix(pressure,temperature,1.0))*X_fluid
		else:
			return pressure*np.exp(self.lnPhi_mix(pressure,temperature,X_fluid))*X_fluid


	def volume(self,P,T,X_fluid):
		""" Calculates the volume of the mixed fluid, by solving Eq (28) of Kerrick and
		Jacobs (1981) using scipy.root_scalar.

		Parameters
		----------
		P   float
			Total pressure of the system, in bars.
		T   float
			Temperature in degC
		X_fluid     float
			Mole fraction of CO2 in the fluid

		Returns
		-------
		float
			Volume of the mixed fluid.
		"""
		if X_fluid != 1.0:
			# x0 = self.volume(P,T,1.0)*X_fluid + self.volume_h(P,T)*(1-X_fluid)
			# print(x0)
			if P >= 20000 and T<800-273.15:
				x0 = (X_fluid*25+(1-X_fluid)*15)
			else:
				x0 = (X_fluid*35+(1-X_fluid)*15)

		else:
			if P >= 20000 and T<800-273.15:
				x0 = 25
			else:
				x0=35
		return root_scalar(self.root_volume,x0=x0,x1=x0*0.9,args=(P,T,X_fluid)).root


	def root_volume(self,v,P,T,X_fluid):
		""" Returns the difference between the lhs and rhs of Eq (28) of Kerrick and Jacobs (1981).
		For use with a root finder to obtain the volume of the mixed fluid.

		Parameters
		----------
		v   float
			Guess for the volume
		P   float
			Total system pressure in bars.
		T   float
			Temperature in degC
		X_fluid     float
			Mole fraction of CO2 in the fluid.

		Returns
		-------
		float
			Difference between lhs and rhs of Eq (28) of Kerrick and Jacobs (1981), in bars.
		"""
		T = T + 273.15
		c = {}
		h = {}

		c['b'] = 58.0
		c['c'] = (28.31 + 0.10721*T - 8.81e-6*T**2)*1e6
		c['d'] = (9380.0 - 8.53*T + 1.189e-3*T**2)*1e6
		c['e'] = (-368654.0 + 715.9*T + 0.1534*T**2)*1e6
		h['b'] = 29.0
		h['c'] = (290.78 - 0.30276*T + 1.4774e-4*T**2)*1e6#3
		h['d'] = (-8374.0 + 19.437*T - 8.148e-3*T**2)*1e6
		h['e'] = (76600.0 - 133.9*T + 0.1071*T**2)*1e6

		if X_fluid == 1:
			bm = c['b']
			cm = c['c']
			c12= c['c']
			dm = c['d']
			d12= c['d']
			em = c['e']
			e12 =c['e']
		else:
			bm = X_fluid*c['b'] + (1-X_fluid)*h['b']
			c12 = (c['c']*h['c'])**0.5
			cm = c['c']*X_fluid**2 + h['c']*(1-X_fluid)**2 + 2*X_fluid*(1-X_fluid)*c12
			d12 = (c['d']*h['d'])**0.5
			dm = c['d']*X_fluid**2 + h['d']*(1-X_fluid)**2 + 2*X_fluid*(1-X_fluid)*d12
			e12 = (c['e']*h['e'])**0.5
			em = c['e']*X_fluid**2 + h['e']*(1-X_fluid)**2 + 2*X_fluid*(1-X_fluid)*e12

		am = cm + dm/v + em/v**2

		y = bm/(4*v)

		pt1 = (83.14 * T * (1 + y + y**2 - y**3)) / (v*(1-y)**3)
		pt2 = - am / (T**0.5 * v * (v+bm))

		return -(P - pt1 - pt2)

	def volume_h(self,P,T):
		""" Calculates the volume of a pure H2O fluid, by solving Eq (14) of
		Kerrick and Jacobs (1981).

		Parameters
		----------
		P   float
			Total pressure in bars.
		T   float
			Temperature in degC.

		Returns
		-------
		Difference between lhs and rhs of Eq (14) of Kerrick and Jacobs (1981), in bars.
		"""
		return root_scalar(self.root_volume_h,x0=15,x1=35,args=(P,T)).root


	def root_volume_h(self,v,P,T):
		""" Returns the difference between the lhs and rhs of Eq (14) of
		Kerrick and Jacobs (1981). For use with a root solver to identify the
		volume of a pure H2O fluid.

		Parameters
		----------
		v   float
			Guess for the volume
		P   float
			Total pressure in bars.
		T   float
			Temperature in degC.

		Returns
		-------
		float
			The difference between the lhs and rhs of Eq (14) of Kerrick and Jacobs (1981),
			in bars.
		"""
		T = T + 273.15
		h = {}
		h['b'] = 29.0
		h['c'] = (290.78 - 0.30276*T + 1.4774e-4*T**2)*1e6#3
		h['d'] = (-8374.0 + 19.437*T - 8.148e-3*T**2)*1e6
		h['e'] = (76600.0 - 133.9*T + 0.1071*T**2)*1e6
		h['a'] = h['c'] + h['d']/v + h['e']/v**2

		y = h['b']/(4*v)

		pt1 = (83.14 * T * (1 + y + y**2 - y**3)) / (v*(1-y)**3)
		pt2 = - h['a'] / (T**0.5 * v * (v+h['b']))

		return -(P - pt1 - pt2)


	def lnPhi_mix(self,P,T,X_fluid):
		""" Calculates the natural log of the fugacity coefficient for CO2 in a
		mixed CO2-H2O fluid. Uses Eq (27) of Kerrick and Jacobs (1981).

		Parameters
		----------
		P   float
			Total pressure in bars.
		T   float
			Temperature in degC
		X_fluid     float
			The mole fraction of CO2 in the fluid.

		Returns
		-------
		float
			The natural log of the fugacity coefficient for CO2 in a mixed fluid.
		"""
		T = T + 273.15
		v = self.volume(P,T-273.15,X_fluid)

		c = {}
		h = {}

		c['b'] = 58.0
		c['c'] = (28.31 + 0.10721*T - 8.81e-6*T**2)*1e6
		c['d'] = (9380.0 - 8.53*T + 1.189e-3*T**2)*1e6
		c['e'] = (-368654.0 + 715.9*T + 0.1534*T**2)*1e6
		h['b'] = 29.0
		h['c'] = (290.78 - 0.30276*T + 1.4774e-4*T**2)*1e6#3
		h['d'] = (-8374.0 + 19.437*T - 8.148e-3*T**2)*1e6
		h['e'] = (76600.0 - 133.9*T + 0.1071*T**2)*1e6

		if X_fluid == 1:
			bm = c['b']
			cm = c['c']
			c12= c['c']
			dm = c['d']
			d12= c['d']
			em = c['e']
			e12 =c['e']
		else:
			bm = X_fluid*c['b'] + (1-X_fluid)*h['b']
			c12 = (c['c']*h['c'])**0.5
			cm = c['c']*X_fluid**2 + h['c']*(1-X_fluid)**2 + 2*X_fluid*(1-X_fluid)*c12
			d12 = (c['d']*h['d'])**0.5
			dm = c['d']*X_fluid**2 + h['d']*(1-X_fluid)**2 + 2*X_fluid*(1-X_fluid)*d12
			e12 = (c['e']*h['e'])**0.5
			em = c['e']*X_fluid**2 + h['e']*(1-X_fluid)**2 + 2*X_fluid*(1-X_fluid)*e12

		am = cm + dm/v + em/v**2

		y = bm/(4*v)

		# Z = (1+y+y**2-y**3)/(1-y)**2 - am/(83.14*T**1.5*(v+bm))
		Z = v*P/(83.14*T)

		lnPhi = 0

		lnPhi += (4*y-3*y**2)/(1-y)**2 + (c['b']/bm * (4*y-2*y**2)/(1-y)**3)
		lnPhi += - (2*c['c']*X_fluid+2*(1-X_fluid)*c12)/(83.14*T**1.5*bm)*np.log((v+bm)/v)
		lnPhi += - cm*c['b']/(83.14*T**1.5*bm*(v+bm))
		lnPhi += cm*c['b']/(83.14*T**1.5*bm**2)*np.log((v+bm)/v)
		lnPhi += - (2*c['d']*X_fluid+2*d12*(1-X_fluid)+dm)/(83.14*T**1.5*bm*v)
		lnPhi += (2*c['d']*X_fluid+2*(1-X_fluid)*d12+dm)/(83.14*T**1.5*bm**2)*np.log((v+bm)/v)
		lnPhi += c['b']*dm/(83.14*T**1.5*v*bm*(v+bm)) + 2*c['b']*dm/(83.14*T**1.5*bm**2*(v+bm))
		lnPhi += - 2*c['b']*dm/(83.14*T**1.5*bm**3)*np.log((v+bm)/v)
		lnPhi += - (2*c['e']*X_fluid + 2*(1-X_fluid)*e12+2*em)/(83.14*T**1.5*2*bm*v**2)
		lnPhi += (2*c['e']*X_fluid+2*e12*(1-X_fluid)+2*em)/(83.14*T**1.5*bm**2*v)
		lnPhi += - (2*c['e']*X_fluid+2*e12*(1-X_fluid)+2*em)/(83.14*T**1.5*bm**3)*np.log((v+bm)/v)
		lnPhi += em*c['b']/(83.14*T**1.5*2*bm*v**2*(v+bm)) - 3*em*c['b']/(83.14*T**1.5*2*bm**2*v*(v+bm))
		lnPhi += 3*em*c['b']/(83.14*T**1.5*bm**4)*np.log((v+bm)/v) - 3*em*c['b']/(83.14*T**1.5*bm**3*(v+bm))
		lnPhi += - np.log(Z)

		return lnPhi



class fugacity_KJ81_h2o(FugacityModel):
	"""Implementation of the Kerrick and Jacobs (1981) EOS for mixed fluids. This class
	will return the properties of the H2O component of the mixed fluid.
	"""
	def __init__(self):
		self.set_calibration_ranges([CalibrationRange('pressure',20000.0,crf_LessThan,'bar','Kerrick and Jacobs (1981) EOS',
													  fail_msg=crmsg_LessThan_fail, pass_msg=crmsg_LessThan_pass, description_msg=crmsg_LessThan_description),
									 CalibrationRange('temperature',1050,crf_LessThan,'oC','Kerrick and Jacobs (1981) EOS',
									 				  fail_msg=crmsg_LessThan_fail, pass_msg=crmsg_LessThan_pass, description_msg=crmsg_LessThan_description)])

	def fugacity(self,pressure,temperature,X_fluid,**kwargs):
		""" Calculates the fugacity of H2O in a mixed CO2-H2O fluid. Above 1050C,
		it assumes H2O and CO2 do not interact, as the equations are not defined
		beyond this point.

		Parameters
		----------
		pressure    float
			Total pressure of the system in bars.
		temperature     float
			Temperature in degC
		X_fluid     float
			Mole fraction of H2O in the fluid.

		Returns
		-------
		float
			fugacity of H2O in bars
		"""
		if X_fluid == 0:
			return 0
		elif temperature >= 1050:
			return pressure*np.exp(self.lnPhi_mix(pressure,temperature,1.0))*X_fluid
		else:
			return pressure*np.exp(self.lnPhi_mix(pressure,temperature,X_fluid))*X_fluid


	def volume(self,P,T,X_fluid):
		""" Calculates the volume of the mixed fluid, by solving Eq (28) of Kerrick and
		Jacobs (1981) using scipy.root_scalar.

		Parameters
		----------
		P   float
			Total pressure of the system, in bars.
		T   float
			Temperature in degC
		X_fluid     float
			Mole fraction of H2O in the fluid

		Returns
		-------
		float
			Volume of the mixed fluid.
		"""
		if X_fluid != 1.0:
			# x0 = self.volume(P,T,1.0)*X_fluid + self.volume_h(P,T)*(1-X_fluid)
			# print(x0)
			if P >= 20000 and T<800-273.15:
				x0 = ((1-X_fluid)*25+X_fluid*15)
			else:
				x0 = ((1-X_fluid)*35+X_fluid*15)

		else:
			if P >= 20000 and T<800-273.15:
				x0 = 10
			else:
				x0=15
		return root_scalar(self.root_volume,x0=x0,x1=x0*0.9,args=(P,T,X_fluid)).root


	def root_volume(self,v,P,T,X_fluid):
		""" Returns the difference between the lhs and rhs of Eq (28) of Kerrick and Jacobs (1981).
		For use with a root finder to obtain the volume of the mixed fluid.

		Parameters
		----------
		v   float
			Guess for the volume
		P   float
			Total system pressure in bars.
		T   float
			Temperature in degC
		X_fluid     float
			Mole fraction of H2O in the fluid.

		Returns
		-------
		float
			Difference between lhs and rhs of Eq (28) of Kerrick and Jacobs (1981), in bars.
		"""
		T = T + 273.15
		c = {}
		h = {}

		c['b'] = 58.0
		c['c'] = (28.31 + 0.10721*T - 8.81e-6*T**2)*1e6
		c['d'] = (9380.0 - 8.53*T + 1.189e-3*T**2)*1e6
		c['e'] = (-368654.0 + 715.9*T + 0.1534*T**2)*1e6
		h['b'] = 29.0
		h['c'] = (290.78 - 0.30276*T + 1.4774e-4*T**2)*1e6#3
		h['d'] = (-8374.0 + 19.437*T - 8.148e-3*T**2)*1e6
		h['e'] = (76600.0 - 133.9*T + 0.1071*T**2)*1e6

		if X_fluid == 1:
			bm = h['b']
			cm = h['c']
			dm = h['d']
			em = h['e']
			c12= h['c']
			d12= h['d']
			e12= h['e']
		else:
			bm = X_fluid*h['b'] + (1-X_fluid)*c['b']
			c12 = (c['c']*h['c'])**0.5
			cm = h['c']*X_fluid**2 + c['c']*(1-X_fluid)**2 + 2*X_fluid*(1-X_fluid)*c12
			d12 = (c['d']*h['d'])**0.5
			dm = h['d']*X_fluid**2 + c['d']*(1-X_fluid)**2 + 2*X_fluid*(1-X_fluid)*d12
			e12 = (c['e']*h['e'])**0.5
			em = h['e']*X_fluid**2 + c['e']*(1-X_fluid)**2 + 2*X_fluid*(1-X_fluid)*e12

		am = cm + dm/v + em/v**2

		y = bm/(4*v)

		pt1 = (83.14 * T * (1 + y + y**2 - y**3)) / (v*(1-y)**3)
		pt2 = - am / (T**0.5 * v * (v+bm))

		return -(P - pt1 - pt2)

	def volume_c(self,P,T):
		""" Calculates the volume of a pure CO2 fluid, by solving Eq (14) of
		Kerrick and Jacobs (1981).

		Parameters
		----------
		P   float
			Total pressure in bars.
		T   float
			Temperature in degC.

		Returns
		-------
		Difference between lhs and rhs of Eq (14) of Kerrick and Jacobs (1981), in bars.
		"""
		return root_scalar(self.root_volume_c,x0=15,x1=35,args=(P,T)).root


	def root_volume_c(self,v,P,T):
		""" Returns the difference between the lhs and rhs of Eq (14) of
		Kerrick and Jacobs (1981). For use with a root solver to identify the
		volume of a pure H2O fluid.

		Parameters
		----------
		v   float
			Guess for the volume
		P   float
			Total pressure in bars.
		T   float
			Temperature in degC.

		Returns
		-------
		float
			The difference between the lhs and rhs of Eq (14) of Kerrick and Jacobs (1981),
			in bars.
		"""
		T = T + 273.15
		c = {}
		c['b'] = 58.0
		c['c'] = (28.31 + 0.10721*T - 8.81e-6*T**2)*1e6
		c['d'] = (9380.0 - 8.53*T + 1.189e-3*T**2)*1e6
		c['e'] = (-368654.0 + 715.9*T + 0.1534*T**2)*1e6
		c['a'] = c['c'] + c['d']/v + c['e']/v**2

		y = c['b']/(4*v)

		pt1 = (83.14 * T * (1 + y + y**2 - y**3)) / (v*(1-y)**3)
		pt2 = - c['a'] / (T**0.5 * v * (v+c['b']))

		return -(P - pt1 - pt2)


	def lnPhi_mix(self,P,T,X_fluid):
		""" Calculates the natural log of the fugacity coefficient for H2O in a
		mixed CO2-H2O fluid. Uses Eq (27) of Kerrick and Jacobs (1981).

		Parameters
		----------
		P   float
			Total pressure in bars.
		T   float
			Temperature in degC
		X_fluid     float
			The mole fraction of H2O in the fluid.

		Returns
		-------
		float
			The natural log of the fugacity coefficient for H2O in a mixed fluid.
		"""
		T = T + 273.15
		v = self.volume(P,T-273.15,X_fluid)

		c = {}
		h = {}

		c['b'] = 58.0
		c['c'] = (28.31 + 0.10721*T - 8.81e-6*T**2)*1e6
		c['d'] = (9380.0 - 8.53*T + 1.189e-3*T**2)*1e6
		c['e'] = (-368654.0 + 715.9*T + 0.1534*T**2)*1e6
		h['b'] = 29.0
		h['c'] = (290.78 - 0.30276*T + 1.4774e-4*T**2)*1e6#3
		h['d'] = (-8374.0 + 19.437*T - 8.148e-3*T**2)*1e6
		h['e'] = (76600.0 - 133.9*T + 0.1071*T**2)*1e6

		if X_fluid == 1:
			bm = h['b']
			cm = h['c']
			dm = h['d']
			em = h['e']
			c12= h['c']
			d12= h['d']
			e12= h['e']
		else:
			bm = X_fluid*h['b'] + (1-X_fluid)*c['b']
			c12 = (c['c']*h['c'])**0.5
			cm = h['c']*X_fluid**2 + c['c']*(1-X_fluid)**2 + 2*X_fluid*(1-X_fluid)*c12
			d12 = (c['d']*h['d'])**0.5
			dm = h['d']*X_fluid**2 + c['d']*(1-X_fluid)**2 + 2*X_fluid*(1-X_fluid)*d12
			e12 = (c['e']*h['e'])**0.5
			em = h['e']*X_fluid**2 + c['e']*(1-X_fluid)**2 + 2*X_fluid*(1-X_fluid)*e12

		am = cm + dm/v + em/v**2

		y = bm/(4*v)

		# Z = (1+y+y**2-y**3)/(1-y)**2 - am/(83.14*T**1.5*(v+bm))
		Z = v*P/(83.14*T)

		lnPhi = 0

		lnPhi += (4*y-3*y**2)/(1-y)**2 + (h['b']/bm * (4*y-2*y**2)/(1-y)**3)
		lnPhi += - (2*h['c']*X_fluid+2*(1-X_fluid)*c12)/(83.14*T**1.5*bm)*np.log((v+bm)/v)
		lnPhi += - cm*h['b']/(83.14*T**1.5*bm*(v+bm))
		lnPhi += cm*h['b']/(83.14*T**1.5*bm**2)*np.log((v+bm)/v)
		lnPhi += - (2*h['d']*X_fluid+2*d12*(1-X_fluid)+dm)/(83.14*T**1.5*bm*v)
		lnPhi += (2*h['d']*X_fluid+2*(1-X_fluid)*d12+dm)/(83.14*T**1.5*bm**2)*np.log((v+bm)/v)
		lnPhi += h['b']*dm/(83.14*T**1.5*v*bm*(v+bm)) + 2*h['b']*dm/(83.14*T**1.5*bm**2*(v+bm))
		lnPhi += - 2*h['b']*dm/(83.14*T**1.5*bm**3)*np.log((v+bm)/v)
		lnPhi += - (2*h['e']*X_fluid + 2*(1-X_fluid)*e12+2*em)/(83.14*T**1.5*2*bm*v**2)
		lnPhi += (2*h['e']*X_fluid+2*e12*(1-X_fluid)+2*em)/(83.14*T**1.5*bm**2*v)
		lnPhi += - (2*h['e']*X_fluid+2*e12*(1-X_fluid)+2*em)/(83.14*T**1.5*bm**3)*np.log((v+bm)/v)
		lnPhi += em*h['b']/(83.14*T**1.5*2*bm*v**2*(v+bm)) - 3*em*h['b']/(83.14*T**1.5*2*bm**2*v*(v+bm))
		lnPhi += 3*em*h['b']/(83.14*T**1.5*bm**4)*np.log((v+bm)/v) - 3*em*h['b']/(83.14*T**1.5*bm**3*(v+bm))
		lnPhi += - np.log(Z)

		return lnPhi




class fugacity_ZD09_co2(FugacityModel):
	""" Implementation of the Zhang and Duan (2009) fugacity model for pure CO2
	fluids."""
	def __init__(self):
		self.set_calibration_ranges([CalibrationRange('pressure',[1,1e5],crf_Between,'bar','Zhang and Duan (2009) EOS',
													  fail_msg=crmsg_Between_fail, pass_msg=crmsg_Between_pass, description_msg=crmsg_Between_description),
									 CalibrationRange('temperature',[200,2300],crf_Between,'oC','Zhang and Duan (2009) EOS',
									 				  fail_msg=crmsg_Between_fail, pass_msg=crmsg_Between_pass, description_msg=crmsg_Between_description)])

	def fugacity(self,pressure,temperature,X_fluid=1.0,**kwargs):
		""" Calculates the fugacity of a pure CO2 fluid, or a mixed fluid assuming
		ideal mixing. Implements eqn (14) of Zhang and Duan (2009).

		Paramters
		---------
		pressure     float
			Pressure in bars
		temperature     float
			Temperature in degC
		X_fluid     float
			Mole fraction of CO2 in the fluid. Default is 1.0.

		Returns
		-------
		float
			Fugacity of CO2, standard state 1 bar.
		"""

		P = pressure/10
		T = temperature + 273.15

		a = np.array([0.0,
					  2.95177298930e-2,
					  -6.33756452413e3,
					  -2.75265428882e5,
					  1.29128089283e-3,
					  -1.45797416153e2,
					  7.65938947237e4,
					  2.58661493537e-6,
					  0.52126532146,
					  -1.39839523753e2,
					  -2.36335007175e-8,
					  5.35026383543e-3,
					  -0.27110649951,
					  2.50387836486e4,
					  0.73226726041,
					  1.5483335997e-2])
		e = 235.0
		s = 3.79

		Pm = 3.0636*P*s**3/e
		Tm = 154*T/e
		Vm = root_scalar(self.Vm,x0=200,x1=100,args=(P,T)).root

		S1 = ((a[1]+a[2]/Tm**2+a[3]/Tm**3)/Vm+
			  (a[4]+a[5]/Tm**2+a[6]/Tm**3)/(2*Vm**2)+
			  (a[7]+a[8]/Tm**2+a[9]/Tm**3)/(4*Vm**4)+
			  (a[10]+a[11]/Tm**2+a[12]/Tm**3)/(5*Vm**5)+
			  (a[13]/(2*a[15]*Tm**3)*(a[14]+1-(a[14]+1+a[15]/Vm**2)*
			   np.exp(-a[15]/Vm**2)))
			 )

		Z = Pm*Vm/(8.314*Tm)

		lnfc = Z - 1 - np.log(Z) + S1

		return P*np.exp(lnfc)*10

	def Vm(self,Vm,P,T):
		""" Function to use for solving for the parameter Vm, defined by eqn (8) of
		Zhang and Duan (2009). Called by scipy.fsolve in the fugacity method.

		Parameters
		----------
		Vm     float
			Guessed value of Vm
		P     float
			Pressure in MPa
		T     float
			Temperature in K

		Returns
		-------
		float
			Difference between (rearranged) LHS and RHS of eqn (8) of Zhang and Duan (2009).
		"""
		Pm = 3.0636*P*3.79**3/235.0
		Tm = 154*T/235.0
		a = np.array([0.0,
					  2.95177298930e-2,
					  -6.33756452413e3,
					  -2.75265428882e5,
					  1.29128089283e-3,
					  -1.45797416153e2,
					  7.65938947237e4,
					  2.58661493537e-6,
					  0.52126532146,
					  -1.39839523753e2,
					  -2.36335007175e-8,
					  5.35026383543e-3,
					  -0.27110649951,
					  2.50387836486e4,
					  0.73226726041,
					  1.5483335997e-2])

		return ((1+(a[1]+a[2]/Tm**2+a[3]/Tm**3)/Vm+
				 (a[4]+a[5]/Tm**2+a[6]/Tm**3)/Vm**2+
				 (a[7]+a[8]/Tm**2+a[9]/Tm**3)/Vm**4)*0.08314*Tm/Pm - Vm
				)

class fugacity_MRK_co2(FugacityModel):
	""" Modified Redlick Kwong fugacity model as used by VolatileCalc. Python implementation by
	D. J. Rasmussen (github.com/DJRgeoscience/VolatileCalcForPython), based on VB code by Newman &
	Lowenstern.
	"""
	def __init__(self):
		self.set_calibration_ranges([])

	def fugacity(self,pressure,temperature,X_fluid=1.0,**kwargs):
		""" Calculates the fugacity of CO2 in a pure or mixed H2O-CO2 fluid (assuming ideal mixing).

		Parameters
		----------
		pressure    float
			Total pressure of the system in bars.
		temperature     float
			Temperature in degC
		X_fluid     float
			Mole fraction of CO2 in the fluid.

		Returns
		-------
		float
			fugacity of CO2 in bars
		"""
		fug = self.MRK(pressure,temperature+273.15)
		return fug*X_fluid

	def FNA(self,TK):
		return (166800000 - 193080 * (TK - 273.15) + 186.4 * (TK - 273.15)**2 - 0.071288 * ((TK - 273.15)**3)) * 1.01325

	def FNB(self,TK):
		return 1.01325 * (73030000 - 71400 * (TK - 273.15) + 21.57 * (TK - 273.15)**2)

	def FNC(self,TK):
		R = 83.14321
		return 1.01325 * (np.exp(-11.071 + 5953 / TK - 2746000 / TK**2 + 464600000 / TK**3) * 0.5 * R * R * TK**2.5 / 1.02668 + 40123800)

	def FNF(self,V,TK,A,B,P):
		R = 83.14321
		return R * TK / (V - B) - A / ((V * V + B * V) * TK**0.5) - P

	def MRK(self,P,TK): #Redlich-Kwong routine to estimate endmember H2O and CO2 fugacities
		R = 83.14321
		B_1 = 14.6
		B_2 = 29.7

		for X_1 in [0,1]:
			B = X_1 * B_1 + (1 - X_1) * B_2
			A = X_1**2 * self.FNA(TK) + 2 * X_1 * (1 - X_1) * self.FNC(TK) + (1 - X_1)**2 * self.FNB(TK)
			Temp2 = B + 5
			Q = 1
			Temp1 = 0
			while abs(Temp2 - Temp1) >= 0.00001:
				Temp1 = Temp2
				F_1 = (self.FNF(Temp1 + 0.01, TK, A, B, P) - self.FNF(Temp1, TK, A, B, P)) / 0.01
				Temp2 = Temp1 - Q * self.FNF(Temp1, TK, A, B, P) / F_1
				F_2 = (self.FNF(Temp2 + 0.01, TK, A, B, P) - self.FNF(Temp2, TK, A, B, P)) / 0.01
				if F_2 * F_1 <= 0:
					Q = Q / 2.
				if abs(Temp2 - Temp1) > 0.00001:
					F_1 = F_2
			V = Temp2
			G_1 = np.log(V / (V - B)) + B_1 / (V - B) - 2 * (X_1 * self.FNA(TK) + (1 - X_1) * self.FNC(TK)) * np.log((V + B) / V) / (R * TK**1.5 * B)
			G_1 = G_1 + (np.log((V + B) / V) - B / (V + B)) * A * B_1 / (R * TK**1.5 * B**2) - np.log(P * V / (R * TK))
			G_1 = np.exp(G_1)
			G_2 = np.log(V / (V - B)) + B_2 / (V - B) - 2 * (X_1 * self.FNC(TK) + (1 - X_1) * self.FNB(TK)) * np.log((V + B) / V) / (R * TK**1.5 * B)
			G_2 = G_2 + (np.log((V + B) / V) - B / (V + B)) * A * B_2 / (R * TK**1.5 * B**2) - np.log(P * V / (R * TK))
			G_2 = np.exp(G_2)
			if X_1 == 0:
				fCO2o = G_2 * P #The fugacity of CO2
				# return fCO2o
			if X_1 == 1:
				fH2Oo = G_1 * P #The fugacity of H2O
				# return fH2Oo
		return fCO2o

class fugacity_MRK_h2o(FugacityModel):
	""" Modified Redlick Kwong fugacity model as used by VolatileCalc. Python implementation by
	D. J. Rasmussen (github.com/DJRgeoscience/VolatileCalcForPython), based on VB code by Newman &
	Lowenstern.
	"""
	def __init__(self):
		self.set_calibration_ranges([])

	def fugacity(self,pressure,temperature,X_fluid=1.0,**kwargs):
		""" Calculates the fugacity of H2O in a pure or mixed H2O-CO2 fluid (assuming ideal mixing).

		Parameters
		----------
		pressure    float
			Total pressure of the system in bars.
		temperature     float
			Temperature in degC
		X_fluid     float
			Mole fraction of H2O in the fluid.

		Returns
		-------
		float
			fugacity of CO2 in bars
		"""
		fug = self.MRK(pressure,temperature+273.15)
		return fug*X_fluid

	def FNA(self,TK):
		return (166800000 - 193080 * (TK - 273.15) + 186.4 * (TK - 273.15)**2 - 0.071288 * ((TK - 273.15)**3)) * 1.01325

	def FNB(self,TK):
		return 1.01325 * (73030000 - 71400 * (TK - 273.15) + 21.57 * (TK - 273.15)**2)

	def FNC(self,TK):
		R = 83.14321
		return 1.01325 * (np.exp(-11.071 + 5953 / TK - 2746000 / TK**2 + 464600000 / TK**3) * 0.5 * R * R * TK**2.5 / 1.02668 + 40123800)

	def FNF(self,V,TK,A,B,P):
		R = 83.14321
		return R * TK / (V - B) - A / ((V * V + B * V) * TK**0.5) - P

	def MRK(self,P,TK): #Redlich-Kwong routine to estimate endmember H2O and CO2 fugacities
		R = 83.14321
		B_1 = 14.6
		B_2 = 29.7

		# X_1 = 1
		for X_1 in [0,1]:
			B = X_1 * B_1 + (1 - X_1) * B_2
			A = X_1**2 * self.FNA(TK) + 2 * X_1 * (1 - X_1) * self.FNC(TK) + (1 - X_1)**2 * self.FNB(TK)
			Temp2 = B + 5
			Q = 1
			Temp1 = 0
			while abs(Temp2 - Temp1) >= 0.00001:
				Temp1 = Temp2
				F_1 = (self.FNF(Temp1 + 0.01, TK, A, B, P) - self.FNF(Temp1, TK, A, B, P)) / 0.01
				Temp2 = Temp1 - Q * self.FNF(Temp1, TK, A, B, P) / F_1
				F_2 = (self.FNF(Temp2 + 0.01, TK, A, B, P) - self.FNF(Temp2, TK, A, B, P)) / 0.01
				if F_2 * F_1 <= 0:
					Q = Q / 2.
				if abs(Temp2 - Temp1) > 0.00001:
					F_1 = F_2
			V = Temp2
			G_1 = np.log(V / (V - B)) + B_1 / (V - B) - 2 * (X_1 * self.FNA(TK) + (1 - X_1) * self.FNC(TK)) * np.log((V + B) / V) / (R * TK**1.5 * B)
			G_1 = G_1 + (np.log((V + B) / V) - B / (V + B)) * A * B_1 / (R * TK**1.5 * B**2) - np.log(P * V / (R * TK))
			G_1 = np.exp(G_1)
			G_2 = np.log(V / (V - B)) + B_2 / (V - B) - 2 * (X_1 * self.FNC(TK) + (1 - X_1) * self.FNB(TK)) * np.log((V + B) / V) / (R * TK**1.5 * B)
			G_2 = G_2 + (np.log((V + B) / V) - B / (V + B)) * A * B_2 / (R * TK**1.5 * B**2) - np.log(P * V / (R * TK))
			G_2 = np.exp(G_2)
			if X_1 == 0:
				fCO2o = G_2 * P #The fugacity of CO2
				# return fCO2o
			if X_1 == 1:
				fH2Oo = G_1 * P #The fugacity of H2O
				# return fH2Oo
		return fH2Oo

class fugacity_HB_co2(FugacityModel):
	"""
	Implementation of the Holloway and Blank (1994) Modified Redlich Kwong EoS for CO2.
	"""
	def __init__(self):
		self.set_calibration_ranges([CalibrationRange('pressure',[1,1e5],crf_Between,'bar','Redlich Kwong EOS',
													  fail_msg=crmsg_Between_fail, pass_msg=crmsg_Between_pass, description_msg=crmsg_Between_description),
									 CalibrationRange('temperature',500.0,crf_GreaterThan,'oC','Redlich Kwong EOS',
									 				  fail_msg=crmsg_GreaterThan_fail, pass_msg=crmsg_GreaterThan_pass, description_msg=crmsg_GreaterThan_description)])
		self.HBmodel = fugacity_HollowayBlank()

	def fugacity(self,pressure,temperature,X_fluid=1.0,**kwargs):
		return self.HBmodel.fugacity(pressure=pressure, temperature=temperature, species='CO2')*X_fluid

class fugacity_HB_h2o(FugacityModel):
	"""
	Implementation of the Holloway and Blank (1994) Modified Redlich Kwong EoS for H2O.
	"""
	def __init__(self):
		self.set_calibration_ranges([CalibrationRange('pressure',[1,1e5],crf_Between,'bar','Redlich Kwong EOS',
													  fail_msg=crmsg_Between_fail, pass_msg=crmsg_Between_pass, description_msg=crmsg_Between_description),
									 CalibrationRange('temperature',500.0,crf_GreaterThan,'oC','Redlich Kwong EOS',
									 				  fail_msg=crmsg_GreaterThan_fail, pass_msg=crmsg_GreaterThan_pass, description_msg=crmsg_GreaterThan_description)])
		self.HBmodel = fugacity_HollowayBlank()

	def fugacity(self,pressure,temperature,X_fluid=1.0,**kwargs):
		return self.HBmodel.fugacity(pressure=pressure, temperature=temperature, species='H2O')*X_fluid

class fugacity_HollowayBlank(FugacityModel):
	"""
	Implementation of the Modified Redlich Kwong presented in Holloway and Blank (1994) Reviews
	in Mineralogy and Geochemistry vol. 30. Originally written in Quickbasic. CO2 calculations
	translated to Matlab by Chelsea Allison and translated to python by K. Iacovino for VESIcal.
	H2O calculations translated to VisualBasic by Gordon M. Moore and translated to python by
	K. Iacovino for VESIcal.

	"""

	def __init__(self):
		self.set_calibration_ranges([CalibrationRange('pressure',[1,1e5],crf_Between,'bar','MRK EOS (Holloway and Blank, 1994)',
													  fail_msg=crmsg_Between_fail, pass_msg=crmsg_Between_pass, description_msg=crmsg_Between_description),
									 CalibrationRange('temperature',500,crf_GreaterThan,'oC','MRK EOS (Holloway and Blank, 1994)',
									 				  fail_msg=crmsg_GreaterThan_fail, pass_msg=crmsg_GreaterThan_pass, description_msg=crmsg_GreaterThan_description)])


	def REDKW(self, BP, A2B):
		"""
		The RK routine. A routine to calculate compressibility factor and fugacity coefficient
		with the Redlich-Kwong equation following Edmister (1968). This solution for supercritical
		fluid.

		Parameters
		----------
		BP: float
			B parameter sum from RKCALC

		A2B: float
			A parameter sum from RKCALC

		Returns
		-------
		float
			XLNFP (fugacity coefficient?)
		"""
		if A2B < 1*10**(-10):
			A2B = 0.001

		#Define constants
		TH = 0.333333
		RR = -A2B*BP**2
		QQ = BP*(A2B-BP-1)
		XN = QQ*TH+RR-0.074074
		XM = QQ-TH
		XNN = XN*XN*0.25
		XMM = XM**3 / 27.0
		ARG = XNN+XMM

		if ARG > 0:
			X = np.sqrt(ARG)
			F = 1
			XN2 = -XN*0.5
			iXMM = XN2+X

			if iXMM < 0:
				F = -1
			XMM = F*((F*iXMM)**TH)
			F = 1
			iXNN = XN2 - X

			if iXNN < 0:
				F = -1
			XNN = F*((F*iXNN)**TH)
			Z = XMM+XNN+TH
			ZBP = Z-BP
			if ZBP < 0.000001:
				ZBP = 0.000001

			BPZ = 1+BP/Z
			FP = Z-1-np.log(ZBP)-A2B*np.log(BPZ)

			if FP < -37 or FP > 37:
				FP = 0.000001

		elif ARG <0:
			COSPHI = np.sqrt(-XNN/XMM)
			if XN > 0:
				COSPHI = -COSPHI

			TANPHI = np.sqrt(1-COSPHI**2)/COSPHI
			PHI = np.arctan(TANPHI)*TH
			FAC = 2*np.sqrt(-XM*TH)

			#sort for largest root
			R1 = np.cos(PHI)
			R2 = np.cos(PHI+2.0944)
			R3 = np.cos(PHI+4.18879)
			RH = R2

			if R1 > R2:
				RH = R1
			if R3 > RH:
				RH = R3

			Z = RH*FAC+TH
			ZBP = Z-BP
			if ZBP < 0.000001:
				ZBP = 0.000001
			BPZ = 1+BP/Z
			FP = Z-1-np.log(ZBP)-A2B*np.log(BPZ)
			if FP < -37 or FP > 37:
				FP = 0.000001
		else:
			FP = 1
			Z = 1
		XLNFP = FP

		return XLNFP

	def Saxena(self, TK, pb):
		"""
		High pressure corresponding states routines from Saxena and Fei (1987) GCA
		vol. 51, 783-791.

		Parameters
		----------
		TK: float
			Temperature in K.

		pb: float
			Pressure in bars.

		Returns
		-------
		float
			XLNF, Natural log of the ratio F(P)/F(4000 bar)
		"""

		#Define integration limit
		PO = 4000

		#Critical temperatures and pressures for CO2
		TR = TK/304.2
		PR = pb/73.9
		PC = 73.9

		#Virial coeficients
		A = 2.0614-2.2351/TR**2 - 0.39411*np.log(TR)
		B = 0.055125/TR + 0.039344/TR**2
		C = -1.8935*10**(-6)/TR - 1.1092*10**(-5)/TR**2 - 2.1892*10**(-5)/TR**3
		D = 5.0527*10**(-11)/TR - 6.3033*10**(-21)/TR**3

		#Calculate molar volume
		Z = A+B*PR+C*PR**2+D*PR**3
		V = Z*83.0117*TK/pb

		#integrate from PO (4000 bars) to P to calculate ln fugacity
		LNF = A*np.log(pb/PO)+(B/PC)*(pb-PO)+(C/(2*PC**2))*(pb**2-PO**2)
		LNF = LNF+(D/(3*PC**3))*(pb**3-PO**3)
		XLNF = LNF

		return XLNF

	def RKCALC(self, temperature, pressure, species):
		"""
		Calculation of pure gas MRK properties following Holloway 1981, 1987

		Parameters
		----------
		temperature: float
			Temperature in degrees K.

		pressure: float
			Pressure in atmospheres.

		Returns
		-------
		float
			Natural log of the fugacity of a pure gas.
		"""
		#Define constants
		R = 82.05736
		RR = 6732.2
		pb = 1.013*pressure
		PBLN = np.log(pb)
		TCEL = temperature-273.15
		RXT = R*temperature
		RT = R*temperature**1.5 * 10**(-6)

		if species == 'CO2':
			#Calculate T-dependent MRK A parameter CO2
			ACO2M = 73.03 - 0.0714*TCEL + 2.157*10**(-5)*TCEL**2

			#Define MRK B parameter for CO2
			BSUM = 29.7

			ASUM = ACO2M / (BSUM*RT)

		elif species == 'H2O':
			#Calculate T-dependent MRK A parameter H2O
			AH2OM = 115.98 - np.double(0.0016295)*temperature - 1.4984*10**(-5)*temperature**2

			#Define MRK B parameter for H2O
			BSUM = 14.5

			ASUM = AH2OM / (BSUM*RT)

		BSUM = pressure*BSUM/RXT
		XLNFP = self.REDKW(BSUM, ASUM)

		#Convert to ln(fugacity)
		PUREG = XLNFP + PBLN
		return PUREG


	def fugacity(self, pressure, temperature, species, **kwargs):
		"""
		Calculates fugacity.

		Parameters
		----------
		temperature: float
			Temperature in degrees C.

		pressure: float
			Pressure in bars.

		species: str
			Choose which species to calculate. Options are 'H2O' and 'CO2'.

		Returns
		-------
		float
			Fugacity coefficient for passed species
		"""

		#convert temp and press to atmospheres and Kelvin
		pressureAtmo = pressure/1.013
		temperatureK = temperature + 273.15
		PO = 4000/1.013

		#Use the MRK below 4,000 bars, Saxena above 4,000 bars
		if pressure > 4000 and species=='CO2':
			iPUREG = self.RKCALC(temperatureK, PO, species)
			XLNF = self.Saxena(temperatureK, pressure)
			PUREG = iPUREG + XLNF
		else:
			PUREG = self.RKCALC(temperatureK, pressureAtmo, species)

		#Convert from ln(fugacity) to fugacity
		stdf = np.exp(PUREG)
		return stdf

class fugacity_RK_co2(FugacityModel):
	"""
	Implementation of the Redlich Kwong EoS for CO2.
	Code derived from http://people.ds.cam.ac.uk/pjb10/thermo/pure.html - Patrick J. Barrie 30 October 2003.
	"""
	def __init__(self):
		self.set_calibration_ranges([CalibrationRange('pressure',[1,1e5],crf_Between,'bar','Redlich Kwong EOS',
													  fail_msg=crmsg_Between_fail, pass_msg=crmsg_Between_pass, description_msg=crmsg_Between_description),
									 CalibrationRange('temperature',[500],crf_GreaterThan,'oC','Redlich Kwong EOS',
									 				  fail_msg=crmsg_GreaterThan_fail, pass_msg=crmsg_GreaterThan_pass, description_msg=crmsg_GreaterThan_description)])
		# self.set_calibration_ranges([cr_Between('pressure',[1.0,1e5],'bar','Redlich Kwong EOS'),
		# 							 cr_GreaterThan('temperature',500,'oC','Redlich Kwong EOS')])
		self.RKmodel = fugacity_RedlichKwong()

	def fugacity(self,pressure,temperature,X_fluid,**kwargs):
		return self.RKmodel.fugacity(pressure, temperature, X_fluid, 'CO2')

class fugacity_RK_h2o(FugacityModel):
	"""
	Implementation of the Redlich Kwong EoS for H2O.
	Code derived from http://people.ds.cam.ac.uk/pjb10/thermo/pure.html - Patrick J. Barrie 30 October 2003.
	"""
	def __init__(self):
		self.set_calibration_ranges([CalibrationRange('pressure',[1,1e5],crf_Between,'bar','Redlich Kwong EOS',
													  fail_msg=crmsg_Between_fail, pass_msg=crmsg_Between_pass, description_msg=crmsg_Between_description),
									 CalibrationRange('temperature',500,crf_GreaterThan,'oC','Redlich Kwong EOS',
									 				  fail_msg=crmsg_GreaterThan_fail, pass_msg=crmsg_GreaterThan_pass, description_msg=crmsg_GreaterThan_description)])
		self.RKmodel = fugacity_RedlichKwong()

	def fugacity(self,pressure,temperature,X_fluid,**kwargs):
		return self.RKmodel.fugacity(pressure, temperature, X_fluid, 'H2O')

class fugacity_RedlichKwong(FugacityModel):
	"""
	Implementation of the Redlich Kwong EoS
	Code derived from http://people.ds.cam.ac.uk/pjb10/thermo/pure.html - Patrick J. Barrie 30 October 2003.
	"""
	def __init__(self):
		self.set_calibration_ranges([CalibrationRange('pressure',[1,1e5],crf_Between,'bar','Redlich Kwong EOS',
													  fail_msg=crmsg_Between_fail, pass_msg=crmsg_Between_pass, description_msg=crmsg_Between_description),
									 CalibrationRange('temperature',500,crf_GreaterThan,'oC','Redlich Kwong EOS',
									 				  fail_msg=crmsg_GreaterThan_fail, pass_msg=crmsg_GreaterThan_pass, description_msg=crmsg_GreaterThan_description)])


	def gamma(self, pressure, temperature, species):
		"""
		Calculates fugacity coefficients.

		Parameters
		----------
		temperature: fload
			Temperature in degrees C.

		pressure: float
			Pressure in bars.

		species: str
			Choose which species to calculate. Options are 'H2O' and 'CO2'.

		Returns
		-------
		float
			Fugacity coefficient for passed species.
		"""

		temperatureK = temperature + 273.15
		R = 8.3145

		fluid_species_names = ['CO2', 'H2O']
		critical_params = {'CO2':{  "cT":   304.15,
							"cP":   73.8659,
							"o":    0.225
								},
							'H2O':{ "cT":   647.25,
							"cP":   221.1925,
							"o":    0.334
								}
							}

		#Calculate a and b parameters (depend only on critical parameters)...
		a = 0.42748 * R**2.0 * critical_params[species]["cT"]**(2.5) / (critical_params[species]["cP"] * 10.0**5)
		b = 0.08664 * R * critical_params[species]["cT"] / (critical_params[species]["cP"] * 10.0**5)
		kappa = 0.0

		#Calculate coefficients in the cubic equation of state...
		#coeffs: (C0, C1, C2, A, B)
		A = a * pressure * 10.0**5 / (np.sqrt(temperatureK) * (R * temperatureK)**2.0)
		B = b * pressure * 10.0**5 / (R * temperatureK)
		C2 = -1.0
		C1 = A - B - B * B
		C0 = -A * B

		#Solve the cubic equation for Z0 - Z2, D...
		Q1 = C2 * C1 / 6.0 - C0 / 2.0 - C2**3.0 / 27.0
		P1 = C2**2.0 / 9.0 - C1 / 3.0
		D = Q1**2.0 - P1**3.0

		if D >= 0:
			kOneThird = 1.0 / 3.0

			absQ1PSqrtD = np.fabs(Q1 + np.sqrt(D))
			temp1 = absQ1PSqrtD**kOneThird
			temp1 *= (Q1 + np.sqrt(D)) / absQ1PSqrtD

			absQ1MSqrtD = np.fabs(Q1 - np.sqrt(D))
			temp2 = absQ1MSqrtD**kOneThird
			temp2 *= (Q1 - np.sqrt(D)) / absQ1MSqrtD

			Z0 = temp1 + temp2 - C2 / 3.0
		else:
			temp1 = Q1**2.0 / (P1**3.0)
			temp2 = np.sqrt(1.0 - temp1) / np.sqrt(temp1)
			temp2 *= Q1 / np.fabs(Q1)

			gamma = np.arctan(temp2)

			if gamma < 0:
				gamma = gamma + np.pi

			Z0 = 2.0 * np.sqrt(P1) * np.cos(gamma/3.0) - C2 / 3.0
			Z1 = 2.0 * np.sqrt(P1) * np.cos((gamma + 2.0 * np.pi) / 3.0) - C2/3.0
			Z2 = 2.0 * np.sqrt(P1) * np.cos((gamma + 4.0 * np.pi) / 3.0) - C2/3.0

			if Z0 < Z1:
				temp0 = Z0
				Z0 = Z1
				Z1 = temp0

			if Z1 < Z2:
				temp0 = Z1
				Z1 = Z2
				Z2 = temp0

			if Z0 < Z1:
				temp0 = Z0
				Z0 = Z1
				Z1 = temp0

		#Calculate Departure Functions
		gamma = np.exp(Z0 - 1.0 - np.log(Z0-B) - A * np.log(1.0+B/Z0)/B)
		Hdep = R * temperatureK * (Z0 - 1.0 - 1.5*A*np.log(1.0+B/Z0)/B)
		Sdep = R * (np.log(Z0-B) - 0.5*A*np.log(1.0+B/Z0)/B)

		return gamma

	def fugacity(self, pressure, temperature, X_fluid=1.0, species='H2O', **kwargs):
		"""
		Calculates the fugacity of H2O in a mixed H2O-CO2 fluid using the universal relationships:
		P_i = f_i/gamma_i = (fpure_i * Xfluid_i) / gamma_i
		See Iacovino (2015) EPSL for further explanation.
		"""

		gammaH2O = self.gamma(pressure, temperature, 'H2O')
		gammaCO2 = self.gamma(pressure, temperature, 'CO2')

		fugacityH2Opure = pressure * gammaH2O
		fugacityCO2pure = pressure * gammaCO2

		if species == 'H2O':
			return fugacityH2Opure * X_fluid
		elif species == 'CO2':
			return fugacityCO2pure * X_fluid
		else:
			raise InputError("Species must be H2O or CO2.")




#---------------ACTVITY MODELS-------------------------------#


class activity_idealsolution(activity_model):
	""" Implements an ideal solution activity model, i.e. it
	will always return the mole fraction.
	"""

	def activity(self,X):
		""" The activity of the component in an ideal solution, i.e., it
		will return the mole fraction.

		Parameters
		----------
		X   float
			The mole fraction of the species in the solution.

		Returns
		-------
		float
			The activity of the species in the solution, i.e., the mole fraction.
		"""
		return X





#------------PURE FLUID MODELS-------------------------------#

class ShishkinaCarbon(Model):
	""" Implementation of the Shishkina et al. (2014) carbon solubility model, as a Model class.
	"""
	def __init__(self):
		self.set_volatile_species(['CO2'])
		self.set_fugacity_model(fugacity_idealgas())
		self.set_activity_model(activity_idealsolution())
		self.set_solubility_dependence(False)
		self.set_calibration_ranges([CalibrationRange('pressure',[500.0,5000.0],crf_Between,'bar','Shishkina et al. carbon',
													  fail_msg=crmsg_Between_fail, pass_msg=crmsg_Between_pass, description_msg=crmsg_Between_description),
									 CalibrationRange('temperature',[1200.0,1300.0],crf_Between,'oC','Shishkina et al. carbon',
									 				  fail_msg=crmsg_Between_fail, pass_msg=crmsg_Between_pass, description_msg=crmsg_Between_description),
									 CalibrationRange('SiO2',[40,57],crf_Between,'wt%','Shishkina et al. carbon',
									 				  fail_msg=crmsg_BC_fail_ShC1, pass_msg=crmsg_Between_pass, description_msg=crmsg_Between_description)])


	def preprocess_sample(self,sample):
		""" Returns sample, unmodified. The Pi* compositional parameter is a ratio of cations,
		therefore the value is not affected by the normalization of the sample. Shishkina et al.
		imply the accuracy of the calculations are little affected whether Fe(tot) or Fe2+ is
		used.

		Parameters
		----------
		sample:        dict or pandas Series
			The major element oxides in wt%.

		Returns
		-------
		dict or pandas Series
			The major element oxides in wt%.

		"""
		if type(sample) != dict and type(sample) != pd.core.series.Series:
			raise InputError("sample must be a dict or a pandas Series.")

		return sample

	def PiStar(self,sample):
		"""Shishkina et al. (2014) Eq (11)

		Calculates the Pi* parameter for use in calculating CO2 solubility.

		Parameters
		----------
		sample:        pandas Series or dict
			Major element oxides in wt%.

		Returns
		-------
		float
			The value of the Pi* compositional parameter.
		"""

		_mols = wtpercentOxides_to_molCations(sample)

		if all(cation in _mols for cation in ['Ca','K','Na','Mg','Fe','Si','Al']) == False:
			raise InputError("To calculate PiStar, values for CaO, K2O, Na2O, MgO, FeO, SiO2, and Al2O3\
								must be provided in sample.")

		_pi = (_mols['Ca'] + 0.8*_mols['K'] + 0.7*_mols['Na'] + 0.4*_mols['Mg'] + 0.4*_mols['Fe'])/\
				(_mols['Si']+_mols['Al'])

		return _pi

	def calculate_dissolved_volatiles(self,pressure,sample,X_fluid=1,**kwargs):
		""" Calculates the dissolved CO2 concentration in wt%, using equation (13) of Shishkina et al. (2014).

		Parameters
		----------
		pressure:    float
			(Total) pressure in bars.
		sample:        dict or pandas Series
			Major element concentrations in wt%. Normalization does not matter.
		X_fluid:    float
			The mol-fraction of the fluid that is CO2. Default is 1, i.e. a pure CO2 fluid.

		Returns
		-------
		float
			The dissolved CO2 concentration in wt%.
		"""

		if X_fluid < 0 or X_fluid > 1:
			raise InputError("X_fluid must have a value between 0 and 1.")
		if pressure < 0:
			raise InputError("pressure must be a positive value.")

		PiStar = self.PiStar(sample)
		fugacity = self.fugacity_model.fugacity(pressure=pressure,X_fluid=X_fluid,**kwargs)

		A = 1.150
		B = 6.71
		C= -1.345

		if fugacity == 0:
			return 0
		else:
			return np.exp(A*np.log(fugacity/10)+B*PiStar+C)/1e4


	def calculate_equilibrium_fluid_comp(self,pressure,sample,**kwargs):
		""" Returns 1.0 if a pure CO2 fluid is saturated.
		Returns 0.0 if a pure CO2 fluid is undersaturated.

		Parameters
		----------
		pressure     float
			The total pressure of the system in bars.
		sample         dict or pandas Series
			Major element oxides in wt%

		Returns
		-------
		float
			1.0 if CO2-fluid saturated, 0.0 otherwise.
		"""
		if self.calculate_saturation_pressure(sample=sample,**kwargs) < pressure:
			return 0.0
		else:
			return 1.0

	def calculate_saturation_pressure(self,sample,**kwargs):
		""" Calculates the pressure at which a pure CO2 fluid is saturated, for the given
		sample composition and CO2 concentration. Calls the scipy.root_scalar routine, which makes
		repeated calls to the calculate_dissolved_volatiles method.

		Parameters
		----------
		sample         dict or pandas Series
			Major elements in wt%, including CO2 (also in wt%).

		Returns
		-------
		float
			Saturation pressure in bar
		"""

		if 'CO2' not in sample:
			raise InputError("sample must contain CO2.")
		if sample['CO2'] < 0:
			raise InputError("CO2 concentration must be greater than 0 wt%.")

		try:
			satP = root_scalar(self.root_saturation_pressure,bracket=[1e-15,1e5],args=(sample,kwargs)).root
		except:
			w.warn("Saturation pressure not found.",RuntimeWarning,stacklevel=2)
			satP = np.nan
		return satP

	def root_saturation_pressure(self,pressure,sample,kwargs):
		""" Function called by scipy.root_scalar when finding the saturation pressure using
		calculate_saturation_pressure.

		Parameters
		----------
		pressure     float
			Pressure guess in bars
		sample         dict or pandas Series
			Major element oxides in wt%, including CO2 (also in wt%).
		kwargs         dictionary
			Additional keyword arguments supplied to calculate_saturation_pressure. Might be required for
			the fugacity or activity models.

		Returns
		-------
		float
			The differece between the dissolved CO2 at the pressure guessed, and the CO2 concentration
			passed in the sample variable.
		"""
		return self.calculate_dissolved_volatiles(pressure=pressure,sample=sample,**kwargs)-sample['CO2']



class ShishkinaWater(Model):
	""" Implementation of the Shishkina et al. (2014) H2O solubility model as a Model class.
	"""
	def __init__(self):
		self.set_volatile_species(['H2O'])
		self.set_fugacity_model(fugacity_idealgas())
		self.set_activity_model(activity_idealsolution())
		self.set_solubility_dependence(False)
		self.set_calibration_ranges([CalibrationRange('pressure',[500.0,5000.0],crf_Between,'bar','Shishkina et al. water',
													  fail_msg=crmsg_Between_fail, pass_msg=crmsg_Between_pass, description_msg=crmsg_Between_description),
									CalibrationRange('temperature',[1050,1400],crf_Between,'oC','Shishkina et al. water',
									 				  fail_msg=crmsg_BC_fail_ShT1, pass_msg=crmsg_Between_pass,
									 				  description_msg=crmsg_Between_description),
									CalibrationRange('SiO2',65,crf_LessThan,'wt%','Shishkina et al. water',
									 				  fail_msg=crmsg_LessThan_fail_ShWaterSi, pass_msg=crmsg_LessThan_pass,
									 				  description_msg=crmsg_LessThan_description),
									CalibrationRange('SiO2',40,crf_GreaterThan,'wt%','Shishkina et al. water',
									 				  fail_msg=crmsg_GreaterThan_fail_ShWaterSi, pass_msg=crmsg_GreaterThan_pass,
									 				  description_msg=crmsg_GreaterThan_description)])


	def preprocess_sample(self,sample):
		""" Returns sample, renormlized so that the major element oxides (excluding volatiles) sum to 100%.
		Normalization must be done this way as the compositional dependence of the solubility takes the
		mole fractions of Na2O and K2O as inputs, presumably assuming no volatiles in the bulk composition.
		Volatile concentrations are left unchanged.

		Parameters
		----------
		sample:         dict or pandas Series
			The major element oxides in wt%.

		Returns
		-------
		dict or pandas Series
			The major element oxides in wt%.

		"""
		return normalize_AdditionalVolatiles(sample)

	def calculate_dissolved_volatiles(self,pressure,sample,X_fluid=1.0,**kwargs):
		"""Calculates the dissolved H2O concentration using Eqn (9) of Shishkina et al. (2014).

		Parameters
		----------
		pressure     float
			Total pressure in bars
		sample         pandas Series or dict
			Major element oxides in wt%. Normalized to zero-volatiles so that the total-alkalis
			mol fraction can be determined accurately.
		X_fluid     float
			The mol fraction of H2O in the fluid

		Returns
		-------
		float
			The H2O concentration in wt%
		"""
		if type(sample) != dict and type(sample) != pd.core.series.Series:
			raise InputError("sample must be a dict or pandas Series.")

		if all(ox in sample for ox in ['Na2O','K2O']) == False:
			raise InputError("Na2O and K2O must be present in sample.")

		if pressure < 0:
			raise InputError("Pressure must be positive.")

		_mols = wtpercentOxides_to_molCations(sample)
		_mol_volatiles = 0
		if 'H' in _mols:
			_mol_volatiles += _mols['H']
		if 'C' in _mols:
			_mol_volatiles += _mols['C']

		total_alkalis = (_mols['Na'] + _mols['K'])/(1-_mol_volatiles)

		fugacity = self.fugacity_model.fugacity(pressure,X_fluid=X_fluid,**kwargs)

		a = 3.36e-7 * (fugacity/10)**3 - 2.33e-4*(fugacity/10)**2 + 0.0711*(fugacity/10) - 1.1309
		b = -1.2e-5*(fugacity/10)**2 + 0.0196*(fugacity/10)+1.1297

		return a*total_alkalis + b


	def calculate_equilibrium_fluid_comp(self,pressure,sample,**kwargs):
		""" Returns 1.0 if a pure H2O fluid is saturated.
		Returns 0.0 if a pure H2O fluid is undersaturated.

		Parameters
		----------
		pressure     float
			The total pressure of the system in bars.
		sample         pandas Series or dict
			Major element oxides in wt%, normalized on the basis of
			no volatiles.

		Returns
		-------
		float
			1.0 if H2O-fluid saturated, 0.0 otherwise.
		"""
		if self.calculate_saturation_pressure(sample=sample,**kwargs) < pressure:
			return 0.0
		else:
			return 1.0

	def calculate_saturation_pressure(self,sample,**kwargs):
		""" Calculates the pressure at which a pure H2O fluid is saturated, for the given
		sample composition and H2O concentration. Calls the scipy.root_scalar routine, which makes
		repeated calls to the calculate_dissolved_volatiles method.

		Parameters
		----------
		sample         pandas Series or dict
			Major elements in wt% (normalized to 100%), including H2O (also in wt%, not included
			in normalization).

		Returns
		-------
		float
			Saturation pressure in bar
		"""
		if 'H2O' not in sample:
			raise InputError("sample must contain H2O")
		if sample['H2O'] < 0:
			raise InputError("H2O concentration must be greater than 0 wt%.")

		if sample['H2O'] < self.calculate_dissolved_volatiles(sample=sample,pressure=0,**kwargs):
			return np.nan

		try:
			satP = root_scalar(self.root_saturation_pressure,bracket=[1e-15,1e5],args=(sample,kwargs)).root
		except:
			w.warn("Saturation pressure not found.",RuntimeWarning,stacklevel=2)
			satP = np.nan
		return satP

	def root_saturation_pressure(self,pressure,sample,kwargs):
		""" Function called by scipy.root_scalar when finding the saturation pressure using
		calculate_saturation_pressure.

		Parameters
		----------
		pressure     float
			Pressure guess in bars
		sample         pandas Series or dict
			Major elements in wt% (normalized to 100%), including H2O (also in wt%, not included
			in normalization).
		kwargs         dictionary
			Additional keyword arguments supplied to calculate_saturation_pressure. Might be required for
			the fugacity or activity models.

		Returns
		-------
		float
			The differece between the dissolved H2O at the pressure guessed, and the H2O concentration
			passed in the sample variable.
		"""
		return self.calculate_dissolved_volatiles(pressure=pressure,sample=sample,**kwargs)-sample['H2O']



class DixonCarbon(Model):
	"""
	Implementation of the Dixon (1997) carbon solubility model, as a Model class.
	"""

	def __init__(self):
		self.set_volatile_species(['CO2'])
		self.set_fugacity_model(fugacity_MRK_co2())
		self.set_activity_model(activity_idealsolution())
		self.set_calibration_ranges([])
		self.set_solubility_dependence(False)

	def preprocess_sample(self,sample):
		""" Returns sample, normalized, keep volatiles unchanged.

		Parameters
		----------
		sample:     pandas Series or dict
			The major element oxides in wt%.

		Returns
		-------
		pandas Series or dict
			The major element oxides in wt%.
		"""
		return normalize_AdditionalVolatiles(sample)

	def calculate_dissolved_volatiles(self,pressure,sample,X_fluid=1.0,**kwargs):
		"""Calculates the dissolved CO2 concentration using Eqn (3) of Dixon (1997).

		Parameters
		----------
		pressure  float
			Total pressure in bars.
		sample      pandas Series or dict
			Major element oxides in wt%.
		X_fluid      float
			The mol fraction of CO2 in the fluid.

		Returns
		-------
		float
			The CO2 concentration in wt%.
		"""

		if X_fluid < 0 or X_fluid > 1:
			raise InputError("X_fluid must have a value between 0 and 1.")
		if pressure < 0:
			raise InputError("Pressure must be positive.")
		if type(sample) != dict and type(sample) != pd.core.series.Series:
			raise InputError("sample must be a dict or pandas Series")
		if 'SiO2' not in sample:
			raise InputError("sample must contain SiO2.")

		if pressure == 0:
			return 0


		Mr = wtpercentOxides_to_formulaWeight(sample,exclude_volatiles=True)
		XCO3 = self.molfrac_molecular(pressure=pressure,sample=sample,X_fluid=X_fluid,**kwargs)
		# return (4400 * XCO3) / (36.6 - 44*XCO3)
		return (4400 * XCO3) / (Mr - 44*XCO3)


	def calculate_equilibrium_fluid_comp(self,pressure,sample,**kwargs):
		""" Returns 1.0 if a pure H2O fluid is saturated.
		Returns 0.0 if a pure H2O fluid is undersaturated.

		Parameters
		----------
		pressure     float
			The total pressure of the system in bars.
		sample         pandas Series or dict
			Major element oxides in wt% (including CO2).

		Returns
		-------
		float
			1.0 if CO2-fluid saturated, 0.0 otherwise.
		"""
		if self.calculate_saturation_pressure(sample=sample,**kwargs) < pressure:
			return 0.0
		else:
			return 1.0

	def calculate_saturation_pressure(self,sample,X_fluid=1.0,**kwargs):
		"""
		Calculates the pressure at which a pure CO2 fluid is saturated, for the given sample
		composition and CO2 concentration. Calls the scipy.root_scalar routine, which makes
		repeated called to the calculate_dissolved_volatiles method.

		Parameters
		----------
		sample         pandas Series or dict
			Major element oxides in wt% (including CO2).
		X_fluid     float
			The mole fraction of CO2 in the fluid. Default is 1.0.

		Returns
		-------
		float
			Calculated saturation pressure in bars.
		"""
		if type(sample) != dict and type(sample) != pd.core.series.Series:
			raise InputError("sample must be a dict or pandas Series.")
		if 'CO2' not in sample:
			raise InputError("sample must contain CO2.")
		if sample['CO2'] < 0:
			raise InputError("Dissolved CO2 concentration must be greater than 0 wt%.")

		try:
			satP = root_scalar(self.root_saturation_pressure,x0=100.0,x1=1000.0,args=(sample,kwargs)).root
		except:
			w.warn("Saturation pressure not found.",RuntimeWarning,stacklevel=2)
			satP = np.nan
		return np.real(satP)

	def molfrac_molecular(self,pressure,sample,X_fluid=1.0,**kwargs):
		"""Calculates the mole fraction of CO3(-2) dissolved when in equilibrium with
		a pure CO2 fluid at 1200C, using Eqn (1) of Dixon (1997).

		Parameters
		----------
		pressure      float
			Total pressure in bars.
		sample         pandas Series or dict
			Major element oxides in wt%.
		X_fluid     float
			Mole fraction of CO2 in the fluid.

		Returns
		-------
		float
			Mole fraction of CO3(2-) dissolved."""

		DeltaVr = 23.14 #cm3 mole-1
		P0 = 1
		R = 83.15
		T0 = 1473.15

		fugacity = self.fugacity_model.fugacity(pressure=pressure,X_fluid=X_fluid,**kwargs)

		XCO3Std = self.XCO3_Std(sample)

		return XCO3Std * fugacity * np.exp(-DeltaVr * (pressure-P0)/(R*T0))

	def XCO3_Std(self,sample):
		""" Calculates the mole fraction of CO3(2-) dissolved when in equilibrium with pure
		CO2 vapour at 1200C and 1 bar, using Eq (8) of Dixon (1997).

		Parameters
		----------
		sample    pandas Series or dict
			Major element oxides in wt%.

		Returns
		-------
		float
			Mole fraction of CO3(2-) dissolved at 1 bar and 1200C.
		"""
		if sample['SiO2'] > 48.9:
			return 3.817e-7
		else:
			return 8.697e-6 - 1.697e-7*sample['SiO2']

	def root_saturation_pressure(self,pressure,sample,kwargs):
		""" The function called by scipy.root_scalar when finding the saturation pressure using
		calculate_saturation_pressure.

		Parameters
		----------
		pressure     float
			Total pressure in bars.
		sample         pandas Series or dict
			Major element oxides in wt% (including CO2).

		Returns
		-------
		float
			The difference between the dissolved CO2 the pressure guessed, and the CO2 concentration
			passed in the sample variable.
		"""
		return self.calculate_dissolved_volatiles(pressure=pressure,sample=sample,**kwargs) - sample['CO2']




class DixonWater(Model):
	"""
	Implementation of the Dixon (1997) water solubility model, as a Model class.
	"""

	def __init__(self):
		self.set_volatile_species(['H2O'])
		self.set_fugacity_model(fugacity_MRK_h2o())
		self.set_activity_model(activity_idealsolution())
		self.set_calibration_ranges([CalibrationRange('pressure',1000,crf_LessThan,'bar','Dixon (1997, Pi-SiO2 simpl.) Water',
													  fail_msg=Dix_1000bar_fail, pass_msg=crmsg_LessThan_pass, description_msg=crmsg_LessThan_description),
									 CalibrationRange('pressure',2000,crf_LessThan,'bar','Dixon (1997, Pi-SiO2 simpl.) Water',
													  fail_msg=Dix_2000bar_fail, pass_msg=crmsg_LessThan_pass, description_msg=crmsg_LessThan_description),                             									 CalibrationRange('pressure',5000,crf_LessThan,'bar','Dixon (1997, Pi-SiO2 simpl.) Water',
													  fail_msg=Dix_5000bar_fail, pass_msg=crmsg_LessThan_pass, description_msg=crmsg_LessThan_description),
									 CalibrationRange('SiO2',49,crf_LessThan,'wt%','Dixon (1997, Pi-SiO2 simpl.) Water',
													  fail_msg=Dix_49_fail, pass_msg=crmsg_LessThan_pass, description_msg=crmsg_LessThan_description),
									 CalibrationRange('SiO2',40,crf_GreaterThan,'wt%','Dixon (1997, Pi-SiO2 simpl.) Water',
													  fail_msg=Dix_40_fail, pass_msg=crmsg_GreaterThan_pass, description_msg=crmsg_GreaterThan_description),
									 CalibrationRange('temperature',[1000,1400],crf_Between,'oC','Dixon (1997, Pi-SiO2 simpl.) Water',
									 				  fail_msg=crmsg_BC_DixonT, pass_msg=crmsg_Between_pass, description_msg=crmsg_Between_description)])
		self.set_solubility_dependence(False)


	def preprocess_sample(self,sample):
		""" Returns sample, normalized, holding volatile concentrations constant.

		Parameters
		----------
		sample:     pandas Series or dict
			The major element oxides in wt%.

		Returns
		-------
		pandas Series or dict
			The major element oxides in wt%.
		"""
		return normalize_AdditionalVolatiles(sample)

	def calculate_dissolved_volatiles(self,pressure,sample,X_fluid=1.0,**kwargs):
		"""Calculates the dissolved H2O concentration using Eqns (5) and (6) of Dixon (1997).

		Parameters
		----------
		pressure  float
			Total pressure in bars.
		sample      pandas Series or dict
			Major element oxides in wt%.
		X_fluid      float
			The mol fraction of H2O in the fluid.

		Returns
		-------
		float
			The H2O concentration in wt%.
		"""
		if type(sample) != dict and type(sample) != pd.core.series.Series:
			raise InputError("sample must be a dict or a pandas Series.")
		if 'SiO2' not in sample:
			raise InputError("sample must contain SiO2.")
		if pressure < 0:
			raise InputError("Pressure must be positive")
		if X_fluid < 0 or X_fluid > 1:
			raise InputError("X_fluid must have a value between 0 and 1.")

		if pressure == 0:
			return 0

		XH2O = self.molfrac_molecular(pressure=pressure,sample=sample,X_fluid=X_fluid,**kwargs)
		XOH = self.XOH(pressure=pressure,sample=sample,X_fluid=X_fluid,**kwargs)

		Mr = wtpercentOxides_to_formulaWeight(sample,exclude_volatiles=True)

		XB = XH2O + 0.5*XOH
		# return 1801.5*XB/(36.6-18.6*XB)
		return 1801.5*XB/(Mr-18.6*XB)


	def calculate_equilibrium_fluid_comp(self,pressure,sample,**kwargs):
		""" Returns 1.0 if a pure H2O fluid is saturated.
		Returns 0.0 if a pure H2O fluid is undersaturated.

		Parameters
		----------
		pressure     float
			The total pressure of the system in bars.
		sample         pandas Series or dict
			Major element oxides in wt% (including H2O).

		Returns
		-------
		float
			1.0 if H2O-fluid saturated, 0.0 otherwise.
		"""
		if self.calculate_saturation_pressure(sample=sample,**kwargs) < pressure:
			return 0.0
		else:
			return 1.0

	def calculate_saturation_pressure(self,sample,X_fluid=1.0,**kwargs):
		"""
		Calculates the pressure at which a pure H2O fluid is saturated, for the given sample
		composition and H2O concentration. Calls the scipy.root_scalar routine, which makes
		repeated called to the calculate_dissolved_volatiles method.

		Parameters
		----------
		sample         pandas Series or dict
			Major element oxides in wt% (including H2O).
		X_fluid     float
			The mole fraction of H2O in the fluid. Default is 1.0.

		Returns
		-------
		float
			Calculated saturation pressure in bars.
		"""
		if 'H2O' not in sample:
			raise InputError("sample must contain H2O")
		if sample['H2O'] < 0:
			raise InputError("H2O concentration must be greater than 0 wt%.")
		try:
			satP = root_scalar(self.root_saturation_pressure,x0=100.0,x1=1000.0,args=(sample,kwargs)).root
		except:
			w.warn("Saturation pressure not found.",RuntimeWarning,stacklevel=2)
			satP = np.nan
		return np.real(satP)

	def molfrac_molecular(self,pressure,sample,X_fluid=1.0,**kwargs):
		"""Calculates the mole fraction of molecular H2O dissolved when in equilibrium with
		a pure H2O fluid at 1200C, using Eqn (2) of Dixon (1997).

		Parameters
		----------
		pressure      float
			Total pressure in bars.
		sample         pandas Series or dict
			Major element oxides in wt%.
		X_fluid     float
			Mole fraction of H2O in the fluid.

		Returns
		-------
		float
			Mole fraction of molecular H2O dissolved.
		"""

		VH2O = 12 #cm3 mole-1
		P0 = 1
		R = 83.15
		T0 = 1473.15

		XH2OStd = self.XH2O_Std(sample)

		fugacity = self.fugacity_model.fugacity(pressure=pressure,X_fluid=X_fluid,**kwargs)

		return XH2OStd * fugacity * np.exp(-VH2O * (pressure-P0)/(R*T0))


	def XH2O_Std(self,sample):
		""" Calculates the mole fraction of molecular H2O dissolved when in equilibrium with pure
		H2O vapour at 1200C and 1 bar, using Eq (9) of Dixon (1997).

		Parameters
		----------
		sample    pandas Series or dict
			Major element oxides in wt%.

		Returns
		-------
		float
			Mole fraction of molecular water dissolved at 1 bar and 1200C.
		"""
		if sample['SiO2'] > 48.9:
			return 3.28e-5
		else:
			return -3.04e-5 + 1.29e-6*sample['SiO2']

	def XOH(self,pressure,sample,X_fluid=1.0,**kwargs):
		"""
		Calculates the mole fraction of hydroxyl groups dissolved by solving Eq (4) of
		Dixon (1997). Calls scipy.root_scalar to find the root of the XOH_root method.

		Parameters
		----------
		pressure     float
			Total pressure in bars.
		sample         pandas Series or dict
			Major element oxides in wt%.
		X_fluid     float
			Mole fraction of H2O in the fluid.

		Returns
		-------
		float
			Mole fraction of hydroxyl groups dissolved.
		"""

		XH2O = self.molfrac_molecular(pressure=pressure,sample=sample,X_fluid=X_fluid,**kwargs)
		if XH2O < 1e-14:
			return 0
		return np.exp(root_scalar(self.XOH_root,x0=np.log(0.5),x1=np.log(0.1),args=(XH2O)).root)

	def XOH_root(self,XOH,XH2O):
		"""
		Method called by scipy.root_scalar when finding the saturation pressure using the
		calculate_saturation_pressure method. Implements Eq (4) of Dixon (1997).

		Parameters
		----------
		XOH         float
			Guess for the mole fraction of hydroxyl groups dissolved in melt.
		XH2O    float
			Mole fraction of molecular water dissolved in melt.

		Returns
		-------
		float
			The difference between the RHS and LHS of Eq (4) of Dixon (1997) for the
			guessed value of XOH.
		"""

		A = 0.403
		B = 15.333
		C = 10.894

		XOH = np.exp(XOH)

		term = (XOH)**2.0/(XH2O*(1.0-XOH-XH2O))

		lhs = - np.log(term)

		rhs = A + B*XOH + C*XH2O

		return rhs - lhs

	def root_saturation_pressure(self,pressure,sample,kwargs):
		""" Function called by scipy.root_scalar when finding the saturation pressure using
		calculate_saturation_pressure.

		Parameters
		----------
		pressure     float
			Pressure guess in bars
		sample         pandas Series or dict
			Major elements in wt% (normalized to 100%), including H2O.
		kwargs         dictionary
			Additional keyword arguments supplied to calculate_saturation_pressure. Might be required for
			the fugacity or activity models.

		Returns
		-------
		float
			The differece between the dissolved H2O at the pressure guessed, and the H2O concentration
			passed in the sample variable.
		"""
		return self.calculate_dissolved_volatiles(pressure=pressure,sample=sample,**kwargs) - sample['H2O']


class IaconoMarzianoWater(Model):
	"""
	Implementation of the Iacono-Marziano et al. (2012) water solubility model, as a Model class. Two
	calibrations are provided- the one incorporating the H2O content as a parameter (hydrous), and the
	one that does not (anhydrous). Specify which should be used when initialising the model, with the
	bool variable hydrous.
	"""

	def __init__(self,hydrous=True):
		"""
		Initialise the model.

		Parameters
		----------
		hydrous     bool
			Whether to use the hydrous parameterization, or not.
		"""
		self.set_volatile_species(['H2O'])
		self.set_fugacity_model(fugacity_idealgas())
		self.set_activity_model(activity_idealsolution())
		self.hydrous = hydrous
		self.set_calibration_ranges([CalibrationRange('temperature',[1000,1400],crf_Between,'oC','IaconoMarzianoWater',
									 				  fail_msg=crmsg_BC_IMT, pass_msg=crmsg_Between_pass, description_msg=crmsg_Between_description),
									 CalibrationRange('pressure',[100,10000],crf_Between,'bars','IaconoMarzianoWater',
									 				  fail_msg=crmsg_BC_IMP, pass_msg=crmsg_Between_pass, description_msg=crmsg_Between_description)])
		self.set_solubility_dependence(False) #Not dependent on CO2 conc, H2O dependence dealt with within model.

	def preprocess_sample(self,sample):
		"""
		Returns sample, normalized to 100 wt%, without changing the wt% of H2O and CO2 if the
		hydrous parameterization is being used (default). If the anhydrous parameterization is
		used, it will normalize without including H2O and CO2.

		Parameters
		----------
		sample     pandas Series or dict
			Major element oxides in wt%.

		Returns
		-------
		pandas Series or dict
			Major element oxides normalized to wt%.
		"""
		if self.hydrous == True:
			return normalize_FixedVolatiles(sample)
		else:
			return normalize_AdditionalVolatiles(sample)

	def calculate_dissolved_volatiles(self,pressure,temperature,sample,X_fluid=1.0,
									  hydrous_coeffs=True,webapp_coeffs=False,**kwargs):
		"""
		Calculates the dissolved H2O concentration, using Eq (13) of Iacono-Marziano et al. (2012).
		If using the hydrous parameterization, it will use the scipy.root_scalar routine to find the
		root of the root_dissolved_volatiles method.

		Parameters
		----------
		pressure    float
			Total pressure in bars.
		temperature     float
			Temperature in C
		sample     pandas Series or dict
			Major element oxides in wt%.
		X_fluid      float
			Mole fraction of H2O in the fluid. Default is 1.0.
		hydrous_coeffs 	bool
			Use the hydrous or anhydrous NBO/O paramterisation (True for hydrous). Default is True.
		webapp_coeffs 	bool
			If True, use the pre-review hydrous coefficients, as implemented in the IM webapp.
			Default is False.

		Returns
		-------
		float
			Dissolved H2O concentration in wt%.
		"""

		temperature = temperature + 273.15 #translate T from C to K

		if type(sample) != dict and type(sample) != pd.core.series.Series:
			raise InputError("sample must be a dict or a pandas Series.")
		if pressure < 0:
			raise InputError("Pressure must be positive.")
		if X_fluid < 0 or X_fluid > 1:
			raise InputError("X_fluid must have a value between 0 and 1.")

		if pressure == 0:
			return 0

		if hydrous_coeffs == True:
			if X_fluid==0:
				return 0
			H2O = root_scalar(self.root_dissolved_volatiles,args=(pressure,temperature,sample,X_fluid,hydrous_coeffs,kwargs),
								x0=1.0,x1=2.0).root
			return H2O
		else:
			a = 0.54
			b = 1.24
			B = -2.95
			C = 0.02

			fugacity = self.fugacity_model.fugacity(pressure=pressure,X_fluid=X_fluid,temperature=temperature-273.15,**kwargs)
			if fugacity == 0:
				return 0
			NBO_O = self.NBO_O(sample=sample,hydrous_coeffs=False)

			H2O = np.exp(a*np.log(fugacity) + b*NBO_O + B + C*pressure/temperature)

			return H2O


	def calculate_equilibrium_fluid_comp(self,pressure,temperature,sample,**kwargs):
		""" Returns 1.0 if a pure H2O fluid is saturated.
		Returns 0.0 if a pure H2O fluid is undersaturated.

		Parameters
		----------
		pressure     float
			The total pressure of the system in bars.
		temperature     float
			The temperature of the system in C.
		sample         pandas Series or dict
			Major element oxides in wt% (including H2O).

		Returns
		-------
		float
			1.0 if H2O-fluid saturated, 0.0 otherwise.
		"""

		if pressure > self.calculate_saturation_pressure(temperature=temperature,sample=sample,**kwargs):
			return 0.0
		else:
			return 1.0

	def calculate_saturation_pressure(self,temperature,sample,**kwargs):
		"""
		Calculates the pressure at which a pure H2O fluid is saturated, for the given sample
		composition and H2O concentration. Calls the scipy.root_scalar routine, which makes
		repeated called to the calculate_dissolved_volatiles method.

		Parameters
		----------
		temperature     float
			The temperature of the system in C.
		sample         pandas Series or dict
			Major element oxides in wt% (including H2O).
		X_fluid     float
			The mole fraction of H2O in the fluid. Default is 1.0.

		Returns
		-------
		float
			Calculated saturation pressure in bars.
		"""

		if type(sample) != dict and type(sample) != pd.core.series.Series:
			raise InputError("sample must be a dict or a pandas Series.")
		if 'H2O' not in sample:
			raise InputError("sample must contain H2O.")
		if sample['H2O'] < 0.0:
			raise InputError("Dissolved H2O must be greater than 0 wt%.")

		try:
			satP = root_scalar(self.root_saturation_pressure,args=(temperature,sample,kwargs),
								bracket=[1e-15,1e5]).root
		except:
			w.warn("Saturation pressure not found.",RuntimeWarning,stacklevel=2)
			satP = np.nan
		return satP

	def root_saturation_pressure(self,pressure,temperature,sample,kwargs):
		""" Function called by scipy.root_scalar when finding the saturation pressure using
		calculate_saturation_pressure.

		Parameters
		----------
		pressure     float
			Pressure guess in bars
		temperature     float
			The temperature of the system in C.
		sample         pandas Series or dict
			Major elements in wt% (normalized to 100%), including H2O.
		kwargs         dictionary
			Additional keyword arguments supplied to calculate_saturation_pressure. Might be required for
			the fugacity or activity models.

		Returns
		-------
		float
			The differece between the dissolved H2O at the pressure guessed, and the H2O concentration
			passed in the sample variable.
		"""

		return sample['H2O'] - self.calculate_dissolved_volatiles(pressure=pressure,temperature=temperature,sample=sample,**kwargs)


	def root_dissolved_volatiles(self,h2o,pressure,temperature,sample,X_fluid,webapp_coeffs,kwargs):
		""" Function called by calculate_dissolved_volatiles method when the hydrous parameterization is
		being used.

		Parameters
		----------
		h2o     float
			Guess for the H2O concentration in wt%.
		pressure     float
			Total pressure in bars.
		temperature     float
			Temperature in K.
		sample         pandas Series or dict
			Major element oxides in wt%.
		X_fluid     float
			Mole fraction of H2O in the fluid.
		kwargs     dictionary
			Keyword arguments

		Returns
		-------
		float
			Difference between H2O guessed and the H2O calculated.
		"""

		if webapp_coeffs == False:
			a = 0.53
			b = 2.35
			B = -3.37
			C = -0.02
		else:
			a = 0.52096846
			b = 2.11575907
			B = -3.24443335
			C = -0.02238884

		sample_copy = sample.copy()

		sample_copy['H2O'] = h2o

		NBO_O = self.NBO_O(sample=sample_copy,hydrous_coeffs=True)
		fugacity = self.fugacity_model.fugacity(pressure=pressure,X_fluid=X_fluid,temperature=temperature,**kwargs)

		return h2o - np.exp(a*np.log(fugacity) + b*NBO_O + B + C*pressure/(temperature+273.15))

	def NBO_O(self,sample,hydrous_coeffs=True):
		"""
		Calculates NBO/O according to Appendix A.1. of Iacono-Marziano et al. (2012). NBO/O
		is calculated on either a hydrous or anhyrous basis, as set when initialising the
		Model class.

		Parameters
		----------
		sample     pandas Series or dict
			Major element oxides in wt% (including H2O if using the hydrous parameterization).

		Returns
		-------
		float
			NBO/O.
		"""
		if type(sample) != dict and type(sample) != pd.core.series.Series:
			raise InputError("sample must be a dict or a pandas Series,")
		if all(ox in sample for ox in ['K2O','Na2O','CaO','MgO','FeO','Al2O3','SiO2','TiO2']) == False:
			raise InputError("sample must contain K2O, Na2O, CaO, MgO, FeO, Al2O3, SiO2, and TiO2.")

		X = wtpercentOxides_to_molOxides(sample)

		if 'Fe2O3' in X:
			Fe2O3 = X['Fe2O3']
		else:
			Fe2O3 = 0

		NBO = 2*(X['K2O']+X['Na2O']+X['CaO']+X['MgO']+X['FeO']+2*Fe2O3-X['Al2O3'])
		O = 2*X['SiO2']+2*X['TiO2']+3*X['Al2O3']+X['MgO']+X['FeO']+2*Fe2O3+X['CaO']+X['Na2O']+X['K2O']

		if hydrous_coeffs == True:
			if 'H2O' not in X:
				raise InputError("sample must contain H2O if using the hydrous parameterization.")
			NBO = NBO + 2*X['H2O']
			O = O + X['H2O']

		return NBO/O


class IaconoMarzianoCarbon(Model):
	"""
	Implementation of the Iacono-Marziano et al. (2012) carbon solubility model, as a Model class. Two
	calibrations are provided- the one incorporating the H2O content as a parameter (hydrous), and the
	one that does not (anhydrous). Specify which should be used when initialising the model, with the
	bool variable hydrous.
	"""

	def __init__(self):
		"""
		Initialise the model.

		Parameters
		----------
		hydrous     bool
			Whether to use the hydrous parameterization, or not.
		"""
		self.set_volatile_species(['CO2'])
		self.set_fugacity_model(fugacity_idealgas())
		self.set_activity_model(activity_idealsolution())
		self.set_calibration_ranges([CalibrationRange('temperature',[1000,1400],crf_Between,'oC','IaconoMarzianoCarbon',
									 				  fail_msg=crmsg_BC_IMT, pass_msg=crmsg_Between_pass, description_msg=crmsg_Between_description),
									 CalibrationRange('pressure',[100,10000],crf_Between,'bars','IaconoMarzianoCarbon',
									 				  fail_msg=crmsg_BC_IMP, pass_msg=crmsg_Between_pass, description_msg=crmsg_Between_description)])
		self.set_solubility_dependence(True)

	def preprocess_sample(self,sample):
		"""
		Returns sample, normalized to 100 wt%, without changing the wt% of H2O and CO2 if the
		hydrous parameterization is being used (default). If the anhydrous parameterization is
		used, it will normalize without including H2O and CO2.

		Parameters
		----------
		sample     pandas Series or dict
			Major element oxides in wt%.

		Returns
		-------
		pandas Series or dict
			Major element oxides normalized to wt%.
		"""
		if self.hydrous == True:
			return normalize_FixedVolatiles(sample)
		else:
			return normalize_AdditionalVolatiles(sample)

	def calculate_dissolved_volatiles(self,pressure,temperature,sample,X_fluid=1,
									  hydrous_coeffs=True, **kwargs):
		"""
		Calculates the dissolved CO2 concentration, using Eq (12) of Iacono-Marziano et al. (2012).
		If using the hydrous parameterization, it will use the scipy.root_scalar routine to find the
		root of the root_dissolved_volatiles method.

		Parameters
		----------
		pressure    float
			Total pressure in bars.
		temperature     float
			Temperature in C
		sample     pandas Series or dict
			Major element oxides in wt%.
		X_fluid      float
			Mole fraction of H2O in the fluid. Default is 1.0.
		hydrous_coeffs 	bool
			Use the hydrous or anhydrous NBO/O paramterisation (True for hydrous). Default is True.

		Returns
		-------
		float
			Dissolved H2O concentration in wt%.
		"""
		temperature = temperature + 273.15 #translate T from C to K

		if type(sample) != dict and type(sample) != pd.core.series.Series:
			raise InputError("sample must be a dict or a pandas Series.")
		if pressure < 0:
			raise InputError("Pressure must be positive.")
		if temperature <= 0:
			raise InputError("Temperature must be greater than 0K.")
		if X_fluid < 0 or X_fluid > 1:
			raise InputError("X_fluid must have a value between 0 and 1.")

		if pressure == 0:
			return 0

		if hydrous_coeffs == True:
			if 'H2O' not in sample:
				raise InputError("sample must contain H2O if using the hydrous parameterization.")
			if sample['H2O'] < 0:
				raise InputError("Dissolved H2O must be positive.")

			im_h2o_model = IaconoMarzianoWater()
			h2o = im_h2o_model.calculate_dissolved_volatiles(pressure=pressure,temperature=temperature-273.15,
														sample=sample,X_fluid=1-X_fluid,**kwargs)
			sample_h2o = sample.copy()
			sample_h2o['H2O'] = h2o

			d = np.array([-16.4,4.4,-17.1,22.8])
			a = 1.0
			b = 17.3
			B = -6.0
			C = 0.12

			NBO_O = self.NBO_O(sample=sample_h2o,hydrous_coeffs=True)

			molarProps = wtpercentOxides_to_molOxides(sample_h2o)

		else:
			d = np.array([2.3,3.8,-16.3,20.1])
			a = 1.0
			b = 15.8
			B = -5.3
			C = 0.14

			NBO_O = self.NBO_O(sample=sample,hydrous_coeffs=False)

			molarProps = wtpercentOxides_to_molOxides(sample)

		fugacity = self.fugacity_model.fugacity(pressure=pressure,X_fluid=X_fluid,temperature=temperature-273.15,**kwargs)

		if fugacity == 0:
			return 0

		if all(ox in molarProps for ox in ['Al2O3','CaO','K2O','Na2O','FeO','MgO','Na2O','K2O']) == False:
			raise InputError("sample must contain Al2O3, CaO, K2O, Na2O, FeO, MgO, Na2O, and K2O.")
		if 'Fe2O3' in molarProps:
			Fe2O3 = molarProps['Fe2O3']
		else:
			Fe2O3 = 0

		x = list()
		if 'H2O' in molarProps:
			x.append(molarProps['H2O'])
		else:
			x.append(0.0)
		x.append(molarProps['Al2O3']/(molarProps['CaO']+molarProps['K2O']+molarProps['Na2O']))
		x.append((molarProps['FeO']+Fe2O3*2+molarProps['MgO']))
		x.append((molarProps['Na2O']+molarProps['K2O']))
		x = np.array(x)

		CO3 = np.exp(np.sum(x*d) + a*np.log(fugacity) + b*NBO_O + B + C*pressure/temperature)

		CO2 = CO3/1e4#/(12+16*3)*(12+16*2)/1e4

		return CO2


	def calculate_equilibrium_fluid_comp(self,pressure,temperature,sample,**kwargs):
		""" Returns 1.0 if a pure CO2 fluid is saturated.
		Returns 0.0 if a pure CO2 fluid is undersaturated.

		Parameters
		----------
		pressure     float
			The total pressure of the system in bars.
		temperature     float
			The temperature of the system in C.
		sample         pandas Series or dict
			Major element oxides in wt% (including H2O).

		Returns
		-------
		float
			1.0 if CO2-fluid saturated, 0.0 otherwise.
		"""

		if pressure > self.calculate_saturation_pressure(temperature=temperature,sample=sample,**kwargs):
			return 0.0
		else:
			return 1.0

	def calculate_saturation_pressure(self,temperature,sample,**kwargs):
		"""
		Calculates the pressure at which a pure CO2 fluid is saturated, for the given sample
		composition and CO2 concentration. Calls the scipy.root_scalar routine, which makes
		repeated called to the calculate_dissolved_volatiles method.

		Parameters
		----------
		temperature     float
			The temperature of the system in C.
		sample         pandas Series or dict
			Major element oxides in wt% (including CO2).

		Returns
		-------
		float
			Calculated saturation pressure in bars.
		"""

		if temperature <= 0:
			raise InputError("Temperature must be greater than 0K.")
		if type(sample) != dict and type(sample) != pd.core.series.Series:
			raise InputError("sample must be a dict or a pandas Series.")
		if 'CO2' not in sample:
			raise InputError("sample must contain CO2")
		if sample['CO2'] < 0:
			raise InputError("Dissolved CO2 must be greater than 0 wt%.")

		try:
			satP = root_scalar(self.root_saturation_pressure,args=(temperature,sample,kwargs),
								bracket=[1e-15,1e5]).root
		except:
			w.warn("Saturation pressure not found.",RuntimeWarning,stacklevel=2)
			satP = np.nan
		return satP

	def root_saturation_pressure(self,pressure,temperature,sample,kwargs):
		""" Function called by scipy.root_scalar when finding the saturation pressure using
		calculate_saturation_pressure.

		Parameters
		----------
		pressure     float
			Pressure guess in bars
		temperature     float
			The temperature of the system in C.
		sample         pandas Series or dict
			Major elements in wt% (normalized to 100%), including CO2.
		kwargs         dictionary
			Additional keyword arguments supplied to calculate_saturation_pressure. Might be required for
			the fugacity or activity models.

		Returns
		-------
		float
			The differece between the dissolved CO2 at the pressure guessed, and the CO2 concentration
			passed in the sample variable.
		"""

		return sample['CO2'] - self.calculate_dissolved_volatiles(pressure=pressure,temperature=temperature,sample=sample,**kwargs)


	def NBO_O(self,sample,hydrous_coeffs=True):
		"""
		Calculates NBO/O according to Appendix A.1. of Iacono-Marziano et al. (2012). NBO/O
		is calculated on either a hydrous or anhyrous basis, as set when initialising the
		Model class.

		Parameters
		----------
		sample     pandas Series or dict
			Major element oxides in wt% (including H2O if using the hydrous parameterization).

		Returns
		-------
		float
			NBO/O.
		"""
		if type(sample) != dict and type(sample) != pd.core.series.Series:
			raise InputError("sample must be a dict or a pandas Series,")
		if all(ox in sample for ox in ['K2O','Na2O','CaO','MgO','FeO','Al2O3','SiO2','TiO2']) == False:
			raise InputError("sample must contain K2O, Na2O, CaO, MgO, FeO, Al2O3, SiO2, and TiO2.")

		X = wtpercentOxides_to_molOxides(sample)

		if 'Fe2O3' in X:
			Fe2O3 = X['Fe2O3']
		else:
			Fe2O3 = 0

		NBO = 2*(X['K2O']+X['Na2O']+X['CaO']+X['MgO']+X['FeO']+2*Fe2O3-X['Al2O3'])
		O = 2*X['SiO2']+2*X['TiO2']+3*X['Al2O3']+X['MgO']+X['FeO']+2*Fe2O3+X['CaO']+X['Na2O']+X['K2O']

		if hydrous_coeffs == True:
			if 'H2O' not in X:
				raise InputError("sample must contain H2O if using the hydrous parameterization.")
			NBO = NBO + 2*X['H2O']
			O = O + X['H2O']

		return NBO/O



class EguchiCarbon(Model):
	"""
	Implementation of the Eguchi and Dasgupta (2018) CO2 solubility model for andesitic melts.
	Uses the Zhang and Duan (2009) CO2 EOS for fugacity calculations, assuming a pure CO2 fluid,
	or ideal mixing for mixed fluids.
	"""

	def __init__(self):
		w.warn("Eguchi model is not working correctly. Do not use any results calculated.")
		self.set_volatile_species(['CO2'])
		self.set_fugacity_model(fugacity_ZD09_co2())
		self.set_activity_model(activity_idealsolution())
		self.set_solubility_dependence(False)
		self.set_calibration_ranges([CalibrationRange('pressure',[500.0,50000.0],crf_Between,'bar','Eguchi & Dasgupta (2018) carbon',
													  fail_msg=crmsg_Between_fail, pass_msg=crmsg_Between_pass, description_msg=crmsg_Between_description),
									 CalibrationRange('temperature',[950.0,1600],crf_Between,'oC','Eguchi & Dasgupta (2018) carbon',
									 				  fail_msg=crmsg_Between_fail, pass_msg=crmsg_Between_pass, description_msg=crmsg_Between_description)])

	def preprocess_sample(self,sample,ferric_total=0.15):
		""" Returns normalized sample composition, with ferric iron. Where a sample
		already contains ferric iron, the composition will be normalized to 100 wt%
		(excluding H2O and CO2). Where a sample contains only FeO, ferric iron will
		be calculated using the ferric/total iron ratio provided.

		Parameters
		----------
		sample     pandas Series or dict
			Major element oxides in wt%.
		ferric_total    float
			Mole ratio of ferric to total iron to be used
			for calculating Fe2O3 and FeO when only FeO is
			provided. Default is 0.15.

		Returns
		-------
		pandas Series or dict
			Normalized major element oxides in wt%.
		"""
		if type(sample) != dict and type(sample) != pd.core.series.Series:
			raise InputError("sample must be a dict or a pandas Series.")
		if 'FeO' not in sample:
			raise InputError("sample must contain FeO.")

		_sample = sample.copy()

		for ox in ['TiO2','P2O5']:
			if ox not in _sample:
				_sample[ox] = 0.0

		if 'Fe2O3' not in _sample:
			Fe_t = _sample['FeO']/oxideMass['FeO']
			Fe3 = ferric_total*Fe_t
			Fe2 = Fe_t - Fe3
			_sample['FeO'] = Fe2*oxideMass['FeO']
			_sample['Fe2O3'] = Fe3*oxideMass['Fe2O3']/2

		return normalize_AdditionalVolatiles(_sample)

	def calculate_dissolved_volatiles(self,pressure,temperature,sample,X_fluid=1.0,**kwargs):
		"""
		Calculates the dissolved (total) CO2 using eqs (9) and (10) of Eguchi and Dasgupta (2018).

		Parameters
		----------
		pressure     float
			Pressure in bars
		temperature     float
			Temperature in C
		sample     pandas Series or dict
			Major element oxides in wt%.
		X_fluid     float
			The mole fraction of CO2 in the fluid.

		Returns
		-------
		float
			Dissolved CO2 concentration.
		"""
		if pressure < 0:
			raise InputError("Pressure must be greater than 0 bar.")

		if pressure == 0:
			return 0

		XCO3 = self.Xi_melt(pressure=pressure,temperature=temperature,sample=sample,X_fluid=X_fluid,species='CO3')
		XCO2 = self.Xi_melt(pressure=pressure,temperature=temperature,sample=sample,X_fluid=X_fluid,species='CO2')

		FW_one = wtpercentOxides_to_formulaWeight(sample)

		CO2_CO2 = ((44.01*XCO2)/(44.01*XCO2+(1-(XCO2+XCO3))*FW_one))*100
		CO2_CO3 = ((44.01*XCO3)/(44.01*XCO3+(1-(XCO2+XCO3))*FW_one))*100

		return CO2_CO2 + CO2_CO3


	def calculate_equilibrium_fluid_comp(self,pressure,temperature,sample,**kwargs):
		""" Returns 1.0 if a pure CO2 fluid is saturated.
		Returns 0.0 if a pure CO2 fluid is undersaturated.

		Parameters
		----------
		pressure     float
			The total pressure of the system in bars.
		temperature     float
			The temperature of the system in C.
		sample         pandas Series or dict
			Major element oxides in wt% (including H2O).

		Returns
		-------
		float
			1.0 if CO2-fluid saturated, 0.0 otherwise.
		"""

		satP = self.calculate_saturation_pressure(temperature=temperature,sample=sample,X_fluid=1.0,**kwargs)
		if pressure < satP:
			return 1.0
		else:
			return 0.0

	def calculate_saturation_pressure(self,temperature,sample,X_fluid=1.0,**kwargs):
		"""
		Calculates the pressure at which a pure CO2 fluid is saturated, for the given sample
		composition and CO2 concentration. Calls the scipy.root_scalar routine, which makes
		repeated called to the calculate_dissolved_volatiles method.

		Parameters
		----------
		temperature     float
			The temperature of the system in C.
		sample         pandas Series or dict
			Major element oxides in wt% (including CO2).
		X_fluid     float
			The mole fraction of H2O in the fluid. Default is 1.0.

		Returns
		-------
		float
			Calculated saturation pressure in bars.
		"""

		if 'CO2' not in sample:
			raise InputError("sample must contain CO2.")
		if sample['CO2'] < 0.0:
			raise InputError("Concentration of CO2 must be greater than 0 wt%.")
		try:
			satP = root_scalar(self.root_saturation_pressure,x0=1000.0,x1=2000.0,
								args=(temperature,sample,X_fluid,kwargs)).root
		except:
			w.warn("Saturation pressure not found.",RuntimeWarning,stacklevel=2)
			satP = np.nan
		return satP

	def root_saturation_pressure(self,pressure,temperature,sample,X_fluid,kwargs):
		""" Function called by scipy.root_scalar when finding the saturation pressure using
		calculate_saturation_pressure.

		Parameters
		----------
		pressure     float
			Pressure guess in bars
		temperature     float
			The temperature of the system in C.
		sample         pandas Series or dict
			Major elements in wt% (normalized to 100%), including CO2.
		kwargs         dictionary
			Additional keyword arguments supplied to calculate_saturation_pressure. Might be required for
			the fugacity or activity models.

		Returns
		-------
		float
			The differece between the dissolved CO2 at the pressure guessed, and the CO2 concentration
			passed in the sample variable.
		"""
		return sample['CO2'] - self.calculate_dissolved_volatiles(pressure=pressure,temperature=temperature,sample=sample,X_fluid=X_fluid,**kwargs)

	def Xi_melt(self,pressure,temperature,sample,species,X_fluid=1.0,**kwargs):
		"""
		Calculates the mole fraction of dissolved molecular CO2 or carbonate CO3(2-), using
		eqn (9) of Eguchi and Dasgupta (2018).

		Parameters
		----------
		pressure    float
			Pressure in bars.
		temperature     float
			Temperature in C.
		sample         pandas Series or dict
			Major element oxides in wt%.
		species        str
			Which species to calculate, molecular CO2 'CO2' or carbonate ion 'CO3'.
		X_fluid     float
			The mole fraction of CO2 in the fluid. Default is 1.0.

		Returns
		-------
		float
			Mole fraction of selected species in the melt
		"""
		temperature = temperature + 273.15 #translate T from C to K

		if all(ox in sample for ox in ['MgO','CaO','FeO','Na2O','K2O','MnO','Al2O3','Fe2O3','SiO2','TiO2','P2O5']) == False:
			raise InputError("sample must contain MgO, CaO, FeO, Na2O, K2O, MnO, Al2O3, Fe2O3, SiO3, TiO2, and P2O5.")
		if X_fluid < 0 or X_fluid > 1:
			raise InputError("X_fluid must have a value between 0 and 1.")
		if pressure < 0:
			raise InputError("Pressure must be positive.")
		if temperature <= 0:
			raise InputError("Temperature must be greater than 0K.")

		if species == 'CO3':
			DH = -1.65e5
			DV = 2.38e-5
			DS = -43.64
			B = 1.47e3
			yNBO = 3.29
			A_CaO = 1.68e5
			A_Na2O = 1.76e5
			A_K2O = 2.11e5
		elif species == 'CO2':
			DH = -9.02e4
			DV = 1.92e-5
			DS = -43.08
			B = 1.12e3
			yNBO = -7.09
			A_CaO = 0
			A_Na2O = 0
			A_K2O = 0
		else:
			raise InputError("species variable must be either 'CO2' or 'CO3'.")
		R = 8.314

		# Calculate NBO term
		cations = wtpercentOxides_to_molSingleO(sample)
		oxides = wtpercentOxides_to_molOxides(sample)



		NM = (cations['Mg'] + cations['Ca'] + cations['Fe'] + cations['Na'] +
			cations['K'] + cations['Mn'])
		Al = cations['Al'] - NM
		if Al > 0:
			Al = NM
		else:
			Al = cations['Al']
		Fe = cations['Fe3'] + Al
		if Al > 0:
			Fe = 0
		if Al < 0 and Fe > 0:
			Fe = - Al
		if Al < 0 and Fe < 0:
			Fe = cations['Fe3']
		Tet = cations['Si'] + cations['Ti'] + cations['P'] + Al + Fe
		NBO = 2 - 4*Tet

		lnfCO2 = np.log(self.fugacity_model.fugacity(pressure=pressure,temperature=temperature-273.15,X_fluid=X_fluid))

		lnXi = ((DH/(R*temperature)-(pressure*1e5*DV)/(R*temperature)+DS/R) +
				(A_CaO*oxides['CaO']+A_Na2O*oxides['Na2O']+A_K2O*oxides['K2O'])/(R*temperature) +
				(B*lnfCO2/temperature) + yNBO*NBO
				)

		return np.exp(lnXi)


class MooreWater(Model):
	"""
	Implementation of the Moore et al. (1998) H2O solubility model for magmas up to 3,000 bars.
	"""

	def __init__(self):
		"""
		Initialize the model.
		"""
		self.set_volatile_species(['H2O'])
		self.set_fugacity_model(fugacity_HB_h2o())
		self.set_activity_model(activity_idealsolution())
		self.set_solubility_dependence(False)
		self.set_calibration_ranges([CalibrationRange('pressure',3000.0,crf_LessThan,'bar','Moore et al. (1998) water',
													  fail_msg=crmsg_LessThan_fail, pass_msg=crmsg_LessThan_pass, description_msg=crmsg_LessThan_description),
									 CalibrationRange('temperature',[700.0,1200],crf_Between,'oC','Moore et al. (1998) water',
									 				  fail_msg=crmsg_Between_fail, pass_msg=crmsg_Between_pass, description_msg=crmsg_Between_description)])


	def preprocess_sample(self, sample):
		"""
		Returns sample with extranneous (non oxide) information removed and any missing oxides given a value of 0.0.
		"""
		for oxide in oxides:
			if oxide in sample.keys():
				pass
			else:
				sample[oxide] = 0.0

		bulk_comp = {oxide:  sample[oxide] for oxide in oxides}
		self.bulk_comp_orig = sample

		return bulk_comp

	def calculate_dissolved_volatiles(self, sample, pressure, temperature, X_fluid=1.0, **kwargs):
		"""
		Parameters
		----------
		sample dict
			Composition of sample in wt% oxides.

		pressure float
			Pressure in bars.

		temperature float
			Temperature in degrees C.

		X_fluid float
			OPTIONAL. Default is 1.0. Mole fraction of H2O in the H2O-CO2 fluid.

		Returns
		-------
		float
			Calculated dissolved H2O concentration in wt%.
		"""

		_sample = sample.copy()
		_sample['H2O'] = 0.0
		_sample['CO2'] = 0.0

		_sample = normalize(_sample)

		fH2O = self.fugacity_model.fugacity(pressure=pressure,temperature=temperature,X_fluid=X_fluid,**kwargs)
		aParam = 2565.0
		bParam_Al2O3 = -1.997
		bParam_FeOt = -0.9275
		bParam_Na2O = 2.736
		cParam = 1.171
		dParam = -14.21

		temperatureK = temperature + 273.15

		sample_molfrac = wtpercentOxides_to_molOxides(_sample)
		FeOtot = sample_molfrac['FeO'] + sample_molfrac['Fe2O3']*0.8998

		b_x_sum = (bParam_Al2O3 * sample_molfrac['Al2O3']) + (bParam_FeOt * FeOtot) + (bParam_Na2O * sample_molfrac['Na2O'])
		two_ln_XH2Omelt = (aParam / temperatureK) + b_x_sum * (pressure/temperatureK) + cParam * np.log(fH2O) + dParam
		ln_XH2Omelt = two_ln_XH2Omelt / 2.0
		XH2Omelt = np.exp(ln_XH2Omelt)
		sample_molfrac['H2O'] = XH2Omelt

		#Normalize mol fractions to sum to 1, while preserving XH2O
		for key, value in sample_molfrac.items():
			if key != 'H2O':
				sample_molfrac.update({key: value/((1/(1-sample_molfrac['H2O'])))})

		sample_wtper = mol_to_wtpercent(sample_molfrac)

		return sample_wtper['H2O']

	def calculate_equilibrium_fluid_comp(self, sample, pressure, temperature, **kwargs):
		"""
		Parameters
		----------
		sample dict
			Composition of sample in wt% oxides.

		pressure float
			Pressure in bars.

		temperature float
			Temperature in degrees C.

		Returns
		-------
		float
			Calculated equilibrium fluid concentration in XH2Ofluid mole fraction.
		"""
		_sample = sample.copy()
		sample_anhy = sample.copy()
		sample_anhy["H2O"] = 0.0
		sample_anhy["CO2"] = 0.0

		aParam = 2565.0
		bParam_Al2O3 = -1.997
		bParam_FeOt = -0.9275
		bParam_Na2O = 2.736
		cParam = 1.171
		dParam = -14.21

		temperatureK = temperature + 273.15

		sample_molfrac_anhy = wtpercentOxides_to_molOxides(sample_anhy)
		sample_molfrac_hy = wtpercentOxides_to_molOxides(_sample)
		FeOtot = sample_molfrac_anhy['FeO'] + sample_molfrac_anhy['Fe2O3']*0.8998

		b_x_sum = (bParam_Al2O3 * sample_molfrac_anhy['Al2O3']) + (bParam_FeOt * FeOtot) + (bParam_Na2O * sample_molfrac_anhy['Na2O'])
		ln_fH2O = (2 * np.log(sample_molfrac_hy['H2O']) - (aParam/temperatureK) - b_x_sum * (pressure/temperatureK) - dParam) / cParam
		fH2O = np.exp(ln_fH2O)
		XH2O_fl = fH2O / pressure

		# SM: I've changed this to return X_H2O only, as otherwise it doesn't conform to other single-volatile
		# models. I'm not sure this is the best solution though.
		# return (XCO2_fl, XH2O_fl)
		return XH2O_fl

	def calculate_saturation_pressure(self,temperature,sample,X_fluid=1.0,**kwargs):
		"""
		Calculates the pressure at which a an H2O-bearing fluid is saturated. Calls the scipy.root_scalar
		routine, which makes repeated called to the calculate_dissolved_volatiles method.

		Parameters
		----------
		sample dict
			Composition of sample in wt% oxides.

		temperature float
			Temperature in degrees C.

		X_fluid float
			OPTIONAL. Default is 1.0. Mole fraction of H2O in the H2O-CO2 fluid.

		Returns
		-------
		float
			Calculated saturation pressure in bars.
		"""
		_sample = sample.copy()

		temperatureK = temperature + 273.15
		if temperatureK <= 0.0:
			raise InputError("Temperature must be greater than 0K.")
		if X_fluid < 0 or X_fluid > 1:
			raise InputError("X_fluid must have a value between 0 and 1.")
		if type(sample) != dict and type(sample) != pd.core.series.Series:
			raise InputError("sample must be a dict or a pandas Series.")
		if 'H2O' not in sample:
			raise InputError("sample must contain H2O.")
		if sample['H2O'] < 0.0:
			raise InputError("Dissolved H2O concentration must be greater than 0 wt%.")

		try:
			satP = root_scalar(self.root_saturation_pressure,args=(temperature,_sample,X_fluid,kwargs),
								x0=100.0,x1=2000.0).root
		except:
			w.warn("Saturation pressure not found.",RuntimeWarning,stacklevel=2)
			satP = np.nan
		return np.real(satP)

	def root_saturation_pressure(self,pressure,temperature,sample,X_fluid,kwargs):
		""" Function called by scipy.root_scalar when finding the saturation pressure using
		calculate_saturation_pressure.

		Parameters
		----------
		pressure     float
			Pressure guess in bars
		temperature     float
			The temperature of the system in C.
		sample         pandas Series or dict
			Major elements in wt% (normalized to 100%), including H2O.
		kwargs         dictionary
			Additional keyword arguments supplied to calculate_saturation_pressure. Might be required for
			the fugacity or activity models.

		Returns
		-------
		float
			The differece between the dissolved H2O at the pressure guessed, and the H2O concentration
			passed in the sample variable.
		"""
		return self.calculate_dissolved_volatiles(pressure=pressure,temperature=temperature,sample=sample,X_fluid=X_fluid,**kwargs) - sample['H2O']

class LiuWater(Model):
	"""
	Implementation of the Liu et al. (2005) H2O solubility model for metaluminous high-silica rhyolitic melts.
	"""

	def __init__(self):
		"""
		Initialize the model.
		"""
		self.set_volatile_species(['H2O'])
		self.set_fugacity_model(fugacity_idealgas())
		self.set_activity_model(activity_idealsolution())
		self.set_solubility_dependence(False)
		self.set_calibration_ranges([CalibrationRange('pressure',[0,5000.0],crf_Between,'bar','Liu et al. (2005) water',
													  fail_msg=crmsg_Between_fail, pass_msg=crmsg_Between_pass, description_msg=crmsg_Between_description),
									 CalibrationRange('temperature',[700.0,1200],crf_Between,'oC','Liu et al. (2005) water',
									 				  fail_msg=crmsg_Between_fail, pass_msg=crmsg_Between_pass, description_msg=crmsg_Between_description),
									 CalibrationRange('SiO2',[LiuWater_SiO2Min,LiuWater_SiO2Max],crf_Between,'wt%','Liu et al. (2005) water',
													  fail_msg=crmsg_BC_fail, pass_msg=crmsg_BC_pass, description_msg=crmsg_Between_description),
									 CalibrationRange('TiO2',[LiuWater_TiO2Min,LiuWater_TiO2Max],crf_Between,'wt%','Liu et al. (2005) water',
													  fail_msg=crmsg_BC_fail, pass_msg=crmsg_BC_pass, description_msg=crmsg_Between_description),
									 CalibrationRange('Al2O3',[LiuWater_Al2O3Min,LiuWater_Al2O3Max],crf_Between,'wt%','Liu et al. (2005) water',
													  fail_msg=crmsg_BC_fail, pass_msg=crmsg_BC_pass, description_msg=crmsg_Between_description),
									 CalibrationRange('FeO',[LiuWater_FeOMin,LiuWater_FeOMax],crf_Between,'wt%','Liu et al. (2005) water',
													  fail_msg=crmsg_BC_fail, pass_msg=crmsg_BC_pass, description_msg=crmsg_Between_description),
									 CalibrationRange('MgO',[LiuWater_MgOMin,LiuWater_MgOMax],crf_Between,'wt%','Liu et al. (2005) water',
													  fail_msg=crmsg_BC_fail, pass_msg=crmsg_BC_pass, description_msg=crmsg_Between_description),
									 CalibrationRange('CaO',[LiuWater_CaOMin,LiuWater_CaOMax],crf_Between,'wt%','Liu et al. (2005) water',
													  fail_msg=crmsg_BC_fail, pass_msg=crmsg_BC_pass, description_msg=crmsg_Between_description),
									 CalibrationRange('Na2O',[LiuWater_Na2OMin,LiuWater_Na2OMax],crf_Between,'wt%','Liu et al. (2005) water',
													  fail_msg=crmsg_BC_fail, pass_msg=crmsg_BC_pass, description_msg=crmsg_Between_description),
									 CalibrationRange('K2O',[LiuWater_K2OMin,LiuWater_K2OMax],crf_Between,'wt%','Liu et al. (2005) water',
													  fail_msg=crmsg_BC_fail, pass_msg=crmsg_BC_pass, description_msg=crmsg_Between_description),
									 CalibrationRange('sample',None,crf_LiuWaterComp,None,None,
									 				  fail_msg=crmsg_LiuWaterComp_fail, pass_msg=crmsg_LiuComp_pass, description_msg=crmsg_LiuComp_description)])


	def preprocess_sample(self, sample):
		"""
		Returns sample with extranneous (non oxide) information removed and any missing oxides given a value of 0.0.
		"""
		for oxide in oxides:
			if oxide in sample.keys():
				pass
			else:
				sample[oxide] = 0.0

		bulk_comp = {oxide:  sample[oxide] for oxide in oxides}
		self.bulk_comp_orig = sample

		return bulk_comp

	def calculate_dissolved_volatiles(self, sample, pressure, temperature, X_fluid=1.0, **kwargs):
		"""
		Parameters
		----------
		sample dict
			Composition of sample in wt% oxides.

		pressure float
			Pressure in bars.

		temperature float
			Temperature in degrees C.

		X_fluid float
			OPTIONAL. Default is 1.0. Mole fraction of H2O in the H2O-CO2 fluid.

		Returns
		-------
		float
			Calculated dissolved H2O concentration in wt%.
		"""
		pressureMPa = pressure / 10.0
		Pw = pressureMPa * X_fluid
		PCO2 = pressureMPa * (1 - X_fluid)

		temperatureK = temperature + 273.15

		H2Ot = ((354.94*Pw**(0.5) + 9.623*Pw - 1.5223*Pw**(1.5)) / temperatureK +
				0.0012439*Pw**(1.5) + PCO2*(-1.084*10**(-4)*Pw**(0.5) - 1.362*10**(-5)*Pw))

		return H2Ot

	def calculate_equilibrium_fluid_comp(self, sample, pressure, temperature, **kwargs):
		"""
		Parameters
		----------
		sample dict
			Composition of sample in wt% oxides.

		pressure float
			Pressure in bars.

		temperature float
			Temperature in degrees C.

		Returns
		-------
		float
			Calculated equilibrium fluid concentration in XH2Ofluid mole fraction.
		"""
		temperatureK = temperature + 273.15
		pressureMPa = pressure / 10.0

		_sample = sample.copy()
		H2Ot = _sample["H2O"]

		#calculate saturation pressure and assert that input P <= SatP
		satP = self.calculate_saturation_pressure(temperature,sample)
		is_saturated = satP - pressure
		if is_saturated >= 0:
			pass
		else:
			w.warn("{:.1f} bars is above the saturation pressure ({:.1f} bars) for this sample. Results from this calculation may be nonsensical.".format(pressure,satP))

		#Use sympy to solve solubility equation for XH2Ofluid
		XH2Ofluid = sympy.symbols('XH2Ofluid') #XH2Ofluid is the variable to solve for

		equation = ((354.94*(XH2Ofluid*pressureMPa)**(0.5) + 9.623*(XH2Ofluid*pressureMPa)
						- 1.5223*(XH2Ofluid*pressureMPa)**(1.5)) / temperatureK
						+ 0.0012439*(XH2Ofluid*pressureMPa)**(1.5)
						+ pressureMPa*(1-XH2Ofluid)*(-1.084*10**(-4)*(XH2Ofluid*pressureMPa)**(0.5)
						- 1.362*10**(-5)*(XH2Ofluid*pressureMPa)) - H2Ot)

		XH2Ofluid = sympy.solve(equation, XH2Ofluid)[0]
		if XH2Ofluid > 1:
			XH2Ofluid = 1
		if XH2Ofluid < 0:
			XH2Ofluid = 0

		return XH2Ofluid

	def calculate_saturation_pressure(self,temperature,sample,X_fluid=1.0,**kwargs):
		"""
		Calculates the pressure at which a an H2O-bearing fluid is saturated. Calls the scipy.root_scalar
		routine, which makes repeated called to the calculate_dissolved_volatiles method.

		Parameters
		----------
		sample dict
			Composition of sample in wt% oxides.

		temperature float
			Temperature in degrees C.

		X_fluid float
			OPTIONAL. Default is 1.0. Mole fraction of H2O in the H2O-CO2 fluid.

		Returns
		-------
		float
			Calculated saturation pressure in bars.
		"""
		_sample = sample.copy()

		temperatureK = temperature + 273.15
		if temperatureK <= 0.0:
			raise InputError("Temperature must be greater than 0K.")
		if X_fluid < 0 or X_fluid > 1:
			raise InputError("X_fluid must have a value between 0 and 1.")
		if type(sample) != dict and type(sample) != pd.core.series.Series:
			raise InputError("sample must be a dict or a pandas Series.")
		if 'H2O' not in sample:
			raise InputError("sample must contain H2O.")
		if sample['H2O'] < 0.0:
			raise InputError("Dissolved H2O concentration must be greater than 0 wt%.")

		try:
			satP = root_scalar(self.root_saturation_pressure,args=(temperature,_sample,X_fluid,kwargs),
								x0=10.0,x1=200.0).root
		except:
			w.warn("Saturation pressure not found.",RuntimeWarning,stacklevel=2)
			satP = np.nan
		return np.real(satP)

	def root_saturation_pressure(self,pressure,temperature,sample,X_fluid,kwargs):
		""" Function called by scipy.root_scalar when finding the saturation pressure using
		calculate_saturation_pressure.

		Parameters
		----------
		pressure     float
			Pressure guess in bars
		temperature     float
			The temperature of the system in C.
		sample         pandas Series or dict
			Major elements in wt% (normalized to 100%), including H2O.
		kwargs         dictionary
			Additional keyword arguments supplied to calculate_saturation_pressure. Might be required for
			the fugacity or activity models.

		Returns
		-------
		float
			The differece between the dissolved H2O at the pressure guessed, and the H2O concentration
			passed in the sample variable.
		"""
		return self.calculate_dissolved_volatiles(pressure=pressure,temperature=temperature,sample=sample,X_fluid=X_fluid,**kwargs) - sample['H2O']

class LiuCarbon(Model):
	"""
	Implementation of the Liu et al. (2005) H2O-CO2 solubility model for metaluminous high-silica rhyolitic melts.
	"""

	def __init__(self):
		"""
		Initialize the model.
		"""
		self.set_volatile_species(['CO2'])
		self.set_fugacity_model(fugacity_idealgas())
		self.set_activity_model(activity_idealsolution())
		self.set_solubility_dependence(False)
		self.set_calibration_ranges([CalibrationRange('pressure',[0,5000.0],crf_Between,'bar','Liu et al. (2005) Carbon',
													  fail_msg=crmsg_Between_fail, pass_msg=crmsg_Between_pass, description_msg=crmsg_Between_description),
									 CalibrationRange('temperature',[700.0,1200],crf_Between,'oC','Liu et al. (2005) Carbon',
									 				  fail_msg=crmsg_Between_fail, pass_msg=crmsg_Between_pass, description_msg=crmsg_Between_description),
									 CalibrationRange('SiO2',[LiuCarbon_SiO2Min,LiuCarbon_SiO2Max],crf_Between,'wt%','Liu et al. (2005) Carbon',
													  fail_msg=crmsg_BC_fail, pass_msg=crmsg_BC_pass, description_msg=crmsg_Between_description),
									 CalibrationRange('TiO2',[LiuCarbon_TiO2Min,LiuCarbon_TiO2Max],crf_Between,'wt%','Liu et al. (2005) Carbon',
													  fail_msg=crmsg_BC_fail, pass_msg=crmsg_BC_pass, description_msg=crmsg_Between_description),
									 CalibrationRange('Al2O3',[LiuCarbon_Al2O3Min,LiuCarbon_Al2O3Max],crf_Between,'wt%','Liu et al. (2005) Carbon',
													  fail_msg=crmsg_BC_fail, pass_msg=crmsg_BC_pass, description_msg=crmsg_Between_description),
									 CalibrationRange('FeO',[LiuCarbon_FeOMin,LiuCarbon_FeOMax],crf_Between,'wt%','Liu et al. (2005) Carbon',
													  fail_msg=crmsg_BC_fail, pass_msg=crmsg_BC_pass, description_msg=crmsg_Between_description),
									 CalibrationRange('MgO',[LiuCarbon_MgOMin,LiuCarbon_MgOMax],crf_Between,'wt%','Liu et al. (2005) Carbon',
													  fail_msg=crmsg_BC_fail, pass_msg=crmsg_BC_pass, description_msg=crmsg_Between_description),
									 CalibrationRange('CaO',[LiuCarbon_CaOMin,LiuCarbon_CaOMax],crf_Between,'wt%','Liu et al. (2005) Carbon',
													  fail_msg=crmsg_BC_fail, pass_msg=crmsg_BC_pass, description_msg=crmsg_Between_description),
									 CalibrationRange('Na2O',[LiuCarbon_Na2OMin,LiuCarbon_Na2OMax],crf_Between,'wt%','Liu et al. (2005) Carbon',
													  fail_msg=crmsg_BC_fail, pass_msg=crmsg_BC_pass, description_msg=crmsg_Between_description),
									 CalibrationRange('K2O',[LiuCarbon_K2OMin,LiuCarbon_K2OMax],crf_Between,'wt%','Liu et al. (2005) Carbon',
													  fail_msg=crmsg_BC_fail, pass_msg=crmsg_BC_pass, description_msg=crmsg_Between_description),
									 CalibrationRange('sample',None,crf_LiuCarbonComp,None,None,
									 				  fail_msg=crmsg_LiuCarbonComp_fail, pass_msg=crmsg_LiuComp_pass, description_msg=crmsg_LiuComp_description)])



	def preprocess_sample(self, sample):
		"""
		Returns sample with extranneous (non oxide) information removed and any missing oxides given a value of 0.0.
		"""
		for oxide in oxides:
			if oxide in sample.keys():
				pass
			else:
				sample[oxide] = 0.0

		bulk_comp = {oxide:  sample[oxide] for oxide in oxides}
		self.bulk_comp_orig = sample

		return bulk_comp

	def calculate_dissolved_volatiles(self, sample, pressure, temperature, X_fluid=1, **kwargs):
		"""
		Parameters
		----------
		sample dict
			Composition of sample in wt% oxides.

		pressure float
			Pressure in bars.

		temperature float
			Temperature in degrees C.

		X_fluid float
			OPTIONAL. Default is 1. Mole fraction of CO2 in the H2O-CO2 fluid.

		Returns
		-------
		float
			Calculated dissolved CO2 concentration in wt%.
		"""
		pressureMPa = pressure / 10.0
		Pw = pressureMPa * (1 - X_fluid)
		PCO2 = pressureMPa * X_fluid #(1 - X_fluid)

		temperatureK = temperature + 273.15

		CO2melt_ppm = (PCO2*(5668 - 55.99*Pw)/temperatureK
					+ PCO2*(0.4133*Pw**(0.5) + 2.041*10**(-3)*Pw**(1.5)))

		CO2melt = CO2melt_ppm / 10000

		return CO2melt

	def calculate_equilibrium_fluid_comp(self, sample, pressure, temperature, **kwargs):
		"""
		Parameters
		----------
		sample dict
			Composition of sample in wt% oxides.

		pressure float
			Pressure in bars.

		temperature float
			Temperature in degrees C.

		Returns
		-------
		float
			Calculated equilibrium fluid concentration in XCO2fluid mole fraction.
		"""
		temperatureK = temperature + 273.15
		pressureMPa = pressure / 10.0

		if temperatureK <= 0.0:
			raise InputError("Temperature must be greater than 0K.")
		if type(sample) != dict and type(sample) != pd.core.series.Series:
			raise InputError("sample must be a dict or a pandas Series.")
		if 'CO2' not in sample:
			raise InputError("sample must contain CO2.")
		if sample['CO2'] < 0.0:
			raise InputError("Dissolved CO2 concentration must be greater than 0 wt%.")

		_sample = sample.copy()
		CO2melt_wt = _sample["CO2"]
		CO2melt_ppm = CO2melt_wt * 10000

		#calculate saturation pressure and assert that input P <= SatP
		satP = self.calculate_saturation_pressure(temperature,sample)
		is_saturated = satP - pressure
		if is_saturated >= 0:
			pass
		else:
			w.warn(str(pressure) + " bars is above the saturation pressure (" + str(satP) + " bars) for this sample. Results from this calculation may be nonsensical.")

		#Use sympy to solve solubility equation for XH2Ofluid
		XCO2fluid = sympy.symbols('XCO2fluid') #XCO2fluid is the variable to solve for

		equation = (((XCO2fluid*pressureMPa)*(5668 - 55.99*(pressureMPa*(1-XCO2fluid)))/temperatureK
					+ (XCO2fluid*pressureMPa)*(0.4133*(pressureMPa*(1-XCO2fluid))**(0.5)
					+ 2.041*10**(-3)*(pressureMPa*(1-XCO2fluid))**(1.5))) - CO2melt_ppm)

		XCO2fluid = sympy.solve(equation, XCO2fluid)[0]
		if XCO2fluid > 1:
			XCO2fluid = 1
		if XCO2fluid < 0:
			XCO2fluid = 0

		return XCO2fluid #1 - XCO2fluid

	def calculate_saturation_pressure(self,temperature,sample,X_fluid=1.0,**kwargs):
		"""
		Calculates the pressure at which a an CO2-bearing fluid is saturated. Calls the scipy.root_scalar
		routine, which makes repeated called to the calculate_dissolved_volatiles method.

		Parameters
		----------
		sample dict
			Composition of sample in wt% oxides.

		temperature float
			Temperature in degrees C.

		X_fluid float
			OPTIONAL. Default is 0. Mole fraction of CO2 in the H2O-CO2 fluid.

		Returns
		-------
		float
			Calculated saturation pressure in bars.
		"""
		_sample = sample.copy()

		temperatureK = temperature + 273.15
		if temperatureK <= 0.0:
			raise InputError("Temperature must be greater than 0K.")
		if X_fluid < 0 or X_fluid > 1:
			raise InputError("X_fluid must have a value between 0 and 1.")
		if type(sample) != dict and type(sample) != pd.core.series.Series:
			raise InputError("sample must be a dict or a pandas Series.")
		if 'CO2' not in sample:
			raise InputError("sample must contain CO2.")
		if sample['CO2'] < 0.0:
			raise InputError("Dissolved CO2 concentration must be greater than 0 wt%.")

		try:
			satP = root_scalar(self.root_saturation_pressure,args=(temperature,_sample,X_fluid,kwargs),
								x0=10.0,x1=2000.0).root
		except:
			w.warn("Saturation pressure not found.",RuntimeWarning,stacklevel=2)
			satP = np.nan
		return np.real(satP)

	def root_saturation_pressure(self,pressure,temperature,sample,X_fluid,kwargs):
		""" Function called by scipy.root_scalar when finding the saturation pressure using
		calculate_saturation_pressure.

		Parameters
		----------
		pressure     float
			Pressure guess in bars
		temperature     float
			The temperature of the system in C.
		sample         pandas Series or dict
			Major elements in wt% (normalized to 100%), including H2O.
		kwargs         dictionary
			Additional keyword arguments supplied to calculate_saturation_pressure. Might be required for
			the fugacity or activity models.

		Returns
		-------
		float
			The differece between the dissolved H2O at the pressure guessed, and the H2O concentration
			passed in the sample variable.
		"""
		return self.calculate_dissolved_volatiles(pressure=pressure,temperature=temperature,sample=sample,X_fluid=X_fluid,**kwargs) - sample['CO2']

class AllisonCarbon(Model):
	"""
	Implementation of the Allison et al. (2019) CO2 solubility model. Which type of fit, and
	which composition must be selected when the Model is initialized. The fit may be either
	thermodynamic or power-law. The composition may be chosen from sunset, sfvf, erebus, vesuvius,
	etna, or stromboli. Default is the power-law fit to sunset.
	"""

	def __init__(self,model_loc='sunset',model_fit='thermodynamic'):
		"""
		Initialize the model.

		Parameters
		----------
		model_fit     str
			Either 'power' for the power-law fits, or 'thermodynamic' for the
			thermodynamic fits.
		model_loc     str
			One of 'sunset', 'sfvf', 'erebus', 'vesuvius', 'etna', 'stromboli'.
		"""

		self.set_volatile_species(['CO2'])
		self.set_fugacity_model(fugacity_HB_co2())
		self.set_activity_model(activity_idealsolution())
		self.set_calibration_ranges([CalibrationRange('temperature',1200,crf_EqualTo,'oC','Allison et al. (2019) carbon',
						 							  fail_msg=crmsg_EqualTo_fail_AllisonTemp, pass_msg=crmsg_EqualTo_pass, description_msg=crmsg_EqualTo_description),
									 CalibrationRange('temperature',[1000,1400],crf_Between,'bar','Allison et al. (2019) carbon',
						 							   fail_msg=crmsg_Between_AllisonTemp, pass_msg=crmsg_Between_pass, description_msg=crmsg_Between_description)])
		self.set_solubility_dependence(False)
		self.model_loc = model_loc
		self.model_fit = model_fit


	def preprocess_sample(self,sample):
		"""
		Returns sample normalized to 100wt%, keeping the concentrations of H2O and CO2 constant.

		Parameters
		----------
		sample     pandas Series or dict
			Major element oxides in wt%.

		Returns
		-------
		pandas Series
			Normalized major element oxides in wt%.
		"""
		return normalize_AdditionalVolatiles(sample)

	def calculate_dissolved_volatiles(self,pressure,temperature=1200,sample=None,X_fluid=1.0,**kwargs):
		"""
		Calclates the dissolved CO2 concentration using (Eqns) 2-7 or 10-11 from Allison et al. (2019).

		Parameters
		----------
		pressure     float
			Pressure in bars.
		temperature     float
			Temperature in C.
		sample         pandas Series, dict or None
			Major element oxides in wt%. Required if using the thermodynamic fits, need not be
			provided if using the power law fits. Default is None.
		X_fluid     float
			The mole fraction of CO2 in the fluid. Default is 1.0.

		Returns
		-------
		float
			Dissolved CO2 concentration in wt%.
		"""
		# temperature = 1200 #temp in degrees C
		temperature = temperature + 273.15 #translate T from C to K

		if pressure < 0.0:
			raise InputError("Pressure must be positive.")
		if X_fluid < 0 or X_fluid > 1:
			raise InputError("X_fluid must have a value between 0 and 1.")

		if self.model_fit not in ['power','thermodynamic']:
			raise InputError("model_fit must be one of 'power', or 'thermodynamic'.")
		if self.model_loc not in ['sunset','sfvf','erebus','vesuvius','etna','stromboli']:
			raise InputError("model_loc must be one of 'sunset', 'sfvf', 'erebus', 'vesuvius', 'etna', or 'stromboli'.")


		if pressure == 0:
			return 0

		if self.model_fit == 'thermodynamic':
			if type(sample) != dict and type(sample) != pd.core.series.Series:
				raise InputError("Thermodynamic fit requires sample to be a dict or a pandas Series.")
			P0 = 1000 # bar
			params = dict({'sunset':[16.4,-14.67],
							'sfvf':[15.02,-14.87],
							'erebus':[15.83,-14.65],
							'vesuvius':[24.42,-14.04],
							'etna':[21.59,-14.28],
							'stromboli':[14.93,-14.68]})
			DV = params[self.model_loc][0]
			lnK0 = params[self.model_loc][1]

			lnK = lnK0 - (pressure-P0)*DV/(10*8.3141*temperature)
			fCO2 = self.fugacity_model.fugacity(pressure=pressure,temperature=temperature-273.15,X_fluid=X_fluid,**kwargs)
			Kf = np.exp(lnK)*fCO2
			XCO3 = Kf/(1-Kf)
			# FWone = wtpercentOxides_to_formulaWeight(sample)#,exclude_volatiles=True)
			FWone = 36.594
			wtCO2 = (44.01*XCO3)/((44.01*XCO3)+(1-XCO3)*FWone)*100

			return wtCO2
		if self.model_fit == 'power':
			params = dict({'stromboli':[1.05,0.883],
							'etna':[2.831,0.797],
							'vesuvius':[4.796,0.754],
							'sfvf':[3.273,0.74],
							'sunset':[4.32,0.728],
							'erebus':[5.145,0.713]})

			fCO2 = self.fugacity_model.fugacity(pressure=pressure,temperature=temperature-273.15,X_fluid=X_fluid,**kwargs)

			return params[self.model_loc][0]*fCO2**params[self.model_loc][1]/1e4



	def calculate_equilibrium_fluid_comp(self,pressure,sample,temperature=1200,**kwargs):
		""" Returns 1.0 if a pure CO2 fluid is saturated.
		Returns 0.0 if a pure CO2 fluid is undersaturated.

		Parameters
		----------
		pressure     float
			The total pressure of the system in bars.
		temperature     float
			The temperature of the system in C.
		sample         pandas Series or dict
			Major element oxides in wt% (including H2O).

		Returns
		-------
		float
			1.0 if CO2-fluid saturated, 0.0 otherwise.
		"""

		satP = self.calculate_saturation_pressure(temperature=temperature,sample=sample,X_fluid=1.0,**kwargs)
		if pressure < satP:
			return 1.0
		else:
			return 0.0

	def calculate_saturation_pressure(self,sample,temperature=1200,X_fluid=1.0,**kwargs):
		"""
		Calculates the pressure at which a pure CO2 fluid is saturated, for the given sample
		composition and CO2 concentration. Calls the scipy.root_scalar routine, which makes
		repeated called to the calculate_dissolved_volatiles method.

		Parameters
		----------
		temperature     float
			The temperature of the system in C.
		sample         pandas Series
			Major element oxides in wt% (including CO2).
		X_fluid     float
			The mole fraction of H2O in the fluid. Default is 1.0.

		Returns
		-------
		float
			Calculated saturation pressure in bars.
		"""

		# temperature = 1200 #temp in degrees C

		if X_fluid < 0 or X_fluid > 1:
			raise InputError("X_fluid must have a value between 0 and 1.")
		if type(sample) != dict and type(sample) != pd.core.series.Series:
			raise InputError("sample must be a dict or a pandas Series.")
		if 'CO2' not in sample:
			raise InputError("sample must contain CO2.")
		if sample['CO2'] < 0.0:
			raise InputError("Dissolved CO2 concentration must be greater than 0 wt%.")

		try:
			satP = root_scalar(self.root_saturation_pressure,args=(temperature,sample,X_fluid,kwargs),
								x0=1000.0,x1=2000.0).root
		except:
			w.warn("Saturation pressure not found.",RuntimeWarning,stacklevel=2)
			satP = np.nan
		return satP

	def root_saturation_pressure(self,pressure,temperature,sample,X_fluid,kwargs):
		""" Function called by scipy.root_scalar when finding the saturation pressure using
		calculate_saturation_pressure.

		Parameters
		----------
		pressure     float
			Pressure guess in bars
		temperature     float
			The temperature of the system in C.
		sample         pandas Series or dict
			Major elements in wt% (normalized to 100%), including CO2.
		kwargs         dictionary
			Additional keyword arguments supplied to calculate_saturation_pressure. Might be required for
			the fugacity or activity models.

		Returns
		-------
		float
			The differece between the dissolved CO2 at the pressure guessed, and the CO2 concentration
			passed in the sample variable.
		"""
		return sample['CO2'] - self.calculate_dissolved_volatiles(pressure=pressure,temperature=temperature,sample=sample,X_fluid=X_fluid,**kwargs)



#------------MIXED FLUID MODELS-------------------------------#
class MixedFluid(Model):
	"""
	Implements the generic framework for mixed fluid solubility. Any set of pure fluid solubility
	models may be specified.
	"""
	def __init__(self,models):
		"""
		Initializes the mixed fluid model.

		Parameters
		----------
		models     dictionary
			Dictionary with names of volatile species as keys, and the model objects as values.
		"""
		self.models = tuple(model for model in models.values())
		self.set_volatile_species(list(models.keys()))

	def preprocess_sample(self,sample):
		""" Returns sample, unmodified.

		Parameters
		----------
		sample         pandas Series or dict
			Major element oxides in wt%.

		Returns
		-------
		pandas Series or dict
			Major element oxides in wt%.
		"""
		if type(sample) != dict and type(sample) != pd.core.series.Series:
			raise InputError("sample must be a dict or a pandas Series.")
		_sample = sample.copy()
		_sample = self.models[0].preprocess_sample(_sample)
		return _sample

	def calculate_dissolved_volatiles(self,pressure,X_fluid,returndict=False,**kwargs):
		"""
		Calculates the dissolved volatile concentrations in wt%, using each model's
		calculate_dissolved_volatiles method. At present the volatile concentrations are
		not propagated through.

		Parameters
		----------
		pressure     float
			The total pressure in bars.
		X_fluid     float, numpy.ndarry, dict, pandas Series
			The mole fraction of each species in the fluid. If the mixed fluid model
			contains only two species (e.g. CO2 and H2O), the value of the first species in
			self.volatile_species may be passed on its own as a float.
		returndict 		bool
			If True, the results will be returned in a dict, otherwise they will be returned
			as a tuple.

		Returns
		-------
		tuple
			Dissolved volatile concentrations of each species in the model, in the order set
			by self.volatile_species.
		"""
		if (type(X_fluid) == float or type(X_fluid) == int) and len(self.volatile_species) == 2:
			X_fluid = (X_fluid,1-X_fluid)
		elif len(X_fluid) != len(self.volatile_species):
			raise InputError("X_fluid must have the same length as the number of volatile species\
			 in the MixedFluids Model class, or it may have length 1 if two species are present\
			 in the MixedFluids Model class.")

		if np.sum(X_fluid) != 1.0:
			raise InputError("X_fluid must sum to 1.0")
		if any(val<0 for val in X_fluid) or any(val>1 for val in X_fluid):
			raise InputError("Each mole fraction in X_fluid must have a value between 0 and 1.")

		if type(X_fluid) == dict or type(X_fluid) == pd.core.series.Series:
			X_fluid = tuple(X_fluid[species] for species in self.volatile_species)

		# If the models don't depend on the concentration of volatiles, themselves.
		if all(model.solubility_dependence == False for model in self.models):
			result = tuple(model.calculate_dissolved_volatiles(pressure=pressure,X_fluid=Xi,**kwargs) for model, Xi in zip(self.models,X_fluid))
		# If one of the models depends on the other volatile concentration
		elif len(self.models) == 2 and self.models[0].solubility_dependence == False and 'sample' in kwargs:
			result0 = self.models[0].calculate_dissolved_volatiles(pressure=pressure,X_fluid=X_fluid[0],**kwargs)
			samplecopy = kwargs['sample'].copy()
			samplecopy[self.volatile_species[0]] = result0
			kwargs['sample'] = samplecopy
			result1 = self.models[1].calculate_dissolved_volatiles(pressure=pressure,X_fluid=X_fluid[1],**kwargs)
			result = (result0,result1)
		elif len(self.models) == 2 and self.models[1].solubility_dependence == False and 'sample' in kwargs:
			result1 = self.models[1].calculate_dissolved_volatiles(pressure=pressure,X_fluid=X_fluid[1],**kwargs)
			samplecopy = kwargs['sample'].copy()
			samplecopy[self.volatile_species[1]] = result1
			kwargs['sample'] = samplecopy
			result0 = self.models[0].calculate_dissolved_volatiles(pressure=pressure,X_fluid=X_fluid[0],**kwargs)
			result = (result0,result1)
		else:
			raise InputError("The solubility dependence of the models is not currently supported by the MixedFluid model.")

		if returndict == True:
			resultsdict = {}
			for i,v in zip(range(len(self.volatile_species)),self.volatile_species):
				resultsdict.update({v+'_liq':result[i]})
			return resultsdict
		else:
			return result


	def calculate_equilibrium_fluid_comp(self,pressure,sample,return_dict=True,**kwargs):
		""" Calculates the composition of the fluid in equilibrium with the dissolved volatile
		concentrations passed. If a fluid phase is undersaturated at the chosen pressure (0,0) will
		be returned. Note, this currently assumes the given H2O and CO2 concentrations are
		the system total, not the total dissolved. If one of the volatile species has a zero or
		negative concentration, the pure fluid model for the other volatile species will be used.

		Parameters
		----------
		pressure     float
			The total pressure in bars.
		sample     pandas Series or dict
			Major element oxides in wt% (including volatiles).
		return_dict     bool
			Set the return type, if true a dict will be returned, if False two floats will be
			returned. Default is True.

		Returns
		-------
		dict or floats
			Mole fractions of the volatile species in the fluid, in the order given by
			self.volatile_species if floats.
		"""
		if len(self.volatile_species) != 2:
			raise InputError("Currently equilibrium fluid compositions can only be calculated when\
			two volatile species are present.")

		dissolved_at_0bar = [self.models[0].calculate_dissolved_volatiles(sample=sample,pressure=0.0,**kwargs),
							 self.models[1].calculate_dissolved_volatiles(sample=sample,pressure=0.0,**kwargs)]

		if sample[self.volatile_species[0]] <= 0.0 or sample[self.volatile_species[0]] <= dissolved_at_0bar[0]:
			Xv0 = 0.0
			Xv1 = self.models[1].calculate_equilibrium_fluid_comp(pressure=pressure,sample=sample,**kwargs)
		elif sample[self.volatile_species[1]] <= 0.0 or sample[self.volatile_species[1]] <= dissolved_at_0bar[1]:
			Xv1 = 0.0
			Xv0 = self.models[0].calculate_equilibrium_fluid_comp(pressure=pressure,sample=sample,**kwargs)
		else:
			satP = self.calculate_saturation_pressure(sample,**kwargs)

			if satP < pressure:
				if return_dict == True:
					return {self.volatile_species[0]:0,self.volatile_species[1]:0}
				else:
					return (0,0)


			molfracs = wtpercentOxides_to_molOxides(sample)
			(Xt0, Xt1) = (molfracs[self.volatile_species[0]],molfracs[self.volatile_species[1]])

			try:
				Xv0 = root_scalar(self.root_for_fluid_comp,bracket=[1e-15,1-1e-15],args=(pressure,Xt0,Xt1,sample,kwargs)).root
				Xv1 = 1 - Xv0
			except:
				try:
					Xv0 = root_scalar(self.root_for_fluid_comp,x0=0.5,x1=0.1,args=(pressure,Xt0,Xt1,sample,kwargs)).root
					Xv1 = 1 - Xv0
				except:
					raise SaturationError("Equilibrium fluid not found. Likely an issue with the numerical solver.")

		if return_dict == True:
			return {self.volatile_species[0]:Xv0,self.volatile_species[1]:Xv1}
		else:
			return Xv0, Xv1


	def calculate_saturation_pressure(self,sample,**kwargs):
		"""
		Calculates the pressure at which a fluid will be saturated, given the dissolved volatile
		concentrations. If one of the volatile species has a zero or negative concentration the
		pure fluid model for the other species will be used. If one of the volatile species has a
		concentration lower than the concentration dissolved at 0 bar, the pure fluid model for the
		other species will be used.

		Parameters
		----------
		sample     pandas Series or dict
			Major element oxides in wt% (including volatiles).

		Returns
		-------
		float
			The saturation pressure in bars.
		"""
		dissolved_at_0bar = [self.models[0].calculate_dissolved_volatiles(sample=sample,pressure=0.0,**kwargs),
							 self.models[1].calculate_dissolved_volatiles(sample=sample,pressure=0.0,**kwargs)]

		if sample[self.volatile_species[0]] <= 0.0 or sample[self.volatile_species[0]] <= dissolved_at_0bar[0]:
			satP = self.models[1].calculate_saturation_pressure(sample=sample,**kwargs)
		elif sample[self.volatile_species[1]] <= 0.0 or sample[self.volatile_species[1]] <= dissolved_at_0bar[1]:
			satP = self.models[0].calculate_saturation_pressure(sample=sample,**kwargs)
		else:
			volatile_concs = np.array(tuple(sample[species] for species in self.volatile_species))

			x0 = 0
			for model in self.models:
				xx0 = model.calculate_saturation_pressure(sample=sample,**kwargs)
				if np.isnan(xx0) == False:
					x0 += xx0

			try:
				satP = root(self.root_saturation_pressure,x0=[x0,0.5],args=(volatile_concs,sample,kwargs)).x[0]
			except:
				w.warn("Saturation pressure not found.",RuntimeWarning,stacklevel=2)
				satP = np.nan

		return satP

	def calculate_isobars_and_isopleths(self,pressure_list,isopleth_list=[0,1],points=51,
										return_dfs=True,extend_to_zero=True,**kwargs):
		"""
		Calculates isobars and isopleths. Isobars can be calculated for any number of pressures. Variables
		required by each of the pure fluid models must be passed, e.g. sample, temperature, etc.

		Parameters
		----------
		pressure_list     list or float
			List of all pressure values at which to calculate isobars, in bars.
		isopleth_list     list
			Default value is None, in which case only isobars will be calculated. List of all
			fluid compositions in mole fraction (of the first species in self.volatile_species) at which
			to calcualte isopleths. Values can range from 0 to 1.
		points     int
			The number of points in each isobar and isopleth. Default value is 101.
		return_dfs     bool
			If True, the results will be returned as two pandas DataFrames, as produced by the MagmaSat
			method. If False the results will be returned as lists of numpy arrays.

		Returns
		-------
		pandas DataFrame object(s) or list(s)
			If isopleth_list is not None, two objects will be returned, one with the isobars and the second with
			the isopleths. If return_dfs is True, two pandas DataFrames will be returned with column names
			'Pressure' or 'XH2O_fl', 'H2O_liq', and 'CO2_liq'. If return_dfs is False, two lists of numpy arrays
			will be returned. Each array is an individual isobar or isopleth, in the order passed via pressure_list
			or isopleth_list. The arrays are the concentrations of H2O and CO2 in the liquid, in the order of the
			species in self.volatile_species.

		"""
		if len(self.volatile_species) != 2 or 'H2O' not in self.volatile_species or 'CO2' not in self.volatile_species:
			raise InputError("calculate_isobars_and_isopleths may only be used with a H2O-CO2 fluid.")

		H2O_id = self.volatile_species.index('H2O')
		CO2_id = self.volatile_species.index('CO2')

		if isinstance(pressure_list, list):
			pass
		elif isinstance(pressure_list, int) or isinstance(pressure_list, float):
			pressure_list = [pressure_list]
		else:
			raise InputError("pressure_list must be a single float (1000.0), int (1000), or list of those [1000, 2000.0, 3000].")

		has_isopleths = True
		if isopleth_list is None:
			has_isopleths = False


		isobars_df = pd.DataFrame(columns=['Pressure','H2O_liq','CO2_liq'])
		isobars = []
		for pressure in pressure_list:
			dissolved = np.zeros([2,points])
			Xv0 = np.linspace(0.0,1.0,points)
			for i in range(points):
				dissolved[:,i] = self.calculate_dissolved_volatiles(pressure=pressure,X_fluid=(Xv0[i],1-Xv0[i]),**kwargs)
				isobars_df = isobars_df.append({'Pressure':pressure,'H2O_liq':dissolved[H2O_id,i],'CO2_liq':dissolved[CO2_id,i]},ignore_index=True)
			isobars.append(dissolved)


		if has_isopleths == True:
			isopleths_df = pd.DataFrame(columns=['XH2O_fl','H2O_liq','CO2_liq'])
			isopleths = []
			for isopleth in isopleth_list:
				dissolved = np.zeros([2,points])
				pmin = np.nanmin(pressure_list)
				pmax = np.nanmax(pressure_list)
				if pmin == pmax:
					pmin = 0.0
				pressure = np.linspace(pmin,pmax,points)
				for i in range(points):
					dissolved[:,i] = self.calculate_dissolved_volatiles(pressure=pressure[i],X_fluid=(isopleth,1-isopleth),**kwargs)
					isopleths_df = isopleths_df.append({'XH2O_fl':[isopleth,1-isopleth][H2O_id],'H2O_liq':dissolved[H2O_id,i],'CO2_liq':dissolved[CO2_id,i]},ignore_index=True)
				isopleths.append(dissolved)

		if return_dfs == True:
			if has_isopleths == True:
				return (isobars_df, isopleths_df)
			else:
				return isobars_df
		else:
			if has_isopleths == True:
				return (isobars, isopleths)
			else:
				return isobars


	def calculate_degassing_path(self,sample,pressure='saturation',fractionate_vapor=0.0,final_pressure=100.0,
									steps=101,return_dfs=True,round_to_zero=True,**kwargs):
		"""
		Calculates the dissolved volatiles in a progressively degassing sample.

		Parameters
		----------
		sample     pandas Series or dict
			Major element oxides in wt% (including volatiles).
		pressure     string, float, int, list, or numpy array
			Defaults to 'saturation', the calculation will begin at the saturation pressure. If a number is passed
			as either a float or int, this will be the starting pressure. If a list of numpy array is passed, the
			pressure values in the list or array will define the degassing path, i.e. final_pressure and steps
			variables will be ignored. Units are bars.
		fractionate_vapor     float
			What proportion of vapor should be removed at each step. If 0.0 (default), the degassing path will
			correspond to closed-system degassing. If 1.0, the degassing path will correspond to open-system
			degassing.
		final_pressure         float
			The final pressure on the degassing path, in bars. Ignored if a list or numpy array is passed as the
			pressure variable. Default is 1 bar.
		steps     int
			The number of steps in the degassing path. Ignored if a list or numpy array are passed as the pressure
			variable.
		return_dfs     bool
			If True, the results will be returned in a pandas DataFrame, if False, two numpy arrays will be returned.
		round_to_zero   bool
			If True, the first entry of FluidProportion_wt will be rounded to zero, rather than being a value
			within numerical error of zero. Default is True.

		Returns
		-------
		pandas DataFrame or numpy arrays
			If return_dfs is True (default), a DataFrame with columns 'Pressure_bars', 'H2O_liq', 'CO2_liq',
			'H2O_fl', 'CO2_fl', and 'FluidProportion_wt', is returned. Dissolved volatiles are in wt%,
			the proportions of volatiles in the fluid are in mole fraction. Otherwise a numpy array containing
			the dissolved volatile concentrations, and a numpy array containing the mole fractions of
			volatiles in the fluid is returned. The columns are in the order of the volatiles in
			self.volatile_species.

		"""

		# if 'model' in kwargs and model=='Liu':
		# 	final_pressure = 1.0

		wtptoxides = sample.copy()
		wtptoxides = normalize_FixedVolatiles(wtptoxides)
		wtm0s, wtm1s = (wtptoxides[self.volatile_species[0]],wtptoxides[self.volatile_species[1]])

		if pressure == 'saturation':
			p0 = self.calculate_saturation_pressure(wtptoxides,**kwargs)
			pressures = np.linspace(p0,final_pressure,steps)
		elif type(pressure) == float or type(pressure) == int:
			pressures = np.linspace(pressure,final_pressure,steps)
		elif type(pressure) == list or type(pressure) == np.ndarray:
			pressures = pressure

		Xv = np.zeros([2,len(pressures)])
		wtm = np.zeros([2,len(pressures)])

		for i in range(len(pressures)):
			try:
				X_fluid = self.calculate_equilibrium_fluid_comp(pressure=pressures[i],sample=wtptoxides,return_dict=False,**kwargs)
				Xv[:,i] = X_fluid
				if X_fluid == (0,0):
					wtm[:,i] = (wtptoxides[self.volatile_species[0]],wtptoxides[self.volatile_species[1]])
				else:
					if X_fluid[0] == 0:
						wtm[0,i] = wtptoxides[self.volatile_species[0]]
						wtm[1,i] = self.calculate_dissolved_volatiles(pressure=pressures[i],sample=wtptoxides,X_fluid=X_fluid,**kwargs)[1]
					elif X_fluid[1] == 0:
						wtm[1,i] = wtptoxides[self.volatile_species[1]]
						wtm[0,i] = self.calculate_dissolved_volatiles(pressure=pressures[i],sample=wtptoxides,X_fluid=X_fluid,**kwargs)[0]
					else:
						wtm[:,i] = self.calculate_dissolved_volatiles(pressure=pressures[i],sample=wtptoxides,X_fluid=X_fluid,**kwargs)

					wtptoxides[self.volatile_species[0]] = wtm[0,i] + (1-fractionate_vapor)*(wtm0s-wtm[0,i])
					wtptoxides[self.volatile_species[1]] = wtm[1,i] + (1-fractionate_vapor)*(wtm1s-wtm[1,i])
				# wtptoxides = normalize_FixedVolatiles(wtptoxides)
			except:
				Xv[:,i] = [np.nan]*np.shape(Xv)[0]
				wtm[:,i] = wtm[:,i-1]

		if return_dfs == True:
			exsolved_degassing_df = pd.DataFrame()
			exsolved_degassing_df['Pressure_bars'] = pressures
			exsolved_degassing_df['H2O_liq'] = wtm[self.volatile_species.index('H2O'),:]
			exsolved_degassing_df['CO2_liq'] = wtm[self.volatile_species.index('CO2'),:]
			exsolved_degassing_df['H2O_fl'] = Xv[self.volatile_species.index('H2O'),:]
			exsolved_degassing_df['CO2_fl'] = Xv[self.volatile_species.index('CO2'),:]
			exsolved_degassing_df['FluidProportion_wt'] = (wtm0s+wtm1s)-exsolved_degassing_df['H2O_liq']-exsolved_degassing_df['CO2_liq']

			if round_to_zero == True and np.round(exsolved_degassing_df.loc[0,'FluidProportion_wt'],2)==0:
				exsolved_degassing_df.loc[0,'FluidProportion_wt'] = 0.0

			return exsolved_degassing_df

		else:
			return (wtm, Xv)


	def root_saturation_pressure(self,x,volatile_concs,sample,kwargs):
		""" Function called by scipy.root when finding the saturation pressure using
		calculate_saturation_pressure.

		Parameters
		----------
		x     numpy array
			The guessed value for the root. x[0] is the pressure (in bars) and x[1] is the
			mole fraction of the first volatile in self.volatile_species.
		volatile_concs     numpy array
			The dissolved volatile concentrations, in the same order as self.volatile_species.
		sample     pandas Series or dict
			Major element oxides in wt% (including volatiles).
		kwargs     dictionary
			Dictionary of keyword arguments, which may be required by the pure-fluid models.

		Returns
		-------
		numpy array
			The difference in the dissolved volatile concentrations, and those predicted with the
			pressure and fluid composition specified by x.
		"""
		if x[1] < 0:
			x[1] = 0
		elif x[1] > 1:
			x[1] = 1
		if x[0] <= 0:
			x[0] = 1e-15
		misfit = np.array(self.calculate_dissolved_volatiles(pressure=x[0],X_fluid=(x[1],1-x[1]),sample=sample,**kwargs)) - volatile_concs
		return misfit


	def root_for_fluid_comp(self,Xv0,pressure,Xt0,Xt1,sample,kwargs):
		""" Function called by scipy.root_scalar when calculating the composition of equilibrium fluid
		in the calculate_equilibrium_fluid_comp method.

		Parameters
		----------
		Xv0     float
			The guessed mole fraction of the first volatile species in self.volatile_species.
		pressure     float
			The total pressure in bars.
		Xt0     float
			The total mole fraction of the first volatile species in self.volatile_species.
		Xt1        float
			The total mole fraction of the second volatile species in self.volatile_species.
		sample     pandas Series
			Major element oxides in wt%
		kwargs     dictionary
			A dictionary of keyword arguments that may be required by the pure fluid models.

		Returns
		-------
		float
			The differene in the LHS and RHS of the mass balance equation. Eq X in manuscript.
		"""

		wtt0 = sample[self.volatile_species[0]]
		wtt1 = sample[self.volatile_species[1]]

		wtm0, wtm1 = self.calculate_dissolved_volatiles(pressure=pressure,X_fluid=(Xv0,1-Xv0),sample=sample,**kwargs)

		Xm0 = Xt0/wtt0*wtm0
		Xm1 = Xt1/wtt1*wtm1

		if self.volatile_species[0] == 'CO2' and Xv0 != Xm0:
			f = (Xt0-Xm0)/(Xv0-Xm0)
			return (1-f)*Xm1 + f*(1-Xv0) - Xt1
		else:
			f = (Xt1-Xm1)/((1-Xv0)-Xm1)
			return (1-f)*Xm0 + f*Xv0 - Xt0

	def check_calibration_range(self,parameters,report_nonexistance=True):
		""" Checks whether the given parameters are within the ranges defined by the
		CalibrationRange objects for each model and its fugacity and activity models. An empty
		string will be returned if all parameters are within the calibration range. If a
		parameter is not within the calibration range, a description of the problem will be
		returned in the string.

		Parameters
		----------
		parameters 	dict
			Dictionary keys are the names of the parameters to be checked, e.g., pressure
			temperature, SiO2, etc. Values are the values of each parameter. A complete set
			need not be given.

		Returns
		-------
		str
			String description of any parameters falling outside of the calibration range.
		"""
		s = ''
		for model in self.models:
			for cr in model.calibration_ranges:
				if cr.check(parameters) == False:
					s += cr.string(parameters,report_nonexistance)
			for cr in model.fugacity_model.calibration_ranges:
				if cr.check(parameters) == False:
					s += cr.string(parameters,report_nonexistance)
			for cr in model.activity_model.calibration_ranges:
				if cr.check(parameters) == False:
					s += cr.string(parameters,report_nonexistance)
		return s

	def get_calibration_range(self):
		""" Returns a string describing the calibration ranges defined by the CalibrationRange
		objects for each model, and its associated fugacity and activity models.

		Returns
		-------
		str
			String description of the calibration range objects."""
		s = ''
		for model in self.models:
			for cr in model.calibration_ranges:
				s += cr.string(None)
			for cr in model.fugacity_model.calibration_ranges:
				s += cr.string(None)
			for cr in model.activity_model.calibration_ranges:
				s += cr.string(None)
		return s


class MagmaSat(Model):
	"""
	An object to instantiate a thermoengine equilibrate class
	"""

	def __init__(self):
		self.melts_version = '1.2.0' #just here so users can see which version is being used

		self.set_volatile_species(['H2O', 'CO2'])
		self.set_calibration_ranges([CalibrationRange('pressure',[0.0,20000.0],crf_Between,'bar','MagmaSat',
													  fail_msg=crmsg_Between_fail, pass_msg=crmsg_Between_pass, description_msg=crmsg_Between_description),
									 CalibrationRange('temperature',[800,1400],crf_Between,'oC','MagmaSat',
									 				  fail_msg=crmsg_Between_fail, pass_msg=crmsg_Between_pass, description_msg=crmsg_Between_description)])

	def preprocess_sample(self,sample):
		"""
		Returns sample with 0.0 values for any oxides not passed.

		Parameters
		----------
		sample: dictionary
			Sample composition in wt% oxides

		Returns
		-------
		tuple of dictionaries
			_sample: any sample information passed
			bulk_comp: only sample data referring to VESIcal oxides in wt%

		"""
		_sample = sample.copy()
		for oxide in oxides:
			if oxide in _sample.keys():
				pass
			else:
				_sample[oxide] = 0.0

		bulk_comp = {oxide:  _sample[oxide] for oxide in oxides}

		self.bulk_comp_orig = sample

		return _sample, bulk_comp

	def check_calibration_range(self,parameters,**kwargs):
		""" Checks whether supplied parameters and calculated results are within the calibration range
		of the model, defined by the CalibrationRange objects. An empty string will be returned if all
		parameters are within the calibration range. If a parameter is not within the calibration range,
		a description of the problem will be returned in the string.

		Parameters
		----------
		parameters 	dict
			Dictionary keys are the names of the parameters to be checked, e.g., pressure
			temperature, SiO2, etc. Values are the values of each parameter. A complete set
			need not be given.

		Returns
		-------
		str
			String description of any parameters falling outside of the calibration range.
		"""
		s = ''
		for cr in self.calibration_ranges:
			if cr.check(parameters) == False:
				s += cr.string(parameters,report_nonexistance=False)
			if 'notsaturated' in kwargs:
				s += "Sample not saturated at these conditions."
		return s

	def get_calibration_range(self):
		""" Returns a string describing the calibration ranges defined by the CalibrationRange
		objects for the model.

		Returns
		-------
		str
			String description of the calibration range objects."""
		s = ''
		for cr in self.calibration_ranges:
			s += cr.string(None)
		return s

	def get_fluid_mass(self, sample, temperature, pressure, H2O, CO2):
		"""An internally used function to calculate fluid mass.

		Parameters
		----------
		sample: dictionary
			Sample composition in wt% oxides

		temperature: float
			Temperature in degrees C.

		pressure: float
			Pressure in bars

		H2O: float
			wt% H2O in the system

		CO2: float
			wt% CO2 in the system

		Returns
		-------
		float
			mass of the fluid in grams
		"""
		pressureMPa = pressure / 10.0

		bulk_comp = {oxide:  sample[oxide] for oxide in oxides}
		bulk_comp["H2O"] = H2O
		bulk_comp["CO2"] = CO2
		feasible = melts.set_bulk_composition(bulk_comp)

		output = melts.equilibrate_tp(temperature, pressureMPa, initialize=True)
		(status, temperature, pressureMPa, xmlout) = output[0]
		fluid_mass = melts.get_mass_of_phase(xmlout, phase_name='Fluid')

		return fluid_mass

	def get_XH2O_fluid(self, sample, temperature, pressure, H2O, CO2):
		"""An internally used function to calculate fluid composition.

		Parameters
		----------
		sample: dictionary
			Sample composition in wt% oxides

		temperature: float
			Temperature in degrees C.

		pressure: float
			Pressure in bars

		H2O: float
			wt% H2O in the system

		CO2: float
			wt% CO2 in the system

		Returns
		-------
		float
			Mole fraction of H2O in the H2O-CO2 fluid

		"""
		sample, bulk_comp = self.preprocess_sample(sample)

		pressureMPa = pressure / 10.0

		bulk_comp["H2O"] = H2O
		bulk_comp["CO2"] = CO2
		feasible = melts.set_bulk_composition(bulk_comp)

		output = melts.equilibrate_tp(temperature, pressureMPa, initialize=True)
		(status, temperature, pressureMPa, xmlout) = output[0]
		fluid_comp = melts.get_composition_of_phase(xmlout, phase_name='Fluid', mode='component')
			#NOTE mode='component' returns endmember component keys with values in mol fraction.

		if "Water" in fluid_comp:
			H2O_fl = fluid_comp["Water"]
		else:
			H2O_fl = 0.0

		# if H2O_fl == 0:
		#   raise SaturationError("Composition not fluid saturated.")

		return H2O_fl

	def calculate_dissolved_volatiles(self, sample, temperature, pressure, X_fluid=1, H2O_guess=0.0, verbose=False, **kwargs):
	#TODO make refinements faster
		"""
		Calculates the amount of H2O and CO2 dissolved in a magma at saturation at the given P/T conditions and fluid composition.
		Fluid composition will be matched to within 0.0001 mole fraction.

		Parameters
		----------
		sample: dict or pandas Series
			Compositional information on one sample in oxides.

		temperature: float or int
			Temperature, in degrees C.

		presure: float or int
			Pressure, in bars.

		X_fluid: float or int
			The default value is 1. The mole fraction of H2O in the H2O-CO2 fluid. X_fluid=1 is a pure H2O fluid. X_fluid=0 is a pure CO2 fluid.

		verbose: bool
			OPTIONAL: Default is False. If set to True, returns H2O and CO2 concentration in the melt, H2O and CO2 concentration in
			the fluid, mass of the fluid in grams, and proportion of fluid in the system in wt%.

		Returns
		-------
		dict
			A dictionary of dissolved volatile concentrations in wt% with keys H2O and CO2.
		"""
		sample, bulk_comp = self.preprocess_sample(sample)

		if isinstance(X_fluid, int) or isinstance(X_fluid, float):
			pass
		else:
			raise InputError("X_fluid must be type int or float")

		if isinstance(H2O_guess, int) or isinstance(H2O_guess, float):
			pass
		else:
			raise InputError("H2O_guess must be type int or float")

		pressureMPa = pressure / 10.0

		if X_fluid != 0 and X_fluid !=1:
			if X_fluid < 0.001 or X_fluid > 0.999:
				raise InputError("X_fluid is calculated to a precision of 0.0001 mole fraction. \
								 Value for X_fluid must be between 0.0001 and 0.9999.")

		H2O_val = H2O_guess
		CO2_val = 0.0
		fluid_mass = 0.0
		while fluid_mass <= 0:
			if X_fluid == 0:
				CO2_val += 0.1
			elif X_fluid >= 0.5:
				H2O_val += 0.2
				CO2_val = (H2O_val / X_fluid) - H2O_val #NOTE this is setting XH2Owt of the system (not of the fluid) to X_fluid
				#TODO this is what needs to be higher for higher XH2O. Slows down computation by a second or two
			else:
				H2O_val += 0.1
				CO2_val = (H2O_val / X_fluid) - H2O_val #NOTE this is setting XH2Owt of the system (not of the fluid) to X_fluid
				#TODO this is what needs to be higher for higher XH2O. Slows down computation by a second or two

			fluid_mass = self.get_fluid_mass(sample, temperature, pressure, H2O_val, CO2_val)

		bulk_comp["H2O"] = H2O_val
		bulk_comp["CO2"] = CO2_val
		feasible = melts.set_bulk_composition(bulk_comp)

		output = melts.equilibrate_tp(temperature, pressureMPa, initialize=True)
		(status, temperature, pressureMPa, xmlout) = output[0]
		liquid_comp = melts.get_composition_of_phase(xmlout, phase_name='Liquid', mode='oxide_wt')
		fluid_comp = melts.get_composition_of_phase(xmlout, phase_name='Fluid', mode='component')

		if "Water" in fluid_comp:
			H2O_fl = fluid_comp["Water"]
		else:
			H2O_fl = 0.0

		XH2O_fluid = H2O_fl

		#------Coarse Check------#
		while XH2O_fluid < X_fluid - 0.1: #too low coarse check
			H2O_val += 0.2
			XH2O_fluid = self.get_XH2O_fluid(sample, temperature, pressure, H2O_val, CO2_val)

		while XH2O_fluid > X_fluid + 0.1: #too high coarse check
			CO2_val += 0.1
			XH2O_fluid = self.get_XH2O_fluid(sample, temperature, pressure, H2O_val, CO2_val)

		#------Refinement 1------#
		while XH2O_fluid < X_fluid - 0.01: #too low refinement 1
			H2O_val += 0.05
			XH2O_fluid = self.get_XH2O_fluid(sample, temperature, pressure, H2O_val, CO2_val)

		while XH2O_fluid > X_fluid + 0.01: #too high refinement 1
			CO2_val += 0.01
			XH2O_fluid = self.get_XH2O_fluid(sample, temperature, pressure, H2O_val, CO2_val)

		#------Refinement 2------#
		while XH2O_fluid < X_fluid - 0.001: #too low refinement 2
			H2O_val += 0.005
			XH2O_fluid = self.get_XH2O_fluid(sample, temperature, pressure, H2O_val, CO2_val)

		while XH2O_fluid > X_fluid + 0.001: #too high refinement 2
			CO2_val += 0.001
			XH2O_fluid = self.get_XH2O_fluid(sample, temperature, pressure, H2O_val, CO2_val)

		#------Final refinement------#
		while XH2O_fluid < X_fluid - 0.0001: #too low final refinement
			H2O_val += 0.001
			XH2O_fluid = self.get_XH2O_fluid(sample, temperature, pressure, H2O_val, CO2_val)

		while XH2O_fluid > X_fluid + 0.0001: #too high final refinement
			CO2_val += 0.0001
			XH2O_fluid = self.get_XH2O_fluid(sample, temperature, pressure, H2O_val, CO2_val)

		#------Get calculated values------#
		bulk_comp["H2O"] = H2O_val
		bulk_comp["CO2"] = CO2_val
		feasible = melts.set_bulk_composition(bulk_comp)

		output = melts.equilibrate_tp(temperature, pressureMPa, initialize=True)
		(status, temperature, pressureMPa, xmlout) = output[0]
		fluid_mass = melts.get_mass_of_phase(xmlout, phase_name='Fluid')
		system_mass = melts.get_mass_of_phase(xmlout, phase_name='System')
		liquid_comp = melts.get_composition_of_phase(xmlout, phase_name='Liquid', mode='oxide_wt')
		fluid_comp = melts.get_composition_of_phase(xmlout, phase_name='Fluid', mode='component')

		if "H2O" in liquid_comp:
			H2O_liq = liquid_comp["H2O"]
		else:
			H2O_liq = 0

		if "CO2" in liquid_comp:
			CO2_liq = liquid_comp["CO2"]
		else:
			CO2_liq = 0

		if "Water" in fluid_comp:
			H2O_fl = fluid_comp["Water"]
		else:
			H2O_fl = 0.0
		if "Carbon Dioxide" in fluid_comp:
			CO2_fl = fluid_comp["Carbon Dioxide"]
		else:
			CO2_fl = 0.0

		XH2O_fluid = H2O_fl

		if verbose == True:
			return {"temperature": temperature, "pressure": pressure,
					"H2O_liq": H2O_liq, "CO2_liq": CO2_liq,
					"XH2O_fl": H2O_fl, "XCO2_fl": CO2_fl,
					"FluidProportion_wt": 100*fluid_mass/system_mass}

		if verbose == False:
			return {"CO2": CO2_liq, "H2O": H2O_liq}

	def calculate_equilibrium_fluid_comp(self, sample, temperature, pressure, verbose=False, **kwargs):
		"""
		Returns H2O and CO2 concentrations in wt% in a fluid in equilibrium with the given sample at the given P/T condition.

		Parameters
		----------
		sample: dict or pandas Series
			Compositional information on one sample in oxides.

		temperature: float or int
			Temperature, in degrees C.

		presure: float or int
			Pressure, in bars.

		verbose: bool
			OPTIONAL: Default is False. If set to True, returns H2O and CO2 concentration in the fluid, mass of the fluid in grams,
			and proportion of fluid in the system in wt%.

		Returns
		-------
		dict
			A dictionary of fluid composition in wt% with keys 'H2O' and 'CO2' is returned.
		"""
		sample, bulk_comp = self.preprocess_sample(sample)

		if isinstance(temperature, float) or isinstance(temperature, int):
			pass
		else:
			raise InputError("temp must be type float or int")

		if isinstance(pressure, float) or isinstance(pressure, int):
			pass
		else:
			raise InputError("presure must be type float or int")

		#Check if only single volatile species is passed. If so, can skip calculations.
		if sample["H2O"] == 0:
			if sample["CO2"] == 0:
				if verbose == False:
					return {'CO2': 0.0, 'H2O': 0.0}
				if verbose == True:
					return {'CO2': 0.0, 'H2O': 0.0, 'FluidMass_grams': 0.0, 'FluidProportion_wt': 0.0}
			else:
				if verbose == False:
					return {'CO2': 1.0, 'H2O': 0.0}
		else:
			if sample["CO2"] == 0:
				if verbose == False:
					return {'CO2': 0.0, 'H2O': 1.0}

		pressureMPa = pressure / 10.0

		feasible = melts.set_bulk_composition(bulk_comp)

		output = melts.equilibrate_tp(temperature, pressureMPa, initialize=True)
		(status, temperature, pressureMPa, xmlout) = output[0]
		fluid_mass = melts.get_mass_of_phase(xmlout, phase_name='Fluid')
		flsystem_wtper = 100 * fluid_mass / (fluid_mass + melts.get_mass_of_phase(xmlout, phase_name='Liquid'))

		if fluid_mass > 0.0:
			fluid_comp = melts.get_composition_of_phase(xmlout, phase_name='Fluid', mode='component')
			fluid_comp_H2O = fluid_comp['Water']
			fluid_comp_CO2 = fluid_comp['Carbon Dioxide']
		else:
			fluid_comp_H2O = 0
			fluid_comp_CO2 = 0

		feasible = melts.set_bulk_composition(bulk_comp) #reset

		if verbose == False:
			return {'CO2': fluid_comp_CO2, 'H2O': fluid_comp_H2O}

		if verbose == True:
			return {'CO2': fluid_comp_CO2, 'H2O': fluid_comp_H2O, 'FluidMass_grams': fluid_mass, 'FluidProportion_wt': flsystem_wtper}

	def calculate_saturation_pressure(self, sample, temperature, verbose=False, **kwargs):
		"""
		Calculates the saturation pressure of a sample composition.

		Parameters
		----------
		sample: dict, pandas Series
			Compositional information on one sample. A single sample can be passed as a dict or pandas Series.

		temperature: flaot or int
			Temperature of the sample in degrees C.

		verbose: bool
			OPTIONAL: Default is False. If set to False, only the saturation pressure is returned. If set to True,
			the saturation pressure, mass of fluid in grams, proportion of fluid in wt%, and H2O and CO2 concentrations
			in the fluid in mole fraction are all returned in a dict.

		Returns
		-------
		float or dict
			If verbose is set to False: Saturation pressure in bars.
			If verbose is set to True: dict of all calculated values.
		"""
		sample, bulk_comp = self.preprocess_sample(sample)
		bulk_comp_orig = sample

		feasible = melts.set_bulk_composition(bulk_comp)
		#Coarse search
		fluid_mass = 0
		pressureMPa = 2000 #NOTE that pressure is in MPa for MagmaSat calculations but reported in bars.
		while fluid_mass <= 0:
			pressureMPa -= 100
			if pressureMPa <= 0:
				break

			output = melts.equilibrate_tp(temperature, pressureMPa, initialize=True)
			(status, temperature, pressureMPa, xmlout) = output[0]
			fluid_mass = melts.get_mass_of_phase(xmlout, phase_name='Fluid')

		pressureMPa+=100

		#Refined search 1
		feasible = melts.set_bulk_composition(bulk_comp)
		fluid_mass = 0
		while fluid_mass <= 0:
			pressureMPa -= 10
			if pressureMPa <= 0:
				break

			output = melts.equilibrate_tp(temperature, pressureMPa, initialize=True)
			(status, temperature, pressureMPa, xmlout) = output[0]
			fluid_mass = melts.get_mass_of_phase(xmlout, phase_name='Fluid')

		pressureMPa += 10

		#Refined search 2
		feasible = melts.set_bulk_composition(bulk_comp)
		fluid_mass = 0
		while fluid_mass <= 0:
			pressureMPa -= 1
			if pressureMPa <= 0:
				break

			output = melts.equilibrate_tp(temperature, pressureMPa, initialize=True)
			(status, temperature, pressureMPa, xmlout) = output[0]
			fluid_mass = melts.get_mass_of_phase(xmlout, phase_name='Fluid')

		if pressureMPa != np.nan:
			satP = pressureMPa*10 #convert pressure to bars
			flmass = fluid_mass
			flsystem_wtper = 100 * fluid_mass / (fluid_mass + melts.get_mass_of_phase(xmlout, phase_name='Liquid'))
			flcomp = melts.get_composition_of_phase(xmlout, phase_name='Fluid', mode='component')
			try:
				flH2O = flcomp['Water']
			except:
				flH2O = 0.0
			try:
				flCO2 = flcomp['Carbon Dioxide']
			except:
				flCO2 = 0.0
		else:
			flmass = np.nan
			flsystem_wtper = np.nan
			flH2O = np.nan
			flCO2 = np.nan
			warnmessage = 'Calculation failed.'

		feasible = melts.set_bulk_composition(bulk_comp_orig) #this needs to be reset always!

		if verbose == False:
			try:
				w.warn(warnmessage)
			except:
				pass
			return satP

		elif verbose == True:
			try:
				w.warn(warnmessage)
			except:
				pass
			return {"SaturationP_bars": satP, "FluidMass_grams": flmass, "FluidProportion_wt": flsystem_wtper,
	 				"XH2O_fl": flH2O, "XCO2_fl": flCO2}

	def calculate_isobars_and_isopleths(self, sample, temperature, pressure_list, isopleth_list=None,
										smooth_isobars=True, smooth_isopleths=True, print_status=True, **kwargs):
		"""
		Calculates isobars and isopleths at a constant temperature for a given sample. Isobars can be calculated
		for any number of pressures. Isobars are calculated using 5 XH2O values (0, 0.25, 0.5, 0.75, 1).

		Parameters
		----------
		sample: dict
			Dictionary with values for sample composition as oxides in wt%.

		temperature: float
			Temperature in degrees C.

		pressure_list: list or float
			List of all pressure values at which to calculate isobars, in bars. If only one value is passed
			it can be as float instead of list.

		isopleth_list: list or float
			OPTIONAL: Default value is None in which case only isobars will be calculated.
			List of all fluid compositions in mole fraction H2O (XH2Ofluid) at which to calcualte isopleths. Values can range from 0-1.
			If only one value is passed it can be as float instead of list.

		smooth_isobars: bool
			OPTIONAL. Default is True. If set to True, polynomials will be fit to the computed isobar points.

		smooth_isopleths: bool
			OPTIONAL. Default is True. If set to True, polynomials will be fit to the computed isopleth points.

		print_status: bool
			OPTIONAL: Default is True. If set to True, progress of the calculations will be printed to the terminal.

		Returns
		-------
		pandas DataFrame objects
			Two pandas DataFrames are returned; the first has isobar data, and the second has isopleth data. Columns in the
			isobar dataframe are 'Pressure', 'H2Omelt', and 'CO2melt', correpsonding to pressure in bars and dissolved H2O
			and CO2 in the liquid in wt%. Columns in the isopleth dataframe are 'Pressure', 'H2Ofl', and 'CO2fl',
			corresponding to pressure in bars and H2O and CO2 concentration in the H2O-CO2 fluid, in wt%.
		"""

		if isinstance(pressure_list, list):
			P_vals = pressure_list
		elif isinstance(pressure_list, int) or isinstance(pressure_list, float):
			P_vals = [pressure_list]
		else:
			raise InputError("pressure_list must be a single float (1000.0), int (1000), or list of those [1000, 2000.0, 3000].")

		if isopleth_list == None:
			has_isopleths = False
		elif isinstance(isopleth_list, list):
			iso_vals = isopleth_list
			has_isopleths = True
		else:
			iso_vals = [isopleth_list]
			has_isopleths = True

		sample, bulk_comp = self.preprocess_sample(sample)

		required_iso_vals = [0, 0.25, 0.5, 0.75, 1]
		all_iso_vals = iso_vals + required_iso_vals
		all_iso_vals = list(dict.fromkeys(all_iso_vals)) #remove duplicates
		all_iso_vals.sort() #sort from smallest to largest

		isobar_data = []
		isopleth_data = []
		for X in iso_vals:
			isopleth_data.append([X, 0.0, 0.0])
		H2O_val = 0.0
		CO2_val = 0.0
		fluid_mass = 0.0
		# Calculate equilibrium phase assemblage for all P/T conditions, check if saturated in fluid...
		for i in P_vals:
			guess = 0.0
			if print_status == True:
				print("Calculating isobar at " + str(i) + " bars")
			for X in all_iso_vals:
				if print_status == True:
					if isopleth_list != None and X in iso_vals:
						print("Calculating isopleth at XH2Ofluid = " + str(X))
					if X not in iso_vals:
						print("Calculating isobar control point at XH2Ofluid = " + str(X))
				saturated_vols = self.calculate_dissolved_volatiles(sample=bulk_comp, temperature=temperature, pressure=i, H2O_guess=guess, X_fluid=X)

				if X in required_iso_vals:
					isobar_data.append([i, saturated_vols['H2O'], saturated_vols['CO2']])
				if X in iso_vals:
					isopleth_data.append([X, saturated_vols['H2O'], saturated_vols['CO2']])

				guess = saturated_vols['H2O']

		if print_status == True:
			print("Done!")

		isobars_df = pd.DataFrame(isobar_data, columns=['Pressure', 'H2O_liq', 'CO2_liq'])
		isopleths_df = pd.DataFrame(isopleth_data, columns=['XH2O_fl', 'H2O_liq', 'CO2_liq'])

		feasible = melts.set_bulk_composition(self.bulk_comp_orig) #reset

		if smooth_isobars == True:
			isobars_smoothed = smooth_isobars_and_isopleths(isobars=isobars_df)
			res_isobars = isobars_smoothed.copy()
		else:
			res_isobars = isobars_df.copy()

		if smooth_isopleths == True:
			isopleths_smoothed = smooth_isobars_and_isopleths(isopleths=isopleths_df)
			res_isopleths = isopleths_smoothed.copy()
		else:
			res_isopleths = isopleths_df.copy()

		return res_isobars, res_isopleths

	def calculate_degassing_path(self, sample, temperature, pressure='saturation', fractionate_vapor=0.0, init_vapor=0.0, steps=50, **kwargs):
		"""
		Calculates degassing path for one sample

		Parameters
		----------
		sample: dict
			Dictionary with values for sample composition as oxides in wt%. If pulling from an uploaded file
			with data for many samples, first call get_sample_oxide_comp() to get the sample desired. Then pass
			the result into this function.

		temperature: float
			Temperature at which to calculate degassing paths, in degrees C.

		pressure: float
			OPTIONAL. The perssure at which to begin the degassing calculations. Default value is 'saturation', which runs the
			calculation with the initial pressure at the saturation pressure. If a pressure greater than the saturation pressure
			is input, the calculation will start at saturation, since this is the first pressure at which any degassing will
			occur.

		fractionate_vapor: float
			OPTIONAL. Proportion of vapor removed at each pressure step.
			Default value is 0.0 (completely closed-system degassing). Specifies the type of calculation performed, either
			closed system (0.0) or open system (1.0) degassing. If any value between <1.0 is chosen, user can also specify the
			'init_vapor' argument (see below). A value in between 0 and 1 will remove that proportion of vapor at each step.
			For example, for a value of 0.2, the calculation will remove 20% of the vapor and retain 80% of the vapor at each
			pressure step.

		init_vapor: float
			OPTIONAL. Default value is 0.0. Specifies the amount of vapor (in wt%) coexisting with the melt before
			degassing.

		steps: int
			OPTIONAL. Default value is 50. Specifies the number of steps in pressure space at which dissolved volatile
			concentrations are calculated.

		Returns
		-------
		pandas DataFrame object

		"""
		_sample, bulk_comp = self.preprocess_sample(sample)
		bulk_comp_orig = sample.copy()

		#MELTS needs to be reloaded here. If an unfeasible composition gets set inside of MELTS, which can
		#happen when running open-system degassing path calcs, the following calls to MELTS will fail. This
		#prevents that from happening.
		melts = equilibrate.MELTSmodel('1.2.0')

		# Suppress phases not required in the melts simulation
		phases = melts.get_phase_names()
		for phase in phases:
			melts.set_phase_inclusion_status({phase: False})
		melts.set_phase_inclusion_status({'Fluid': True, 'Liquid': True})

		bulk_comp = normalize(bulk_comp)
		feasible = melts.set_bulk_composition(bulk_comp)

		# Get saturation pressure
		data = self.calculate_saturation_pressure(sample=_sample, temperature=temperature, verbose=True)

		if pressure == 'saturation' or pressure >= data["SaturationP_bars"]:
			SatP_MPa = data["SaturationP_bars"] / 10.0
		else:
			SatP_MPa = pressure / 10.0

		#convert number of steps to step size
		MPa_step = SatP_MPa / steps
		if MPa_step < 1:
			MPa_step = 1

		P_array = np.arange(1.0, SatP_MPa, MPa_step)
		P_array = -np.sort(-P_array)
		fl_wtper = data["FluidProportion_wt"]

		while fl_wtper <= init_vapor:
			output = melts.equilibrate_tp(temperature, SatP_MPa, initialize=True)
			(status, temperature, p, xmlout) = output[0]
			fl_mass = melts.get_mass_of_phase(xmlout, phase_name='Fluid')
			liq_mass = melts.get_mass_of_phase(xmlout, phase_name='Liquid')
			fl_comp = melts.get_composition_of_phase(xmlout, phase_name='Fluid')
			fl_wtper = 100 * fl_mass / (fl_mass+liq_mass)
			try:
				bulk_comp["H2O"] += fl_comp["H2O"]*0.0005
			except:
				bulk_comp["H2O"] = bulk_comp["H2O"] * 1.1
			try:
				bulk_comp["CO2"] += fl_comp["CO2"]*0.0005
			except:
				bulk_comp["CO2"] = bulk_comp["CO2"] * 1.1
			bulk_comp = normalize(bulk_comp)
			feasible = melts.set_bulk_composition(bulk_comp)

		pressure = []
		H2Oliq = []
		CO2liq = []
		H2Ofl = []
		CO2fl = []
		fluid_wtper = []
		for i in P_array:
			fl_mass = 0.0
			feasible = melts.set_bulk_composition(bulk_comp)
			output = melts.equilibrate_tp(temperature, i, initialize=True)
			(status, temperature, p, xmlout) = output[0]
			liq_comp = melts.get_composition_of_phase(xmlout, phase_name='Liquid')
			fl_comp = melts.get_composition_of_phase(xmlout, phase_name='Fluid', mode='component')
			liq_mass = melts.get_mass_of_phase(xmlout, phase_name='Liquid')
			fl_mass = melts.get_mass_of_phase(xmlout, phase_name='Fluid')
			fl_wtper = 100 * fl_mass / (fl_mass+liq_mass)

			if fl_mass > 0:
				pressure.append(p * 10.0)
				try:
					H2Oliq.append(liq_comp["H2O"])
				except:
					H2Oliq.append(0)
				try:
					CO2liq.append(liq_comp["CO2"])
				except:
					CO2liq.append(0)
				try:
					H2Ofl.append(fl_comp["Water"])
				except:
					H2Ofl.append(0)
				try:
					CO2fl.append(fl_comp["Carbon Dioxide"])
				except:
					CO2fl.append(0)
				fluid_wtper.append(fl_wtper)

				try:
					bulk_comp["H2O"] = liq_comp["H2O"] + (bulk_comp["H2O"] - liq_comp["H2O"]) * (1.0-fractionate_vapor)
				except:
					bulk_comp["H2O"] = 0
				try:
					bulk_comp["CO2"] = liq_comp["CO2"] + (bulk_comp["CO2"] - liq_comp["CO2"]) * (1.0-fractionate_vapor)
				except:
					bulk_comp["CO2"] = 0
			bulk_comp = normalize(bulk_comp)

		feasible = melts.set_bulk_composition(bulk_comp_orig) #this needs to be reset always!
		open_degassing_df = pd.DataFrame(list(zip(pressure, H2Oliq, CO2liq, H2Ofl, CO2fl, fluid_wtper)),
									columns =['Pressure_bars', 'H2O_liq', 'CO2_liq', 'XH2O_fl', 'XCO2_fl', 'FluidProportion_wt'])

		open_degassing_df = open_degassing_df[open_degassing_df.CO2_liq > 0.000001]
		open_degassing_df = open_degassing_df[open_degassing_df.H2O_liq > 0.000001]

		return open_degassing_df

#-----------MAGMASAT PLOTTING FUNCTIONS-----------#
def smooth_isobars_and_isopleths(isobars=None, isopleths=None):
	"""
	Takes in a dataframe with calculated isobar and isopleth information (e.g., output from calculate_isobars_and_isopleths)
	and smooths the data for plotting.

	Parameters
	----------
	isobars: pandas DataFrame
		OPTIONAL. DataFrame object containing isobar information as calculated by calculate_isobars_and_isopleths.

	isopleths: pandas DataFrame
		OPTIONAL. DataFrame object containing isopleth information as calculated by calculate_isobars_and_isopleths.

	Returns
	-------
	pandas DataFrame
		DataFrame with x and y values for all isobars and all isopleths. Useful if a user wishes to do custom plotting
		with isobar and isopleth data rather than using the built-in `plot_isobars_and_isopleths()` function.
	"""
	np.seterr(divide='ignore', invalid='ignore') #turn off numpy warning
	w.filterwarnings("ignore", message="Polyfit may be poorly conditioned")

	if isobars is not None:
		P_vals = isobars.Pressure.unique()
		isobars_lists = isobars.values.tolist()
		# add zero values to volatiles list
		isobars_lists.append([0.0, 0.0, 0.0, 0.0])

		isobars_pressure = []
		isobars_H2O_liq = []
		isobars_CO2_liq = []
		# do some data smoothing
		for pressure in P_vals:
			Pxs = [item[1] for item in isobars_lists if item[0] == pressure]
			Pys = [item[2] for item in isobars_lists if item[0] == pressure]

			try:
				## calcualte polynomial
				Pz = np.polyfit(Pxs, Pys, 3)
				Pf = np.poly1d(Pz)

				## calculate new x's and y's
				Px_new = np.linspace(Pxs[0], Pxs[-1], 50)
				Py_new = Pf(Px_new)

				# Save x's and y's
				Px_new_list = list(Px_new)
				isobars_H2O_liq += Px_new_list

				Py_new_list = list(Py_new)
				isobars_CO2_liq += Py_new_list

				pressure_vals_for_list = [pressure]*len(Px_new)
				isobars_pressure += pressure_vals_for_list

			except:
				Px_list = list(Pxs)
				isobars_H2O_liq += Px_list

				Py_list = list(Pys)
				isobars_CO2_liq += Py_list

				pressure_vals_for_list = [pressure]*len(Pxs)
				isobars_pressure += pressure_vals_for_list

		isobar_df = pd.DataFrame({"Pressure": isobars_pressure,
							  "H2O_liq": isobars_H2O_liq,
							  "CO2_liq": isobars_CO2_liq})

	if isopleths is not None:
		XH2O_vals = isopleths.XH2O_fl.unique()
		isopleths_lists = isopleths.values.tolist()

		isopleths_XH2O_fl = []
		isopleths_H2O_liq = []
		isopleths_CO2_liq = []
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

				# Save x's and y's
				Xx_new_list = list(Xx_new)
				isopleths_H2O_liq += Xx_new_list

				Xy_new_list = list(Xy_new)
				isopleths_CO2_liq += Xy_new_list

				XH2Ofl_vals_for_list = [Xfl]*len(Xx_new)
				isopleths_XH2O_fl += XH2Ofl_vals_for_list

			except:
				Xx_list = list(Xxs)
				isopleths_H2O_liq += Xx_list

				Xy_list = list(Xys)
				isopleths_CO2_liq += Xy_list

				XH2Ofl_vals_for_list = [Xfl]*len(Xxs)
				isopleths_XH2O_fl += XH2Ofl_vals_for_list

		isopleth_df = pd.DataFrame({"XH2O_fl": isopleths_XH2O_fl,
							  "H2O_liq": isopleths_H2O_liq,
							  "CO2_liq": isopleths_CO2_liq})

	np.seterr(divide='warn', invalid='warn') #turn numpy warning back on
	w.filterwarnings("always", message="Polyfit may be poorly conditioned")

	if isobars is not None:
		if isopleths is not None:
			return isobar_df, isopleth_df
		else:
			return isobar_df
	else:
		if isopleths is not None:
			return isopleth_df


def plot(isobars=None, isopleths=None, degassing_paths=None, custom_H2O=None, custom_CO2=None,
		 isobar_labels=None, isopleth_labels=None, degassing_path_labels=None, custom_labels=None,
		 custom_colors="VESIcal", custom_symbols=None, markersize=10,
		 extend_isobars_to_zero=True, smooth_isobars=False, smooth_isopleths=False, **kwargs):
	"""
	Custom automatic plotting of model calculations in VESIcal.
	Isobars, isopleths, and degassing paths can be plotted. Labels can be specified for each.
	Any combination of isobars, isopleths, and degassing paths can be plotted.

	Parameters
	----------
	isobars: pandas DataFrame or list
		OPTIONAL. DataFrame object containing isobar information as calculated by calculate_isobars_and_isopleths. Or a list
		of DataFrame objects.

	isopleths: pandas DataFrame or list
		OPTIONAL. DataFrame object containing isopleth information as calculated by calculate_isobars_and_isopleths. Or a list
		of DataFrame objects.

	degassing_paths: list
		OPTIONAL. List of DataFrames with degassing information as generated by calculate_degassing_path().

	custom_H2O: list
		OPTIONAL. List of groups of H2O values to plot as points. For example myfile.data['H2O'] is one group of H2O values.
		Must be passed with custom_CO2 and must be same length as custom_CO2.

	custom_CO2: list
		OPTIONAL. List of groups of CO2 values to plot as points.For example myfile.data['CO2'] is one group of CO2 values.
		Must be passed with custom_H2O and must be same length as custom_H2O.

	isobar_labels: list
		OPTIONAL. Labels for the plot legend. Default is None, in which case each plotted line will be given the generic
		legend name of "Isobars n", with n referring to the nth isobars passed. Isobar pressure is given in parentheses.
		The user can pass their own labels as a list of strings. If more than one set of isobars is passed, the labels should
		refer to each set of isobars, not each pressure.

	isopleth_labels: list
		OPTIONAL. Labels for the plot legend. Default is None, in which case each plotted isopleth will be given the generic
		legend name of "Isopleth n", with n referring to the nth isopleths passed. Isopleth XH2O values are given in
		parentheses. The user can pass their own labels as a list of strings. If more than one set of isopleths is passed,
		the labels should refer to each set of isopleths, not each XH2O value.

	degassing_path_labels: list
		OPTIONAL. Labels for the plot legend. Default is None, in which case each plotted line will be given the generic
		legend name of "Pathn", with n referring to the nth degassing path passed. The user can pass their own labels
		as a list of strings.

	custom_labels: list
		OPTIONAL. Labels for the plot legend. Default is None, in which case each group of custom points will be given the
		generic legend name of "Customn", with n referring to the nth degassing path passed. The user can pass their own labels
		as a list of strings.

	custom_colors: list
		OPTIONAL. Default value is "VESIcal", which uses VESIcal's color ramp. A list of color values readable by matplotlib
		can be passed here if custom symbol colors are desired. The length of this list must match that of custom_H2O and
		custom_CO2.

	custom_symbols: list
		OPTIONAL. Default value is None, in which case data are plotted as filled circles.. A list of symbol tyles readable
		by matplotlib can be passed here if custom symbol types are desired. The length of this list must match that of
		custom_H2O and custom_CO2.

	markersize: int
		OPTIONAL. Default value is 10. Same as markersize kwarg in matplotlib. Any numeric value passed here will set the
		marker size for (custom_H2O, custom_CO2) points.

	extend_isobars_to_zero: bool
		OPTIONAL. If True (default), isobars will be extended to zero, even if there is a finite solubility at zero partial pressure.

	smooth_isobars: bool
		OPTIONAL. Default is False. If set to True, isobar data will be fit to a polynomial and plotted. If False, the raw input data
		will be plotted.

	smooth_isopleths: bool
		OPTIONAL. Default is False. If set to True, isopleth data will be fit to a polynomial and plotted. If False, the raw input data
		will be plotted.

	Returns
	-------
	matplotlib object
		Plot with x-axis as H2O wt% in the melt and y-axis as CO2 wt% in the melt. Isobars, or lines of
		constant pressure at which the sample magma composition is saturated, and isopleths, or lines of constant
		fluid composition at which the sample magma composition is saturated, are plotted if passed. Degassing
		paths, or the concentration of dissolved H2O and CO2 in a melt equilibrated along a path of decreasing
		pressure, is plotted if passed.
	"""
	np.seterr(divide='ignore', invalid='ignore') #turn off numpy warning
	w.filterwarnings("ignore", message="Polyfit may be poorly conditioned")

	if custom_H2O is not None:
		if custom_CO2 is None:
			raise InputError("If x data is passed, y data must also be passed.")
		else:
			if len(custom_H2O) == len(custom_CO2):
				pass
			else:
				raise InputError("x and y data must be same length")
	if custom_CO2 is not None:
		if custom_H2O is None:
			raise InputError("If y data is passed, x data must also be passed.")

	if custom_colors == "VESIcal":
		use_colors = color_list
	elif isinstance(custom_colors, list):
		use_colors = custom_colors
	else:
		raise InputError("Argument custom_colors must be type list. Just passing one item? Try putting square brackets, [], around it.")

	plt.figure(figsize=(12,8))
	if 'custom_x' in kwargs:
		plt.xlabel(kwargs['xlabel'])
		plt.ylabel(kwargs['ylabel'])
	else:
		plt.xlabel('H$_2$O wt%')
		plt.ylabel('CO$_2$ wt%')

	labels = []

	if isobars is not None:
		if isinstance(isobars, pd.DataFrame):
			isobars = [isobars]

		for i in range(len(isobars)):
			P_vals = isobars[i].Pressure.unique()
			isobars_lists = isobars[i].values.tolist()

			# add zero values to volatiles list
			isobars_lists.append([0.0, 0.0, 0.0, 0.0])

			P_iter = 0
			for pressure in P_vals:
				P_iter += 1
				Pxs = [item[1] for item in isobars_lists if item[0] == pressure]
				Pys = [item[2] for item in isobars_lists if item[0] == pressure]

				if len(isobars) > 1:
					if P_iter == 1:
						P_list = [int(i) for i in P_vals]
						if isinstance(isobar_labels, list):
							labels.append(str(isobar_labels[i]) + ' (' + ', '.join(map(str, P_list)) + " bars)")
						else:
							labels.append('Isobars ' + str(i+1) + ' (' + ', '.join(map(str, P_list)) + " bars)")
					else:
						labels.append('_nolegend_')
				if smooth_isobars == True:
				# do some data smoothing
					try:
						## calcualte polynomial
						Pz = np.polyfit(Pxs, Pys, 3)
						Pf = np.poly1d(Pz)

						## calculate new x's and y's
						Px_new = np.linspace(Pxs[0], Pxs[-1], 50)
						Py_new = Pf(Px_new)

						if extend_isobars_to_zero == True and Px_new[0]*Py_new[0] != 0.0:
							if Px_new[0] > Py_new[0]:
								Px_newer = np.zeros(np.shape(Px_new)[0]+1)
								Px_newer[0] = 0
								Px_newer[1:] = Px_new
								Px_new = Px_newer

								Py_newer = np.zeros(np.shape(Py_new)[0]+1)
								Py_newer[0] = Py_new[0]
								Py_newer[1:] = Py_new
								Py_new = Py_newer
							else:
								Px_newer = np.zeros(np.shape(Px_new)[0]+1)
								Px_newer[0] = Px_new[0]
								Px_newer[1:] = Px_new
								Px_new = Px_newer

								Py_newer = np.zeros(np.shape(Py_new)[0]+1)
								Py_newer[0] = 0
								Py_newer[1:] = Py_new
								Py_new = Py_newer

						if extend_isobars_to_zero == True and Px_new[-1]*Py_new[-1] != 0.0:
							if Px_new[-1] < Py_new[-1]:
								Px_newer = np.zeros(np.shape(Px_new)[0]+1)
								Px_newer[-1] = 0
								Px_newer[:-1] = Px_new
								Px_new = Px_newer

								Py_newer = np.zeros(np.shape(Py_new)[0]+1)
								Py_newer[-1] = Py_new[-1]
								Py_newer[:-1] = Py_new
								Py_new = Py_newer
							else:
								Px_newer = np.zeros(np.shape(Px_new)[0]+1)
								Px_newer[-1] = Px_new[-1]
								Px_newer[:-1] = Px_new
								Px_new = Px_newer

								Py_newer = np.zeros(np.shape(Py_new)[0]+1)
								Py_newer[-1] = 0
								Py_newer[:-1] = Py_new
								Py_new = Py_newer

						# Plot some stuff
						if len(isobars) > 1:
							plt.plot(Px_new, Py_new, color=color_list[i])
						else:
							plt.plot(Px_new, Py_new)
					except:
						if len(isobars) > 1:
							plt.plot(Pxs, Pys, color=color_list[i])
						else:
							plt.plot(Pxs, Pys)

				elif smooth_isobars == False:
					if extend_isobars_to_zero == True and Pxs[0]*Pys[0] != 0.0:
						if Pxs[0] > Pys[0]:
							Px_newer = np.zeros(np.shape(Pxs)[0]+1)
							Px_newer[0] = 0
							Px_newer[1:] = Pxs
							Pxs = Px_newer

							Py_newer = np.zeros(np.shape(Pys)[0]+1)
							Py_newer[0] = Pys[0]
							Py_newer[1:] = Pys
							Pys = Py_newer
						else:
							Px_newer = np.zeros(np.shape(Pxs)[0]+1)
							Px_newer[0] = Pxs[0]
							Px_newer[1:] = Pxs
							Pxs = Px_newer

							Py_newer = np.zeros(np.shape(Pys)[0]+1)
							Py_newer[0] = 0
							Py_newer[1:] = Pys
							Pys = Py_newer

					if extend_isobars_to_zero == True and Pxs[-1]*Pys[-1] != 0.0:
						if Pxs[-1] < Pys[-1]:
							Px_newer = np.zeros(np.shape(Pxs)[0]+1)
							Px_newer[-1] = 0
							Px_newer[:-1] = Pxs
							Pxs = Px_newer

							Py_newer = np.zeros(np.shape(Pys)[0]+1)
							Py_newer[-1] = Pys[-1]
							Py_newer[:-1] = Pys
							Pys = Py_newer
						else:
							Px_newer = np.zeros(np.shape(Pxs)[0]+1)
							Px_newer[-1] = Pxs[-1]
							Px_newer[:-1] = Pxs
							Pxs = Px_newer

							Py_newer = np.zeros(np.shape(Pys)[0]+1)
							Py_newer[-1] = 0
							Py_newer[:-1] = Pys
							Pys = Py_newer
					if len(isobars) > 1:
						plt.plot(Pxs, Pys, color=color_list[i])
					else:
						plt.plot(Pxs, Pys)

			if len(isobars) == 1:
				labels = [str(P_val) + " bars" for P_val in P_vals]

	if isopleths is not None:
		if isinstance(isopleths, pd.DataFrame):
			isopleths = [isopleths]

		for i in range(len(isopleths)):
			XH2O_vals = isopleths[i].XH2O_fl.unique()
			isopleths_lists = isopleths[i].values.tolist()

			H_iter = 0
			for Xfl in XH2O_vals:
				H_iter += 1
				Xxs = [item[1] for item in isopleths_lists if item[0] == Xfl]
				Xys = [item[2] for item in isopleths_lists if item[0] == Xfl]

				if len(isopleths) > 1:
					if H_iter == 1:
						H_list = [i for i in XH2O_vals]
						if isinstance(isopleth_labels, list):
							labels.append(str(isopleth_labels[i]) + ' (' + ', '.join(map(str, H_list)) + " XH2Ofluid)")
						else:
							labels.append('Isopleths ' + str(i+1) + ' (' + ', '.join(map(str, H_list)) + " XH2Ofluid)")
					else:
						labels.append('_nolegend_')
				if smooth_isopleths == True:
				# do some data smoothing
					try:
						## calcualte polynomial
						Xz = np.polyfit(Xxs, Xys, 2)
						Xf = np.poly1d(Xz)

						## calculate new x's and y's
						Xx_new = np.linspace(Xxs[0], Xxs[-1], 50)
						Xy_new = Xf(Xx_new)

						# Plot some stuff
						if len(isopleths) == 1:
							plt.plot(Xx_new, Xy_new, ls='dashed', color='k')
						else:
							plt.plot(Xx_new, Xy_new, ls='dashed', color=color_list[i])
					except:
						if len(isopleths) == 1:
							plt.plot(Xxs, Xys, ls='dashed', color='k')
						else:
							plt.plot(Xxs, Xys, ls='dashed', color=color_list[i])

				elif smooth_isopleths == False:
					if len(isopleths) == 1:
							plt.plot(Xxs, Xys, ls='dashed', color='k')
					else:
						plt.plot(Xxs, Xys, ls='dashed', color=color_list[i])

			if len(isopleths) == 1:
				H_list = [i for i in XH2O_vals]
				iso_label_iter = 0
				for i in XH2O_vals:
					iso_label_iter += 1
					if iso_label_iter == 1:
						labels.append('Isopleths (' + ', '.join(map(str, H_list)) + " XH2Ofluid)")
					else:
						labels.append('_nolegend_')

	if degassing_paths is not None:
		if isinstance(degassing_paths, pd.DataFrame):
			degassing_paths = [degassing_paths]

		degassing_colors = color_list.copy()
		#degassing_colors.reverse()
		iterno = 0
		for i in range(len(degassing_paths)):
			if degassing_path_labels == None:
				iterno += 1
				labels.append('Path%s' %iterno)
				plt.plot(degassing_paths[i]["H2O_liq"], degassing_paths[i]["CO2_liq"], ls='dotted', color=degassing_colors[i])
			else:
				labels.append(degassing_path_labels[iterno])
				plt.plot(degassing_paths[i]["H2O_liq"], degassing_paths[i]["CO2_liq"], ls='dotted', color=degassing_colors[i])
				iterno += 1

		for i in range(len(degassing_paths)):
			plt.plot(degassing_paths[i]["H2O_liq"].max(), degassing_paths[i]["CO2_liq"].max(), 'o', color=degassing_colors[i])
			labels.append('_nolegend_')

	if custom_H2O is not None and custom_CO2 is not None:
		if isinstance(custom_H2O, pd.DataFrame):
			custom_H2O = [custom_H2O]
		if isinstance(custom_CO2, pd.DataFrame):
			custom_CO2 = [custom_CO2]

		if custom_symbols == None:
			use_marker = ['o'] * len(custom_H2O)
		else:
			use_marker = custom_symbols

		iterno = 0
		for i in range(len(custom_H2O)):
			if custom_labels == None:
				iterno +=1
				labels.append('Custom%s' %iterno)
				plt.plot(custom_H2O[i], custom_CO2[i], use_marker[i], color=use_colors[i], markersize=markersize)
			else:
				labels.append(custom_labels[iterno])
				plt.plot(custom_H2O[i], custom_CO2[i], use_marker[i], color=use_colors[i], markersize=markersize)
				iterno += 1

	if 'custom_x' in kwargs:
			custom_x = kwargs['custom_x']
			custom_y = kwargs['custom_y']
			xlabel = kwargs['xlabel']
			ylabel = kwargs['ylabel']

			if isinstance(custom_x, pd.core.series.Series):
				custom_x = [list(custom_x.values)]
			if isinstance(custom_y, pd.core.series.Series):
				custom_y = [list(custom_y.values)]

			if custom_symbols == None:
				use_marker = ['o'] * len(custom_x)
			else:
				use_marker = custom_symbols

			iterno = 0
			for i in range(len(custom_x)):
				if custom_labels == None:
					iterno +=1
					labels.append('Custom%s' %iterno)
					plt.plot(custom_x[i], custom_y[i], use_marker[i], color=use_colors[i], markersize=markersize)
				else:
					labels.append(custom_labels[iterno])
					plt.plot(custom_x[i], custom_y[i], use_marker[i], color=use_colors[i], markersize=markersize)
					iterno += 1


	plt.legend(labels, bbox_to_anchor=(1.01,1), loc='upper left')

	if 'custom_x' not in kwargs:
		plt.xlim(left=0)
		plt.ylim(bottom=0)

	np.seterr(divide='warn', invalid='warn') #turn numpy warning back on
	w.filterwarnings("always", message="Polyfit may be poorly conditioned")

	return plt.show()

def scatterplot(custom_x, custom_y, xlabel=None, ylabel=None, **kwargs):
	"""
	Custom x-y plotting using VESIcal's built-in plot() function, built on Matplotlib's plot and scatter functions.

	Parameters
	----------
	custom_x: list
		List of groups of x-values to plot as points or lines

	custom_y: list
		List of groups of y-values to plot as points or lines

	xlabel: str
		OPTIONAL. What to display along the x-axis.

	ylabel: str
		OPTIONAL. What to display along the y-axis.

	kwargs:
		Can take in any key word agruments that can be passed to `plot()`.

	Returns
	-------
	matplotlib object
		X-y plot with custom x and y axis values and labels.
	"""

	if isinstance(custom_x, list) and isinstance(custom_y, list):
		if len(custom_x) != len(custom_y):
			raise InputError("X and y lists must be same length")

	if xlabel is not None:
		if isinstance(xlabel, str):
			pass
		else:
			raise InputError("xlabel must be string")

	if ylabel is not None:
		if isinstance(ylabel, str):
			pass
		else:
			raise InputError("ylabel must be string")

	return plot(custom_x=custom_x, custom_y=custom_y, xlabel=xlabel, ylabel=ylabel, **kwargs)

#====== Define some standard model options =======================================================#

default_models = {'Shishkina':                MixedFluid({'H2O':ShishkinaWater(),'CO2':ShishkinaCarbon()}),
				  'Dixon':                    MixedFluid({'H2O':DixonWater(),'CO2':DixonCarbon()}),
				  'IaconoMarziano':           MixedFluid({'H2O':IaconoMarzianoWater(),'CO2':IaconoMarzianoCarbon()}),
				  'Liu':					  MixedFluid({'H2O':LiuWater(),'CO2':LiuCarbon()}),
				  'ShishkinaCarbon':          ShishkinaCarbon(),
				  'ShishkinaWater':           ShishkinaWater(),
				  'DixonCarbon':              DixonCarbon(),
				  'DixonWater':               DixonWater(),
				  'IaconoMarzianoCarbon':     IaconoMarzianoCarbon(),
				  'IaconoMarzianoWater':      IaconoMarzianoWater(),
				  #'EguchiCarbon':             EguchiCarbon(),
				  'AllisonCarbon':			  AllisonCarbon(model_loc='vesuvius'),
				  'AllisonCarbon_sunset':     AllisonCarbon(model_loc='sunset'),
				  'AllisonCarbon_sfvf':       AllisonCarbon(model_loc='sfvf'),
				  'AllisonCarbon_erebus':     AllisonCarbon(model_loc='erebus'),
				  'AllisonCarbon_vesuvius':   AllisonCarbon(model_loc='vesuvius'),
				  'AllisonCarbon_etna':       AllisonCarbon(model_loc='etna'),
				  'AllisonCarbon_stromboli':  AllisonCarbon(model_loc='stromboli'),
				  'MooreWater':               MooreWater(),
				  'LiuWater':				  LiuWater(),
				  'LiuCarbon':				  LiuCarbon()
}

# Return homogenised calibration range checks for mixed fluid models:
default_models['Liu'].models[1].set_calibration_ranges([])

_crs_to_update = default_models['Liu'].models[0].calibration_ranges
for _cr in _crs_to_update:
	_cr.model_name = 'Liu et al. (2005)'
default_models['Liu'].models[0].set_calibration_ranges([CalibrationRange('pressure',[0,5000.0],crf_Between,'bar','Liu et al. (2005)',
													  fail_msg=crmsg_Between_fail, pass_msg=crmsg_Between_pass, description_msg=crmsg_Between_description),
									 CalibrationRange('temperature',[700.0,1200],crf_Between,'oC','Liu et al. (2005)',
									 				  fail_msg=crmsg_Between_fail, pass_msg=crmsg_Between_pass, description_msg=crmsg_Between_description),
									 CalibrationRange('SiO2',[Liu_SiO2Min,Liu_SiO2Max],crf_Between,'wt%','Liu et al. (2005)',
													  fail_msg=crmsg_BC_fail, pass_msg=crmsg_BC_pass, description_msg=crmsg_Between_description),
									 CalibrationRange('TiO2',[Liu_TiO2Min,Liu_TiO2Max],crf_Between,'wt%','Liu et al. (2005)',
													  fail_msg=crmsg_BC_fail, pass_msg=crmsg_BC_pass, description_msg=crmsg_Between_description),
									 CalibrationRange('Al2O3',[Liu_Al2O3Min,Liu_Al2O3Max],crf_Between,'wt%','Liu et al. (2005)',
													  fail_msg=crmsg_BC_fail, pass_msg=crmsg_BC_pass, description_msg=crmsg_Between_description),
									 CalibrationRange('FeO',[Liu_FeOMin,Liu_FeOMax],crf_Between,'wt%','Liu et al. (2005)',
													  fail_msg=crmsg_BC_fail, pass_msg=crmsg_BC_pass, description_msg=crmsg_Between_description),
									 CalibrationRange('MgO',[Liu_MgOMin,Liu_MgOMax],crf_Between,'wt%','Liu et al. (2005)',
													  fail_msg=crmsg_BC_fail, pass_msg=crmsg_BC_pass, description_msg=crmsg_Between_description),
									 CalibrationRange('CaO',[Liu_CaOMin,Liu_CaOMax],crf_Between,'wt%','Liu et al. (2005)',
													  fail_msg=crmsg_BC_fail, pass_msg=crmsg_BC_pass, description_msg=crmsg_Between_description),
									 CalibrationRange('Na2O',[Liu_Na2OMin,Liu_Na2OMax],crf_Between,'wt%','Liu et al. (2005)',
													  fail_msg=crmsg_BC_fail, pass_msg=crmsg_BC_pass, description_msg=crmsg_Between_description),
									 CalibrationRange('K2O',[Liu_K2OMin,Liu_K2OMax],crf_Between,'wt%','Liu et al. (2005)',
													  fail_msg=crmsg_BC_fail, pass_msg=crmsg_BC_pass, description_msg=crmsg_Between_description),
									 CalibrationRange('sample',None,crf_LiuComp,None,None,
									 				  fail_msg=crmsg_LiuComp_fail, pass_msg=crmsg_LiuComp_pass, description_msg=crmsg_LiuComp_description)])




_crs_to_update = default_models['AllisonCarbon_sunset'].calibration_ranges
default_models['AllisonCarbon_sunset'].set_calibration_ranges(_crs_to_update+
				 [CalibrationRange('pressure',7000,crf_LessThan,'bar','Allison et al. (2019) sunset carbon',
													  fail_msg=crmsg_LessThan_fail_Allison, pass_msg=crmsg_LessThan_pass),
				 CalibrationRange('pressure',4071,crf_GreaterThan,'bar','Allison et al. (2019) sunset carbon',
													  fail_msg=crmsg_GreaterThan_fail_Allison, pass_msg=crmsg_GreaterThan_pass),
				 CalibrationRange('H2O',0.5,crf_LessThan,'wt%','Allison et al. (2019) sunset carbon',
													  fail_msg=crmsg_LessThan_fail_AllisonH2O, pass_msg=crmsg_LessThan_pass),
				 CalibrationRange('SiO2',[Allison_SiO2_sunset[0],Allison_SiO2_sunset[1]],crf_Between,'wt%','Allison et al. (2019) sunset carbon',
													  fail_msg=crmsg_BC_fail, pass_msg=crmsg_BC_pass, description_msg=crmsg_Between_description),
				 CalibrationRange('TiO2',[Allison_TiO2_sunset[0],Allison_TiO2_sunset[1]],crf_Between,'wt%','Allison et al. (2019) sunset carbon',
													  fail_msg=crmsg_BC_fail, pass_msg=crmsg_BC_pass, description_msg=crmsg_Between_description),
				 CalibrationRange('Al2O3',[Allison_Al2O3_sunset[0],Allison_Al2O3_sunset[1]],crf_Between,'wt%','Allison et al. (2019) sunset carbon',
													  fail_msg=crmsg_BC_fail, pass_msg=crmsg_BC_pass, description_msg=crmsg_Between_description),
				 CalibrationRange('FeO',[Allison_FeO_sunset[0],Allison_FeO_sunset[1]],crf_Between,'wt%','Allison et al. (2019) sunset carbon',
													  fail_msg=crmsg_BC_fail, pass_msg=crmsg_BC_pass, description_msg=crmsg_Between_description),
				 CalibrationRange('MgO',[Allison_MgO_sunset[0],Allison_MgO_sunset[1]],crf_Between,'wt%','Allison et al. (2019) sunset carbon',
													  fail_msg=crmsg_BC_fail, pass_msg=crmsg_BC_pass, description_msg=crmsg_Between_description),
				 CalibrationRange('CaO',[Allison_CaO_sunset[0],Allison_CaO_sunset[1]],crf_Between,'wt%','Allison et al. (2019) sunset carbon',
													  fail_msg=crmsg_BC_fail, pass_msg=crmsg_BC_pass, description_msg=crmsg_Between_description),
				 CalibrationRange('Na2O',[Allison_Na2O_sunset[0],Allison_Na2O_sunset[1]],crf_Between,'wt%','Allison et al. (2019) sunset carbon',
													  fail_msg=crmsg_BC_fail, pass_msg=crmsg_BC_pass, description_msg=crmsg_Between_description),
				 CalibrationRange('K2O',[Allison_K2O_sunset[0],Allison_K2O_sunset[1]],crf_Between,'wt%','Allison et al. (2019) sunset carbon',
													  fail_msg=crmsg_BC_fail, pass_msg=crmsg_BC_pass, description_msg=crmsg_Between_description),
				CalibrationRange('sample',None,crf_AllisonComp_sunset,None,None,
									 				  fail_msg=crmsg_AllisonComp_fail, pass_msg=crmsg_AllisonComp_pass)])


_crs_to_update = default_models['AllisonCarbon_sfvf'].calibration_ranges
default_models['AllisonCarbon_sfvf'].set_calibration_ranges(_crs_to_update+
				 [CalibrationRange('pressure',7000,crf_LessThan,'bar','Allison et al. (2019) sfvf carbon',
													  fail_msg=crmsg_LessThan_fail_Allison, pass_msg=crmsg_LessThan_pass),
				 CalibrationRange('pressure',4133,crf_GreaterThan,'bar','Allison et al. (2019) sfvf carbon',
													  fail_msg=crmsg_GreaterThan_fail_Allison, pass_msg=crmsg_GreaterThan_pass),
				 CalibrationRange('H2O',0.5,crf_LessThan,'wt%','Allison et al. (2019) sfvf carbon',
													  fail_msg=crmsg_LessThan_fail_AllisonH2O, pass_msg=crmsg_LessThan_pass),
				 CalibrationRange('SiO2',[Allison_SiO2_sfvf[0],Allison_SiO2_sfvf[1]],crf_Between,'wt%','Allison et al. (2019) sfvf carbon',
													  fail_msg=crmsg_BC_fail, pass_msg=crmsg_BC_pass, description_msg=crmsg_Between_description),
				 CalibrationRange('TiO2',[Allison_TiO2_sfvf[0],Allison_TiO2_sfvf[1]],crf_Between,'wt%','Allison et al. (2019) sfvf carbon',
													  fail_msg=crmsg_BC_fail, pass_msg=crmsg_BC_pass, description_msg=crmsg_Between_description),
				 CalibrationRange('Al2O3',[Allison_Al2O3_sfvf[0],Allison_Al2O3_sfvf[1]],crf_Between,'wt%','Allison et al. (2019) sfvf carbon',
													  fail_msg=crmsg_BC_fail, pass_msg=crmsg_BC_pass, description_msg=crmsg_Between_description),
				 CalibrationRange('FeO',[Allison_FeO_sfvf[0],Allison_FeO_sfvf[1]],crf_Between,'wt%','Allison et al. (2019) sfvf carbon',
													  fail_msg=crmsg_BC_fail, pass_msg=crmsg_BC_pass, description_msg=crmsg_Between_description),
				 CalibrationRange('MgO',[Allison_MgO_sfvf[0],Allison_MgO_sfvf[1]],crf_Between,'wt%','Allison et al. (2019) sfvf carbon',
													  fail_msg=crmsg_BC_fail, pass_msg=crmsg_BC_pass, description_msg=crmsg_Between_description),
				 CalibrationRange('CaO',[Allison_CaO_sfvf[0],Allison_CaO_sfvf[1]],crf_Between,'wt%','Allison et al. (2019) sfvf carbon',
													  fail_msg=crmsg_BC_fail, pass_msg=crmsg_BC_pass, description_msg=crmsg_Between_description),
				 CalibrationRange('Na2O',[Allison_Na2O_sfvf[0],Allison_Na2O_sfvf[1]],crf_Between,'wt%','Allison et al. (2019) sfvf carbon',
													  fail_msg=crmsg_BC_fail, pass_msg=crmsg_BC_pass, description_msg=crmsg_Between_description),
				 CalibrationRange('K2O',[Allison_K2O_sfvf[0],Allison_K2O_sfvf[1]],crf_Between,'wt%','Allison et al. (2019) sfvf carbon',
													  fail_msg=crmsg_BC_fail, pass_msg=crmsg_BC_pass, description_msg=crmsg_Between_description),
				CalibrationRange('sample',None,crf_AllisonComp_sfvf,None,None,
									 				  fail_msg=crmsg_AllisonComp_fail, pass_msg=crmsg_AllisonComp_pass)])

_crs_to_update = default_models['AllisonCarbon_erebus'].calibration_ranges
default_models['AllisonCarbon_erebus'].set_calibration_ranges(_crs_to_update+
				 [CalibrationRange('pressure',7000,crf_LessThan,'bar','Allison et al. (2019) erebus carbon',
													  fail_msg=crmsg_LessThan_fail_Allison, pass_msg=crmsg_LessThan_pass),
				 CalibrationRange('pressure',4078,crf_GreaterThan,'bar','Allison et al. (2019) erebus carbon',
													  fail_msg=crmsg_GreaterThan_fail_Allison, pass_msg=crmsg_GreaterThan_pass),
				 CalibrationRange('H2O',0.5,crf_LessThan,'wt%','Allison et al. (2019) erebus carbon',
													  fail_msg=crmsg_LessThan_fail_AllisonH2O, pass_msg=crmsg_LessThan_pass),
				 CalibrationRange('SiO2',[Allison_SiO2_erebus[0],Allison_SiO2_erebus[1]],crf_Between,'wt%','Allison et al. (2019) erebus carbon',
													  fail_msg=crmsg_BC_fail, pass_msg=crmsg_BC_pass, description_msg=crmsg_Between_description),
				 CalibrationRange('TiO2',[Allison_TiO2_erebus[0],Allison_TiO2_erebus[1]],crf_Between,'wt%','Allison et al. (2019) erebus carbon',
													  fail_msg=crmsg_BC_fail, pass_msg=crmsg_BC_pass, description_msg=crmsg_Between_description),
				 CalibrationRange('Al2O3',[Allison_Al2O3_erebus[0],Allison_Al2O3_erebus[1]],crf_Between,'wt%','Allison et al. (2019) erebus carbon',
													  fail_msg=crmsg_BC_fail, pass_msg=crmsg_BC_pass, description_msg=crmsg_Between_description),
				 CalibrationRange('FeO',[Allison_FeO_erebus[0],Allison_FeO_erebus[1]],crf_Between,'wt%','Allison et al. (2019) erebus carbon',
													  fail_msg=crmsg_BC_fail, pass_msg=crmsg_BC_pass, description_msg=crmsg_Between_description),
				 CalibrationRange('MgO',[Allison_MgO_erebus[0],Allison_MgO_erebus[1]],crf_Between,'wt%','Allison et al. (2019) erebus carbon',
													  fail_msg=crmsg_BC_fail, pass_msg=crmsg_BC_pass, description_msg=crmsg_Between_description),
				 CalibrationRange('CaO',[Allison_CaO_erebus[0],Allison_CaO_erebus[1]],crf_Between,'wt%','Allison et al. (2019) erebus carbon',
													  fail_msg=crmsg_BC_fail, pass_msg=crmsg_BC_pass, description_msg=crmsg_Between_description),
				 CalibrationRange('Na2O',[Allison_Na2O_erebus[0],Allison_Na2O_erebus[1]],crf_Between,'wt%','Allison et al. (2019) erebus carbon',
													  fail_msg=crmsg_BC_fail, pass_msg=crmsg_BC_pass, description_msg=crmsg_Between_description),
				 CalibrationRange('K2O',[Allison_K2O_erebus[0],Allison_K2O_erebus[1]],crf_Between,'wt%','Allison et al. (2019) erebus carbon',
													  fail_msg=crmsg_BC_fail, pass_msg=crmsg_BC_pass, description_msg=crmsg_Between_description),
				CalibrationRange('sample',None,crf_AllisonComp_erebus,None,None,
									 				  fail_msg=crmsg_AllisonComp_fail, pass_msg=crmsg_AllisonComp_pass)])


_crs_to_update = default_models['AllisonCarbon_etna'].calibration_ranges
default_models['AllisonCarbon_etna'].set_calibration_ranges(_crs_to_update+
				 [CalibrationRange('pressure',7000,crf_LessThan,'bar','Allison et al. (2019) etna carbon',
													  fail_msg=crmsg_LessThan_fail_Allison, pass_msg=crmsg_LessThan_pass),
				 CalibrationRange('pressure',485,crf_GreaterThan,'bar','Allison et al. (2019) etna carbon',
													  fail_msg=crmsg_GreaterThan_fail_Allison, pass_msg=crmsg_GreaterThan_pass),
				 CalibrationRange('H2O',0.5,crf_LessThan,'wt%','Allison et al. (2019) etna carbon',
													  fail_msg=crmsg_LessThan_fail_AllisonH2O, pass_msg=crmsg_LessThan_pass),
				 CalibrationRange('SiO2',[Allison_SiO2_etna[0],Allison_SiO2_etna[1]],crf_Between,'wt%','Allison et al. (2019) etna carbon',
													  fail_msg=crmsg_BC_fail, pass_msg=crmsg_BC_pass, description_msg=crmsg_Between_description),
				 CalibrationRange('TiO2',[Allison_TiO2_etna[0],Allison_TiO2_etna[1]],crf_Between,'wt%','Allison et al. (2019) etna carbon',
													  fail_msg=crmsg_BC_fail, pass_msg=crmsg_BC_pass, description_msg=crmsg_Between_description),
				 CalibrationRange('Al2O3',[Allison_Al2O3_etna[0],Allison_Al2O3_etna[1]],crf_Between,'wt%','Allison et al. (2019) etna carbon',
													  fail_msg=crmsg_BC_fail, pass_msg=crmsg_BC_pass, description_msg=crmsg_Between_description),
				 CalibrationRange('FeO',[Allison_FeO_etna[0],Allison_FeO_etna[1]],crf_Between,'wt%','Allison et al. (2019) etna carbon',
													  fail_msg=crmsg_BC_fail, pass_msg=crmsg_BC_pass, description_msg=crmsg_Between_description),
				 CalibrationRange('MgO',[Allison_MgO_etna[0],Allison_MgO_etna[1]],crf_Between,'wt%','Allison et al. (2019) etna carbon',
													  fail_msg=crmsg_BC_fail, pass_msg=crmsg_BC_pass, description_msg=crmsg_Between_description),
				 CalibrationRange('CaO',[Allison_CaO_etna[0],Allison_CaO_etna[1]],crf_Between,'wt%','Allison et al. (2019) etna carbon',
													  fail_msg=crmsg_BC_fail, pass_msg=crmsg_BC_pass, description_msg=crmsg_Between_description),
				 CalibrationRange('Na2O',[Allison_Na2O_etna[0],Allison_Na2O_etna[1]],crf_Between,'wt%','Allison et al. (2019) etna carbon',
													  fail_msg=crmsg_BC_fail, pass_msg=crmsg_BC_pass, description_msg=crmsg_Between_description),
				 CalibrationRange('K2O',[Allison_K2O_etna[0],Allison_K2O_etna[1]],crf_Between,'wt%','Allison et al. (2019) etna carbon',
													  fail_msg=crmsg_BC_fail, pass_msg=crmsg_BC_pass, description_msg=crmsg_Between_description),
				CalibrationRange('sample',None,crf_AllisonComp_etna,None,None,
									 				  fail_msg=crmsg_AllisonComp_fail, pass_msg=crmsg_AllisonComp_pass)])

_crs_to_update = default_models['AllisonCarbon_vesuvius'].calibration_ranges
default_models['AllisonCarbon_vesuvius'].set_calibration_ranges(_crs_to_update+
				 [CalibrationRange('pressure',7000,crf_LessThan,'bar','Allison et al. (2019) vesuvius carbon',
													  fail_msg=crmsg_LessThan_fail_Allison, pass_msg=crmsg_LessThan_pass),
				 CalibrationRange('pressure',269,crf_GreaterThan,'bar','Allison et al. (2019) vesuvius carbon',
													  fail_msg=crmsg_GreaterThan_fail_Allison, pass_msg=crmsg_GreaterThan_pass),
				 CalibrationRange('H2O',0.5,crf_LessThan,'wt%','Allison et al. (2019) vesuvius carbon',
													  fail_msg=crmsg_LessThan_fail_AllisonH2O, pass_msg=crmsg_LessThan_pass),
				 CalibrationRange('SiO2',[Allison_SiO2_vesuvius[0],Allison_SiO2_vesuvius[1]],crf_Between,'wt%','Allison et al. (2019) vesuvius carbon',
													  fail_msg=crmsg_BC_fail, pass_msg=crmsg_BC_pass, description_msg=crmsg_Between_description),
				 CalibrationRange('TiO2',[Allison_TiO2_vesuvius[0],Allison_TiO2_vesuvius[1]],crf_Between,'wt%','Allison et al. (2019) vesuvius carbon',
													  fail_msg=crmsg_BC_fail, pass_msg=crmsg_BC_pass, description_msg=crmsg_Between_description),
				 CalibrationRange('Al2O3',[Allison_Al2O3_vesuvius[0],Allison_Al2O3_vesuvius[1]],crf_Between,'wt%','Allison et al. (2019) vesuvius carbon',
													  fail_msg=crmsg_BC_fail, pass_msg=crmsg_BC_pass, description_msg=crmsg_Between_description),
				 CalibrationRange('FeO',[Allison_FeO_vesuvius[0],Allison_FeO_vesuvius[1]],crf_Between,'wt%','Allison et al. (2019) vesuvius carbon',
													  fail_msg=crmsg_BC_fail, pass_msg=crmsg_BC_pass, description_msg=crmsg_Between_description),
				 CalibrationRange('MgO',[Allison_MgO_vesuvius[0],Allison_MgO_vesuvius[1]],crf_Between,'wt%','Allison et al. (2019) vesuvius carbon',
													  fail_msg=crmsg_BC_fail, pass_msg=crmsg_BC_pass, description_msg=crmsg_Between_description),
				 CalibrationRange('CaO',[Allison_CaO_vesuvius[0],Allison_CaO_vesuvius[1]],crf_Between,'wt%','Allison et al. (2019) vesuvius carbon',
													  fail_msg=crmsg_BC_fail, pass_msg=crmsg_BC_pass, description_msg=crmsg_Between_description),
				 CalibrationRange('Na2O',[Allison_Na2O_vesuvius[0],Allison_Na2O_vesuvius[1]],crf_Between,'wt%','Allison et al. (2019) vesuvius carbon',
													  fail_msg=crmsg_BC_fail, pass_msg=crmsg_BC_pass, description_msg=crmsg_Between_description),
				 CalibrationRange('K2O',[Allison_K2O_vesuvius[0],Allison_K2O_vesuvius[1]],crf_Between,'wt%','Allison et al. (2019) vesuvius carbon',
													  fail_msg=crmsg_BC_fail, pass_msg=crmsg_BC_pass, description_msg=crmsg_Between_description),
				CalibrationRange('sample',None,crf_AllisonComp_vesuvius,None,None,
									 				  fail_msg=crmsg_AllisonComp_fail, pass_msg=crmsg_AllisonComp_pass)])

_crs_to_update = default_models['AllisonCarbon_stromboli'].calibration_ranges
default_models['AllisonCarbon_stromboli'].set_calibration_ranges(_crs_to_update+
				 [CalibrationRange('pressure',7000,crf_LessThan,'bar','Allison et al. (2019) stromboli carbon',
													  fail_msg=crmsg_LessThan_fail_Allison, pass_msg=crmsg_LessThan_pass),
				 CalibrationRange('pressure',524,crf_GreaterThan,'bar','Allison et al. (2019) stromboli carbon',
													  fail_msg=crmsg_GreaterThan_fail_Allison, pass_msg=crmsg_GreaterThan_pass),
				 CalibrationRange('H2O',0.5,crf_LessThan,'wt%','Allison et al. (2019) stromboli carbon',
													  fail_msg=crmsg_LessThan_fail_AllisonH2O, pass_msg=crmsg_LessThan_pass),
				 CalibrationRange('SiO2',[Allison_SiO2_stromboli[0],Allison_SiO2_stromboli[1]],crf_Between,'wt%','Allison et al. (2019) stromboli carbon',
													  fail_msg=crmsg_BC_fail, pass_msg=crmsg_BC_pass, description_msg=crmsg_Between_description),
				 CalibrationRange('TiO2',[Allison_TiO2_stromboli[0],Allison_TiO2_stromboli[1]],crf_Between,'wt%','Allison et al. (2019) stromboli carbon',
													  fail_msg=crmsg_BC_fail, pass_msg=crmsg_BC_pass, description_msg=crmsg_Between_description),
				 CalibrationRange('Al2O3',[Allison_Al2O3_stromboli[0],Allison_Al2O3_stromboli[1]],crf_Between,'wt%','Allison et al. (2019) stromboli carbon',
													  fail_msg=crmsg_BC_fail, pass_msg=crmsg_BC_pass, description_msg=crmsg_Between_description),
				 CalibrationRange('FeO',[Allison_FeO_stromboli[0],Allison_FeO_stromboli[1]],crf_Between,'wt%','Allison et al. (2019) stromboli carbon',
													  fail_msg=crmsg_BC_fail, pass_msg=crmsg_BC_pass, description_msg=crmsg_Between_description),
				 CalibrationRange('MgO',[Allison_MgO_stromboli[0],Allison_MgO_stromboli[1]],crf_Between,'wt%','Allison et al. (2019) stromboli carbon',
													  fail_msg=crmsg_BC_fail, pass_msg=crmsg_BC_pass, description_msg=crmsg_Between_description),
				 CalibrationRange('CaO',[Allison_CaO_stromboli[0],Allison_CaO_stromboli[1]],crf_Between,'wt%','Allison et al. (2019) stromboli carbon',
													  fail_msg=crmsg_BC_fail, pass_msg=crmsg_BC_pass, description_msg=crmsg_Between_description),
				 CalibrationRange('Na2O',[Allison_Na2O_stromboli[0],Allison_Na2O_stromboli[1]],crf_Between,'wt%','Allison et al. (2019) stromboli carbon',
													  fail_msg=crmsg_BC_fail, pass_msg=crmsg_BC_pass, description_msg=crmsg_Between_description),
				 CalibrationRange('K2O',[Allison_K2O_stromboli[0],Allison_K2O_stromboli[1]],crf_Between,'wt%','Allison et al. (2019) stromboli carbon',
													  fail_msg=crmsg_BC_fail, pass_msg=crmsg_BC_pass, description_msg=crmsg_Between_description),
				CalibrationRange('sample',None,crf_AllisonComp_stromboli,None,None,
									 				  fail_msg=crmsg_AllisonComp_fail, pass_msg=crmsg_AllisonComp_pass)])

def get_models(models='all'):
	"""
	Returns model names as a list
	Parameters
	----------
	models:	str
		OPTIONAL. Default value is 'all' in which case all keys in defaule_models are returned.
		If 'mixed' is passed, only the MixedFluid model names are returned.
	"""
	if models == 'all':
		return list(default_models.keys())
	if models == 'mixed':
		return ['Shishkina', 'Dixon', 'IaconoMarziano', 'Liu']

class calculate_dissolved_volatiles(Calculate):
	""" Calculates the dissolved volatile concentration using a chosen model (default is MagmaSat).
	Using this interface will preprocess the sample, run the calculation, and then check
	the calibration ranges. All parameters required by the chosen model must be passed.

	Parameters
	----------
	sample:     dict or pandas Series
		The major element oxides in wt%.
	pressure:   float
		Total pressure in bars.
	model:  string or Model object
		Model to be used. If using one of the default models, this can be
		the string corresponding to the model in the default_models dict.
	silence_warnings 	bool
		If set to True, no warnings will be raised automatically when calibration checks fail.
	preprocess_sample 	bool
		If True (default), the sample will be preprocessed according to the preprocessing operations within
		the models. If you obtain unexpected results, try setting to False.

	Returns
	-------
	Calculate object
		Calculate object, access results by fetching the result property. Dissolved
		volatile concentrations (in wt%), in order (CO2, H2O, if using a mixed fluid
		default model).
	"""
	def calculate(self,sample,pressure,**kwargs):
		dissolved = self.model.calculate_dissolved_volatiles(pressure=pressure,sample=sample,returndict=True,**kwargs)
		return dissolved

	def check_calibration_range(self,sample,pressure,**kwargs):
		parameters = kwargs
		parameters['sample'] = sample
		parameters.update(dict(sample))
		parameters['pressure'] = pressure
		if len(self.model.volatile_species) == 1:
			volspec = self.model.volatile_species[0]
			volconc = self.result
			parameters.update({volspec:volconc})
		else:
			 parameters.update(self.result)

		calib_check = self.model.check_calibration_range(parameters)
		return calib_check

class calculate_equilibrium_fluid_comp(Calculate):
	""" Calculates the equilibrium fluid composition using a chosen model (default is MagmaSat).
	Using this interface will preprocess the sample, run the calculation, and then check
	the calibration ranges. All parameters required by the chosen model must be passed.

	Parameters
	----------
	sample:     dict or pandas Series
		The major element oxides in wt%.
	pressure:   float or None
		Total pressure in bars. If None, the saturation pressure will be used.
	model:  string or Model object
		Model to be used. If using one of the default models, this can be
		the string corresponding to the model in the default_models dict.
	silence_warnings 	bool
		If set to True, no warnings will be raised automatically when calibration checks fail.
	preprocess_sample 	bool
		If True (default), the sample will be preprocessed according to the preprocessing operations within
		the models. If you obtain unexpected results, try setting to False.

	Returns
	-------
	Calculate object
		Calculate object, access result by fetching the result property. Mole fractions
		of each volatile species, in order (CO2, then H2O, if using a mixed-fluid default
		model).
	"""
	def calculate(self,sample,pressure=None,**kwargs):
		if type(pressure) == type(None):
			pressure = float(self.model.calculate_saturation_pressure(sample=sample,verbose=False,**kwargs))

		fluid_comp = self.model.calculate_equilibrium_fluid_comp(pressure=pressure,sample=sample,**kwargs)
		return fluid_comp

	def check_calibration_range(self,sample,pressure=None,**kwargs):
		if type(pressure) == type(None):
			pressure = self.model.calculate_saturation_pressure(sample=sample,**kwargs)
		parameters = kwargs
		parameters.update(dict(sample))
		parameters['sample'] = sample
		parameters['pressure'] = pressure
		if len(self.model.volatile_species) == 1:
			volspec = self.model.volatile_species
			volconc = {volspec[0]:self.result}
			parameters.update(volconc)
		elif type(self.model.volatile_species) == list:
			 parameters.update(self.result)

		calib_check = self.model.check_calibration_range(parameters)
		return calib_check


class calculate_isobars_and_isopleths(Calculate):
	""" Calculates isobars and isopleths using a chosen model (default is MagmaSat).
	Using this interface will preprocess the sample, run the calculation, and then check
	the calibration ranges. All parameters required by the chosen model must be passed.

	Parameters
	----------
	sample:     dict or pandas Series
		The major element oxides in wt%.
	pressure_list:   list
		List of all pressure values at which to calculate isobars, in bars.
	isopleth_list:   list
		OPTIONAL: Default value is None, in which case only isobars will be calculated. List of all
		fluid compositions in mole fraction (of the first species in self.volatile_species) at which
		to calcualte isopleths. Values can range from 0 to 1.
	points:     int
		The number of points in each isobar and isopleth. Default value is 101.
	model:  string or Model object
		Model to be used. If using one of the default models, this can be
		the string corresponding to the model in the default_models dict.
	silence_warnings 	bool
		If set to True, no warnings will be raised automatically when calibration checks fail.
	preprocess_sample 	bool
		If True (default), the sample will be preprocessed according to the preprocessing operations within
		the models. If you obtain unexpected results, try setting to False.

	Returns
	-------
	Calculate object
		Calculate object, access results by fetching the result property.
		If isopleth_list is not None, two objects will be returned, one with the isobars and the second with
		the isopleths. If return_dfs is True, two pandas DataFrames will be returned with column names
		'Pressure' or 'XH2O_fl', 'H2O_liq', and 'CO2_liq'. If return_dfs is False, two lists of numpy arrays
		will be returned. Each array is an individual isobar or isopleth, in the order passed via pressure_list
		or isopleth_list. The arrays are the concentrations of H2O and CO2 in the liquid, in the order of the
		species in self.volatile_species.
	"""
	def calculate(self,sample,pressure_list,isopleth_list=[0,1],points=101,**kwargs):
		check = getattr(self.model, "calculate_isobars_and_isopleths", None)
		if callable(check):
			samplenorm = sample.copy()
			samplenorm = normalize_AdditionalVolatiles(samplenorm)
			isobars, isopleths = self.model.calculate_isobars_and_isopleths(sample=samplenorm,pressure_list=pressure_list,isopleth_list=isopleth_list,points=points,**kwargs)
			return isobars, isopleths
		else:
			raise InputError("This model does not have a calculate_isobars_and_isopleths method built in, most likely because it is a pure fluid model.")

	def check_calibration_range(self,sample,pressure_list,**kwargs):
		parameters = kwargs
		parameters.update(dict(sample))
		parameters['sample'] = sample
		s = ''
		s += self.model.check_calibration_range(parameters)
		parameters = {}
		if isinstance(pressure_list, list):
			pass
		else:
			pressure_list = [pressure_list]
		for pressure in pressure_list:
			parameters['pressure'] = pressure
			s += self.model.check_calibration_range(parameters,report_nonexistance=False)
		return s


class calculate_saturation_pressure(Calculate):
	"""
	Calculates the pressure at which a fluid will be saturated, given the dissolved volatile
	concentrations. Using this interface will preprocess the sample, run the calculation, and then check
	the calibration ranges. All parameters required by the chosen model must be passed.

	Parameters
	----------
	sample     pandas Series or dict
		Major element oxides in wt% (including volatiles).
	model:  string or Model object
		Model to be used. If using one of the default models, this can be
		the string corresponding to the model in the default_models dict.
	silence_warnings 	bool
		If set to True, no warnings will be raised automatically when calibration checks fail.
	preprocess_sample 	bool
		If True (default), the sample will be preprocessed according to the preprocessing operations within
		the models. If you obtain unexpected results, try setting to False.

	Returns
	-------
	Calculate object
		Calculate object, access results by fetching the result property.
		The saturation pressure in bars as a float.
	"""
	def calculate(self,sample,**kwargs):
		satP = self.model.calculate_saturation_pressure(sample=sample,**kwargs)
		return satP

	def check_calibration_range(self,sample,**kwargs):
		parameters = kwargs
		if isinstance(self.result, dict): #handles cases where verbose=True
			parameters['pressure'] = next(iter(self.result.values()))
		else:
			parameters['pressure'] = self.result
		parameters.update(dict(sample))
		parameters['sample'] = sample
		s = self.model.check_calibration_range(parameters)
		return s

class calculate_degassing_path(Calculate):
	"""
	Calculates the dissolved volatiles in a progressively degassing sample.

	Parameters
	----------
	sample     pandas Series or dict
		Major element oxides in wt% (including volatiles).
	pressure     string, float, int, list, or numpy array
		Defaults to 'saturation', the calculation will begin at the saturation pressure. If a number is passed
		as either a float or int, this will be the starting pressure. If a list of numpy array is passed, the
		pressure values in the list or array will define the degassing path, i.e. final_pressure and steps
		variables will be ignored. Units are bars.
	fractionate_vapor     float
		What proportion of vapor should be removed at each step. If 0.0 (default), the degassing path will
		correspond to closed-system degassing. If 1.0, the degassing path will correspond to open-system
		degassing.
	final_pressure         float
		The final pressure on the degassing path, in bars. Ignored if a list or numpy array is passed as the
		pressure variable. Default is 1 bar.
	steps     int
		The number of steps in the degassing path. Ignored if a list or numpy array are passed as the pressure
		variable.
	model:  string or Model object
		Model to be used. If using one of the default models, this can be
		the string corresponding to the model in the default_models dict.
	silence_warnings 	bool
		If set to True, no warnings will be raised automatically when calibration checks fail.
	preprocess_sample 	bool
		If True (default), the sample will be preprocessed according to the preprocessing operations within
		the models. If you obtain unexpected results, try setting to False.

	Returns
	-------
	Calculate object
		Calculate object, access results by fetching the result property.
		A DataFrame with columns 'Pressure', 'H2O_liq', 'CO2_liq',
		'H2O_fl', 'CO2_fl', and 'FluidProportion_wt', is returned. Dissolved volatiles are in wt%,
		the proportions of volatiles in the fluid are in mole fraction.
	"""
	def calculate(self,sample,pressure='saturation',fractionate_vapor=0.0,
				  final_pressure=100.0,**kwargs):
		check = getattr(self.model, "calculate_degassing_path", None)
		if callable(check):
			data = self.model.calculate_degassing_path(sample=sample, pressure=pressure, fractionate_vapor=fractionate_vapor,
													   final_pressure=final_pressure, **kwargs)
			return data
		else:
			raise InputError("This model does not have a calculate_isobars_and_isopleths method built in, most likely because it is a pure fluid model.")

	def check_calibration_range(self,sample,**kwargs):
		parameters = kwargs
		parameters.update(sample)
		parameters['sample'] = sample
		s = self.model.check_calibration_range(parameters)
		parameters = {}
		parameters['pressure'] = np.nanmax(self.result.Pressure_bars)
		# if istype(kwargs.get("pressure"), float) or istype(kwargs.get("pressure"), int):
		# 	parameters['pressure'] = kwargs.get("pressure")
		s += self.model.check_calibration_range(parameters,report_nonexistance=False)
		return s

#-------Define custom plotting tools for checking calibrations-------#
#----------------------------------------------------------#
#    			  TAS PLOT PYTHON SCRIPT        	       #
#														   #
#  COPYRIGHT:  (C) 2015 John A Stevenson / @volcan01010    #
#                       Joaquin Cortés					   #
#  WEBSITE: http://all-geo.org/volcan01010				   #
#----------------------------------------------------------#
def add_LeMaitre_fields(plot_axes, fontsize=12, color=(0.6, 0.6, 0.6)):
	"""Add fields for geochemical classifications from LeMaitre et al (2002)
	to pre-existing axes.  If necessary, the axes object can be retrieved via
	plt.gca() command. e.g.

	ax1 = plt.gca()
	add_LeMaitre_fields(ax1)
	ax1.plot(silica, total_alkalis, 'o')

	Fontsize and color options can be used to change from the defaults.

	It may be necessary to follow the command with plt.draw() to update
	the plot.

	Le Maitre RW (2002) Igneous rocks : IUGS classification and glossary of
		terms : recommendations of the International Union of Geological
		Sciences Subcommission on the Systematics of igneous rocks, 2nd ed.
		Cambridge University Press, Cambridge
	"""
	from collections import namedtuple
	# Prepare the field information
	FieldLine = namedtuple('FieldLine', 'x1 y1 x2 y2')
	lines = (FieldLine(x1=41, y1=0, x2=41, y2=7),
			 FieldLine(x1=41, y1=7, x2=52.5, y2=14),
			 FieldLine(x1=45, y1=0, x2=45, y2=5),
			 FieldLine(x1=41, y1=3, x2=45, y2=3),
			 FieldLine(x1=45, y1=5, x2=61, y2=13.5),
			 FieldLine(x1=45, y1=5, x2=52, y2=5),
			 FieldLine(x1=52, y1=5, x2=69, y2=8),
			 FieldLine(x1=49.4, y1=7.3, x2=52, y2=5),
			 FieldLine(x1=52, y1=5, x2=52, y2=0),
			 FieldLine(x1=48.4, y1=11.5, x2=53, y2=9.3),
			 FieldLine(x1=53, y1=9.3, x2=57, y2=5.9),
			 FieldLine(x1=57, y1=5.9, x2=57, y2=0),
			 FieldLine(x1=52.5, y1=14, x2=57.6, y2=11.7),
			 FieldLine(x1=57.6, y1=11.7, x2=63, y2=7),
			 FieldLine(x1=63, y1=7, x2=63, y2=0),
			 FieldLine(x1=69, y1=12, x2=69, y2=8),
			 FieldLine(x1=45, y1=9.4, x2=49.4, y2=7.3),
			 FieldLine(x1=69, y1=8, x2=77, y2=0))

	FieldName = namedtuple('FieldName', 'name x y rotation')
	names = (FieldName('Picro\nbasalt', 43, 2, 0),
			 FieldName('Basalt', 48.5, 2, 0),
			 FieldName('Basaltic\nandesite', 54.5, 2, 0),
			 FieldName('Andesite', 60, 2, 0),
			 FieldName('Dacite', 68.5, 2, 0),
			 FieldName('Rhyolite', 76, 9, 0),
			 FieldName('Trachyte\n(Q < 20%)\n\nTrachydacite\n(Q > 20%)',
					   64.5, 11.5, 0),
			 FieldName('Basaltic\ntrachyandesite', 53, 8, -20),
			 FieldName('Trachy-\nbasalt', 49, 6.2, 0),
			 FieldName('Trachyandesite', 57.2, 9, 0),
			 FieldName('Phonotephrite', 49, 9.6, 0),
			 FieldName('Tephriphonolite', 53.0, 11.8, 0),
			 FieldName('Phonolite', 57.5, 13.5, 0),
			 FieldName('Tephrite\n(Ol < 10%)', 45, 8, 0),
			 FieldName('Foidite', 44, 11.5, 0),
			 FieldName('Basanite\n(Ol > 10%)', 43.5, 6.5, 0))

	# Plot the lines and fields
	for line in lines:
		plot_axes.plot([line.x1, line.x2], [line.y1, line.y2],
					   '-', color=color, zorder=0)
	for name in names:
		plot_axes.text(name.x, name.y, name.name, color=color, size=fontsize,
				 horizontalalignment='center', verticalalignment='top',
				 rotation=name.rotation, zorder=0)

def calib_plot(user_data=None, model='all', plot_type='TAS', zoom=None, save_fig=False, **kwargs):
	"""
	Plots user data and calibration set of any or all models on any x-y plot or a total alkalis vs silica (TAS) diagram.
	TAS diagram boundaries provided by tasplot python module, copyright John A Stevenson.

	Parameters
	----------
	user_data: ExcelFile object, pandas DataFrame, pandas Series, or dict
		OPTIONAL. Default value is None, in which case only the model calibration set is plotted.
		User provided sample data describing the oxide composition of one or more samples. Multiple samples
		can be passed as an ExcelFile object or pandas DataFrame. A single sample can be passed as a pandas
		Series.

	model: str or list
		OPTIONAL. Default value is 'all', in which case all model calibration datasets will be plotted.
		'Mixed' can be used to plot all mixed fluid models.
		String of the name of the model calibration dataset to plot (e.g., 'Shishkina'). Multiple models
		can be plotted by passing them as strings within a list (e.g., ['Shishkina', 'Dixon']).

	plot_type: str
		OPTIONAL. Default value is 'TAS', which returns a total alkali vs silica (TAS) diagram. Any two oxides can
		be plotted as an x-y plot by setting plot_type='xy' and specifying x- and y-axis oxides, e.g., x='SiO2', y='Al2O3'

	zoom: str or list
		OPTIONAL. Default value is None in which case axes will be set to the default of 35<x<100 wt% and 0<y<25 wt% for
		TAS type plots and the best values to show the data for xy type plots. Can pass "user_data" to plot the figure
		where the x and y axes are scaled down to zoom in and only show the region surrounding the user_data. A list of
		tuples may be passed to manually specify x and y limits. Pass in data as  [(x_min, x_max), (y_min, y_max)].
		For example, the default limits here would be passed in as [(35,100), (0,25)].

	save_fig: False or str
		OPTIONAL. Default value is False, in which case the figure will not be saved. If a string is passed,
		the figure will be saved with the string as the filename. The string must include the file extension.

	Returns
	-------
	matplotlib object
	"""
	sys.path.insert(0, 'Calibration/')
	import calibrations

	#Get x and y axis limits, if user passed them
	if zoom == None:
		user_xmin = 35
		user_xmax = 100
		user_ymin = 0
		user_ymax = 25
	elif zoom == 'user_data':
		if isinstance(user_data, ExcelFile) or isinstance(user_data, pd.DataFrame):
			print("'user_data' type zoom for more than one sample is not implemented yet.")
			user_xmin = 35
			user_xmax = 100
			user_ymin = 0
			user_ymax = 25
		elif isinstance(user_data, pd.core.series.Series) or isinstance(user_data, dict):
			user_xmin = user_data['SiO2'] - 5
			user_xmax = user_data['SiO2'] + 5
			user_ymin = user_data['Na2O'] + user_data['K2O'] - 2
			if user_ymin <0:
				user_ymin = 0
			user_ymax = user_data['Na2O'] + user_data['K2O'] + 2
	elif isinstance(zoom, list):
		user_xmin, user_xmax = zoom[0]
		user_ymin, user_ymax = zoom[1]

	#Create the figure
	fig, ax1 = plt.subplots(figsize = (17,8))
	font = {'family': 'sans-serif',
				'color':  'black',
				'weight': 'normal',
				'size': 20,
				}

	#TAS figure
	if plot_type == 'TAS':
		ax1.set_xlim([user_xmin, user_xmax]) # adjust x limits here if you want to focus on a specific part of compostional space
		ax1.set_ylim([user_ymin, user_ymax]) # adjust y limits here
		plt.xlabel('SiO$_2$, wt%', fontdict=font, labelpad = 15)
		plt.ylabel('Na$_2$O+K$_2$O, wt%', fontdict=font, labelpad = 15)
		if zoom == None:
			add_LeMaitre_fields(ax1)
	elif plot_type == 'xy':
		if 'x' in kwargs and 'y' in kwargs:
			x = kwargs['x']
			y = kwargs['y']
			if zoom != None:
				ax1.set_xlim([user_xmin, user_xmax])
				ax1.set_ylim([user_ymin, user_ymax])
			plt.xlabel(str(x)+", wt%", fontdict=font, labelpad = 15)
			plt.ylabel(str(y)+", wt%", fontdict=font, labelpad = 15)
		else:
			raise InputError("If plot_type is 'xy', then x and y values must be passed as strings. For example, x='SiO2', y='Al2O3'.")

	#Plot Calibration Data
	if model == 'all':
		model = ['MagmaSat',
				'Shishkina',
			   'Dixon',
			   'IaconoMarziano',
			   'Liu',
			   #'EguchiCarbon',
			   'AllisonCarbon',
			   'MooreWater']
	if model == 'mixed':
		model = ['MagmaSat',
				 'Shishkina',
				 'Dixon',
				 'IaconoMarziano',
				 'Liu']

	if isinstance(model, str):
		model = [model]

	if isinstance(model, list):
		for modelname in model:
			calibdata = calibrations.return_calibration(modelname)
			if isinstance(calibdata, str):
				w.warn(calibdata)
			else:
				if 'CO2' in calibdata.keys():
					if plot_type == 'TAS':
						try:
							plt.scatter(calibdata['CO2']['SiO2'], calibdata['CO2']['Na2O'] + calibdata['CO2']['K2O'],
										marker='o', facecolors=calibdata['facecolor'], edgecolors='k', label=str(modelname)+" CO2")
						except:
							plt.scatter(calibdata['CO2']['SiO2'], calibdata['CO2']['Na2O+K2O'],
									marker='o', facecolors=calibdata['facecolor'], edgecolors='k', label=str(modelname)+" CO2")
					if plot_type == 'xy':
						try:
							plt.scatter(calibdata['CO2'][x], calibdata['CO2'][y],
									marker='o', facecolors=calibdata['facecolor'], edgecolors='k', label=str(modelname)+" CO2")
						except:
							w.warn("The requested oxides were not found in the calibration dataset for " + str(modelname) + ".")

				if 'H2O' in calibdata.keys():
					if plot_type == 'TAS':
						try:
							plt.scatter(calibdata['H2O']['SiO2'], calibdata['H2O']['Na2O'] + calibdata['H2O']['K2O'],
										marker='s', facecolors=calibdata['facecolor'], edgecolors='k', label=str(modelname)+" H2O")
						except:
							plt.scatter(calibdata['H2O']['SiO2'], calibdata['H2O']['Na2O+K2O'],
										marker='s', facecolors=calibdata['facecolor'], edgecolors='k', label=str(modelname)+" H2O")
					if plot_type == 'xy':
						try:
							plt.scatter(calibdata['H2O'][x], calibdata['H2O'][y],
										marker='s', facecolors=calibdata['facecolor'], edgecolors='k', label=str(modelname)+" H2O")
						except:
							w.warn("The requested oxides were not found in the calibration dataset for " + str(modelname) + ".")
				if 'Mixed' in calibdata.keys():
					if plot_type == 'TAS':
						try:
							plt.scatter(calibdata['Mixed']['SiO2'], calibdata['Mixed']['Na2O'] + calibdata['Mixed']['K2O'],
										marker='d', facecolors=calibdata['facecolor'], edgecolors='k', label=str(modelname)+" Mixed")
						except:
							plt.scatter(calibdata['Mixed']['SiO2'], calibdata['Mixed']['Na2O+K2O'],
										marker='d', facecolors=calibdata['facecolor'], edgecolors='k', label=str(modelname)+" Mixed")
					if plot_type == 'xy':
						try:
							plt.scatter(calibdata['Mixed'][x], calibdata['Mixed'][y],
										marker='d', facecolors=calibdata['facecolor'], edgecolors='k', label=str(modelname)+" Mixed")
						except:
							w.warn("The requested oxides were not found in the calibration dataset for " + str(modelname) + ".")
	else:
		raise InputError("model must be of type str or list")

	#Plot user data
	if user_data is None:
		pass
	else:
		if isinstance(user_data, ExcelFile):
			user_data = user_data.data
		if plot_type == 'TAS':
			_sample = user_data.copy()
			try:
				_sample["TotalAlkalis"] = _sample["Na2O"] + _sample["K2O"]
			except:
				InputError("Na2O and K2O data must be in user_data")
			plt.scatter(_sample['SiO2'], _sample['TotalAlkalis'],
						s=150, edgecolors='w', facecolors='red', marker='P',
						label = 'User Data')
		if plot_type == 'xy':
			_sample = user_data.copy()
			plt.scatter(_sample[x], _sample[y],
						s=150, edgecolors='w', facecolors='red', marker='P',
						label = 'User Data')

	plt.legend(bbox_to_anchor=(1.04,1), loc="upper left")
	fig.tight_layout()
	if isinstance(save_fig, str):
		fig.savefig(save_fig)

	return plt.show()

def test_ExcelFile(filename=None):
	"""
	A test routine that takes in an excel file and runs multiple calculations on that file.

	Parameters
	----------
	filename: string
		OPTIONAL. Name with extension of file to be tested. If no file is passed, data will
		be generated automatially.
	"""
	pd.set_option('display.max_rows', None, 'display.max_columns', None)

	print("\n================================\n= MAGMASATPLUS BATCH TESTING ROUTINE =\n================================")

	print("\n This routine will check that key methods run using typical values of variables or given user data. \
The routine does not check that the results are correct (though this may be obvious from the outputs),\
nor does it check every possible iteration of methods and input types. It will check that an update hasn't \
COMPLETELY broken the module.")

	#SET UP THE DATA
	if filename == None:
		fakedata = pd.DataFrame({'Label': ['Samp1', 'Samp2', 'Samp3'],
								'SiO2': [47.95, 69.02, 55.4],
								'TiO2': [1.67, 0.78, 1.01],
								'Al2O3': [17.32, 15.2, 10.1],
								'FeO': [10.24, 4.2, 8.9],
								'Fe2O3': [0.1, 0.2, 0.5],
								'MgO': [5.76, 0.3, 3.0],
								'CaO': [10.93, 12.99, 10.9],
								'Na2O': [3.45, 5.67, 4.22],
								'K2O': [1.99, 3.2, 3.2],
								'P2O5': [0.51, 0.2, 0.5],
								'MnO': [0.1, 0.15, 0.1],
								'CO2': [0.8, 0.2, 0.3],
								'H2O': [4.0, 6.0, 2.0]})
		fakedata = fakedata.set_index('Label')
		myfile = ExcelFile(filename=None, dataframe=fakedata)
	else:
		myfile = ExcelFile(filename)

	test_temperature = 1000
	test_pressure = 2000
	test_X_fluid = 1

	#CALCULTE SATURATION PRESSURE
	print("\n Saturation Pressures:")
	print(" =====================")
	satPs = myfile.calculate_saturation_pressure(temperature=test_temperature)
	print(satPs)

	#CALCULATE DISSOLVED VOLATILES
	print("\n Dissolved Volatile Concentrations:")
	print(" ==================================")
	dissolved = myfile.calculate_dissolved_volatiles(temperature=test_temperature, pressure=test_pressure, X_fluid=test_X_fluid)
	print(dissolved)

	#CALCULATE EQUILIBRIUM FLUID COMPOSITIONS
	print("\n Equilibrium Fluid Compositions:")
	print(" ===============================")
	eqfluid = myfile.calculate_equilibrium_fluid_comp(temperature=test_temperature, pressure=test_pressure)
	print(eqfluid)

	print("\nTesting routine complete.\n")



def test():
	"""
	This is a set of tests for the module, firstly to ensure all of the functions that should run, do run.
	Secondly, the output can be used to check the module reproduces the results in each of the manuscripts
	describing the models.
	"""

	test_sample = {'SiO2':47.95,
				 'TiO2':1.67,
				 'Al2O3':17.32,
				 'FeO':10.24,
				 'Fe2O3':0.1,
				 'MgO':5.76,
				 'CaO':10.93,
				 'Na2O':3.45,
				 'K2O':1.99,
				 'P2O5':0.51,
				 'MnO':0.1,
				 'CO2':0.8,
				 'H2O':4.0}

	test_pressure = 2000.0
	test_temperature = 1473.15
	test_pressure_list = [1000.0,2000.0,5000.0]
	test_isopleth_list = [0.0,0.5,1.0]


	print("\n================================\n= MAGMASATPLUS TESTING ROUTINE =\n================================")

	print("\n This routine will check that key methods run using typical values of variables. \
The routine does not check that the results are correct (though this may be obvious from the outputs),\
nor does it check every possible iteration of methods and input types. It will check that an update hasn't \
COMPLETELY broken the module.")

	for i, model_name in zip(range(len(default_models)),list(default_models.keys())):
		print("\nTesting model {:d} of {:d}: {:s}".format(i+1,len(default_models),model_name))

		### calculate_dissolved_volatiles
		model = default_models[model_name]
		print("Model contains "+" ".join(model.volatile_species))

		print("Testing calculate_dissolved_volatiles method...")

		if len(model.volatile_species) == 1:
			dissolved = model.calculate_dissolved_volatiles(pressure=test_pressure,temperature=test_temperature,
															sample=test_sample)
		else:
			X_fluid = 1.0/len(model.volatile_species)
			X_fluid = tuple(X_fluid for x in range(len(model.volatile_species)))
			print("Setting X_fluid to "+str(X_fluid))
			dissolved = model.calculate_dissolved_volatiles(pressure=test_pressure,temperature=test_temperature,
															sample=test_sample,X_fluid=X_fluid)

		if len(model.volatile_species) == 1:
			print("  {:s} solubility at {:.0f} bars and {:.0f} K is {:.3f} wt%".format(model.volatile_species[0],
																					 test_pressure,test_temperature,
																					 dissolved))
		else:
			for i,volatile in zip(range(len(model.volatile_species)),model.volatile_species):
				print("  {:s} solubility at {:.0f} bars and {:.0f} K is {:.3f} wt%".format(volatile,
																						 test_pressure,test_temperature,
																						 dissolved[i]))

		print("Testing calculate_dissolved_volatiles class interface...")
		if len(model.volatile_species) == 1:
			result = calculate_dissolved_volatiles(sample=test_sample,pressure=test_pressure,
												   temperature=test_temperature,model=model_name)
		else:
			result = calculate_dissolved_volatiles(sample=test_sample,pressure=test_pressure,
												   temperature=test_temperature,model=model_name,
												   X_fluid=X_fluid)

		if len(model.volatile_species) == 1:
			print("  {:s} solubility at {:.0f} bars and {:.0f} K is {:.3f} wt%".format(model.volatile_species[0],
																					 test_pressure,test_temperature,
																					 result.result))
		else:
			for volatile in model.volatile_species:
				print("  {:s} solubility at {:.0f} bars and {:.0f} K is {:.3f} wt%".format(volatile,
																						 test_pressure,test_temperature,
																						 result.result[volatile+'_liq']))



		### calculate_saturation_pressure
		print("Testing calculate_saturation_pressure method...")
		satP = model.calculate_saturation_pressure(sample=test_sample,temperature=test_temperature)
		if len(model.volatile_species) == 1:
			print("  A concentration of {:.2f} wt% {:s} is saturated at {:.0f} bars at {:.0f} K".format(test_sample[model.volatile_species[0]],
																									   model.volatile_species[0],satP,
																									   test_temperature))
		else:
			concstr = ""
			for volatile in model.volatile_species:
				concstr += "{:.2f}".format(test_sample[volatile])
				concstr += " wt% "
				concstr += volatile
				concstr += ", "
			print("  Concentrations of "+concstr[:-1]+" are saturated at {:.0f} bars at {:.0f} K".format(satP,test_temperature))

		print("Testing calculate_saturation_pressure class interface...")
		satP = calculate_saturation_pressure(model=model_name,sample=test_sample,temperature=test_temperature).result
		if len(model.volatile_species) == 1:
			print("  A concentration of {:.2f} wt% {:s} is saturated at {:.0f} bars at {:.0f} K".format(test_sample[model.volatile_species[0]],
																									   model.volatile_species[0],satP,
																									   test_temperature))
		else:
			concstr = ""
			for volatile in model.volatile_species:
				concstr += "{:.2f}".format(test_sample[volatile])
				concstr += " wt% "
				concstr += volatile
				concstr += ", "
			print("  Concentrations of "+concstr[:-1]+" are saturated at {:.0f} bars at {:.0f} K".format(satP,test_temperature))


		### calculate_equilibrium_fluid_comp
		print("Testing calculate_equilibrium_fluid_comp method...")
		fluid = model.calculate_equilibrium_fluid_comp(sample=test_sample,temperature=test_temperature,pressure=test_pressure)

		if len(model.volatile_species) == 1:
			print("  A mole fraction of {:.2f} of {:s} is present in the fluid.".format(fluid,model.volatile_species[0]))
		else:
			fluidstr = ""
			for i,volatile in zip(range(len(model.volatile_species)),model.volatile_species):
				fluidstr += "{:.2f}".format(fluid[model.volatile_species[i]])
				fluidstr += " "
				fluidstr += volatile
				fluidstr += ", "
			print("  Mole fractions of "+fluidstr +"are present in the fluid.")
			if np.sum(list(fluid.values())) != 0.0 and np.sum(list(fluid.values())) != 1.0:
				print("  WARNING: MOLE FRACTIONS DO NOT SUM TO 1.0")

		print("Testing calculate_equilibrium_fluid_comp class interface...")
		fluid = model.calculate_equilibrium_fluid_comp(model=model_name,sample=test_sample,
													   temperature=test_temperature,pressure=test_pressure)
		if len(model.volatile_species) == 1:
			print("  A mole fraction of {:.2f} of {:s} is present in the fluid.".format(fluid,model.volatile_species[0]))
		else:
			fluidstr = ""
			for i,volatile in zip(range(len(model.volatile_species)),model.volatile_species):
				fluidstr += "{:.2f}".format(fluid[model.volatile_species[i]])
				fluidstr += " "
				fluidstr += volatile
				fluidstr += ", "
			print("  Mole fractions of "+fluidstr +"are present in the fluid.")
			if np.sum(list(fluid.values())) != 0.0 and np.sum(list(fluid.values())) != 1.0:
				print("  WARNING: MOLE FRACTIONS DO NOT SUM TO 1.0")

		### calculate_isobars_and_isopleths
		if len(model.volatile_species) > 1:
			print("Testing calculate_isobars_and_isopleths method...")
			isobars, isopleths = model.calculate_isobars_and_isopleths(pressure_list=test_pressure_list,
																	   isopleth_list=test_isopleth_list,
																	   sample=test_sample,
																	   temperature=test_temperature)
			print("Isobars:")
			print(isobars)
			print("\nIsopleths:")
			print(isopleths)

			print("Testing calculate_isobars_and_isopleths class interface...")
			isobars, isopleths = calculate_isobars_and_isopleths(model=model_name,pressure_list=test_pressure_list,
																 isopleth_list=test_isopleth_list,
																 sample=test_sample,
																 temperature=test_temperature).result
			print("Isobars:")
			print(isobars)
			print("\nIsopleths:")
			print(isopleths)

		### calculate_degassing_path
		if len(model.volatile_species) > 1:
			print("Testing calculate_degassing_path method...")
			degassing = model.calculate_degassing_path(sample=test_sample,temperature=test_temperature)
			print("  Degassing path:")
			print(degassing)

			print("Testing calculate_degassing_path class interface...")
			degassing = calculate_degassing_path(model=model_name,sample=test_sample,
														temperature=test_temperature).result
			print("  Degassing path:")
			print(degassing)


	print("\nTesting routine complete.\n")




if __name__ == '__main__':
	test()

"""
						,,,                                     .*****
					   ,***,*                                  ,* *****
					  ,***,,,,*                              ,* ,*******
					 .****,,,.,**                          ,,*.,*********
					 *****,,,,.****                     ,,,,*.,********* *
					 *,**,,,,,,,******                ,,,,,* ,,********,.*                                     .
					.*,*,, //,,.,*******            ,,,,,**.//,*******,,.,                                   //
					,**,, //// ,,******** ,,,,,,,, **,,,*** ((,,,,,,,,,,,*                                  //*####
					.,* , ///((,,,******. ,,,,,,, *********(((,.#### ,,.,*                                   ///#######
					 ,,** .((((,******** ,,,*,,, *********.(((#####/.,,**    .,,                              ///#*#######
					 ,,**(((((.******* ,,**** **.*********.#(####(( .,**.******                                ///##/#######
					  ,,***(.( ****** ,*******#*********** #*###( ,,** ******,                              /////(###########
					,,, ,**,,,,,,. ,,,**,**.### ****.************.  .******,                        .,,,,,****, ##############
					.,,,,,,.,,,,,,,,,* *,#### ************# *************,                    ,,,,********,***/#######(#######
					  ,,,,,, /(   (#,*##### *****###/ #### ************                 ,,,,,,**********,.**.#######(#########(
						,,,,, ((((((# ### * * #..######## *****************.##        ,,,,,****************,####. (## #######
					  ,,,,,,,, *(((( *###.****,#  ##### *********** #### ###         ,,,,,*********************,######*#######
					 ./ ,,,,,,,   ,****########*****,       .**########### ,         ,,,*********************** (##./##.##/#*
					   /////*,  ,,  *****,*,******   ,,.,((( ########### ,****      .,,.*****************,######  //#/////,
						//////  /(,. *****###.****  *((,*&((########## *****       .****,************ ###########(/////*  //
						  /////,*#((/,(((*########* *((#(*##########********** .,   *****,*******    .###########//////,
							 ///// ,* (((((####### %#  ##########,**********,************* **  **** /##### ///////*/ /
								./////((((((################# ###***************************** ,,,,,  ///////*  ////
							  ///// .////(((############ .(((###**#.******,********************     .    //////.
								* //// ////(######## ///(#########*******,****.*****************
								 *//////.////### */////(#########.*****,,****** ****************
							   /////////////////////##/##########* .* ,,*******,,************ **.
								   ,///////////////################(,,,********,*  ************ *,*
										, ////////################,,,,*********  ..,******,** **.***
										, //////##/############# ,,,,,******* *****.,,*.,,, *********,,
									   ,,,//////(#/############  ,,,,,************ ((( ,, .***********
									   ,,///  .  ///###########((,,,,,*********.* //   ,  ,** *  *,
									   ,,,,,,,,,,  //#########,(( ,,,,******** *       ,, ,,**   *,
									   ,,,,,,,,,.   *//###( #,/// ,,,,******* *        .,,, ,,***..
									   .,,,***,,      ///# *,//*,,.,,,****** ,           ,,,,,,.,,,,
									   ,,,*****         *  // ,,,. ,,,***** *              ,,,,,,,,*.
									   , ,*****              ,, ,, *,,******                 ,,,,,,**
									   ,,,.***               ,,,,, ** *.****                  , ,,,,,
									   ,,,****              ,,,,,  ********                    ,,,,,,
									   ,,****               ,,,,,  *******                      ,,,,,
									   ,,***.             .,,,,    ******                       ,,,,*
									   ,****                       *****,                       ,,**
									   .***                       ,*****                        , ,
																  *****
																 ,****


						   ***** *      **        ***** **          *******            *****  *                       ***
						******  *    *****     ******  **** *     *       ***       ******  *                          ***
					   **   *  *       *****  **   *  * ****     *         **      **   *  *                            **
					  *    *  **       * **  *    *  *   **      **        *      *    *  *                             **
						  *  ***      *          *  *             ***                 *  *                              **
						 **   **      *         ** **            ** ***              ** **         ****       ****      **
						 **   **      *         ** **             *** ***            ** **        * ***  *   * ***  *   **
						 **   **     *          ** ******           *** ***        **** **       *   ****   *   ****    **
						 **   **     *          ** *****              *** ***     * *** **      **         **    **     **
						 **   **     *          ** **                   ** ***       ** **      **         **    **     **
						  **  **    *           *  **                    ** **  **   ** **      **         **    **     **
						   ** *     *              *                      * *  ***   *  *       **         **    **     **
							***     *          ****         *   ***        *    ***    *        ***     *  **    **     **
							 *******          *  ***********   *  *********      ******          *******    ***** **    *** *
							   ***           *     ******     *     *****          ***            *****      ***   **    ***
											 *                *
											  **               **

logo courtesy Twai (http://www.twitter.com/_twai)
"""

def WhatDoesTheFoxSay():
	print("https://www.youtube.com/watch?v=jofNR_WkoCE")
