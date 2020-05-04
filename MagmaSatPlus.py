# Python 3.5
# Script written by Kayla Iacovino (kayla.iacovino@nasa.gov)
# VERSION 0.1- MARCH 2020

import pandas as pd
import numpy as np
from thermoengine import equilibrate
import matplotlib.pyplot as plt
from abc import ABC, abstractmethod
from scipy.optimize import root_scalar
from scipy.optimize import root
from scipy.optimize import minimize

print('Write code to deal with volatile undersaturation, i.e. low pressures have an upper bound on volatile solubility, so will cause errors when calculating satP')
print("Build in a function for re-normalising wtpt after reporting H2O")
print("Build in some knowledge of how the parameterisations want to treat Fe, H2O and CO2 in bulk comp.")

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
			 'NiO': 'Ni', 'CoO': 'Co', 'Fe2O3': 'Fe3', 'H2O': 'H', 'CO2': 'C'}
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
msp_fontdict = {'family': 'serif',
				 'color': 'darkblue',
				 'weight': 'normal',
				 'size': 18,}

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
   table = pd.DataFrame([v for v in myDict.values()], columns = ['value'],
						 index = [k for k in myDict.keys()])

   return table

#----------DEFINE SOME BASIC DATA TRANSFORMATION METHODS-----------#

def mol_to_wtpercent(dataframe):
	"""
	Takes in a pandas DataFrame containing multi-sample input and returns a pandas DataFrame object
	with oxide values converted from mole percent to wt percent.

	Parameters
	----------
	dataframe: pandas DataFrame object
		Variable name referring to the pandas DataFrame object that contains user-imported data
	"""
	data = dataframe

	for key, value in oxideMass.items():
		data.loc[:, key] *= value

	data["MPOSum"] = sum([data[oxide] for oxide in oxides])

	for oxide in oxides:
		data.loc[:, oxide] /= data['MPOSum']
		data.loc[:, oxide] *= 100
	del data['MPOSum']

	return data

def wtpercentOxides_to_molCations(oxides):
	"""Takes in a pandas Series containing major element oxides in wt%, and converts it
	to molar proportions of cations (normalised to 1).

	Parameters
	----------
	oxides 		pandas Series or dictionary
		Major element oxides in wt%.

	Returns
	-------
	pandas Series
		Molar proportions of cations, normalised to 1.
	"""
	molCations = {}
	if type(oxides) == dict:
		_oxides = pd.Series(oxides.copy())
	elif type(oxides) != pd.core.series.Series:
		return InputError("The composition input must be a pandas Series or dictionary.")
	else:
		_oxides = oxides.copy()
	for ox in _oxides.index:
		cation = oxides_to_cations[ox]
		molCations[cation] = CationNum[ox]*_oxides[ox]/oxideMass[ox]
	molCations = pd.Series(molCations)

	molCations = molCations/molCations.sum()

	return molCations

def wtpercentOxides_to_molOxides(oxides):
	""" Takes in a pandas Series containing major element oxides in wt%, and converts it
	to molar proportions (normalised to 1).

	Parameters
	----------
	oxides 		pandas Series or dictionary
		Major element oxides in wt%

	Returns
	-------
	pandas Series
		Molar proportions of major element oxides, normalised to 1.
	"""
	molOxides = {}
	if type(oxides) == dict:
		_oxides = pd.Series(oxides.copy())
	elif type(oxides) != pd.core.series.Series:
		return InputError("The composition input must be a pandas Series or dictionary.")
	else:
		_oxides = oxides.copy()
	for ox in _oxides.index:
		molOxides[ox] = _oxides[ox]/oxideMass[ox]
	molOxides = pd.Series(molOxides)

	molOxides = molOxides/molOxides.sum()

	return molOxides

def wtpercentOxides_to_molSingleO(oxides):
	""" Takes in a pandas Series containing major element oxides in wt%, and constructs
	the chemical formula, on a single oxygen basis.

	Parameters
	----------
	oxides 		pandas Series or dictionary
		Major element oxides in wt%

	Returns
	-------
	pandas Series
		The chemical formula of the composition, on a single oxygen basis. Each element is
		a separate entry in the Series.
	"""
	molCations = {}
	if type(oxides) == dict:
		_oxides = pd.Series(oxides.copy())
	elif type(oxides) != pd.core.series.Series:
		return InputError("The composition input must be a pandas Series or dictionary.")
	else:
		_oxides = oxides
	for ox in _oxides.index:
		cation = oxides_to_cations[ox]
		molCations[cation] = CationNum[ox]/OxygenNum[ox]*oxides[ox]/oxideMass[ox]
	molCations = pd.Series(molCations)

	molCations = molCations/molCations.sum()

	return molCations

def wtpercentOxides_to_formulaWeight(sample):
	""" Converts major element oxides in wt% to the formula weight (on a 1 oxygen basis).
	Parameters
	----------
	sample 	pandas Series
		Major element oxides in wt%.

	Returns
	-------
	float
		The formula weight of the composition, on a one oxygen basis.
	"""
	if type(sample) == dict:
		_sample = pd.Series(sample.copy())
	elif type(sample) != pd.core.series.Series:
		return InputError("The composition input must be a pandas Series or dictionary.")
	else:
		_sample = sample.copy()
	cations = wtpercentOxides_to_molSingleO(_sample)
	FW = 15.999
	for cation in list(cations.keys()):
		FW += cations[cation]*CationMass[cations_to_oxides[cation]]
	return FW

#----------DEFINE SOME NORMALIZATION METHODS-----------#	

def normalize(sample):
	"""Normalizes an input composition to 100%. This is the 'standard' normalization routine.

	Parameters
	----------
	sample:	pandas Series, dictionary, pandas DataFrame, or ExcelFile object
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
	sample:	pandas Series, dictionary, pandas DataFrame, or ExcelFile object
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
		normalized = pd.Series({})
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
		return InputError("The composition input must be a pandas Series or dictionary for single sample \
							or a pandas DataFrame or ExcelFile object for multi-sample.")


def normalize_AdditionalVolatiles(sample):
	"""Normalises major element oxide wt% to 100%, assuming it is volatile-free. If
	H2O or CO2 are passed to the function, their un-normalized values will be retained
	in addition to the normalized non-volatile oxides, summing to >100%.

	Parameters
	----------
	sample:	pandas Series, dictionary, pandas DataFrame, or ExcelFile object
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
		normalized = pd.Series({})
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
		return InputError("The composition input must be a pandas Series or dictionary for single sample \
							or a pandas DataFrame or ExcelFile object for multi-sample.")


#------------DEFINE MAJOR CLASSES-------------------#
class ExcelFile(object):
	"""An excel file with sample names and oxide compositions

	Attributes
	----------
		input_type: str
			String defining whether the oxide composition is given in wt percent ("wtpercent", which is the default),
			mole percent ("molpercent"), or mole fraction ("molfrac).
	"""

	def __init__(self, filename, input_type='wtpercent'):
		"""Return an ExcelFile object whoes parameters are defined here."""
		self.input_type = input_type

		data = pd.read_excel(filename)

		try:
			data = data.set_index('Label')
		except:
			raise InputError(
				"Imported file must contain a column of sample names with the column name \'Label\'") #TODO test

		for oxide in oxides:
			if oxide in data.columns:
				pass
			else:
				data[oxide] = 0.0

		# TODO test all input types produce correct values
		if input_type == "wtpercent":
			pass

		if input_type == "molpercent":
			data = mol_to_wtpercent(data)

		if input_type == "molfrac":
			data = mol_to_wtpercent(data)

		self.data = data

	def preprocess_sample(self,sample):
		#--------------MELTS preamble---------------#
		# instantiate thermoengine equilibrate MELTS instance
		melts = equilibrate.MELTSmodel('1.2.0')

		# Suppress phases not required in the melts simulation
		self.oxides = melts.get_oxide_names()
		self.phases = melts.get_phase_names()

		for phase in self.phases:
			melts.set_phase_inclusion_status({phase: False})
		melts.set_phase_inclusion_status({'Fluid': True, 'Liquid': True})
		self.melts = melts
		#-------------------------------------------#
		oxides = self.oxides
		for oxide in oxides:
			if oxide in self.data.columns:
				pass
			else:
				self.data[oxide] = 0.0

		return sample

	def get_sample_oxide_comp(self, sample, norm='none'):
		"""
		Returns oxide composition of a single sample from a user-imported excel file as a dictionary

		Parameters
		----------
		sample: string
			Name of the desired sample

		norm_style: string
			OPTIONAL. Default value is 'standard'. This specifies the style of normalization applied to the sample.

			'standard' normalizes the entire input composition (including any volatiles) to 100%.

			'fixedvolatiles' normalizes oxides to 100%, including volatiles. The volatile
			wt% will remain fixed, whilst the other major element oxides are reduced proportionally
			so that the total is 100 wt%.

			'anhydrous' normalizes oxides to 100%, assuming it is volatile-free. If
			H2O or CO2 are passed to the function, their un-normalized values will be retained
			in addition to the normalized non-volatile oxides, summing to >100%.

			'none' returns the value-for-value un-normalized composition.

		Returns
		-------
		dictionary
			Composition of the sample as oxides
		"""
		data = self.data
		my_sample = pd.DataFrame(data.loc[sample])
		sample_dict = (my_sample.to_dict()[sample])
		sample_oxides = {}
		for item, value in sample_dict.items():
			if item in oxides:
				sample_oxides.update({item: value})

		if norm == 'standard':
			return normalize(sample_oxides)
		if norm == 'fixedvolatiles':
			return normalize_FixedVolatiles(sample_oxides).to_dict()
		if norm == 'anhydrous':
			return normalize_AdditionalVolatiles(sample_oxides).to_dict()
		if norm == 'none':
			return sample_oxides

	def save_excel_file(self):
		#TODO write this code!
		return("You haven't written this code yet!")

	def calculate_saturation_pressure(self, temperature, print_status=False): #TODO fix weird printing
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
			OPTIONAL: Default is False. If set to True, progress of the calculations will be printed to the terminal.

		Returns
		-------
		pandas DataFrame object
			Values returned are saturation pressure in bars, the mass of fluid present, and the composition of the
			fluid present. 
		"""

		oxides = self.oxides
		data = preprocess_sample(self.data)

		if isinstance(temperature, str):
			file_has_temp = True
		elif isinstance(temperature, float) or isinstance(temperature, int):
			file_has_temp = False
		else:
			raise InputError("temperature must be type str or float or int")


		# Do the melts equilibrations
		bulk_comp = {}
		startingP = []
		startingP_ref = []
		satP = []
		flmass = []
		flH2O = []
		flCO2 = []
		flsystem_wtper = []
		iterno = 0
		for index, row in data.iterrows():
			if iterno == 0:
				bulk_comp_orig = {oxide:  row[oxide] for oxide in oxides}

			bulk_comp = {oxide:  row[oxide] for oxide in oxides}
			feasible = melts.set_bulk_composition(bulk_comp)

			if file_has_temp == True:
				temperature = row[temperature]

			fluid_mass = 0.0
			pressureMPa = 2000.0
			while fluid_mass <= 0.0:
				pressureMPa -= 100.0

				output = melts.equilibrate_tp(temperature, pressureMPa, initialize=True)
				(status, temperature, pressureMPa, xmlout) = output[0]
				fluid_mass = melts.get_mass_of_phase(xmlout, phase_name='Fluid')

				if pressureMPa <= 0:
					break

			startingP.append(pressureMPa+100.0)
			iterno += 1

		data["StartingP"] = startingP

		for index, row in data.iterrows():
			bulk_comp = {oxide:  row[oxide] for oxide in oxides}
			feasible = melts.set_bulk_composition(bulk_comp)

			if file_has_temp == True:
				temperature = row[temperature]

			fluid_mass = 0.0
			pressureMPa = row["StartingP"]
			while fluid_mass <= 0.0:
				pressureMPa -= 10.0

				output = melts.equilibrate_tp(temperature, pressureMPa, initialize=True)
				(status, temperature, pressureMPa, xmlout) = output[0]
				fluid_mass = melts.get_mass_of_phase(xmlout, phase_name='Fluid')

				if pressureMPa <= 0:
					break

			startingP_ref.append(pressureMPa+10.0)

		data["StartingP_ref"] = startingP_ref

		for index, row in data.iterrows():
			bulk_comp = {oxide:  row[oxide] for oxide in oxides}
			feasible = melts.set_bulk_composition(bulk_comp)

			if file_has_temp == True:
				temperature = row[temperature]

			fluid_mass = 0.0
			pressureMPa = row["StartingP_ref"]
			while fluid_mass <= 0.0:
				pressureMPa -= 1.0

				output = melts.equilibrate_tp(temperature, pressureMPa, initialize=True)
				(status, temperature, pressureMPa, xmlout) = output[0]
				fluid_mass = melts.get_mass_of_phase(xmlout, phase_name='Fluid')

				if pressureMPa <= 0:
					break

			satP.append(pressureMPa * 10)
			flmass.append(fluid_mass)
			flsystem_wtper.append(100 * fluid_mass / (fluid_mass +
								  melts.get_mass_of_phase(xmlout, phase_name='Liquid')))

			flcomp = melts.get_composition_of_phase(xmlout, phase_name='Fluid')
			flH2O.append(flcomp["H2O"])
			flCO2.append(flcomp["CO2"])

			if print_status == True:
				print(index)
				print("Pressure (bars) = " + str(pressureMPa * 10))
				print("Fluid mass = " + str(fluid_mass))
				print("\n")

		data["SaturationPressure_bars"] = satP
		data["FluidMassAtSaturation_grams"] = flmass
		data["H2Ofluid_wtper"] = flH2O
		data["CO2fluid_wtper"] = flCO2
		data["FluidSystem_wtper"] = flsystem_wtper
		del data["StartingP"]
		del data["StartingP_ref"]

		feasible = melts.set_bulk_composition(bulk_comp_orig) #this needs to be reset always!

		return data

	def calculate_equilibrium_fluid_comp(self, temperature, pressure): 
		"""
		Returns H2O and CO2 concentrations in wt% in a fluid in equilibrium with the given sample(s) at the given P/T condition.

		Parameters
		----------
		sample: ExcelFile object
			Compositional information on samples in oxides. 

		temperature: float, int, or str
			Temperature, in degrees C. Can be passed as float, in which case the
			passed value is used as the temperature for all samples. Alternatively, temperature information for each individual
			sample may already be present in the passed ExcelFile object. If so, pass the str value corresponding to the column
			title in the passed ExcelFile object.

		presure: float, int, or str
			Pressure, in bars. Can be passed as float or int, in which case the
			passed value is used as the pressure for all samples. Alternatively, pressure information for each individual
			sample may already be present in the passed ExcelFile object. If so, pass the str value corresponding to the column
			title in the passed ExcelFile object.

		Returns
		-------
		pandas DataFrame 
			Original data passed plus newly calculated values are returned.
		"""
		data = self.preprocess_sample(self.data)
		oxides = self.oxides
		melts = self.melts

		if isinstance(temperature, str):
			file_has_temp = True
		elif isinstance(temperature, float) or isinstance(temperature, int):
			file_has_temp = False
		else:
			raise InputError("temp must be type str or float or int")

		if isinstance(pressure, str):
			file_has_press = True
		elif isinstance(pressure, float) or isinstance(pressure, int):
			file_has_press = False
		else:
			raise InputError("pressure must be type str or float or int")

		fluid_comp_H2O = []
		fluid_comp_CO2 = []
		for index, row in data.iterrows():
			bulk_comp = {oxide:  row[oxide] for oxide in oxides}
			feasible = melts.set_bulk_composition(bulk_comp)

			if file_has_temp == True:
				temperature = row[temperature]
			if file_has_press == True:
				pressure = row[pressure]
			
			pressureMPa = pressure / 10.0

			output = melts.equilibrate_tp(temperature, pressureMPa, initialize=True)
			(status, temperature, pressureMPa, xmlout) = output[0]
			fluid_mass = melts.get_mass_of_phase(xmlout, phase_name='Fluid')

			if fluid_mass > 0.0:
				fluid_comp = melts.get_composition_of_phase(xmlout, phase_name='Fluid')
				fluid_comp_H2O.append(fluid_comp['H2O'])
				fluid_comp_CO2.append(fluid_comp['CO2'])
			else:
				fluid_comp_H2O.append(0)
				fluid_comp_CO2.append(0)

		data["H2Ofluid_wtper"] = fluid_comp_H2O
		data["CO2fluid_wtper"] = fluid_comp_CO2

		return data

class Model(object):
	"""The model object implements a volatile solubility model. It is composed
	of the methods needed to evaluate calculate_dissolved_volatiles,
	calculate_equilibrium_fluid_comp, and calculate_saturation_pressure. The
	fugacity and activity models for the volatiles species must be specified,
	defaulting to ideal.
	"""

	def __init__(self):
		self.set_volatile_species(None)
		self.set_fugacity_model(fugacity_idealgas())
		self.set_activity_model(activity_idealsolution())

	def set_volatile_species(self,volatile_species):
		self.volatile_species = volatile_species

	def set_fugacity_model(self,fugacity_model):
		self.fugacity_model = fugacity_model

	def set_activity_model(self,activity_model):
		self.activity_model = activity_model

	@abstractmethod
	def calculate_dissolved_volatiles(self,**kwargs):
		"""
		WRITE STUFF
		"""

	@abstractmethod
	def calculate_equilibrium_fluid_comp(self,**kwargs):
		"""
		WRITE STUFF
		"""

	@abstractmethod
	def calculate_saturation_pressure(self,**kwargs):
		"""
		WRITE STUFF
		"""

	@abstractmethod
	def preprocess_sample(self,**kwargs):
		"""
		WRITE stuff
		"""

	@abstractmethod
	def check_calibration_range(self,**kwargs):
		"""
		WRITE stuff
		"""


class FugacityModel(object):
	""" The fugacity model object is for implementations of fugacity models
	for individual volatile species, though it may depend on the mole
	fraction of other volatile species. It contains all the methods required
	to calculate the fugacity at a given pressure and mole fraction.
	"""

	@abstractmethod
	def fugacity(self,pressure,**kwargs):
		"""
		"""

	@abstractmethod
	def check_calibration_range(self,parameters):
		"""
		"""



class activity_model(object):
	""" The activity model object is for implementing activity models
	for volatile species in melts. It contains all the methods required to
	evaluate the activity.
	"""

	@abstractmethod
	def activity(self,X,**kwargs):
		"""
		"""

	@abstractmethod
	def check_calibration_range(self,parameters):
		"""
		"""



class Calculate(object):
	""" The Calculate object is a template for implementing user-friendly methods for
	running calculations using the volatile solubility models. All Calculate methods
	have a common workflow- sample is read in, preprocessed, the calculation is performed,
	the calibration range is checked, and the results stored.
	"""
	def __init__(self,sample,model='MagmaSat',**kwargs):
		if type(model) == str:
			self.model = default_models[model]
		else:
			self.model = model

		self.sample = sample.copy()
		self.sample = self.model.preprocess_sample(self.sample)

		self.result = self.calculate(sample=self.sample,**kwargs)
		self.calib_check = self.check_calibration_range(sample=self.sample,**kwargs)


	@abstractmethod
	def calculate(self):
		""" """

	@abstractmethod
	def check_calibration_range(self):
		""" """

#-------------FUGACITY MODELS--------------------------------#

class fugacity_idealgas(FugacityModel):
	""" An instance of FugacityModel for an ideal gas.
	"""

	def fugacity(self,pressure,X_fluid=1.0,**kwargs):
		""" Returns the fugacity of an ideal gas, i.e., the partial pressure.

		Parameters
		----------
		pressure 	float
			Total pressure of the system, in bars.
		X_fluid 	float
			The mole fraction of the species in the vapour phase.

		Returns
		-------
		float
			Fugacity (partial pressure) in bars
		"""
		return pressure*X_fluid

	def check_calibration_range(self,parameters):
		""" Checks that the parameters are within the calbrated range of
		the model, for use with the Calculate methods. Since this is a
		statement of ideality, everything is within its calibrated range.

		Parameters
		----------
		parameters 	dictionary
			Parameter names (keys) and their values, for checking.

		Returns
		-------
		dictionary
			Dictionary with parameter names as keys. The values are booleans
			representing whether the parameter is in the calibrated range, or not.
			Always True for this instance.

		"""
		results = {}
		for param in list(parameters.keys()):
			results[param] = True
		return results

class fugacity_KJ81_co2(FugacityModel):
	""" Implementation of the Kerrick and Jacobs (1981) EOS for mixed fluids. This class
	will return the properties of the CO2 component of the mixed fluid.
	"""

	def fugacity(self,pressure,temperature,X_fluid,**kwargs):
		""" Calculates the fugacity of CO2 in a mixed CO2-H2O fluid. Above 1050C,
		it assumes H2O and CO2 do not interact, as the equations are not defined
		beyond this point.

		Parameters
		----------
		pressure 	float
			Total pressure of the system in bars.
		temperature 	float
			Temperature in K
		X_fluid 	float
			Mole fraction of CO2 in the fluid.

		Returns
		-------
		float
			fugacity of CO2 in bars
		"""
		if X_fluid == 0:
			return 0
		elif temperature >= 1050.0 + 273.15:
			return pressure*np.exp(self.lnPhi_mix(pressure,temperature,1.0))*X_fluid
		else:
			return pressure*np.exp(self.lnPhi_mix(pressure,temperature,X_fluid))*X_fluid


	def volume(self,P,T,X_fluid):
		""" Calculates the volume of the mixed fluid, by solving Eq (28) of Kerrick and
		Jacobs (1981) using scipy.root_scalar.

		Parameters
		----------
		P 	float
			Total pressure of the system, in bars.
		T 	float
			Temperature in K
		X_fluid 	float
			Mole fraction of CO2 in the fluid

		Returns
		-------
		float
			Volume of the mixed fluid.
		"""
		if X_fluid != 1.0:
			# x0 = self.volume(P,T,1.0)*X_fluid + self.volume_h(P,T)*(1-X_fluid)
			# print(x0)
			if P >= 20000 and T<800:
				x0 = (X_fluid*25+(1-X_fluid)*15)
			else:
				x0 = (X_fluid*35+(1-X_fluid)*15)

		else:
			if P >= 20000 and T<800:
				x0 = 25
			else:
				x0=35
		return root_scalar(self.root_volume,x0=x0,x1=x0*0.9,args=(P,T,X_fluid)).root


	def root_volume(self,v,P,T,X_fluid):
		""" Returns the difference between the lhs and rhs of Eq (28) of Kerrick and Jacobs (1981).
		For use with a root finder to obtain the volume of the mixed fluid.

		Parameters
		----------
		v 	float
			Guess for the volume
		P 	float
			Total system pressure in bars.
		T 	float
			Temperature in K
		X_fluid 	float
			Mole fraction of CO2 in the fluid.

		Returns
		-------
		float
			Difference between lhs and rhs of Eq (28) of Kerrick and Jacobs (1981), in bars.
		"""
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
		P 	float
			Total pressure in bars.
		T 	float
			Temperature in K.

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
		v 	float
			Guess for the volume
		P 	float
			Total pressure in bars.
		T 	float
			Temperature in K.

		Returns
		-------
		float
			The difference between the lhs and rhs of Eq (14) of Kerrick and Jacobs (1981),
			in bars.
		"""
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
		P 	float
			Total pressure in bars.
		T 	float
			Temperature in K
		X_fluid 	float
			The mole fraction of CO2 in the fluid.

		Returns
		-------
		float
			The natural log of the fugacity coefficient for CO2 in a mixed fluid.
		"""
		v = self.volume(P,T,X_fluid)

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

	def check_calibration_range(self,parameters,**kwargs):
		"""Checks whether supplied parameters and calculated results are within the calibration range
		of the model. Designed for use with the Calculate methods.

		Parameters supported currently are pressure and temperature.

		Parameters
		----------
		parameters 		dictionary
			Parameters to check calibration range for, the parameter name should be given as the key, and
			its value as the value.

		Returns
		-------
		dictionary
			Dictionary with parameter names as keys. The values are booleans
			representing whether the parameter is in the calibrated range, or not.
		"""
		results = {}
		if 'pressure' in parameters.keys():
			results['pressure'] = (parameters['pressure']<20000)
		if 'temperature' in parameters.keys():
			results['temperature'] = (parameters['temperature']<1323)
		return results


class fugacity_KJ81_h2o(FugacityModel):
	"""Implementation of the Kerrick and Jacobs (1981) EOS for mixed fluids. This class
	will return the properties of the H2O component of the mixed fluid.
	"""

	def fugacity(self,pressure,temperature,X_fluid,**kwargs):
		""" Calculates the fugacity of H2O in a mixed CO2-H2O fluid. Above 1050C,
		it assumes H2O and CO2 do not interact, as the equations are not defined
		beyond this point.

		Parameters
		----------
		pressure 	float
			Total pressure of the system in bars.
		temperature 	float
			Temperature in K
		X_fluid 	float
			Mole fraction of H2O in the fluid.

		Returns
		-------
		float
			fugacity of H2O in bars
		"""
		if X_fluid == 0:
			return 0
		elif temperature >= 1050+273.15:
			return pressure*np.exp(self.lnPhi_mix(pressure,temperature,1.0))*X_fluid
		else:
			return pressure*np.exp(self.lnPhi_mix(pressure,temperature,X_fluid))*X_fluid


	def volume(self,P,T,X_fluid):
		""" Calculates the volume of the mixed fluid, by solving Eq (28) of Kerrick and
		Jacobs (1981) using scipy.root_scalar.

		Parameters
		----------
		P 	float
			Total pressure of the system, in bars.
		T 	float
			Temperature in K
		X_fluid 	float
			Mole fraction of H2O in the fluid

		Returns
		-------
		float
			Volume of the mixed fluid.
		"""
		if X_fluid != 1.0:
			# x0 = self.volume(P,T,1.0)*X_fluid + self.volume_h(P,T)*(1-X_fluid)
			# print(x0)
			if P >= 20000 and T<800:
				x0 = ((1-X_fluid)*25+X_fluid*15)
			else:
				x0 = ((1-X_fluid)*35+X_fluid*15)

		else:
			if P >= 20000 and T<800:
				x0 = 10
			else:
				x0=15
		return root_scalar(self.root_volume,x0=x0,x1=x0*0.9,args=(P,T,X_fluid)).root


	def root_volume(self,v,P,T,X_fluid):
		""" Returns the difference between the lhs and rhs of Eq (28) of Kerrick and Jacobs (1981).
		For use with a root finder to obtain the volume of the mixed fluid.

		Parameters
		----------
		v 	float
			Guess for the volume
		P 	float
			Total system pressure in bars.
		T 	float
			Temperature in K
		X_fluid 	float
			Mole fraction of H2O in the fluid.

		Returns
		-------
		float
			Difference between lhs and rhs of Eq (28) of Kerrick and Jacobs (1981), in bars.
		"""
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
		P 	float
			Total pressure in bars.
		T 	float
			Temperature in K.

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
		v 	float
			Guess for the volume
		P 	float
			Total pressure in bars.
		T 	float
			Temperature in K.

		Returns
		-------
		float
			The difference between the lhs and rhs of Eq (14) of Kerrick and Jacobs (1981),
			in bars.
		"""
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
		P 	float
			Total pressure in bars.
		T 	float
			Temperature in K
		X_fluid 	float
			The mole fraction of H2O in the fluid.

		Returns
		-------
		float
			The natural log of the fugacity coefficient for H2O in a mixed fluid.
		"""
		v = self.volume(P,T,X_fluid)

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

	def check_calibration_range(self,parameters,**kwargs):
		"""Checks whether supplied parameters and calculated results are within the calibration range
		of the model. Designed for use with the Calculate methods.

		Parameters supported currently are pressure and temperature.

		Parameters
		----------
		parameters 		dictionary
			Parameters to check calibration range for, the parameter name should be given as the key, and
			its value as the value.

		Returns
		-------
		dictionary
			Dictionary with parameter names as keys. The values are booleans
			representing whether the parameter is in the calibrated range, or not.
		"""
		results = {}
		if 'pressure' in parameters.keys():
			results['pressure'] = (parameters['pressure']<20000)
		if 'temperature' in parameters.keys():
			results['temperature'] = (parameters['temperature']<1323)
		return results

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
		X 	float
			The mole fraction of the species in the solution.

		Returns
		-------
		float
			The activity of the species in the solution, i.e., the mole fraction.
		"""
		return X

	def check_calibration_range(self,parameters):
		""" Checks that the parameters are within the calbrated range of
		the model, for use with the Calculate methods. Since this is a
		statement of ideality, everything is within its calibrated range.

		Parameters
		----------
		parameters 	dictionary
			Parameter names (keys) and their values, for checking.

		Returns
		-------
		dictionary
			Dictionary with parameter names as keys. The values are booleans
			representing whether the parameter is in the calibrated range, or not.
			Always True for this instance.

		"""
		results = {}
		for param in list(parameters.keys()):
			results[param] = True
		return results



#------------PURE FLUID MODELS-------------------------------#

class ShishkinaCarbon(Model):
	""" Implementation of the Shishkina et al. (2014) carbon solubility model, as a Model class.
	"""
	def __init__(self):
		self.set_volatile_species(['CO2'])
		self.set_fugacity_model(fugacity_idealgas())
		self.set_activity_model(activity_idealsolution())

	def preprocess_sample(self,sample):
		""" Returns sample, unmodified. The Pi* compositional parameter is a ratio of cations,
		therefore the value is not affected by the normalization of the sample. Shishkina et al.
		imply the accuracy of the calculations are little affected whether Fe(tot) or Fe2+ is
		used.

		Parameters
		----------
		sample:		 pandas Series
			The major element oxides in wt%.

		Returns
		-------
		pandas Series
			The major element oxides in wt%.

		"""
		return sample

	def PiStar(self,sample):
		"""Shishkina et al. (2014) Eq (11)

		Calculates the Pi* parameter for use in calculating CO2 solubility.

		Parameters
		----------
		sample:		pandas Series
			Major element oxides in wt%.

		Returns
		-------
		float
			The value of the Pi* compositional parameter.
		"""
		_mols = wtpercentOxides_to_molCations(sample)
		_pi = (_mols['Ca'] + 0.8*_mols['K'] + 0.7*_mols['Na'] + 0.4*_mols['Mg'] + 0.4*_mols['Fe'])/\
				(_mols['Si']+_mols['Al'])
		return _pi

	def calculate_dissolved_volatiles(self,pressure,sample,X_fluid=1,**kwargs):
		""" Calculates the dissolved CO2 concentration in wt%, using equation (13) of Shishkina et al. (2014).

		Parameters
		----------
		pressure:	float
			(Total) pressure in bars.
		sample:		pandas Series
			Major element concentrations in wt%. Normalization does not matter.
		X_fluid:	float
			The mol-fraction of the fluid that is CO2. Default is 1, i.e. a pure CO2 fluid.

		Returns
		-------
		float
			The dissolved CO2 concentration in wt%.
		"""

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
		pressure 	float
			The total pressure of the system in bars.
		sample 		pandas Series
			Major element oxides in wt%

		Returns
		-------
		float
			1.0 if CO2-fluid saturated, 0.0 otherwise.
		"""
		if self.calculate_saturation_pressure(self,sample=sample,**kwargs) > pressure:
			return 0.0
		else:
			return 1.0

	def calculate_saturation_pressure(self,sample,**kwargs):
		""" Calculates the pressure at which a pure CO2 fluid is saturated, for the given
		sample composition and CO2 concentration. Calls the scipy.root_scalar routine, which makes
		repeated calls to the calculate_dissolved_volatiles method.

		Parameters
		----------
		sample 		pandas Series
			Major elements in wt%, including CO2 (also in wt%).

		Returns
		-------
		float
			Saturation pressure in bar
		"""
		satP = root_scalar(self.root_saturation_pressure,x0=1000.0,x1=2000.0,args=(sample,kwargs)).root
		return satP

	def root_saturation_pressure(self,pressure,sample,kwargs):
		""" Function called by scipy.root_scalar when finding the saturation pressure using
		calculate_saturation_pressure.

		Parameters
		----------
		pressure 	float
			Pressure guess in bars
		sample 		pandas Series
			Major element oxides in wt%, including CO2 (also in wt%).
		kwargs 		dictionary
			Additional keyword arguments supplied to calculate_saturation_pressure. Might be required for
			the fugacity or activity models.

		Returns
		-------
		float
			The differece between the dissolved CO2 at the pressure guessed, and the CO2 concentration
			passed in the sample variable.
		"""
		return self.calculate_dissolved_volatiles(pressure=pressure,sample=sample,**kwargs)-sample['CO2']

	def check_calibration_range(self,parameters,**kwargs):
		""" Checks whether supplied parameters and calculated results are within the calibration range
		of the model. Designed for use with the Calculate methods. Calls the check_calibration_range
		functions for the fugacity and activity models.

		Parameters supported currently are pressure and temperature.

		Parameters
		----------
		parameters 		dictionary
			Parameters to check calibration range for, the parameter name should be given as the key, and
			its value as the value.

		Returns
		-------
		dictionary
			Dictionary with parameter names as keys. The values are dictionarys, which have the model component
			as the keys, and bool values, indicating whether the parameter is within the calibration range.
		"""
		results = {}
		if 'pressure' in parameters.keys():
			pressure_results = {'Shishkina Model': (parameters['pressure']>500.0) and (parameters['pressure']<5000.0),
								'Fugacity Model': self.fugacity_model.check_calibration_range(parameters)['pressure'],
								'Activity Model': self.activity_model.check_calibration_range(parameters)['pressure']}
			results['pressure'] = pressure_results

		if 'temperature' in parameters.keys():
			temperature_results = {'Shishkina Model': (parameters['temperature']>1200+273.15) and (parameters['temperature']<1250+273.15),
									'Fugacity Model': self.fugacity_model.check_calibration_range(parameters)['temperature'],
									'Activity Model': self.activity_model.check_calibration_range(parameters)['temperature']}
			results['temperature'] = temperature_results

		return results

class ShishkinaWater(Model):
	""" Implementation of the Shishkina et al. (2014) H2O solubility model as a Model class.
	"""
	def __init__(self):
		self.set_volatile_species(['H2O'])
		self.set_fugacity_model(fugacity_idealgas())
		self.set_activity_model(activity_idealsolution())

	def preprocess_sample(self,sample):
		""" Returns sample, renormlized so that the major element oxides (excluding volatiles) sum to 100%.
		Normalization must be done this way as the compositional dependence of the solubility takes the
		mole fractions of Na2O and K2O as inputs, presumably assuming no volatiles in the bulk composition.
		Volatile concentrations are left unchanged.

		Parameters
		----------
		sample:		 pandas Series
			The major element oxides in wt%.

		Returns
		-------
		pandas Series
			The major element oxides in wt%.

		"""
		return normalize_AdditionalVolatiles(sample)

	def calculate_dissolved_volatiles(self,pressure,sample,X_fluid=1.0,**kwargs):
		"""Calculates the dissolved H2O concentration using Eqn (9) of Shishkina et al. (2014).

		Parameters
		----------
		pressure 	float
			Total pressure in bars
		sample 		pandas Series
			Major element oxides in wt%. Normalized to zero-volatiles so that the total-alkalis
			mol fraction can be determined accurately.
		X_fluid 	float
			The mol fraction of H2O in the fluid

		Returns
		-------
		float
			The H2O concentration in wt%
		"""
		_mols = wtpercentOxides_to_molCations(sample)
		total_alkalis = _mols['Na'] + _mols['K']

		fugacity = self.fugacity_model.fugacity(pressure,X_fluid=X_fluid,**kwargs)

		a = 3.36e-7 * (fugacity/10)**3 - 2.33e-4*(fugacity/10)**2 + 0.0711*(fugacity/10) - 1.1309
		b = -1.2e-5*(fugacity/10)**2 + 0.0196*(fugacity/10)+1.1297

		return a*total_alkalis + b


	def calculate_equilibrium_fluid_comp(self,pressure,sample,**kwargs):
		""" Returns 1.0 if a pure H2O fluid is saturated.
		Returns 0.0 if a pure H2O fluid is undersaturated.

		Parameters
		----------
		pressure 	float
			The total pressure of the system in bars.
		sample 		pandas Series
			Major element oxides in wt%, normalized on the basis of
			no volatiles.

		Returns
		-------
		float
			1.0 if H2O-fluid saturated, 0.0 otherwise.
		"""
		if self.calculate_saturation_pressure(self,sample=sample,**kwargs) > pressure:
			return 0.0
		else:
			return 1.0

	def calculate_saturation_pressure(self,sample,**kwargs):
		""" Calculates the pressure at which a pure H2O fluid is saturated, for the given
		sample composition and H2O concentration. Calls the scipy.root_scalar routine, which makes
		repeated calls to the calculate_dissolved_volatiles method.

		Parameters
		----------
		sample 		pandas Series
			Major elements in wt% (normalized to 100%), including H2O (also in wt%, not included
			in normalization).

		Returns
		-------
		float
			Saturation pressure in bar
		"""
		satP = root_scalar(self.root_saturation_pressure,x0=1000.0,x1=2000.0,args=(sample,kwargs)).root
		return satP

	def root_saturation_pressure(self,pressure,sample,kwargs):
		""" Function called by scipy.root_scalar when finding the saturation pressure using
		calculate_saturation_pressure.

		Parameters
		----------
		pressure 	float
			Pressure guess in bars
		sample 		pandas Series
			Major elements in wt% (normalized to 100%), including H2O (also in wt%, not included
			in normalization).
		kwargs 		dictionary
			Additional keyword arguments supplied to calculate_saturation_pressure. Might be required for
			the fugacity or activity models.

		Returns
		-------
		float
			The differece between the dissolved H2O at the pressure guessed, and the H2O concentration
			passed in the sample variable.
		"""
		return self.calculate_dissolved_volatiles(pressure=pressure,sample=sample,**kwargs)-sample['H2O']

	def check_calibration_range(self,parameters,**kwargs):
		""" Checks whether supplied parameters and calculated results are within the calibration range
		of the model. Designed for use with the Calculate methods. Calls the check_calibration_range
		functions for the fugacity and activity models.

		Parameters supported currently are pressure and temperature.

		Parameters
		----------
		parameters 		dictionary
			Parameters to check calibration range for, the parameter name should be given as the key, and
			its value as the value.

		Returns
		-------
		dictionary
			Dictionary with parameter names as keys. The values are dictionarys, which have the model component
			as the keys, and bool values, indicating whether the parameter is within the calibration range.
		"""
		results = {}
		if 'pressure' in parameters.keys():
			pressure_results = {'Shishkina Model': (parameters['pressure']>500.0) and (parameters['pressure']<5000.0),
								'Fugacity Model': self.fugacity_model.check_calibration_range(parameters)['pressure'],
								'Activity Model': self.activity_model.check_calibration_range(parameters)['pressure']}
			results['pressure'] = pressure_results

		if 'temperature' in parameters.keys():
			temperature_results = {'Shishkina Model': (parameters['temperature']>1200+273.15) and (parameters['temperature']<1250+273.15),
									'Fugacity Model': self.fugacity_model.check_calibration_range(parameters)['temperature'],
									'Activity Model': self.activity_model.check_calibration_range(parameters)['temperature']}
			results['temperature'] = temperature_results

		return results


class DixonCarbon(Model):
	"""
	"""

	def __init__(self):
		self.set_volatile_species('CO2')
		self.set_fugacity_model(fugacity_KJ81_co2())
		self.set_activity_model(activity_idealsolution())

	def preprocess_sample(self,sample):
		return sample

	def calculate_dissolved_volatiles(self,pressure,sample,X_fluid=1.0,**kwargs):
		"""Dixon (1997) Eqn (3).

		Calculates the concentration (in ppm) of CO2 in the melt given the mole
		fraction of carbonate ion.

		Parameters
		----------
		XCO3: float
			Mole fraction of carbonate ion in the melt."""


		XCO3 = self.molfrac_molecular(pressure=pressure,sample=sample,X_fluid=X_fluid,**kwargs)
		return (4400 * XCO3) / (36.6 - 44*XCO3)


	def calculate_equilibrium_fluid_comp(self,pressure,sample,**kwargs):
		"""
		WRITE STUFF
		"""
		if self.calculate_saturation_pressure(sample=sample,**kwargs) > pressure:
			return 0.0
		else:
			return 1.0

	def calculate_saturation_pressure(self,sample,X_fluid=1.0,**kwargs):
		"""
		WRITE STUFF
		"""
		satP = root_scalar(self.root_saturation_pressure,x0=1000.0,x1=2000.0,args=(sample,kwargs)).root
		return np.real(satP)

	def molfrac_molecular(self,pressure,sample,X_fluid=1.0,**kwargs):
		"""Dixon (1997) Eqn (1).

		Calculates the mole fraction of carbonate ion dissolved in a tholeiitic basalt in
		equilibrium with pure CO2 fluid at 1200C, at a chosen pressure.

		Parameters
		----------
		P:  float
			Pressure in bar
		XCO3Std:    float
			Mole fraction of carbonate dissolved in tholeiitic basalt in equilibrium with
			pure CO2 fluid at 1200C and 1 bar. Default value is 3.8e-7, as suggested in paper.
			If MajorElements variable not set to string, this will be replaced.
		fCO2:   float
			Fugacity of CO2 at the pressure of interest. If value set to 0.0, fCO2 will be
			set to that calculated using the Pitzer & Sterner equations by the CHO software.
			This is the default option.
		MajorElements:  series
			If set to string XH2OStd variable passed will be used (or default value).
			Series of major element oxide concentrations in wt%, requiring SiO2 as minimum.
			Used to calculated XCO3std from compositional parametrisation."""

		DeltaVr = 23 #cm3 mole-1
		P0 = 1
		R = 83.15
		T0 = 1473.15

		fugacity = self.fugacity_model.fugacity(pressure=pressure,X_fluid=X_fluid,**kwargs)

		XCO3Std = self.XCO3_Std(sample)

		return XCO3Std * fugacity * np.exp(-DeltaVr * (pressure-P0)/(R*T0))

	def XCO3_Std(self,sample):
		"""Dixon (1997) Eq (8).

		The compositional parameterisation for the mole fraction of carbonate ions
		dissolved in a melt in equilibrium with CO2 vapour at 1200C and 1 bar.

		Parameters
		----------
		MajorElements:    dict or series (or str)
			A dictionary or series containing 'SiO2' as a label for the SiO2 content
			of the melt. If the variable is a str, it will return the default XCO3_std
			value."""
		return 8.7e-6 - 1.7e-7*sample['SiO2']

	def root_saturation_pressure(self,pressure,sample,kwargs):
		return self.calculate_dissolved_volatiles(pressure=pressure,sample=sample,**kwargs) - sample['CO2']

	def check_calibration_range(self,parameters,**kwargs):
		return 0



class DixonWater(Model):
	"""
	"""

	def __init__(self):
		self.set_volatile_species('H2O')
		self.set_fugacity_model(fugacity_KJ81_h2o())
		self.set_activity_model(activity_idealsolution())

	def preprocess_sample(self,sample):
		return sample

	def calculate_dissolved_volatiles(self,pressure,sample,X_fluid=1.0,**kwargs):
		"""Dixon (1997) Eq (5) and Eq (6).

		Calculates total H2O (wt%) dissolved in the melt, given the molar fractions
		of molecular water and hydroxyl groups dissolved.

		Parameters
		----------
		XH2O: float
			Mole fraction of molecular water dissolved in the melt.
		XOH: float
			Mole fraction of hydroxyl groups dissolved in the melt."""
		XH2O = self.molfrac_molecular(pressure=pressure,sample=sample,X_fluid=X_fluid,**kwargs)
		XOH = self.XOH(pressure=pressure,sample=sample,X_fluid=X_fluid,**kwargs)

		XB = XH2O + 0.5*XOH
		return 1801.5*XB/(36.6-18.6*XB)


	def calculate_equilibrium_fluid_comp(self,pressure,sample,**kwargs):
		"""
		WRITE STUFF
		"""
		if self.calculate_saturation_pressure(sample=sample,**kwargs) > pressure:
			return 0.0
		else:
			return 1.0

	def calculate_saturation_pressure(self,sample,X_fluid=1.0,**kwargs):
		"""
		WRITE STUFF
		"""
		satP = root_scalar(self.root_saturation_pressure,x0=1000.0,x1=2000.0,args=(sample,kwargs)).root

		return np.real(satP)

	def molfrac_molecular(self,pressure,sample,X_fluid=1.0,**kwargs):
		"""Dixon (1997) Eqn (2).

		Calculates the mole fraction of molecular H2O dissolved in a tholeiitic basalt in
		equilibrium with pure H2O fluid at 1200C, at a chosen pressure.

		Parameters
		----------
		P: float
			Pressure in bar
		XH2OStd: float
			Mole fraction of cmolecular H2O in tholeiitic basalt in equilibrium with
			pure H2O fluid at 1200C and 1 bar. Default value is 3.28e-5, as suggested in paper.
			If MajorElements variable not set to string, this will be replaced.
		fH2O: float
			Fugacity of H2O at the pressure of interest. If value set to 0.0, fH2O will be
			set to that calculated using the Pitzer & Sterner equations by the CHO software.
			This is the default option.
		MajorElements:  series
			If set to string XH2OStd variable passed will be used (or default value).
			Series of major element oxide concentrations in wt%, requiring SiO2 as minimum.
			Used to calculated XH2Ostd from compositional parametrisation."""
		VH2O = 12 #cm3 mole-1
		P0 = 1
		R = 83.15
		T0 = 1473.15

		XH2OStd = self.XH2O_Std(sample)

		fugacity = self.fugacity_model.fugacity(pressure=pressure,X_fluid=X_fluid,**kwargs)

		return XH2OStd * fugacity * np.exp(-VH2O * (pressure-P0)/(R*T0))

	def XH2O_Std(self,sample):
		"""Dixon (1997) Eq (9).

		The compositional parameterisation for the mole fraction of molecular water
		dissolved in a melt in equilibrium with H2O vapour at 1200C and 1 bar.

		Parameters
		----------
		MajorElements:    dict or series (or str)
			A dictionary or series containing 'SiO2' as a label for the SiO2 content
			of the melt. If the variable is a str, it will return the default XH2O_std
			value."""
		return -3.04e-5 + 1.29e-6*sample['SiO2']

	def XOH(self,pressure,sample,X_fluid=1.0,**kwargs):
		"""Solves Dixon (1997) Eq (4).

		Finds the mole fraction of hydroxyl groups dissolved in the melt, given the
		mole fraction of molecular water. Uses the Newton method to find the root of
		Eq (4).

		Often struggles to find root when XH2O -> 0.

		Parameters
		----------
		XH2O: float
			Mole fraction of molecular water dissolved in melt."""

		XH2O = self.molfrac_molecular(pressure=pressure,sample=sample,X_fluid=X_fluid,**kwargs)
		if XH2O < 1e-14:
			return 0
		else:
			return root_scalar(self.XOH_root,bracket=(1e-14,1),args=(XH2O)).root

	def XOH_root(self,XOH,XH2O):
		"""Dixon (1997) Eq (4).

		Returns the difference between the two sides of equation (4) which relates
		the mole fraction of hydroxyl groups to the mole fraction of molecular water.
		For use in finding the root of the equation.

		It is a regular solution model (Silver and Stolper, 1989).

		Parameters
		----------
		XOH: float
			Mole fraction of hydroxyl groups dissolved in melt.
		XH2O: float
			Mole fraction of molecular water dissolved in melt."""
		A = 0.403
		B = 15.333
		C = 10.894

		lhs = - np.log(XOH**2.0/(XH2O*(1.0-XH2O)))
		rhs = A + B*XOH + C*XH2O

		return rhs - lhs

	def root_saturation_pressure(self,pressure,sample,kwargs):
		return self.calculate_dissolved_volatiles(pressure=pressure,sample=sample,**kwargs) - sample['H2O']

	def check_calibration_range(self,**kwargs):
		return 0

class IaconoMarzianoWater(Model):
	"""
	"""

	def __init__(self):
		self.set_volatile_species('H2O')
		self.set_fugacity_model(fugacity_idealgas())
		self.set_activity_model(activity_idealsolution())

	def preprocess_sample(self,sample):
		return sample

	def calculate_dissolved_volatiles(self,pressure,temperature,sample,X_fluid=1,hydrous=True,**kwargs):
		"""
		Iacono-Marziano (2012) Eq (8). Calculates the concentration of CO2 in
		ppm dissolved in the melt at vapour saturation, for specified pressure,
		temperature and partial pressure of CO2. Modifies the output of the equation
		to give CO2 rather than CO3.

		Parameters
		----------
		P: float
			Pressure in bar
		T: float
			Temperature in K
		MajorElements: series
			Series labelled with major element oxides in wt%. All oxides passed will be
			included in calculating molar proportions.
		P_CO2: float
			Partial pressure of CO2 in bar
		hydrous: bool
			Use the hydrous or anhydrous calibration
		"""

		if hydrous == True:
			if X_fluid==0:
				return 0
			H2O = root_scalar(self.root_dissolved_volatiles,args=(pressure,temperature,sample,X_fluid,kwargs),
								x0=1.0,x1=2.0).root
			return H2O/(100+H2O)*100
		else:
			a = 0.54
			b = 1.24
			B = -2.95
			C = 0.02

			fugacity = self.fugacity_model.fugacity(pressure=pressure,X_fluid=X_fluid,temperature=temperature,**kwargs)
			if fugacity == 0:
				return 0
			NBO_O = self.NBO_O(sample=sample,hydrous=False)

			H2O = np.exp(a*np.log(fugacity) + b*NBO_O + B + C*pressure/temperature)

			return H2O


	def calculate_equilibrium_fluid_comp(self,pressure,temperature,sample,hydrous=True,**kwargs):
		"""
		WRITE STUFF
		"""
		if pressure > self.calculate_saturation_pressure(temperature=temperature,sample=sample,hydrous=hydrous,**kwargs):
			return 0.0
		else:
			return 1.0

	def calculate_saturation_pressure(self,temperature,sample,hydrous=True,**kwargs):
		"""
		WRITE STUFF
		"""
		return root_scalar(self.root_saturation_pressure,args=(temperature,sample,hydrous,kwargs),
							x0=1000.0,x1=2000.0).root

	def root_saturation_pressure(self,pressure,temperature,sample,hydrous,kwargs):
		return sample['H2O'] - self.calculate_dissolved_volatiles(pressure=pressure,temperature=temperature,sample=sample,hydrous=hydrous,**kwargs)


	def root_dissolved_volatiles(self,h2o,pressure,temperature,sample,X_fluid,kwargs):
		a = 0.53
		b = 2.35
		B = -3.37
		C = -0.02

		sample_h2o = sample.copy()
		sample_h2o['H2O'] = h2o
		NBO_O = self.NBO_O(sample=sample_h2o,hydrous=True)
		fugacity = self.fugacity_model.fugacity(pressure=pressure,X_fluid=X_fluid,temperature=temperature,**kwargs)

		return h2o - np.exp(a*np.log(fugacity) + b*NBO_O + B + C*pressure/temperature)

	def NBO_O(self,sample,hydrous):
		X = wtpercentOxides_to_molOxides(sample)

		NBO = 2*(X['K2O']+X['Na2O']+X['CaO']+X['MgO']+X['FeO']-X['Al2O3'])
		O = 2*X['SiO2']+2*X['TiO2']+3*X['Al2O3']+X['MgO']+X['FeO']+X['CaO']+X['Na2O']+X['K2O']

		if hydrous == True:
			NBO = NBO + 2*X['H2O']
			O = O + X['H2O']

		return NBO/O

	def check_calibration_range(self,**kwargs):
		return 0

class IaconoMarzianoCarbon(Model):
	"""
	"""

	def __init__(self):
		self.set_volatile_species('CO2')
		self.set_fugacity_model(fugacity_idealgas())
		self.set_activity_model(activity_idealsolution())

	def preprocess_sample(self,sample):
		return sample

	def calculate_dissolved_volatiles(self,pressure,temperature,sample,X_fluid=1,hydrous=True,**kwargs):
		"""
		Iacono-Marziano (2012) Eq (8). Calculates the concentration of CO2 in
		ppm dissolved in the melt at vapour saturation, for specified pressure,
		temperature and partial pressure of CO2. Modifies the output of the equation
		to give CO2 rather than CO3.

		Parameters
		----------
		P: float
			Pressure in bar
		T: float
			Temperature in K
		MajorElements: series
			Series labelled with major element oxides in wt%. All oxides passed will be
			included in calculating molar proportions.
		P_CO2: float
			Partial pressure of CO2 in bar
		hydrous: bool
			Use the hydrous or anhydrous calibration
		"""

		if hydrous == True:
			im_h2o_model = IaconoMarzianoWater()
			h2o = im_h2o_model.calculate_dissolved_volatiles(pressure=pressure,temperature=temperature,
														sample=sample,X_fluid=1-X_fluid,hydrous=True,**kwargs)
			sample_h2o = sample.copy()
			sample_h2o['H2O'] = h2o

			d = np.array([-16.4,4.4,-17.1,22.8])
			a = 1.0
			b = 17.3
			B = -6.0
			C = 0.12

			NBO_O = self.NBO_O(sample=sample_h2o,hydrous=True)

			molarProps = wtpercentOxides_to_molOxides(sample_h2o)

		else:
			d = np.array([2.3,3.8,-16.3,20.1])
			a = 1.0
			b = 15.8
			B = -5.3
			C = 0.14

			NBO_O = self.NBO_O(sample=sample,hydrous=False)

			molarProps = wtpercentOxides_to_molOxides(sample)

		fugacity = self.fugacity_model.fugacity(pressure=pressure,X_fluid=X_fluid,temperature=temperature,**kwargs)

		if fugacity == 0:
			return 0


		x = list()
		x.append(molarProps['H2O'])
		x.append(molarProps['Al2O3']/(molarProps['CaO']+molarProps['K2O']+molarProps['Na2O']))
		x.append((molarProps['FeO']+molarProps['MgO']))
		x.append((molarProps['Na2O']+molarProps['K2O']))
		x = np.array(x)

		CO3 = np.exp(np.sum(x*d) + a*np.log(fugacity) + b*NBO_O + B + C*pressure/temperature)

		CO2 = CO3/1e4#/(12+16*3)*(12+16*2)/1e4

		return CO2


	def calculate_equilibrium_fluid_comp(self,pressure,temperature,sample,hydrous=True,**kwargs):
		"""
		WRITE STUFF
		"""
		if pressure > self.calculate_saturation_pressure(temperature=temperature,sample=sample,hydrous=hydrous,**kwargs):
			return 0.0
		else:
			return 1.0

	def calculate_saturation_pressure(self,temperature,sample,hydrous=True,**kwargs):
		"""
		WRITE STUFF
		"""
		return root_scalar(self.root_saturation_pressure,args=(temperature,sample,hydrous,kwargs),
							x0=1000.0,x1=2000.0).root

	def root_saturation_pressure(self,pressure,temperature,sample,hydrous,kwargs):
		return sample['CO2'] - self.calculate_dissolved_volatiles(pressure=pressure,temperature=temperature,sample=sample,hydrous=hydrous,**kwargs)


	def NBO_O(self,sample,hydrous):
		X = wtpercentOxides_to_molOxides(sample)

		NBO = 2*(X['K2O']+X['Na2O']+X['CaO']+X['MgO']+X['FeO']-X['Al2O3'])
		O = 2*X['SiO2']+2*X['TiO2']+3*X['Al2O3']+X['MgO']+X['FeO']+X['CaO']+X['Na2O']+X['K2O']

		if hydrous == True:
			NBO = NBO + 2*X['H2O']
			O = O + X['H2O']

		return NBO/O

	def check_calibration_range(self,**kwargs):
		return 0

class EguchiCarbon(Model):
	"""
	"""

	def __init__(self):
		self.set_volatile_species('CO2')
		self.set_fugacity_model(fugacity_KJ81_co2())
		print("Warning: Eguchi and Dasgupta model should use the Zhang and Duan EOS.")
		self.set_activity_model(activity_idealsolution())

	def preprocess_sample(self,sample):
		return sample

	def calculate_dissolved_volatiles(self,pressure,temperature,sample,X_fluid=1.0,**kwargs):
		"""
		WRITE STUFF
		"""
		XCO3 = self.Xi_melt(pressure=pressure,temperature=temperature,sample=sample,X_fluid=X_fluid,species='CO3')
		XCO2 = self.Xi_melt(pressure=pressure,temperature=temperature,sample=sample,X_fluid=X_fluid,species='CO2')

		FW_one = wtpercentOxides_to_formulaWeight(sample)

		CO2_CO2 = ((44.01*XCO2)/(44.01*XCO2+(1-(XCO2+XCO3))*FW_one))*100
		CO2_CO3 = ((44.01*XCO3)/(44.01*XCO3+(1-(XCO2+XCO3))*FW_one))*100

		return CO2_CO2 + CO2_CO3


	def calculate_equilibrium_fluid_comp(self,pressure,temperature,sample,**kwargs):
		"""
		WRITE STUFF
		"""
		satP = self.calculate_saturation_pressure(temperature=temperature,sample=sample,X_fluid=1.0,**kwargs)
		if pressure > satP:
			return 1.0
		else:
			return 0.0

	def calculate_saturation_pressure(self,temperature,sample,X_fluid=1.0,**kwargs):
		"""
		WRITE STUFF
		"""
		return root_scalar(self.root_saturation_pressure,x0=1000.0,x1=2000.0,args=(temperature,sample,X_fluid,kwargs)).root

	def root_saturation_pressure(self,pressure,temperature,sample,X_fluid,kwargs):
		return sample['CO2'] - self.calculate_dissolved_volatiles(pressure=pressure,temperature=temperature,sample=sample,X_fluid=X_fluid,**kwargs)

	def Xi_melt(self,pressure,temperature,sample,species,X_fluid=1.0,**kwargs):
		"""
		Equation (9) in Eguchi and Dasgupta (2018). Calculates the mole fraction
		of dissolved molecular CO2 or carbonate CO3(2-).

		Parameters
		----------
		P (float):  Pressure (bar)
		T (float):  Temperature (K)
		MajorElements (pandas Series):  Major Element Oxides.
		species (str):  Which species to calculate, molecular CO2 or carbonate ion.

		Returns
		-------
		Xi (float):     Mole fraction of selected species in the melt

		"""

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
			print("WRITE CODE TO RAISE AN ERROR!")
		R = 8.314


		# Calculate NBO term
		cations = wtpercentOxides_to_molSingleO(sample)
		oxides = wtpercentOxides_to_molOxides(sample)

		for cation in ['Mg','Ca','Fe','Na','K','Mn','Fe3']:
			if cation not in list(cations.keys()):
				cations[cation] = 0

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

		lnfCO2 = np.log(self.fugacity_model.fugacity(pressure=pressure,temperature=temperature,X_fluid=X_fluid))

		lnXi = ((DH/(R*temperature)-(pressure*1e5*DV)/(R*temperature)+DS/R) +
				(A_CaO*oxides['CaO']+A_Na2O*oxides['Na2O']+A_K2O*oxides['K2O'])/(R*temperature) +
				(B*lnfCO2/temperature) + yNBO*NBO
				)

		return np.exp(lnXi)

	def check_calibration_range(self,**kwargs):
		return 0

class AllisonCarbon(Model):
	"""
	"""

	def __init__(self,model_fit='power',model_loc='sunset'):
		self.set_volatile_species('CO2')
		self.set_fugacity_model(fugacity_KJ81_co2())
		self.set_activity_model(activity_idealsolution())
		self.model_fit = model_fit
		self.model_loc = model_loc

	def preprocess_sample(self,sample):
		return sample

	def calculate_dissolved_volatiles(self,pressure,temperature,sample=None,X_fluid=1.0,**kwargs):
		"""
		WRITE STUFF
		"""
		if self.model_fit == 'thermodynamic':
			if type(sample) == type(None):
				print("Generate an error...")
			P0 = 1000 # bar
			params = dict({'sunset':[16.4,-14.67],
							'sfvf':[15.02,-14.87],
							'erebus':[15.83,-14.65],
							'vesuvius':[24.42,-14.04],
							'etna':[21.59,-14.28],
							'stromboli':[14.93,-14.68]})
			DV = params[self.model_loc][0]
			lnK0 = params[self.model_loc][1]

			lnK = lnK0 - (pressure-P0)*DV/(8.3141*temperature)
			fCO2 = self.fugacity_model.fugacity(pressure=pressure,temperature=temperature,X_fluid=X_fluid,**kwargs)
			print(fCO2)
			Kf = np.exp(lnK)*fCO2
			XCO3 = Kf/(1-Kf)
			FWone = wtpercentOxides_to_formulaWeight(sample)
			wtCO2 = (44.01*XCO3)/((44.01*XCO3)+(1-XCO3)*FWone)*100

			return wtCO2
		if self.model_fit == 'power':
			params = dict({'stromboli':[1.05,0.883],
							'etna':[2.831,0.797],
							'vesuvius':[4.796,0.754],
							'sfvf':[3.273,0.74],
							'sunset':[4.32,0.728],
							'erebus':[5.145,0.713]})

			fCO2 = self.fugacity_model.fugacity(pressure=pressure,temperature=temperature,X_fluid=X_fluid,**kwargs)

			return params[self.model_loc][0]*fCO2**params[self.model_loc][1]/1e4



	def calculate_equilibrium_fluid_comp(self,pressure,temperature,sample,**kwargs):
		"""
		WRITE STUFF
		"""
		satP = self.calculate_saturation_pressure(temperature=temperature,sample=sample,X_fluid=1.0,**kwargs)
		if pressure > satP:
			return 1.0
		else:
			return 0.0

	def calculate_saturation_pressure(self,temperature,sample,X_fluid=1.0,**kwargs):
		"""
		WRITE STUFF
		"""
		return root_scalar(self.root_saturation_pressure,args=(temperature,sample,X_fluid,kwargs),x0=1000.0,x1=2000.0).root

	def root_saturation_pressure(self,pressure,temperature,sample,X_fluid,kwargs):
		return sample['CO2'] - self.calculate_dissolved_volatiles(pressure=pressure,temperature=temperature,sample=sample,X_fluid=X_fluid,**kwargs)

	def check_calibration_range(self,**kwargs):
		return 0


#------------MIXED FLUID MODELS-------------------------------#
class MixedFluids(Model):
	""" HELLO!
	"""
	def __init__(self,models):
		print('Write code to check models input makes sense.')
		self.models = tuple(model for model in models.values())
		self.set_volatile_species(list(models.keys()))

	def preprocess_sample(self,sample):
		return sample

	def calculate_dissolved_volatiles(self,pressure,X_fluid,**kwargs):
		if len(X_fluid) != len(self.volatile_species):
			print('YOU NEED TO WRITE THE CODE TO RAISE AN ERROR! 0')

		if type(X_fluid) == dict or type(X_fluid) == pd.core.series.Series:
			X_fluid = tuple(X_fluid[species] for species in self.volatile_species)
		elif type(X_fluid) != tuple and type(X_fluid) != list:
			print('YOU NEED TO WRITE THE CODE TO RAISE AN ERROR! 1')

		return tuple(model.calculate_dissolved_volatiles(pressure=pressure,X_fluid=Xi,**kwargs) for model, Xi in zip(self.models,X_fluid))

	def calculate_equilibrium_fluid_comp(self,pressure,sample,**kwargs):
		if len(self.volatile_species) != 2:
			print('WRITE CODE TO RAISE AN ERROR- CAN ONLY HANDLE TWO VOLATILE SPECIES')

		satP = self.calculate_saturation_pressure(sample,**kwargs)

		if satP < pressure:
			return (0,0)

		sample_mod = sample.copy()
		sample_mod = normalize_FixedVolatiles(sample_mod)

		molfracs = wtpercentOxides_to_molOxides(sample_mod)
		(Xt0, Xt1) = (molfracs[self.volatile_species[0]],molfracs[self.volatile_species[1]])


		Xv0 = root_scalar(self.root_for_fluid_comp,args=(pressure,Xt0,Xt1,sample,kwargs),bracket=(0,1)).root
		Xv1 = 1-Xv0

		return (Xv0,Xv1)

	def calculate_saturation_pressure(self,sample,**kwargs):
		sample_mod = normalize_FixedVolatiles(sample)
		volatile_concs = np.array(tuple(sample_mod[species] for species in self.volatile_species))

		x0 = 0
		for model in self.models:
			x0 = x0 + model.calculate_saturation_pressure(sample=sample_mod,**kwargs)

		satP = root(self.root_saturation_pressure,x0=[x0,0.5],args=(volatile_concs,sample_mod,kwargs)).x[0]

		return satP

	def calculate_isobars_and_isopleths(self,pressure_list,isopleth_list=[0,1],points=101,**kwargs):
		isobars = []
		for pressure in pressure_list:
			dissolved = np.zeros([2,points])
			Xv0 = np.linspace(0.0,1.0,points)
			for i in range(points):
				dissolved[:,i] = self.calculate_dissolved_volatiles(pressure=pressure,X_fluid=(Xv0[i],1-Xv0[i]),**kwargs)
			isobars.append(dissolved)

		isopleths = []
		for isopleth in isopleth_list:
			dissolved = np.zeros([2,points])
			pressure = np.linspace(np.nanmin(pressure_list),np.nanmax(pressure_list),points)
			for i in range(points):
				dissolved[:,i] = self.calculate_dissolved_volatiles(pressure=pressure[i],X_fluid=(isopleth,1-isopleth),**kwargs)
			isopleths.append(dissolved)

		return (isobars,isopleths)

	def calculate_degassing_paths(self,pressures,sample,fractionate_vapor=1.0,**kwargs):

		Xv = np.zeros([2,len(pressures)])
		wtm = np.zeros([2,len(pressures)])

		wtptoxides = sample.copy()
		wtptoxides = normalize_FixedVolatiles(wtptoxides)
		wtm0s, wtm1s = (wtptoxides[self.volatile_species[0]],wtptoxides[self.volatile_species[1]])

		for i in range(len(pressures)):
			X_fluid = self.calculate_equilibrium_fluid_comp(pressure=pressures[i],sample=wtptoxides,**kwargs)
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
				wtptoxides = normalize_FixedVolatiles(wtptoxides)

		return wtm, Xv


	def root_saturation_pressure(self,x,volatile_concs,sample,kwargs):
		misfit = np.array(self.calculate_dissolved_volatiles(pressure=x[0],X_fluid=(x[1],1-x[1]),sample=sample,**kwargs)) - volatile_concs
		return misfit


	def root_for_fluid_comp(self,Xv0,pressure,Xt0,Xt1,sample,kwargs):
		# print("Need to fix volatile normalization issue.")
		wtm0, wtm1 = self.calculate_dissolved_volatiles(pressure=pressure,X_fluid=(Xv0,1-Xv0),sample=sample,**kwargs)
		sample_mod = sample.copy()
		sample_mod[self.volatile_species[0]] = wtm0
		sample_mod[self.volatile_species[1]] = wtm1
		sample_mod = normalize_FixedVolatiles(sample_mod)
		cations = wtpercentOxides_to_molOxides(sample_mod)
		Xm0 = cations[self.volatile_species[0]]
		Xm1 = cations[self.volatile_species[1]]
		if Xv0 == 0:
			return Xt0 - Xm0*(Xt1-1)/(Xm1-1)
		elif Xv0 == 1:
			return -(Xt1 - Xm1*(Xt0-1)/(Xm0-1))
		return (Xt0-Xm0)/(Xv0-Xm0) - (Xt1-Xm1)/(1-Xv0-Xm1)

	def check_calibration_range(self,parameters,**kwargs):
		results = {}
		for species, model in zip(self.volatile_species,self.models):
			results[species] = model.check_calibration_range(parameters=parameters)
		return results

class MagmaSat(Model):
	"""
	An object to instantiate a thermoengine equilibrate class
	"""

	def __init__(self):
		self.melts_version = '1.2.0'

	def preprocess_sample(self,sample):
		#--------------MELTS preamble---------------#
		# instantiate thermoengine equilibrate MELTS instance
		melts = equilibrate.MELTSmodel('1.2.0')

		# Suppress phases not required in the melts simulation
		self.oxides = melts.get_oxide_names()
		self.phases = melts.get_phase_names()

		for phase in self.phases:
			melts.set_phase_inclusion_status({phase: False})
		melts.set_phase_inclusion_status({'Fluid': True, 'Liquid': True})
		self.melts = melts
		#-------------------------------------------#
		oxides = self.oxides
		for oxide in oxides:
			if oxide in sample.keys():
				pass
			else:
				sample[oxide] = 0.0

		return sample

	def check_calibration_range(self,**kwargs):
		return 0

	def calculate_dissolved_volatiles(self, sample, temperature, pressure, XH2Ofluid):
		"""
		Calculates the amount of H2O and CO2 dissolved in a magma at the given P/T conditions and fluid composition.

		Parameters
		----------
		sample: dict or pandas Series
			Compositional information on one sample in oxides.

		temperature: float or int
			Temperature, in degrees C. 

		presure: float or int
			Pressure, in bars. 

		XH2Ofluid: float or int
			The mole fraction of H2O in the H2O-CO2 fluid. XH2Ofluid=1 is a pure H2O fluid. XH2O=0 is a pure CO2 fluid.

		Returns
		-------
		dict
			A dictionary of dissolved volatile concentrations in wt% with keys H2O and CO2.
		"""
		melts = self.melts
		oxides = self.oxides

		pressureMPa = pressure / 10.0

		bulk_comp = {oxide:  sample[oxide] for oxide in oxides}
		feasible = melts.set_bulk_composition(bulk_comp)

		output = melts.equilibrate_tp(temperature, pressureMPa, initialize=True)

		pass #TODO Write this code...

	def calculate_equilibrium_fluid_comp(self, sample, temperature, pressure, verbose=False): #TODO fix weird printing
		"""
		Returns H2O and CO2 concentrations in wt% in a fluid in equilibrium with the given sample at the given P/T condition.

		Parameters
		----------
		sample: dict or pandas Series
			Compositional information on one sample in oxides.

		temperature: float or int
			Temperature, in degrees C. 

		presure: float or int
			Pressure, in bars. #TODO check units

		verbose: bool
			OPTIONAL: Default is False. If set to True, returns H2O and CO2 concentration in the fluid, mass of the fluid in grams, 
			and proportion of fluid in the system in wt%.

		Returns
		-------
		dict 
			A dictionary of fluid composition in wt% with keys 'H2O' and 'CO2' is returned. #TODO make list?
		"""
		melts = self.melts
		oxides = self.oxides

		if isinstance(temperature, float) or isinstance(temperature, int):
			pass
		else:
			raise InputError("temp must be type float or int")

		if isinstance(pressure, float) or isinstance(pressure, int):
			pass
		else:
			raise InputError("presure must be type float or int")

		pressureMPa = pressure / 10.0

		bulk_comp = {oxide:  sample[oxide] for oxide in oxides}
		feasible = melts.set_bulk_composition(bulk_comp)

		output = melts.equilibrate_tp(temperature, pressureMPa, initialize=True)
		(status, temperature, pressureMPa, xmlout) = output[0]
		fluid_mass = melts.get_mass_of_phase(xmlout, phase_name='Fluid')
		flsystem_wtper = 100 * fluid_mass / (fluid_mass + melts.get_mass_of_phase(xmlout, phase_name='Liquid'))

		if fluid_mass > 0.0:
			fluid_comp = melts.get_composition_of_phase(xmlout, phase_name='Fluid')
			fluid_comp_H2O = fluid_comp['H2O']
			fluid_comp_CO2 = fluid_comp['CO2']
		else:
			fluid_comp_H2O = 0
			fluid_comp_CO2 = 0

		feasible = melts.set_bulk_composition(bulk_comp) #reset

		if verbose == False:
			return {'H2O': fluid_comp_H2O, 'CO2': fluid_comp_CO2}

		if verbose == True:
			return {'H2O': fluid_comp_H2O, 'CO2': fluid_comp_CO2, 'FluidMass_grams': fluid_mass, 'FluidProportion_wtper': flsystem_wtper}

	def calculate_isobars_and_isopleths(self, sample, temperature, print_status=False, pressure_min='', pressure_max='', pressure_int='', pressure_list=''):
		"""
		Calculates isobars and isopleths at a constant temperature for a given sample. Isobars can be calculated
		for any number of pressures. Pressures can be passed as min, max, interval (100.0, 500.0, 100.0 would result
		in pressures of 100.0, 200.0, 300.0, 400.0, and 500.0 MPa). Alternatively pressures can be passed as a list of all
		desired pressures ([100.0, 200.0, 250.0, 300.0] would calculate isobars for each of those pressures in MPa).

		Parameters
		----------
		sample: dict
			Dictionary or pandas Series. with values for sample composition as oxides in wt%.

		temperature: float
			Temperature in degrees C.

		pressure_min: float
			OPTIONAL. If passed, also requires pressure_max and pressure_int be passed. If passed, do not pass
			pressure_list. Minimum pressure	value in MPa.

		pressure_max: float
			OPTIONAL. If passed, also requires pressure_min and pressure_int be passed. If passed, do not pass
			pressure_list. Maximum pressure value in MPa.

		pressure_int: float
			OPTIONAL: If passed, also requires pressure_min and pressure_max be passed. If passed, do not pass
			pressure_list. Interval between pressure values in MPa.

		pressure_list: list
			OPTIONAL: If passed, do not pass pressure_min, pressure_max, or pressure_int. List of all pressure
			values in MPa.

		print_status: bool
			OPTIONAL: Default is False. If set to True, progress of the calculations will be printed to the terminal. 

		Returns
		-------
		pandas DataFrame object
			DataFrame containing calcualted isobar and isopleth information for the passed melt composition. Column titles
			are 'Pressure', 'H2Omelt', 'CO2melt', 'H2Ofl', and 'CO2fl'.
		"""
		#--------------Preamble required for every MagmaSat method within Modeller---------------#
		# instantiate thermoengine equilibrate MELTS instance
		melts = equilibrate.MELTSmodel(self.model_version)
		bulk_comp_orig = sample #for reset

		# Suppress phases not required in the melts simulation
		self.oxides = melts.get_oxide_names()
		self.phases = melts.get_phase_names()

		for phase in self.phases:
			melts.set_phase_inclusion_status({phase: False})
		melts.set_phase_inclusion_status({'Fluid': True, 'Liquid': True})
		#---------------------------------------------------------------------------------------#

		if isinstance(pressure_min, float) and isinstance(pressure_list, list):
			raise InputError(
				"Enter pressure either as min, max, int OR as list. Not both.")
		if isinstance(pressure_max, float) and isinstance(pressure_list, list):
			raise InputError(
				"Enter pressure either as min, max, int OR as list. Not both.")
		if isinstance(pressure_int, float) and isinstance(pressure_list, list):
			raise InputError(
				"Enter pressure either as min, max, int OR as list. Not both.")

		if isinstance(pressure_min, float):
			P_vals = np.arange(pressure_min, pressure_max+pressure_int, pressure_int)

		if isinstance(pressure_list, list):
			P_vals = pressure_list

		bulk_comp = sample
		phases = self.phases
		oxides = self.oxides
		phases = melts.get_phase_names()

		volatiles_at_saturation = []
		H2O_val = 0
		CO2_val = 0
		fluid_mass = 0.0

		# Calculate equilibrium phase assemblage for all P/T conditions, check if saturated in fluid...
		for i in P_vals:
			if print_status == True:
				print("Calculating isobars at " + str(i) + " MPa")

			for j in np.arange(0, 15.5, 0.5):
				bulk_comp["H2O"] = j
				while fluid_mass <= 0.0:
					bulk_comp["CO2"] = CO2_val

					melts.set_bulk_composition(bulk_comp)

					output = melts.equilibrate_tp(temperature, i, initialize=True)
					(status, temperature, i, xmlout) = output[0]
					fluid_mass = melts.get_mass_of_phase(xmlout, phase_name='Fluid')

					CO2_val = CO2_val + 0.1

				if fluid_mass > 0.0:
					liquid_comp = melts.get_composition_of_phase(
						xmlout, phase_name='Liquid', mode='oxide_wt')
					fluid_comp = melts.get_composition_of_phase(
						xmlout, phase_name='Fluid')

					if "H2O" in liquid_comp:
								H2O_liq = liquid_comp["H2O"]
					else:
						H2O_liq = 0

					if "CO2" in liquid_comp:
						CO2_liq = liquid_comp["CO2"]
					else:
						CO2_liq = 0

					if "H2O" in fluid_comp:
						H2O_fl = fluid_comp["H2O"]
					else:
						H2O_fl = 0.0
					if "CO2" in fluid_comp:
						CO2_fl = fluid_comp["CO2"]
					else:
						CO2_fl = 0.0
					volatiles_at_saturation.append(
						[i, H2O_liq, CO2_liq, H2O_fl, CO2_fl])
					CO2_val = 0.0
					fluid_mass = 0.0

		if print_status == True:
			print("Done!")

		isobars_df = pd.DataFrame(volatiles_at_saturation, columns=[
								  'Pressure', 'H2Omelt', 'CO2melt', 'H2Ofl', 'CO2fl'])

		feasible = melts.set_bulk_composition(bulk_comp_orig) #reset

		return isobars_df

	def plot_isobars_and_isopleths(self, isobars_df):
		"""
		Takes in a dataframe with calculated isobar and isopleth information (e.g., output from calculate_isobars_and_isopleths)
		and plots data as isobars (lines of constant pressure) and isopleths (lines of constant fluid composition). These lines
		represent the saturation pressures of the melt composition used to calculate the isobar and isopleth information.

		Parameters
		----------
		isobars_df: pandas DataFrame
			DataFrame object containing isobar and isopleth information as calculated by calculate_isobars_and_isopleths.

		Returns
		-------
		matplotlib object
			Plot with x-axis as H2O wt% in the melt and y-axis as CO2 wt% in the melt. Isobars, or lines of
			constant pressure at which the sample magma composition is saturated, and isopleths, or lines of constant
			fluid composition at which the sample magma composition is saturated, are plotted.
		"""
		#--------------Preamble required for every MagmaSat method within Modeller---------------#
		# instantiate thermoengine equilibrate MELTS instance
		melts = equilibrate.MELTSmodel(self.model_version)

		# Suppress phases not required in the melts simulation
		self.oxides = melts.get_oxide_names()
		self.phases = melts.get_phase_names()

		for phase in self.phases:
			melts.set_phase_inclusion_status({phase: False})
		melts.set_phase_inclusion_status({'Fluid': True, 'Liquid': True})
		#---------------------------------------------------------------------------------------#

		P_vals = isobars_df.Pressure.unique()
		isobars_lists = isobars_df.values.tolist()

		# make a list of isopleth values to plot
		iso_step = 20.0
		isopleth_vals = np.arange(0+iso_step, 100.0, iso_step)

		# add zero values to volatiles list
		isobars_lists.append([0.0, 0.0, 0.0, 0.0])

		# draw the figure
		fig, ax1 = plt.subplots()

		# turn on interactive plotting
		plt.ion()

		plt.xlabel('H2O wt%')
		plt.ylabel('CO2 wt%')

		# Plot some stuff
		for pressure in P_vals:
			ax1.plot([item[1] for item in isobars_lists if item[0] == pressure],
					 [item[2] for item in isobars_lists if item[0] == pressure])

		for val in isopleth_vals:
			val_min = val-1.0
			val_max = val+1.0
			x_vals_iso = [item[1]
				for item in isobars_lists if val_min <= item[3] <= val_max]
			x_vals_iso.append(0)
			x_vals_iso = sorted(x_vals_iso)
			x_vals_iso = np.array(x_vals_iso)
			y_vals_iso = [item[2]
				for item in isobars_lists if val_min <= item[3] <= val_max]
			y_vals_iso.append(0)
			y_vals_iso = sorted(y_vals_iso)
			y_vals_iso = np.array(y_vals_iso)

			ax1.plot(x_vals_iso, y_vals_iso, ls='dashed', color='k')

		labels = P_vals
		ax1.legend(labels)

		return ax1

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
			in the fluid are all returned in a dict.

		Returns
		-------
		float or dict
			If verbose is set to False: Saturation pressure in bars.
			If verbose is set to True: dict of all calculated values.
		"""
		melts = self.melts
		bulk_comp_orig = sample

		bulk_comp = {oxide:  sample[oxide] for oxide in oxides}
		feasible = melts.set_bulk_composition(bulk_comp)

		#Coarse search
		fluid_mass = 0.0
		pressureMPa = 2000.0 #NOTE that pressure is in MPa for MagmaSat calculations but reported in bars.
		while fluid_mass <= 0.0:
			pressureMPa -= 100.0

			output = melts.equilibrate_tp(temperature, pressureMPa, initialize=True)
			(status, temperature, pressureMPa, xmlout) = output[0]
			fluid_mass = melts.get_mass_of_phase(xmlout, phase_name='Fluid')

			if pressureMPa <= 0:
				break

		startingP = pressureMPa+100.0 

		#Refined search 1
		feasible = melts.set_bulk_composition(bulk_comp)
		fluid_mass = 0.0
		pressureMPa = startingP
		while fluid_mass <= 0.0:
			pressureMPa -= 10.0

			output = melts.equilibrate_tp(temperature, pressureMPa, initialize=True)
			(status, temperature, pressureMPa, xmlout) = output[0]
			fluid_mass = melts.get_mass_of_phase(xmlout, phase_name='Fluid')

			if pressureMPa <= 0:
				break

		startingP_ref = pressureMPa +10.0

		#Refined search 2
		feasible = melts.set_bulk_composition(bulk_comp)

		fluid_mass = 0.0
		pressureMPa = startingP_ref
		while fluid_mass <= 0.0:
			pressureMPa -= 1.0

			output = melts.equilibrate_tp(temperature, pressureMPa, initialize=True)
			(status, temperature, pressureMPa, xmlout) = output[0]
			fluid_mass = melts.get_mass_of_phase(xmlout, phase_name='Fluid')

			if pressureMPa <= 0:
				break

		satP = pressureMPa*10 #convert pressure to bars
		flmass = fluid_mass
		flsystem_wtper = 100 * fluid_mass / (fluid_mass + melts.get_mass_of_phase(xmlout, phase_name='Liquid'))
		flcomp = melts.get_composition_of_phase(xmlout, phase_name='Fluid')
		flH2O = flcomp["H2O"]
		flCO2 = flcomp["CO2"]

		feasible = melts.set_bulk_composition(bulk_comp_orig) #this needs to be reset always!

		if verbose == False:
			return satP
		elif verbose == True:
			return {"SaturationP_bars": satP, "FluidMass_grams": flmass, "FluidProportion_wtper": flsystem_wtper,
					"H2Ofluid_wtper": flH2O, "CO2fluid_wtper": flCO2}

	def calculate_degassing_paths(self, sample, temperature, system='closed', init_vapor='None'):
		"""
		Calculates degassing path for one sample

		Parameters
		----------
		sample: dict
			Dictionary with values for sample composition as oxides in wt%. If pulling from an uploaded file
			with data for many samples, first call get_sample_oxide_comp() to get the sample desired. Then pass
			the result into this function.

		temp: float
			Temperature at which to calculate degassing paths, in degrees C.

		system: str
			OPTIONAL. Default value is 'closed'. Specifies the type of calculation performed, either closed system or closed
			system degassing. If closed is chosen, user can also specify the 'init_vapor' argument (see below).
			Possible inputs are 'open' and 'closed'.

		init_vapor: float
			OPTIONAL. Default value is 0.0. Specifies the amount of vapor (in wt%) coexisting with the melt before
			degassing.

		Returns
		-------
		pandas DataFrame object

		"""
		#--------------Preamble required for every MagmaSat method within Modeller---------------#
		# instantiate thermoengine equilibrate MELTS instance
		melts = equilibrate.MELTSmodel(self.model_version)

		# Suppress phases not required in the melts simulation
		self.oxides = melts.get_oxide_names()
		self.phases = melts.get_phase_names()

		for phase in self.phases:
			melts.set_phase_inclusion_status({phase: False})
		melts.set_phase_inclusion_status({'Fluid': True, 'Liquid': True})
		#---------------------------------------------------------------------------------------#

		sample = normalize(sample)
		bulk_comp_orig = sample

		feasible = melts.set_bulk_composition(sample)

		# Get saturation pressure
		data = self.calculate_saturation_pressure(sample, temperature)

		if system == 'closed':
			if init_vapor == 'None':
				P_array = np.arange(1.0, data['SaturationPressure_MPa']+10.0, 10)
				P_array = -np.sort(-P_array)
				output = melts.equilibrate_tp(temperature, P_array)

				pressure = []
				H2Oliq = []
				CO2liq = []
				H2Ofl = []
				CO2fl = []
				fluid_wtper = []
				for i in range(len(output)):
					(status, temperature, p, xmlout) = output[i]
					liq_comp = melts.get_composition_of_phase(xmlout, phase_name='Liquid')
					fl_comp = melts.get_composition_of_phase(xmlout, phase_name='Fluid')
					liq_mass = melts.get_mass_of_phase(xmlout, phase_name='Liquid')
					fl_mass = melts.get_mass_of_phase(xmlout, phase_name='Fluid')
					fl_wtper = 100 * fl_mass / (fl_mass+liq_mass)

					pressure.append(p)
					try:
						H2Oliq.append(liq_comp["H2O"])
					except:
						H2Oliq.append(0)
					try:
						CO2liq.append(liq_comp["CO2"])
					except:
						CO2liq.append(0)
					try:
						H2Ofl.append(fl_comp["H2O"])
					except:
						H2Ofl.append(0)
					try:
						CO2fl.append(fl_comp["CO2"])
					except:
						CO2fl.append(0)
					fluid_wtper.append(fl_wtper)

					try:
						sample["H2O"] = liq_comp["H2O"]
					except:
						sample["H2O"] = 0
					try:
						sample["CO2"] = liq_comp["CO2"]
					except:
						sample["CO2"] = 0
					fluid_wtper.append(fl_wtper)

				sample = bulk_comp_orig
				feasible = melts.set_bulk_composition(sample)
				closed_degassing_df = pd.DataFrame(list(zip(pressure, H2Oliq, CO2liq, H2Ofl, CO2fl, fluid_wtper)),
											columns =['pressure', 'H2Oliq', 'CO2liq', 'H2Ofl', 'CO2fl', 'fluid_wtper'])

				return closed_degassing_df

			else:
				P_array = np.arange(1.0, data['SaturationPressure_MPa'], 10)
				P_array = -np.sort(-P_array)
				fl_wtper = data["FluidSystem_wtper"]

				while fl_wtper <= init_vapor:
					output = melts.equilibrate_tp(temperature, data["SaturationPressure_MPa"])
					(status, temperature, p, xmlout) = output[0]
					fl_mass = melts.get_mass_of_phase(xmlout, phase_name='Fluid')
					liq_mass = melts.get_mass_of_phase(xmlout, phase_name='Liquid')
					fl_comp = melts.get_composition_of_phase(xmlout, phase_name='Fluid')
					fl_wtper = 100 * fl_mass / (fl_mass+liq_mass)
					sample["H2O"] += fl_comp["H2O"]*0.0005
					sample["CO2"] += fl_comp["CO2"]*0.0005
					sample = normalize(sample)
					feasible = melts.set_bulk_composition(sample)

				output = melts.equilibrate_tp(temperature, P_array)

				pressure = []
				H2Oliq = []
				CO2liq = []
				H2Ofl = []
				CO2fl = []
				fluid_wtper = []
				for i in range(len(output)):
					(status, temperature, p, xmlout) = output[i]
					liq_comp = melts.get_composition_of_phase(xmlout, phase_name='Liquid')
					fl_comp = melts.get_composition_of_phase(xmlout, phase_name='Fluid')
					liq_mass = melts.get_mass_of_phase(xmlout, phase_name='Liquid')
					fl_mass = melts.get_mass_of_phase(xmlout, phase_name='Fluid')
					fl_wtper = 100 * fl_mass / (fl_mass+liq_mass)

					pressure.append(p)
					try:
						H2Oliq.append(liq_comp["H2O"])
					except:
						H2Oliq.append(0)
					try:
						CO2liq.append(liq_comp["CO2"])
					except:
						CO2liq.append(0)
					try:
						H2Ofl.append(fl_comp["H2O"])
					except:
						H2Ofl.append(0)
					try:
						CO2fl.append(fl_comp["CO2"])
					except:
						CO2fl.append(0)
					fluid_wtper.append(fl_wtper)

					try:
						sample["H2O"] = liq_comp["H2O"]
					except:
						sample["H2O"] = 0
					try:
						sample["CO2"] = liq_comp["CO2"]
					except:
						sample["CO2"] = 0
					fluid_wtper.append(fl_wtper)

				sample = bulk_comp_orig
				feasible = melts.set_bulk_composition(sample)
				fl_wtper = data["FluidSystem_wtper"]
				exsolved_degassing_df = pd.DataFrame(list(zip(pressure, H2Oliq, CO2liq, H2Ofl, CO2fl, fluid_wtper)),
											columns =['pressure', 'H2Oliq', 'CO2liq', 'H2Ofl', 'CO2fl', 'fluid_wtper'])

				return exsolved_degassing_df

		elif system == 'open':
			P_array = np.arange(1.0, data['SaturationPressure_MPa']+10.0, 10)
			P_array = -np.sort(-P_array)

			pressure = []
			H2Oliq = []
			CO2liq = []
			H2Ofl = []
			CO2fl = []
			fluid_wtper = []
			for i in P_array:
				fl_mass = 0.0
				output = melts.equilibrate_tp(temperature, i)
				(status, temperature, p, xmlout) = output[0]
				liq_comp = melts.get_composition_of_phase(xmlout, phase_name='Liquid')
				fl_comp = melts.get_composition_of_phase(xmlout, phase_name='Fluid')
				liq_mass = melts.get_mass_of_phase(xmlout, phase_name='Liquid')
				fl_mass = melts.get_mass_of_phase(xmlout, phase_name='Fluid')
				fl_wtper = 100 * fl_mass / (fl_mass+liq_mass)

				if fl_mass > 0:
					pressure.append(p)
					try:
						H2Oliq.append(liq_comp["H2O"])
					except:
						H2Oliq.append(0)
					try:
						CO2liq.append(liq_comp["CO2"])
					except:
						CO2liq.append(0)
					try:
						H2Ofl.append(fl_comp["H2O"])
					except:
						H2Ofl.append(0)
					try:
						CO2fl.append(fl_comp["CO2"])
					except:
						CO2fl.append(0)
					fluid_wtper.append(fl_wtper)

					try:
						sample["H2O"] = liq_comp["H2O"]
					except:
						sample["H2O"] = 0
					try:
						sample["CO2"] = liq_comp["CO2"]
					except:
						sample["CO2"] = 0
					sample = normalize(sample)

					feasible = melts.set_bulk_composition(sample)

			sample = bulk_comp_orig
			feasible = melts.set_bulk_composition(sample) #this needs to be reset always!
			open_degassing_df = pd.DataFrame(list(zip(pressure, H2Oliq, CO2liq, H2Ofl, CO2fl, fluid_wtper)),
										columns =['pressure', 'H2Oliq', 'CO2liq', 'H2Ofl', 'CO2fl', 'fluid_wtper'])

			return open_degassing_df

	


#====== Define some standard model options =======================================================#

default_models = {'Shishkina':				MixedFluids({'CO2':ShishkinaCarbon(),'H2O':ShishkinaWater()}),
				  'Dixon':					MixedFluids({'CO2':DixonCarbon(),'H2O':DixonWater()}),
				  'IaconoMarziano': 		MixedFluids({'CO2':IaconoMarzianoCarbon(),'H2O':IaconoMarzianoWater()}),
				  'ShishkinaCarbon':		ShishkinaCarbon(),
				  'ShishkinaWater': 		ShishkinaWater(),
				  'DixonCarbon':			DixonCarbon(),
				  'DixonWater':				DixonWater(),
				  'IaconoMarzianoCarbon':	IaconoMarzianoCarbon(),
				  'IaconoMarzianoWater':	IaconoMarzianoWater(),
				  'EguchiCarbon':			EguchiCarbon(),
				  'AllisonCarbon':			AllisonCarbon(),
				  'MagmaSat':				MagmaSat()
}

class calculate_dissolved_volatiles(Calculate):
	""" This will be used as the user interface to doing the calculation. This provides a generic
	way of accessing the relevant functions, and here will be when calibration checking can be done.
	"""
	def calculate(self,sample,pressure,**kwargs):
		dissolved = self.model.calculate_dissolved_volatiles(pressure=pressure,sample=sample,**kwargs)
		return dissolved

	def check_calibration_range(self,sample,pressure,**kwargs):
		calib_check = self.model.check_calibration_range({'pressure':pressure,
														'dissolved_volatiles':self.result,
														'sample':sample})
		return calib_check


class calculate_equilibrium_fluid_comp(Calculate):
	""" This will be used as the user interface to doing the calculation. This provides a generic
	way of accessing the relevant functions, and here will be when calibration checking can be done.
	"""
	def calculate(self,sample,pressure,**kwargs):
		fluid_comp = self.model.calculate_equilibrium_fluid_comp(pressure=pressure,sample=sample,**kwargs)
		return fluid_comp
	def check_calibration_range(self,pressure,**kwargs):
		return 0


class calculate_isobars_and_isopleths(Calculate):
	"""
	"""
	def calculate(self,sample,pressure_list,isopleth_list=[0,1],points=101,**kwargs):
		check = getattr(self.model, "calculate_isobars_and_isopleths", None)
		if callable(check):
			isobars, isopleths = self.model.calculate_isobars_and_isopleths(sample=sample,pressure_list=pressure_list,isopleth_list=isopleth_list,points=points,**kwargs)
			return isobars, isopleths
		else:
			print("Write code to report error!!!!!")

	def check_calibration_range(self,**kwargs):
		return 0


class calculate_saturation_pressure(Calculate):
	"""
	"""
	def calculate(self,sample,**kwargs):
		satP = self.model.calculate_saturation_pressure(sample=sample,**kwargs)
		return satP

	def check_calibration_range(self,**kwargs):
		return 0

class calculate_degassing_paths(Calculate):
	"""
	"""
	def calculate(self,pressures,sample,fractionate_vapor=1.0,**kwargs):
		check = getattr(self.model, "calculate_degassing_paths", None)
		if callable(check):
			melt, fluid = self.model.calculate_degassing_paths(sample=sample,pressures=pressures,
															fractionate_vapor=fractionate_vapor,**kwargs)
			return [melt,fluid]
		else:
			print("Write code to report error!!!!!")

	def check_calibration_range(self,**kwargs):
		return 0




# class fugacity_ZD09_co2(FugacityModel):
# 	def __init__(self):
# 		self.a = np.array([0.0,
# 	                  	   2.95177298930e-2,
# 	         			   -6.33756452413e3,
# 						   -2.75265428882e5,
# 						   1.29128089283e-3,
# 						   -1.45797416153e2,
# 						   7.65938947237e4,
# 						   2.58661493537e-6,
# 						   0.52126532146,
# 						   -1.39839523753e2,
# 						   -2.36335007175e-8,
# 						   5.35026383543e-3,
# 						   -0.27110649951,
# 						   2.50387836486e4,
# 						   0.73226726041,
# 						   1.5483335997e-2])
# 		self.e = {'H2O':510.0,'CO2':235.0}
# 		self.s = {'H2O':2.88,'CO2':3.79}
# 		self.k1 = 0.85
# 		self.k2 = 1.02
#
# 	def fugacity(self,pressure,temperature,X_fluid,**kwargs):
# 		phi = np.exp(self.logPhi(P=pressure,T=temperature,x=X_fluid))
# 		return pressure*X_fluid*phi
#
# 	def volume(self,P,T,x):
# 		""" write stuff"""
# 		return root_scalar(self.root_volumez,args=(P,T,x),x0=20,method='newton',
# 							fprime=self.root_dpdvm).root
#
# 	def root_volume(self,v,P,T,x):
# 		e = x**2*self.e['CO2'] + (1-x)**2*self.e['H2O'] + 2*self.k1*x*(1-x)*(self.e['H2O']*self.e['CO2'])**0.5
# 		s = x**2*self.s['CO2'] + (1-x)**2*self.s['H2O'] + self.k2*x*(1-x)*(self.s['CO2']+self.s['H2O'])
# 		a = self.a
#
# 		vm = v/1000*(3.691/s)**3
#
# 		Pm = 3.0636*P*s**3/e
# 		Tm = 154*T/e
#
# 		Z = 1 + (a[1]+a[2]/Tm**2 + a[3]/Tm**3)/vm + (a[4]+a[5]/Tm**2+a[6]/Tm**3)/vm**2 \
# 			+ (a[7]+a[8]/Tm**2+a[9]/Tm**3)/vm**4 + (a[10]+a[11]/Tm**2+a[12]/Tm**3)/vm**5 \
# 			+ a[13]/(Tm**3*vm**2)*(a[14]+a[15]/vm**2)*np.exp(-a[15]/vm**2)
#
# 		p = Z*0.08314*Tm/vm
#
# 		return p - P
#
# 	def root_volumez(self,v,P,T,x):
# 		e = x**2*self.e['CO2'] + (1-x)**2*self.e['H2O'] + 2*self.k1*x*(1-x)*(self.e['H2O']*self.e['CO2'])**0.5
# 		s = x**2*self.s['CO2'] + (1-x)**2*self.s['H2O'] + self.k2*x*(1-x)*(self.s['CO2']+self.s['H2O'])
# 		a = self.a
#
# 		vm = v/1000*(3.691/s)**3
#
# 		Pm = 3.0636*P*s**3/e
# 		Tm = 154*T/e
#
# 		Z = 1 + (a[1]+a[2]/Tm**2 + a[3]/Tm**3)/vm + (a[4]+a[5]/Tm**2+a[6]/Tm**3)/vm**2 \
# 			+ (a[7]+a[8]/Tm**2+a[9]/Tm**3)/vm**4 + (a[10]+a[11]/Tm**2+a[12]/Tm**3)/vm**5 \
# 			+ a[13]/(Tm**3*vm**2)*(a[14]+a[15]/vm**2)*np.exp(-a[15]/vm**2)
#
#
#
# 		return Z - Pm*vm/0.08314/Tm
#
# 	def volume_est(self,P,T,x):
#
# 		delv = 1.0
# 		vPrevious = 1.0
# 		delvPrevious = 1.0
# 		vLowest = 0.017
#
# 		e = x**2*self.e['CO2'] + (1-x)**2*self.e['H2O'] + 2*self.k1*x*(1-x)*(self.e['H2O']*self.e['CO2'])**0.5
# 		s = x**2*self.s['CO2'] + (1-x)**2*self.s['H2O'] + self.k2*x*(1-x)*(self.s['CO2']+self.s['H2O'])
# 		a = self.a
#
# 		v = 0.08314467*T/P
# 		Pm = 3.063*s**3*P/e
# 		Tm = 154*T/e
# 		vm = v*(3.691/s)**3
#
# 		z = self.root_z(vm,T,x)
#
# 		iter = 0
# 		while (z < 0.0) and (iter < 200):
# 			dzdvm = self.root_dzdm(vm,T,x)
# 			vm += (0.1-z)/dzdvm
# 			z = self.root_z(vm,T,x)
# 			iter += 1
#
# 		iter = 0
# 		while iter < 200:
# 			delv = z*0.08314467*Tm/Pm - vm
# 			if ((iter > 1) and (delv*delvprevious < 0.0)) or (np.abs(delv)<vm*1000*e):
# 				break
# 			vPrevious = vm
# 			delvPrevious = delv
# 			vm = (z*0.08314467*Tm/Pm + vm)/2.0
# 			increment = vm - vPrevious
#
# 			while self.root_z(vm,T,x) < 0.0:
# 				increment = increment/2
# 				if np.abs(increment) < 100*e:
# 					return False
# 				vm = vPrevious + increment
# 			if vm < vLowest:
# 				vm = vLowest
# 			iter += 1
# 			if iter == 200:
# 				if np.abs(vm-vLowest) <= 100*e:
# 					vLowest = 0.9*vLowest
# 					iter = 0
#
# 		if np.abs(delv) > vm*100*e:
# 			if delv < 0.0:
# 				dx = vPrevious - vm
# 				rtb = vm
# 			else:
# 				dx = vm-vPrevious
# 				rtb = vPrevious
# 			iter = 0
# 			while iter<200:
# 				dx = dx*0.5
# 				vm = rtb + dx
# 				z = self.root_z(vm,T,x)
# 				delv = z*0.08314467*Tm/Pm - vm
# 				if delv <= 0.0:
# 					rtb = vm
# 				if (np.abs(delv)<100*e) or delv == 0.0:
# 					break
# 				iter += 1
# 			if (iter == 200) or (np.abs(dx) > 100.0*e):
# 				return False
# 		elif iter == 200:
# 			return False
#
# 		return 1000.0*vm*(s/3.691)**3
#
#
#
# 	def root_z(self,vm,T,x):
# 		a = self.a
# 		e = x**2*self.e['CO2'] + (1-x)**2*self.e['H2O'] + 2*self.k1*x*(1-x)*(self.e['H2O']*self.e['CO2'])**0.5
# 		Tm = 154*T/e
#
# 		z = 1 + (a[1]+a[2]/Tm**2 + a[3]/Tm**3)/vm + (a[4]+a[5]/Tm**2+a[6]/Tm**3)/vm**2 \
# 			+ (a[7]+a[8]/Tm**2+a[9]/Tm**3)/vm**4 + (a[10]+a[11]/Tm**2+a[12]/Tm**3)/vm**5 \
# 			+ a[13]/(Tm**3*vm**2)*(a[14]+a[15]/vm**2)*np.exp(-a[15]/vm**2)
# 		return z
#
# 	def root_dzdm(self,vm,T,x):
# 		a = self.a
# 		e = x**2*self.e['CO2'] + (1-x)**2*self.e['H2O'] + 2*self.k1*x*(1-x)*(self.e['H2O']*self.e['CO2'])**0.5
# 		Tm = 154*T/e
#
# 		bv = (a[1]+a[2]/Tm**2 + a[3]/Tm**3)
# 		cv = (a[4]+a[5]/Tm**2+a[6]/Tm**3)
# 		dv = (a[7]+a[8]/Tm**2+a[9]/Tm**3)
# 		ev = (a[10]+a[11]/Tm**2+a[12]/Tm**3)
# 		fv = a[13]*a[14]/Tm**3
# 		gv = a[13]*a[15]/Tm**3
# 		gammav = a[15]
#
# 		dzdvm = - bv/vm**2 - 2.0*cv/vm**3+ - 4.0*dv/vm**5 - 5.0*ev/vm**6 \
# 				- 2.0*(fv/vm**3) * np.exp(-gammav/vm**2) + 2.0*(fv/vm**2) * (gammav/vm**3) * np.exp(-gammav/vm**2) \
# 				- 4.0*(gv/vm**5) * np.exp(-gammav/vm**2) + 2.0*(gv/vm**4) * (gammav/vm**3) * np.exp(-gammav/vm**2)
#
# 		return dzdvm
#
# 	def root_dpdvm(self,vm,P,T,x):
# 		a = self.a
# 		e = x**2*self.e['CO2'] + (1-x)**2*self.e['H2O'] + 2*self.k1*x*(1-x)*(self.e['H2O']*self.e['CO2'])**0.5
# 		s = x**2*self.s['CO2'] + (1-x)**2*self.s['H2O'] + self.k2*x*(1-x)*(self.s['CO2']+self.s['H2O'])
# 		Tm = 154*T/e
# 		Pm = 3.063*s**3*P/e
#
# 		bv = (a[1]+a[2]/Tm**2 + a[3]/Tm**3)
# 		cv = (a[4]+a[5]/Tm**2+a[6]/Tm**3)
# 		dv = (a[7]+a[8]/Tm**2+a[9]/Tm**3)
# 		ev = (a[10]+a[11]/Tm**2+a[12]/Tm**3)
# 		fv = a[13]*a[14]/Tm**3
# 		gv = a[13]*a[15]/Tm**3
# 		gammav = a[15]
#
# 		dzdvm = - bv/vm**2 - 2.0*cv/vm**3+ - 4.0*dv/vm**5 - 5.0*ev/vm**6 \
# 				- 2.0*(fv/vm**3) * np.exp(-gammav/vm**2) + 2.0*(fv/vm**2) * (gammav/vm**3) * np.exp(-gammav/vm**2) \
# 				- 4.0*(gv/vm**5) * np.exp(-gammav/vm**2) + 2.0*(gv/vm**4) * (gammav/vm**3) * np.exp(-gammav/vm**2) \
# 				- Pm/0.08314/Tm
#
# 		return dzdvm
#
#
#
# 	def logPhi(self,P,T,x):
# 		e = x**2*self.e['CO2'] + (1-x)**2*self.e['H2O'] + 2*self.k1*x*(1-x)*(self.e['H2O']*self.e['CO2'])**0.5
# 		s = x**2*self.s['CO2'] + (1-x)**2*self.s['H2O'] + self.k2*x*(1-x)*(self.s['CO2']+self.s['H2O'])
# 		a = self.a
#
# 		v = self.volume(P,T,x)
#
# 		# vm = v/1000*(3.691/s)**3
# 		vm = v*(3.691/s)**3
#
#
# 		Pm = 3.0636*P*s**3/e
# 		Tm = 154*T/e
#
# 		Z = 1 + (a[1]+a[2]/Tm**2 + a[3]/Tm**3)/vm + (a[4]+a[5]/Tm**2+a[6]/Tm**3)/vm**2 \
# 			+ (a[7]+a[8]/Tm**2+a[9]/Tm**3)/vm**4 + (a[10]+a[11]/Tm**2+a[12]/Tm**3)/vm**5 \
# 			+ a[13]/(Tm**3*vm**2)*(a[14]+a[15]/vm**2)*np.exp(-a[15]/vm**2)
#
# 		S1 = (a[1] + a[2]/Tm**2 + a[3]/Tm**3)/vm + (a[4] + a[5]/Tm**2 + a[6]/Tm**3)/(2*vm**2) \
# 			+ (a[7] + a[8]/Tm**2 + a[9]/Tm**3)/(4*vm**4) + (a[10] + a[11]/Tm**2 + a[12]/Tm**3)/(5*vm**5) \
# 			+ a[13]/(2*a[15]*Tm**3)*(a[14]+1-(a[14]+1+a[15]/vm**2)*np.exp(-a[15]/vm**2))
#
# 		S2 = (2*a[2]/Tm**2 + 3*a[3]/Tm**3)/vm + (2*a[5]/Tm**2 + 3*a[6]/Tm**3)/(2*vm**2) \
# 			+ (2*a[8]/Tm**2 + 3*a[9]/Tm**3)/(4*vm**4) + (2*a[11]/Tm**2 + 3*a[12]/Tm**3)/(5*vm**5) \
# 			+ 3*a[13]/(2*a[15]*Tm**3)*(a[14]+1-(a[14]+1+a[15]/vm**2)*np.exp(-a[15]/vm**2))
#
# 		logPhi = Z - 1 - np.log(Z) + S1 + 2*S2* \
# 			(1-(self.k1*self.e['CO2'] + self.k1*(1-x)*(self.e['CO2']*self.e['H2O'])**0.5)) \
# 			+ 6*(1-Z)*(1-(self.k2*x*self.s['CO2'] + self.k2*(1-x)*(self.s['CO2']+self.s['H2O'])*0.5))
#
# 		return logPhi
