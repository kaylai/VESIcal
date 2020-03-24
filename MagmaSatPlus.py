# Python 3.5
# Script written by Kayla Iacovino (kayla.iacovino@nasa.gov)
# VERSION 0.1- MARCH 2020

import pandas as pd

#----------DEFINE SOME CONSTANTS-------------#
oxides = ['SiO2', 'TiO2', 'Al2O3', 'Fe2O3', 'Cr2O3', 'FeO', 'MnO', 'MgO', 'NiO', 'CoO', 'CaO', 'Na2O', 'K2O', 'P2O5', 
		  'H2O', 'CO2']
oxideMass = {'SiO2':28.085+32,'MgO':24.305+16,'FeO':55.845+16,'CaO':40.078+16,'Al2O3':2*26.982+16*3,'Na2O':22.99*2+16,
             'K2O':39.098*2+16,'MnO':54.938+16,'TiO2':47.867+32,'P2O5':2*30.974+5*16,'Cr2O3':51.996*2+3*16,
             'NiO':58.693+16, 'CoO':28.01+16, 'Fe2O3':55.845*2+16*3,
             'H2O':18.02, 'CO2':44.01}
CationNum = {'SiO2':1,'MgO':1,'FeO':1,'CaO':1,'Al2O3':2,'Na2O':2,
             'K2O':2,'MnO':1,'TiO2':1,'P2O5':2,'Cr2O3':2,
             'NiO':1,'CoO':1,'Fe2O3':2,'H2O':2, 'CO2':1}
CationCharge = {'SiO2':4,'MgO':2,'FeO':2,'CaO':2,'Al2O3':3,'Na2O':1,
             'K2O':1,'MnO':2,'TiO2':4,'P2O5':5,'Cr2O3':3,
             'NiO':2,'CoO':2,'Fe2O3':3,'H2O':1, 'CO2':4}
CationMass = {'SiO2':28.085,'MgO':24.305,'FeO':55.845,'CaO':40.078,'Al2O3':26.982,'Na2O':22.990,
             'K2O':39.098,'MnO':54.938,'TiO2':47.867,'P2O5':30.974,'Cr2O3':51.996,
             'NiO':58.693,'CoO':28.01,'Fe2O3':55.845,'H2O':2, 'CO2':12.01}

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

    def __init__(self, expression, message):
        self.expression = expression
        self.message = message

#----------DEFINE SOME BASIC METHODS-----------#
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
		data.loc[:,key] *= value 

	data["MPOSum"] = sum([data[oxide] for oxide in oxides])

	for oxide in oxides:
		data.loc[:,oxide] /= data['MPOSum']
		data.loc[:,oxide] *= 100
	del data['MPOSum']

	return data


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
			raise InputError("Imported file must contain a column of sample names with the column name \'Label\'")

		for oxide in oxides:
			if oxide in data.columns:
				pass
			else:
				data[oxide] = 0.0

		#TODO test all input types produce correct values
		if input_type == "wtpercent":
			pass

		if input_type == "molpercent":
			data = mol_to_wtpercent(data)

		if input_type == "molfrac":
			data = mol_to_wtpercent(data)

		self.data = data

	def get_sample_oxide_comp(self, sample='default'):
		"""
		Returns oxide composition of a single sample from a user-imported excel file as a dictionary

		Parameters
		----------
		sample: string
			Name of the desired sample

		Returns
		-------
		dictionary
			Composition of the sample as oxides
		"""
		sample = sample
		data = self.data
		my_sample = pd.DataFrame(data.loc[sample])
		sample_dict = (my_sample.to_dict()[sample])
		sample_oxides = {}
		for item, value in sample_dict.items():
			if item in oxides:
				sample_oxides.update({item:value})
		return sample_oxides

#------------DEFINE SOME PLOTTING METHODS-----------#
def plot_isobars_and_isopleths(sample, temp, pressure_min=None, pressure_max=None, pressure_int=None, pressure_list=None):
	"""
	Plots isobars and isopleths at a constant temperature for a given sample. Isobars can be calculated
	for any number of pressures. Pressures can be passed as min, max, interval (100.0, 500.0, 100.0 would result
	in pressures of 100.0, 200.0, 300.0, 400.0, and 500.0 MPa). Alternatively pressures can be passed as a list of all 
	desired pressures ([100.0, 200.0, 250.0, 300.0] would calculate isobars for each of those pressures in MPa).

	Parameters
	----------
	sample: dict
		Dictionary with values for sample composition as oxides in wt%.

	temp: float
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

	Returns
	-------
	matplotlib object
		Plot with x-axis as H2O wt% in the melt and y-axis as CO2 wt% in the melt. Isobars, or lines of
		constant pressure at which the sample magma composition is saturated, and isopleths, or lines of constant
		fluid composition at which the sample magma composition is saturated, are plotted.
	"""
	if pressure_min is not None or pressure_max is not None or pressure_int is not None and pressure_list is not None:
		raise InputError('Enter pressure either as min, max, int OR as list. Not both.')

	if pressure_min is not None:
		P_vals = np.arange(pressure_min, pressure_max+pressure_int, pressure_int)

	if pressure_list is not None:
		P_vals = pressure_list
















