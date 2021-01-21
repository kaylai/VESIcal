# ---------- DEFINE SOME CONSTANTS ------------- #
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
