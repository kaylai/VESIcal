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

class sample(object):
	""" WORK IN PROGRESS.
	The sample class stores compositional information for samples, and contains methods
	for normalization and other compositional calculations.
	"""

	def __init__(self, composition, type='oxide_wtpt', default_normalization='standard', default_type='oxides',
				 default_units='wtpt'):
		""" Initialises the sample class.

		The composition is stored as wtpt. If the composition
		is provided as wtpt, no normalization will be applied. If the composition is supplied as
		mols, the composition will be normalized to 100 wt%.

		Parameters
		----------
		composition 	dict or pandas.Series
			The composition of the sample in the format specified by the composition_type
			parameter. Defulat is oxides in wtpt.

		type 	str
			Specifies the units and type of compositional information passed in the
			composition parameter. Choose from 'oxide_wtpt', 'oxide_mols', 'cation_mols'.

		default_normalization: 	None or str
			The type of normalization to apply to the data by default. One of:
				- None (no normalization)
				- 'standard' (default): Normalizes an input composition to 100%.
				- 'fixedvolatiles': Normalizes major element oxides to 100 wt%, including volatiles.
				The volatile wt% will remain fixed, whilst the other major element oxides are reduced
				proportionally so that the total is 100 wt%.
				- 'additionalvolatiles': Normalises major element oxide wt% to 100%, assuming it is
				volatile-free. If H2O or CO2 are passed to the function, their un-normalized values will
				be retained in addition to the normalized non-volatile oxides, summing to >100%.

		default_type 	str
			The type of composition to return by default, one of:
			- oxides (default)
			- cations (used only with mole fractions)

		default_units 	str
			The units in which the composition should be returned by default, one of:
			- wtpt (weight percent, default)
			- mols (mole fraction, summing to 1)
		"""
		if type != 'oxide_wtpt':
			w.warn("Presently the sample class can only be initialised with oxides in wtpt.",RuntimeWarning,stacklevel=2)

		self._composition = composition

		self.set_default_normalization(default_normalization)
		self.set_default_type(default_type)


	def set_default_normalization(self, default_normalization):
		""" Set the default type of normalization to use with the get_composition() method.

		Parameters
		----------
		default_normalization:	str
			The type of normalization to apply to the data. One of:
				- 'none' (no normalization)
				- 'standard' (default): Normalizes an input composition to 100%.
				- 'fixedvolatiles': Normalizes major element oxides to 100 wt%, including volatiles.
				The volatile wt% will remain fixed, whilst the other major element oxides are reduced
				proportionally so that the total is 100 wt%.
				- 'additionalvolatiles': Normalises major element oxide wt% to 100%, assuming it is
				volatile-free. If H2O or CO2 are passed to the function, their un-normalized values will
				be retained in addition to the normalized non-volatile oxides, summing to >100%.
		"""
		self.default_normalization = default_normalization

	def set_default_type(self, default_type):
		""" Set the default type of composition to return when using the get_composition() method.

		Parameters
		----------
		default_type 	str
			The type of composition to return, one of:
			- oxides (default)
			- cations (used only with mole fractions)
		"""
		self.default_type = default_type


	def get_composition(self,normalization=None,type=None):
		""" Returns the composition in the format requested, normalized as requested.

		Parameters
		----------
		normalization: 	NoneType or str
			The type of normalization to apply to the data. One of:
				- 'none' (no normalization)
				- 'standard' (default): Normalizes an input composition to 100%.
				- 'fixedvolatiles': Normalizes major element oxides to 100 wt%, including volatiles.
				The volatile wt% will remain fixed, whilst the other major element oxides are reduced
				proportionally so that the total is 100 wt%.
				- 'additionalvolatiles': Normalises major element oxide wt% to 100%, assuming it is
				volatile-free. If H2O or CO2 are passed to the function, their un-normalized values will
				be retained in addition to the normalized non-volatile oxides, summing to >100%.
			If NoneType is passed the default normalization option will be used (self.default_normalization).

		type: 	NoneType or str
			The type of composition to return, one of:
			- wtpt_oxides (default)
			- mol_oxides
			- mol_cations
			- mol_singleO
			If NoneType is passed the default type option will be used (self.default_type).

		Returns
		-------
		dict
			The sample composition, as specified.
		"""

		# Fetch the default return types if not specified in function call
		if normalization == None:
			normalization = self.default_normalization
		if type == None:
			type = self.default_type

		# Do requested normalization
		if normalization == 'none':
			normed = self._composition
		elif normalization == 'standard':
			normed = self._normalize_Standard(self._composition)
		elif normalization == 'fixedvolatiles':
			normed = self._normalize_FixedVolatiles(self._composition)
		elif normalization == 'additionalvolatiles':
			normed = self._normalize_AdditionalVolatiles(self._composition)
		else:
			raise InputError("The normalization method must be one of 'none', 'standard', 'fixedvolatiles',\
			 or 'additionalvolatiles'.")

		# Get the requested type of composition
		if type == 'wtpt_oxides':
			return normed
		elif type == 'mol_oxides':
			return self._wtpercentOxides_to_molOxides(self._composition)
		elif type == 'mol_cations':
			return self._wtpercentOxides_to_molCations(self._composition)
		elif type == 'mol_singleO':
			return self._wtpercentOxides_to_molSingleO(self._composition)
		else:
			raise InputError("The type must be one of 'wtpt_oxides', 'mol_oxides', 'mol_cations', \
			or 'mol_singleO'.")


	def get_formulaweight(self):
		""" Converts major element oxides in wt% to the formula weight (on a 1 oxygen basis).

		Returns
		-------
		float
			The formula weight of the composition, on a one oxygen basis.
		"""

		cations = self.get_composition(type='mol_singleO')

		if type(cations) != dict:
			cations = dict(cations)

		FW = 15.999
		for cation in list(cations.keys()):
			FW += cations[cation]*CationMass[cations_to_oxides[cation]]

		return FW

	def _normalize_Standard(self, composition):
		"""
		Normalizes the given composition to 100 wt%, including volatiles. This method
		is intended only to be called by the get_composition() method.

		Parameters
		----------
		composition: 	pandas.Series
			A rock composition with oxide names as keys and wt% concentrations as values.

		Returns
		-------
		pandas.Series
			Normalized oxides in wt%.
		"""
		comp = composition.copy()
		return {k: 100.0 * v / sum(comp.values()) for k, v in comp.items()}

	def _normalize_FixedVolatiles(self, composition):
		"""
		Normalizes major element oxides to 100 wt%, including volatiles. The volatile
		wt% will remain fixed, whilst the other major element oxides are reduced proportionally
		so that the total is 100 wt%.

		Intended to be called only by the get_composition() method.

		Parameters
		----------
		composition:    pandas Series
			Major element oxides in wt%

		Returns
		-------
		pandas Series
			Normalized major element oxides.
		"""
		comp = composition.copy()
		normalized = pd.Series({},dtype=float)
		volatiles = 0
		if 'CO2' in list(comp.index):
			volatiles += comp['CO2']
		if 'H2O' in list(comp.index):
			volatiles += comp['H2O']

		for ox in list(comp.index):
			if ox != 'H2O' and ox != 'CO2':
				normalized[ox] = comp[ox]

		normalized = normalized/np.sum(normalized)*(100-volatiles)

		if 'CO2' in list(comp.index):
			normalized['CO2'] = comp['CO2']
		if 'H2O' in list(sample.index):
			normalized['H2O'] = comp['H2O']

		return normalized

	def _normalize_AdditionalVolatiles(self, composition):
		"""
		Normalises major element oxide wt% to 100%, assuming it is volatile-free. If
		H2O or CO2 are passed to the function, their un-normalized values will be retained
		in addition to the normalized non-volatile oxides, summing to >100%.

		Intended to be called only by the get_composition() method.

		Parameters
		----------
		sample:    pandas.Series
			Major element oxides in wt%

		Returns
		-------
		pandas.Series
			Normalized major element oxides.
		"""
		comp = composition.copy()
		normalized = pd.Series({}, dtype=float)
		for ox in list(comp.index):
			if ox != 'H2O' and ox != 'CO2':
				normalized[ox] = comp[ox]

		normalized = normalized/np.sum(normalized)*100
		if 'H2O' in comp.index:
			normalized['H2O'] = comp['H2O']
		if 'CO2' in comp.index:
			normalized['CO2'] = comp['CO2']

		return normalized

	def _wtpercentOxides_to_molOxides(self, composition):
		"""
		Converts a wt% oxide composition to mol oxides, normalised to 1 mol.

		Intended to be called only by the get_composition() method.

		Parameters
		----------
		composition:	pandas.Series
			Major element oxides in wt%

		Returns
		-------
		pandas.Series
			Molar proportions of major element oxides, normalised to 1.
		"""
		molOxides = {}
		comp = composition.copy()
		oxideslist = list(comp.index)

		for ox in oxideslist:
			molOxides[ox] = comp[ox]/oxideMass[ox]

		molOxides = pd.Series(molOxides)
		molOxides = molOxides/molOxides.sum()

		return molOxides

	def _wtpercentOxides_to_molCations(self, composition):
		"""
		Converts a wt% oxide composition to molar proportions of cations (normalised to 1).

		Intended to be called only by the get_composition() method.

		Parameters
		----------
		composition		pandas.Series
			Major element oxides in wt%.

		Returns
		-------
		pandas.Series
			Molar proportions of cations, normalised to 1.
		"""
		molCations = {}
		comp = composition.copy()
		oxideslist = list(comp.index)

		for ox in oxideslist:
			cation = oxides_to_cations[ox]
			molCations[cation] = CationNum[ox]*comp[ox]/oxideMass[ox]

		molCations = pd.Series(molCations)
		molCations = molCations/molCations.sum()

		return molCations

	def _wtpercentOxides_to_molSingleO(self, composition):
		"""
		Constructs the chemical formula, on a single oxygen basis, from wt% oxides.

		Intended to be called only by the get_composition() method.

		Parameters
		----------
		composition		pandas.Series
			Major element oxides in wt%

		Returns
		-------
		pandas.Series
			The chemical formula of the composition, on a single oxygen basis. Each element is
			a separate entry in the Series.
		"""
		molCations = {}
		comp = composition.copy()

		oxideslist = list(comp.index)

		total_O = 0.0
		for ox in oxideslist:
			cation = oxides_to_cations[ox]
			molCations[cation] = CationNum[ox]*comp[ox]/oxideMass[ox]
			total_O += OxygenNum[ox]*comp[ox]/oxideMass[ox]

		molCations = pd.Series(molCations)
		molCations = molCations/total_O

		return molCations

	def _molOxides_to_wtpercentOxides(self, composition):
		""""""
		return composition

	def _molOxides_to_molCations(self, composition):
		""""""
		return composition

	def _molCations_to_wtpercentOxides(self, composition):
		""""""
		return composition

	def _molCations_to_molOxides(self, composition):
		""""""
		return composition
