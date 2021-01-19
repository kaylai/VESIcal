
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

		self.composition = composition

		self.set_default_normalization(default_normalization)
		self.set_default_type(default_type)
		self.set_default_units(default_units)


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

	def set_default_units(self, default_units):
		""" Set the default units for reporting the sample composition when using the get_composition() method.

		Parameters
		----------
		default_units 	str
			The units in which the composition should be returned, one of:
			- wtpt (weight percent, default)
			- mols (mole fraction, summing to 1)
		"""
		self.default_units = default_units

	def get_composition(self,normalization=None,type=None,units=None):
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
			- oxides (default)
			- cations (used only with mole fractions)
			If NoneType is passed the default type option will be used (self.default_type).

		units: 	str
			The units in which the composition should be returned, one of:
			- wtpt (weight percent, default)
			- mols (mole fraction, summing to 1)
			If NoneType is passed the default units will be used (self.default_units).

		Returns
		-------
		dict
			The sample composition, as specified.
		"""

		if normalization == None:
			normalization = self.default_normalization

		if type == None:
			type = self.default_type

		if units == None:
			units = self.default_units

		if normalization != None or type != 'oxides' or units != 'wtpt':
			w.warn("The method presently only returns wtpt oxides.",RuntimeWarning,stacklevel=2)

		return self.composition

	def get_formulaweight(self):
		""""""
		return np.nan

	def normalize_Standard(self, composition):
		""""""
		return composition

	def normalize_FixedVolatiles(self, composition):
		""""""
		return composition

	def normalize_AdditionalVolatiles(self, composition):
		""""""
		return composition

	def wtpercentOxides_to_molOxides(self, composition):
		""""""
		return composition

	def wtpercentOxides_to_molCations(self, composition):
		""""""
		return composition

	def molOxides_to_wtpercentOxides(self, composition):
		""""""
		return composition

	def molOxides_to_molCations(self, composition):
		""""""
		return composition

	def molCations_to_wtpercentOxides(self, composition):
		""""""
		return composition

	def molCations_to_molOxides(self, composition):
		""""""
		return composition
