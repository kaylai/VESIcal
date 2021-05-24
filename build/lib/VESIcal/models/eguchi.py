# Code retained in case it is useful in the future. Note that the workflow has not been maintained to be
# consistent with the rest of VESIcal. Trying to use this without updates (e.g. employing the new Sample
# class) will result in errors.
# SM Feb 2021.

# class EguchiCarbon(Model):
# 	"""
# 	Implementation of the Eguchi and Dasgupta (2018) CO2 solubility model for andesitic melts.
# 	Uses the Zhang and Duan (2009) CO2 EOS for fugacity calculations, assuming a pure CO2 fluid,
# 	or ideal mixing for mixed fluids.
# 	"""
#
# 	def __init__(self):
# 		w.warn("Eguchi model is not working correctly. Do not use any results calculated.")
# 		self.set_volatile_species(['CO2'])
# 		self.set_fugacity_model(fugacity_ZD09_co2())
# 		self.set_activity_model(activity_idealsolution())
# 		self.set_solubility_dependence(False)
# 		self.set_calibration_ranges([CalibrationRange('pressure',[500.0,50000.0],crf_Between,'bar','Eguchi & Dasgupta (2018) carbon',
# 													  fail_msg=crmsg_Between_fail, pass_msg=crmsg_Between_pass, description_msg=crmsg_Between_description),
# 									 CalibrationRange('temperature',[950.0,1600],crf_Between,'oC','Eguchi & Dasgupta (2018) carbon',
# 									 				  fail_msg=crmsg_Between_fail, pass_msg=crmsg_Between_pass, description_msg=crmsg_Between_description)])
#
# 	def preprocess_sample(self,sample,ferric_total=0.15):
# 		""" Returns normalized sample composition, with ferric iron. Where a sample
# 		already contains ferric iron, the composition will be normalized to 100 wt%
# 		(excluding H2O and CO2). Where a sample contains only FeO, ferric iron will
# 		be calculated using the ferric/total iron ratio provided.
#
# 		Parameters
# 		----------
# 		sample     pandas Series or dict
# 			Major element oxides in wt%.
# 		ferric_total    float
# 			Mole ratio of ferric to total iron to be used
# 			for calculating Fe2O3 and FeO when only FeO is
# 			provided. Default is 0.15.
#
# 		Returns
# 		-------
# 		pandas Series or dict
# 			Normalized major element oxides in wt%.
# 		"""
# 		if type(sample) != dict and type(sample) != pd.core.series.Series:
# 			raise core.InputError("sample must be a dict or a pandas Series.")
# 		if 'FeO' not in sample:
# 			raise core.InputError("sample must contain FeO.")
#
# 		_sample = sample.copy()
#
# 		for ox in ['TiO2','P2O5']:
# 			if ox not in _sample:
# 				_sample[ox] = 0.0
#
# 		if 'Fe2O3' not in _sample:
# 			Fe_t = _sample['FeO']/oxideMass['FeO']
# 			Fe3 = ferric_total*Fe_t
# 			Fe2 = Fe_t - Fe3
# 			_sample['FeO'] = Fe2*oxideMass['FeO']
# 			_sample['Fe2O3'] = Fe3*oxideMass['Fe2O3']/2
#
# 		return normalize_AdditionalVolatiles(_sample)
#
# 	def calculate_dissolved_volatiles(self,pressure,temperature,sample,X_fluid=1.0,**kwargs):
# 		"""
# 		Calculates the dissolved (total) CO2 using eqs (9) and (10) of Eguchi and Dasgupta (2018).
#
# 		Parameters
# 		----------
# 		pressure     float
# 			Pressure in bars
# 		temperature     float
# 			Temperature in C
# 		sample     pandas Series or dict
# 			Major element oxides in wt%.
# 		X_fluid     float
# 			The mole fraction of CO2 in the fluid.
#
# 		Returns
# 		-------
# 		float
# 			Dissolved CO2 concentration.
# 		"""
# 		if pressure < 0:
# 			raise core.InputError("Pressure must be greater than 0 bar.")
#
# 		if pressure == 0:
# 			return 0
#
# 		XCO3 = self.Xi_melt(pressure=pressure,temperature=temperature,sample=sample,X_fluid=X_fluid,species='CO3')
# 		XCO2 = self.Xi_melt(pressure=pressure,temperature=temperature,sample=sample,X_fluid=X_fluid,species='CO2')
#
# 		FW_one = wtpercentOxides_to_formulaWeight(sample)
#
# 		CO2_CO2 = ((44.01*XCO2)/(44.01*XCO2+(1-(XCO2+XCO3))*FW_one))*100
# 		CO2_CO3 = ((44.01*XCO3)/(44.01*XCO3+(1-(XCO2+XCO3))*FW_one))*100
#
# 		return CO2_CO2 + CO2_CO3
#
#
# 	def calculate_equilibrium_fluid_comp(self,pressure,temperature,sample,**kwargs):
# 		""" Returns 1.0 if a pure CO2 fluid is saturated.
# 		Returns 0.0 if a pure CO2 fluid is undersaturated.
#
# 		Parameters
# 		----------
# 		pressure     float
# 			The total pressure of the system in bars.
# 		temperature     float
# 			The temperature of the system in C.
# 		sample         pandas Series or dict
# 			Major element oxides in wt% (including H2O).
#
# 		Returns
# 		-------
# 		float
# 			1.0 if CO2-fluid saturated, 0.0 otherwise.
# 		"""
#
# 		satP = self.calculate_saturation_pressure(temperature=temperature,sample=sample,X_fluid=1.0,**kwargs)
# 		if pressure < satP:
# 			return 1.0
# 		else:
# 			return 0.0
#
# 	def calculate_saturation_pressure(self,temperature,sample,X_fluid=1.0,**kwargs):
# 		"""
# 		Calculates the pressure at which a pure CO2 fluid is saturated, for the given sample
# 		composition and CO2 concentration. Calls the scipy.root_scalar routine, which makes
# 		repeated called to the calculate_dissolved_volatiles method.
#
# 		Parameters
# 		----------
# 		temperature     float
# 			The temperature of the system in C.
# 		sample         pandas Series or dict
# 			Major element oxides in wt% (including CO2).
# 		X_fluid     float
# 			The mole fraction of H2O in the fluid. Default is 1.0.
#
# 		Returns
# 		-------
# 		float
# 			Calculated saturation pressure in bars.
# 		"""
#
# 		if 'CO2' not in sample:
# 			raise core.InputError("sample must contain CO2.")
# 		if sample['CO2'] < 0.0:
# 			raise core.InputError("Concentration of CO2 must be greater than 0 wt%.")
# 		try:
# 			satP = root_scalar(self.root_saturation_pressure,x0=1000.0,x1=2000.0,
# 								args=(temperature,sample,X_fluid,kwargs)).root
# 		except:
# 			w.warn("Saturation pressure not found.",RuntimeWarning,stacklevel=2)
# 			satP = np.nan
# 		return satP
#
# 	def root_saturation_pressure(self,pressure,temperature,sample,X_fluid,kwargs):
# 		""" Function called by scipy.root_scalar when finding the saturation pressure using
# 		calculate_saturation_pressure.
#
# 		Parameters
# 		----------
# 		pressure     float
# 			Pressure guess in bars
# 		temperature     float
# 			The temperature of the system in C.
# 		sample         pandas Series or dict
# 			Major elements in wt% (normalized to 100%), including CO2.
# 		kwargs         dictionary
# 			Additional keyword arguments supplied to calculate_saturation_pressure. Might be required for
# 			the fugacity or activity models.
#
# 		Returns
# 		-------
# 		float
# 			The differece between the dissolved CO2 at the pressure guessed, and the CO2 concentration
# 			passed in the sample variable.
# 		"""
# 		return sample['CO2'] - self.calculate_dissolved_volatiles(pressure=pressure,temperature=temperature,sample=sample,X_fluid=X_fluid,**kwargs)
#
# 	def Xi_melt(self,pressure,temperature,sample,species,X_fluid=1.0,**kwargs):
# 		"""
# 		Calculates the mole fraction of dissolved molecular CO2 or carbonate CO3(2-), using
# 		eqn (9) of Eguchi and Dasgupta (2018).
#
# 		Parameters
# 		----------
# 		pressure    float
# 			Pressure in bars.
# 		temperature     float
# 			Temperature in C.
# 		sample         pandas Series or dict
# 			Major element oxides in wt%.
# 		species        str
# 			Which species to calculate, molecular CO2 'CO2' or carbonate ion 'CO3'.
# 		X_fluid     float
# 			The mole fraction of CO2 in the fluid. Default is 1.0.
#
# 		Returns
# 		-------
# 		float
# 			Mole fraction of selected species in the melt
# 		"""
# 		temperature = temperature + 273.15 #translate T from C to K
#
# 		if all(ox in sample for ox in ['MgO','CaO','FeO','Na2O','K2O','MnO','Al2O3','Fe2O3','SiO2','TiO2','P2O5']) == False:
# 			raise core.InputError("sample must contain MgO, CaO, FeO, Na2O, K2O, MnO, Al2O3, Fe2O3, SiO3, TiO2, and P2O5.")
# 		if X_fluid < 0 or X_fluid > 1:
# 			raise core.InputError("X_fluid must have a value between 0 and 1.")
# 		if pressure < 0:
# 			raise core.InputError("Pressure must be positive.")
# 		if temperature <= 0:
# 			raise core.InputError("Temperature must be greater than 0K.")
#
# 		if species == 'CO3':
# 			DH = -1.65e5
# 			DV = 2.38e-5
# 			DS = -43.64
# 			B = 1.47e3
# 			yNBO = 3.29
# 			A_CaO = 1.68e5
# 			A_Na2O = 1.76e5
# 			A_K2O = 2.11e5
# 		elif species == 'CO2':
# 			DH = -9.02e4
# 			DV = 1.92e-5
# 			DS = -43.08
# 			B = 1.12e3
# 			yNBO = -7.09
# 			A_CaO = 0
# 			A_Na2O = 0
# 			A_K2O = 0
# 		else:
# 			raise core.InputError("species variable must be either 'CO2' or 'CO3'.")
# 		R = 8.314
#
# 		# Calculate NBO term
# 		cations = wtpercentOxides_to_molSingleO(sample)
# 		oxides = wtpercentOxides_to_molOxides(sample)
#
#
#
# 		NM = (cations['Mg'] + cations['Ca'] + cations['Fe'] + cations['Na'] +
# 			cations['K'] + cations['Mn'])
# 		Al = cations['Al'] - NM
# 		if Al > 0:
# 			Al = NM
# 		else:
# 			Al = cations['Al']
# 		Fe = cations['Fe3'] + Al
# 		if Al > 0:
# 			Fe = 0
# 		if Al < 0 and Fe > 0:
# 			Fe = - Al
# 		if Al < 0 and Fe < 0:
# 			Fe = cations['Fe3']
# 		Tet = cations['Si'] + cations['Ti'] + cations['P'] + Al + Fe
# 		NBO = 2 - 4*Tet
#
# 		lnfCO2 = np.log(self.fugacity_model.fugacity(pressure=pressure,temperature=temperature-273.15,X_fluid=X_fluid))
#
# 		lnXi = ((DH/(R*temperature)-(pressure*1e5*DV)/(R*temperature)+DS/R) +
# 				(A_CaO*oxides['CaO']+A_Na2O*oxides['Na2O']+A_K2O*oxides['K2O'])/(R*temperature) +
# 				(B*lnfCO2/temperature) + yNBO*NBO
# 				)
#
# 		return np.exp(lnXi)
