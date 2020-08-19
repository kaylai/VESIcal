import VESIcal as v 

#-----------------IMPORTS-----------------#
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from abc import ABC, abstractmethod
from scipy.optimize import root_scalar
from scipy.optimize import root
from scipy.optimize import minimize
import sys
import sympy

oxides = ['SiO2', 'TiO2', 'Al2O3', 'Fe2O3', 'Cr2O3', 'FeO', 'MnO', 'MgO', 'NiO', 'CoO', 'CaO', 'Na2O', 'K2O', 'P2O5',
		  'H2O', 'CO2']

#--------------MELTS preamble---------------#
from thermoengine import equilibrate
# instantiate thermoengine equilibrate MELTS instance
melts = equilibrate.MELTSmodel('1.2.0')

# Suppress phases not required in the melts simulation
phases = melts.get_phase_names()
for phase in phases:
	melts.set_phase_inclusion_status({phase: False})
melts.set_phase_inclusion_status({'Fluid': True, 'Liquid': True})



myfile = v.ExcelFile('testDataSets/example_data.xlsx')
samp = myfile.get_sample_oxide_comp('10*')

sample = samp 
temperature = 1200
pressure = 'saturation'
fractionate_vapor=0.0
init_vapor = 0.0




sample = v.normalize(sample)
bulk_comp_orig = sample

bulk_comp = {oxide:  sample[oxide] for oxide in oxides}
feasible = melts.set_bulk_composition(bulk_comp)

# Get saturation pressure
data = v.calculate_saturation_pressure(sample=sample, temperature=temperature, verbose=True).result

if pressure == 'saturation' or pressure >= data["SaturationP_bars"]:
	SatP_MPa = data["SaturationP_bars"] / 10
else:
	SatP_MPa = pressure / 10

#If pressure is low, use smaller P steps
if SatP_MPa >= 50:
	MPa_step = 10
elif SatP_MPa < 50:
	MPa_step = 1

P_array = np.arange(SatP_MPa, 1, -10)
fl_wtper = data["FluidProportion_wt"]

if fractionate_vapor == 0: #closed-system
	while fl_wtper <= init_vapor:
		output = melts.equilibrate_tp(temperature, SatP_MPa)
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

	output = melts.equilibrate_tp(temperature, P_array)

	pressure_list = []
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

		pressure_list.append(p * 10)
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
			bulk_comp["H2O"] = liq_comp["H2O"]
		except:
			bulk_comp["H2O"] = 0
		try:
			bulk_comp["CO2"] = liq_comp["CO2"]
		except:
			bulk_comp["CO2"] = 0
		fluid_wtper.append(fl_wtper)

	feasible = melts.set_bulk_composition(bulk_comp_orig)
	fl_wtper = data["FluidProportion_wt"]
	exsolved_degassing_df = pd.DataFrame(list(zip(pressure_list, H2Oliq, CO2liq, H2Ofl, CO2fl, fluid_wtper)),
								columns =['Pressure_bars', 'H2O_liq', 'CO2_liq', 'H2O_fl', 'CO2_fl', 'FluidProportion_wt'])

	print(exsolved_degassing_df)
else:
	pressure = []
	H2Oliq = []
	CO2liq = []
	H2Ofl = []
	CO2fl = []
	fluid_wtper = []
	for i in P_array:
		fl_mass = 0
		feasible = melts.set_bulk_composition(bulk_comp)
		output = melts.equilibrate_tp(temperature, i, initialize=True)
		(status, temperature, p, xmlout) = output[0]
		liq_comp = melts.get_composition_of_phase(xmlout, phase_name='Liquid')
		fl_comp = melts.get_composition_of_phase(xmlout, phase_name='Fluid', mode='component')
		liq_mass = melts.get_mass_of_phase(xmlout, phase_name='Liquid')
		fl_mass = melts.get_mass_of_phase(xmlout, phase_name='Fluid')
		fl_wtper = 100 * fl_mass / (fl_mass+liq_mass)

		if fl_mass > 0:
			pressure.append(p * 10)
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
				bulk_comp["H2O"] = liq_comp["H2O"] + (bulk_comp["H2O"] - liq_comp["H2O"]) * (1-fractionate_vapor)
			except:
				bulk_comp["H2O"] = 0
			try:
				bulk_comp["CO2"] = liq_comp["CO2"] + (bulk_comp["CO2"] - liq_comp["CO2"]) * (1-fractionate_vapor)
			except:
				bulk_comp["CO2"] = 0
			#bulk_comp = normalize(bulk_comp)

	feasible = melts.set_bulk_composition(bulk_comp_orig) #this needs to be reset always!
	open_degassing_df = pd.DataFrame(list(zip(pressure, H2Oliq, CO2liq, H2Ofl, CO2fl, fluid_wtper)),
								columns =['Pressure_bars', 'H2O_liq', 'CO2_liq', 'XH2O_fl', 'XCO2_fl', 'FluidProportion_wt'])

	print(open_degassing_df)