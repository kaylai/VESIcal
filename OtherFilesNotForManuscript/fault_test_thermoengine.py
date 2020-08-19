import numpy as np

#--------------MELTS preamble---------------#
from thermoengine import equilibrate
# instantiate thermoengine equilibrate MELTS instance
melts = equilibrate.MELTSmodel('1.2.0')

# Suppress phases not required in the melts simulation
phases = melts.get_phase_names()
for phase in phases:
	melts.set_phase_inclusion_status({phase: False})
melts.set_phase_inclusion_status({'Fluid': True, 'Liquid': True})

bulk_comp = {'SiO2':50.63,
             'TiO2':0.95,
             'Al2O3':19.81,
             'FeO':9.90,
             'MgO':3.81,
             'MnO':0.12,
             'CaO':11.51,
             'Na2O':2.65,
             'K2O':0.42,
             'P2O5':0.10,
             'CO2':0.02426744,
             'H2O':1.6}

temperature = 950



bulk_comp_orig = bulk_comp

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

print({"SaturationP_bars": satP, "FluidMass_grams": flmass, "FluidProportion_wt": flsystem_wtper,
	 				"XH2O_fl": flH2O, "XCO2_fl": flCO2})