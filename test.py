import MagmaSatPlus as m
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

wtpt_oxides = pd.Series({'SiO2':50.17,
                             'TiO2':0.92,
                             'Al2O3':18.28,
                             'FeO':9.37,
                             'MnO':0.17,
                             'MgO':7.00,
                             'CaO':11.37,
                             'Na2O':2.33,
                             'K2O':0.23,
                             'P2O5':0.15})


c_model = m.ShishkinaCarbon()
h_model = m.ShishkinaWater()

model = m.MixedFluids({'H2O':h_model,'CO2':c_model})

print('Calculating the concentrations in the melt...')
dissolved = model.calculate_dissolved_volatiles(pressure=5000.0,X_fluid=(0.5,0.5),sample=wtpt_oxides)

for string in ('%s %f wtpt' % (species, conc) for (species, conc) in zip(model.volatile_species,dissolved)):
    print(string)

print('\nCalculating the mole fractions in the melt...')
molfrac = model.molfrac_molecular(pressure=5000.0,X_fluid=(0.5,0.5),sample=wtpt_oxides)

for string in ('%s %f' % (species, conc) for (species, conc) in zip(model.volatile_species,molfrac)):
    print(string)


wtpt_oxides['H2O'] = 6.0
wtpt_oxides['CO2'] = 0.2

print('\nCalculating the composition of the fluid...')
fluid_comp = model.calculate_equilibrium_fluid_comp(pressure=7000.0,sample=wtpt_oxides)
for string in ('%s %f' % (species, conc) for (species, conc) in zip(model.volatile_species,fluid_comp)):
    print(string)


print('\nCalculating Saturation Pressure...')
print('Water only...')
print('Saturation pressure for {:.1f} wt% H2O is {:.1f} bar.'.format(wtpt_oxides['H2O'],h_model.calculate_saturation_pressure(wtpt_oxides)))
print('Carbon only...')
print('Saturation pressure for {:.1f} wt% CO2 is {:.1f} bar.'.format(wtpt_oxides['CO2'],c_model.calculate_saturation_pressure(wtpt_oxides)))
print('Mixed fluid...')
# print(model.calculate_saturation_pressure(wtpt_oxides))
print('Saturation pressure is {:.1f} bar.'.format(model.calculate_saturation_pressure(wtpt_oxides)))

isobars = model.calculate_isobars(pressure_list=[1000,2000,3000,4000,5000,6000,7000,8000],sample=wtpt_oxides)
f,a = plt.subplots()
for isobar in isobars:
    a.plot(isobar[0,:],isobar[1,:])
plt.show()
