import MagmaSatPlus as m
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# wtpt_oxides = pd.Series({'SiO2':50.17,
#                              'TiO2':0.92,
#                              'Al2O3':18.28,
#                              'FeO':9.37,
#                              'MnO':0.17,
#                              'MgO':7.00,
#                              'CaO':11.37,
#                              'Na2O':2.33,
#                              'K2O':0.23,
#                              'P2O5':0.15,
#                              'CO2':0.02,
#                              'H2O':2.0})

wtpt_oxides = pd.Series({'SiO2':47.95,
                         'TiO2':1.67,
                         'Al2O3':17.32,
                         'FeO':10.24,
                         'MgO':5.76,
                         'CaO':10.93,
                         'Na2O':3.45,
                         'K2O':1.99,
                         'P2O5':0.51,
                         'CO2':0.8,
                         'H2O':4.0})

#Mono Craters, CA from Liu et al (2005)
mono = pd.Series({'SiO2':77.19,
                         'TiO2':0.06,
                         'Al2O3':12.80,
                         'FeO':0.94,
                         'MgO':0.03,
                         'CaO':0.53,
                         'Na2O':3.98,
                         'K2O':4.65,
                         'CO2':0.05,
                         'H2O':0.26})

mono_incVols = pd.Series({'SiO2':77.19,
                         'TiO2':0.06,
                         'Al2O3':12.80,
                         'FeO':0.94,
                         'MgO':0.03,
                         'CaO':0.53,
                         'Na2O':3.98,
                         'K2O':4.65,
                         'CO2':0.55,
                         'H2O':2.26})

#Experimental sample from Tamic et al. (2001)
Ech6 = pd.Series({'SiO2':77.04,
                         'TiO2':0.11,
                         'Al2O3':12.76,
                         'FeO':0.68,
                         'MnO': 0.07,
                         'MgO':0.08,
                         'CaO':0.58,
                         'Na2O':4.07,
                         'K2O':4.79,
                         'CO2':0.066,
                         'H2O':3.84})



print(m.calculate_dissolved_volatiles(pressure=5000.0,X_fluid=(0.5,0.5),sample=wtpt_oxides,temperature=1473.15,model='Shishkina').result)
print(m.calculate_dissolved_volatiles(pressure=5000.0,X_fluid=(0.5,0.5),sample=wtpt_oxides,temperature=1473.15,model='Shishkina').calib_check)


isobars, isopleths = m.calculate_isobars_and_isopleths(pressure_list=[500.0,1000,2000,3000,4000,5000],#,6000,7000,8000],
                                                               isopleth_list=[0.0,0.2,0.4,0.6,0.8,1.0],
                                                               sample=wtpt_oxides,temperature=1473.15,points=101,
                                                               model='Shishkina').result

degas_melt, degas_fluid = m.calculate_degassing_paths(pressures=np.linspace(6000.0,1.0,50),temperature=1473.15,sample=wtpt_oxides,fractionate_vapor=0.0,model='Shishkina').result
degas_melto, degas_fluid = m.calculate_degassing_paths(pressures=np.linspace(6000.0,1.0,50),temperature=1473.15,sample=wtpt_oxides,fractionate_vapor=1.0,model='Shishkina').result


f,a = plt.subplots()
for isobar in isobars:
    a.plot(isobar[1,:],isobar[0,:],c='k')
for isopleth in isopleths:
    a.plot(isopleth[1,:],isopleth[0,:],c='k')
a.plot(degas_melt[1,:],degas_melt[0,:],c='r')
a.plot(degas_melto[1,:],degas_melto[0,:],c='b')
plt.show()



# h_model = m.DixonWater()
# print(h_model.calculate_dissolved_volatiles(pressure=1000.0,temperature=1473.15,sample=wtpt_oxides,X_fluid=0.0))

c_model = m.IaconoMarzianoCarbon()
h_model = m.IaconoMarzianoWater()

print(c_model.calculate_dissolved_volatiles(sample=wtpt_oxides,pressure=4000.0,temperature=973.15))
print(c_model.calculate_equilibrium_fluid_comp(sample=wtpt_oxides,pressure=1000.0,temperature=973.15))
print(c_model.calculate_saturation_pressure(sample=wtpt_oxides,temperature=973.15))

print(h_model.calculate_dissolved_volatiles(sample=wtpt_oxides,pressure=1000.0,temperature=973.15))
print(h_model.calculate_equilibrium_fluid_comp(sample=wtpt_oxides,pressure=3000.0,temperature=973.15))
print(h_model.calculate_saturation_pressure(sample=wtpt_oxides,temperature=973.15))
#
#
#
# c_model = m.ShishkinaCarbon()
# h_model = m.ShishkinaWater()

# for c_model, h_model in zip([m.ShishkinaCarbon(),m.DixonCarbon(),m.IaconoMarzianoCarbon(),m.EguchiCarbon(),m.AllisonCarbon()],
#                             [m.ShishkinaWater(),m.DixonWater(),m.IaconoMarzianoWater(),m.DixonWater(),m.DixonWater()]):
# for c_model, h_model in zip([m.ShishkinaCarbon()],[m.ShishkinaWater()]):
#     model = m.MixedFluids({'H2O':h_model,'CO2':c_model})
#
#     print('Calculating the concentrations in the melt...')
#     dissolved = model.calculate_dissolved_volatiles(pressure=5000.0,temperature=973.15,X_fluid=(0.5,0.5),sample=wtpt_oxides)
#
#     for string in ('%s %f wtpt' % (species, conc) for (species, conc) in zip(model.volatile_species,dissolved)):
#         print(string)
#
#     # print('\nCalculating the mole fractions in the melt...')
#     # molfrac = model.molfrac_molecular(pressure=5000.0,temperature=973.15,X_fluid=(0.5,0.5),sample=wtpt_oxides)
#     #
#     # for string in ('%s %f' % (species, conc) for (species, conc) in zip(model.volatile_species,molfrac)):
#     #     print(string)
#
#
#     wtpt_oxides['H2O'] = 2.0
#     wtpt_oxides['CO2'] = 0.8
#
#     print('\nCalculating the composition of the fluid...')
#     fluid_comp = model.calculate_equilibrium_fluid_comp(pressure=1000.0,temperature=973.15,sample=wtpt_oxides)
#     for string in ('%s %f' % (species, conc) for (species, conc) in zip(model.volatile_species,fluid_comp)):
#         print(string)
#
#     print('\nCalculating Saturation Pressure...')
#     print('Water only...')
#     print('Saturation pressure for {:.1f} wt% H2O is {:.1f} bar.'.format(wtpt_oxides['H2O'],h_model.calculate_saturation_pressure(sample=wtpt_oxides,temperature=973.15)))
#     print('Carbon only...')
#     print('Saturation pressure for {:.1f} wt% CO2 is {:.1f} bar.'.format(wtpt_oxides['CO2'],c_model.calculate_saturation_pressure(sample=wtpt_oxides,temperature=973.15)))
#     print('Mixed fluid...')
#     # print(model.calculate_saturation_pressure(wtpt_oxides))
#     print('Saturation pressure is {:.1f} bar.'.format(model.calculate_saturation_pressure(sample=wtpt_oxides,temperature=973.15)))
#
#     isobars, isopleths = model.calculate_isobars_and_isopleths(pressure_list=[500.0,1000,2000,3000,4000,5000],#,6000,7000,8000],
#                                                                isopleth_list=[0.0,0.2,0.4,0.6,0.8,1.0],
#                                                                sample=wtpt_oxides,temperature=1473.15,points=101)
#     wtpt_oxides = m.normalize_FixedVolatiles(wtpt_oxides)
#     satP = model.calculate_saturation_pressure(sample=wtpt_oxides,temperature=1473.15)
#     degas_melt, degas_fluid = model.calculate_degassing_paths(pressures=np.linspace(satP,1.0,50),temperature=1473.15,sample=wtpt_oxides,fractionate_vapor=0.0)
#     degas_melto, degas_fluid = model.calculate_degassing_paths(pressures=np.linspace(satP,1.0,50),temperature=1473.15,sample=wtpt_oxides,fractionate_vapor=1.0)
#     f,a = plt.subplots()
#     for isobar in isobars:
#         a.plot(isobar[0,:],isobar[1,:],c='k')
#     for isopleth in isopleths:
#         a.plot(isopleth[0,:],isopleth[1,:],c='k')
#     a.plot(degas_melt[0,:],degas_melt[1,:],c='r')
#     a.plot(degas_melto[0,:],degas_melto[1,:],c='b')
#     plt.show()
