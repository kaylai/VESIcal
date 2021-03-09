import VESIcal as v
import pickle

mysample = v.Sample({'SiO2':  77.3,
 'TiO2':   0.08,
 'Al2O3': 12.6,
 'Fe2O3':  0.207,
 'Cr2O3':  0.0,
 'FeO':    0.473,
 'MnO':    0.0,
 'MgO':    0.03,
 'NiO':    0.0,
 'CoO':    0.0,
 'CaO':    0.43,
 'Na2O':   3.98,
 'K2O':    4.88,
 'P2O5':   0.0,
 'H2O':    6.5,
 'CO2':    0.05})

# satP = v.calculate_saturation_pressure(mysample, temperature=800).result
# print("SatP")
# print(satP)
# print("\n")

# dissvol = v.calculate_dissolved_volatiles(mysample, temperature=800, pressure=1000).result
# print("Dissolved volatiles")
# print(dissvol)
# print("\n")

# eqfluid = v.calculate_equilibrium_fluid_comp(mysample, temperature=800, pressure=1000).result
# print("Equilibrium Fluid")
# print(eqfluid)
# print("\n")

# try:
# 	isobars = pickle.load(open("isobars.p", "rb"))
# 	isopleths = pickle.load(open("isopleths.p", "rb"))
# except:
# 	isobars, isopleths = v.calculate_isobars_and_isopleths(mysample, temperature=800, pressure_list=[100, 500, 1000], isopleth_list=[.25]).result
# 	pickle.dump(isobars, open("isobars.p", "wb"))
# 	pickle.dump(isopleths, open("isopleths.p", "wb"))

# fig, ax = v.plot(isobars=isobars, isopleths=isopleths)
# v.show()


try:
	dpath = pickle.load(open("dpath.p", "rb"))
except:
	dpath = v.calculate_degassing_path(mysample, temperature=800).result
	pickle.dump(dpath, open("dpath.p", "wb"))

fig, ax = v.plot(degassing_paths=dpath)
v.show()