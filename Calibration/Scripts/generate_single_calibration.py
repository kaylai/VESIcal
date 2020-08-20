import pandas as pd
import numpy as np
df_Liu_H2O = pd.read_excel('Liu_2005.xlsx', sheet_name='H2O')
df_Liu_CO2 = pd.read_excel('Liu_2005.xlsx', sheet_name='CO2')
df_Liu_H2O.fillna("No Data")
df_Liu_CO2.fillna("No Data")

import math
f= open("hard_coded_Liu_calib.py","w+")
oxides = ['SiO2', 'TiO2', 'Al2O3', 'FeOt', 'MgO', 'CaO', 'Na2O', 'K2O']

f.write("\n")
f.write("df_Liu_H2O = pd.DataFrame({ \n")
for oxide in oxides:
	iterno = 1
	f.write("'" + (str(oxide)+"': ["))
	for index, row in df_Liu_H2O.iterrows():
		try:
			value = float(row[oxide])
			if str(value) == 'nan':
				f.write("'No Data'")
			else:
				f.write(str(value))
		except:
			f.write("'" + str(row[oxide]) + "'")
		if iterno < len(df_Liu_H2O.index):
			f.write(",")
		iterno += 1
	f.write("], \n")
f.write(" }) \n")

f.write("\n")
f.write("df_Liu_CO2 = pd.DataFrame({ \n")
for oxide in oxides:
	iterno = 1
	f.write("'" + (str(oxide)+"': ["))
	for index, row in df_Liu_CO2.iterrows():
		try:
			value = float(row[oxide])
			if str(value) == 'nan':
				f.write("'No Data'")
			else:
				f.write(str(value))
		except:
			f.write("'" + str(row[oxide]) + "'")
		if iterno < len(df_Liu_CO2.index):
			f.write(",")
		iterno += 1
	f.write("], \n")
f.write(" }) \n")