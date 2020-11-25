import math
import pandas as pd
import sys
sys.path.insert(0, '../')

#Load in datasets
filename = '../Solubility_Datasets.xlsx'
#CO2 only
df_Eguchi_CO2= pd.read_excel(filename, sheet_name='Eguchi_CO2', index_col=0)
df_Allison_CO2= pd.read_excel(filename, sheet_name='Allison_CO2', index_col=0)
df_Dixon_CO2= pd.read_excel(filename, sheet_name='Dixon_CO2', index_col=0)
df_MagmaSat_CO2= pd.read_excel(filename, sheet_name='MagmaSat_CO2', index_col=0)
df_Shishkina_CO2=pd.read_excel(filename, sheet_name='Shishkina_PureCO2', index_col=0)

#H2O Only
df_Iacono_H2O= pd.read_excel(filename, sheet_name='Iacono_H2O', index_col=0)
df_Shishkina_H2O=pd.read_excel(filename, sheet_name='Shishkina_H2O', index_col=0)
df_MagmaSat_H2O= pd.read_excel(filename, sheet_name='MagmSat_H2OExt', index_col=0)
df_Dixon_H2O=pd.read_excel(filename, sheet_name='Dixon_H2O', index_col=0)
df_Moore_H2O=pd.read_excel(filename, sheet_name='Moore_H2O', index_col=0)
df_Liu_H2O=pd.read_excel(filename, sheet_name='Liu_H2O', index_col=0)

#Mixed CO2-H2O
df_Iacono_CO2H2O= pd.read_excel(filename, sheet_name='Iacono_H2O-CO2', index_col=0)
df_MagmaSat_CO2H2O= pd.read_excel(filename, sheet_name='MagmaSat_CO2H2O', index_col=0)
df_Shishkina_CO2H2O = pd.read_excel(filename, sheet_name='Shishkina_MixedCO2', index_col=0)
df_Liu_CO2H2O = pd.read_excel(filename, sheet_name='Liu_CO2H2O', index_col=0)

#Create calibration file
f= open("hard_coded_calibrations.py","w+")
list_of_models=[df_Eguchi_CO2,
                df_Allison_CO2,
                df_Dixon_CO2,
                df_Dixon_H2O,
                df_Iacono_H2O,
                df_Iacono_CO2H2O,
                df_Liu_H2O,
                df_Liu_CO2H2O,
                df_MagmaSat_CO2,
                df_MagmaSat_H2O,
                df_MagmaSat_CO2H2O,
                df_Moore_H2O,
                df_Shishkina_CO2,
                df_Shishkina_H2O,
                df_Shishkina_CO2H2O
                ]

list_of_modelnames=['df_Eguchi_CO2',
                    'df_Allison_CO2',
                    'df_Dixon_CO2',
                    'df_Dixon_H2O',
                    'df_Iacono_H2O',
                    'df_Iacono_CO2H2O',
                    'df_Liu_CO2H2O',
                    'df_Liu_H2O',
                    'df_MagmaSat_CO2',
                    'df_MagmaSat_H2O',
                    'df_MagmaSat_CO2H2O',
                    'df_Moore_H2O',
                    'df_Shishkina_CO2',
                    'df_Shishkina_H2O',
                    'df_Shishkina_CO2H2O'
                    ]

for i in range(len(list_of_models)):
    current = list_of_models[i]
    oxides = ['SiO2', 'TiO2', 'Al2O3', 'Fe2O3', 'Cr2O3', 'FeO', 'MnO', 'MgO', 'NiO', 'CoO', 'CaO', 'Na2O', 'K2O', 'P2O5',
              'H2O', 'CO2', 'Na2O+K2O']

    f.write("\n")
    f.write(str(list_of_modelnames[i])+" = pd.DataFrame({ \n")
    for oxide in oxides:
        iterno = 1
        if oxide in current.columns.values:
            f.write("'" + (str(oxide)+"': ["))
            for index, row in current.iterrows():
                if math.isnan(row[oxide]):
                    f.write("float('nan')")
                else:
                    f.write(str(row[oxide]))
                if iterno < len(current.index):
                    f.write(",")
                iterno += 1
            f.write("], \n")
    f.write(" }) \n")