import VESIcal as v
import anvil.server
import anvil.media
import pandas as pd
from anvil.tables import app_tables

oxides = ['SiO2', 'TiO2', 'Al2O3', 'Fe2O3', 'Cr2O3', 'FeO', 'MnO', 'MgO', 'NiO', 'CoO', 'CaO', 'Na2O', 'K2O', 'P2O5',
		  'H2O', 'CO2']

anvil.server.connect("KYKUMVNXIZ23YPRPDDGL5DM6-E4XIIJ7V23RSDNS5")
app_tables.user_data.delete_all_rows()

@anvil.server.callable
def get_user_data():
	return app_tables.user_data.search()

@anvil.server.callable
def import_ExcelFile(file):
	with anvil.media.TempFile(file) as filename:
		myfile = v.ExcelFile(filename)

	#this isn't working :(
	myfile.data = myfile.data.reset_index().set_index('Label', drop=False) #turn Label column into normal column so that it's not lost

	for d in myfile.data.to_dict(orient="records"):
		# d is now a dict of {columnname -> value} for this row
		# We use Python's **kwargs syntax to pass the whole dict as
		# keyword arguments
		app_tables.user_data.add_row(**d)

	return("Your file has been loaded!")

@anvil.server.callable
def anvil_calculate_saturation_pressure(model, **kwargs):
	user_data = app_tables.user_data.search()

	columnnames = ['Label','SiO2', 'TiO2', 'Al2O3', 'Fe2O3', 'Cr2O3', 'FeO', 'MnO', 'MgO', 'NiO', 'CoO', 'CaO', 'Na2O', 'K2O', 'P2O5',
		  'H2O', 'CO2']

	myfile_data = pd.DataFrame({oxide: row[oxide] for oxide in oxides} for row in user_data)

	iterno = 0
	for index, row in myfile_data.iterrows():
		iterno+=1
		result = v.calculate_saturation_pressure(sample=row, temperature=1000, model=model).result
		app_tables.output_data.add_row(Name="Sample No."+str(iterno),
										SaturationP_bars=result)

	return("Done!")

@anvil.server.callable
def pull_data():
	output_data = app_tables.output_data.search()

	results_df = pd.DataFrame({'Name': row['Name'], 'SaturationP_bars_VESIcal': row['SaturationP_bars']} for row in output_data)

	with pd.ExcelWriter('test_output.xlsx') as writer:
			return results_df.to_excel(writer)

anvil.server.wait_forever()