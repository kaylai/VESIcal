#!/opt/conda/bin/python

import VESIcal as v
import anvil.server
import anvil.media
import pandas as pd
from anvil.tables import app_tables

oxides = v.oxides

anvil.server.connect("KYKUMVNXIZ23YPRPDDGL5DM6-E4XIIJ7V23RSDNS5")
app_tables.user_data.delete_all_rows()

gMyFile = None

def add_df_to_anvil_table(dataframe, anvil_table):
    for d in dataframe.to_dict(orient="records"):
          # d is now a dict of {columnname -> value} for this row
          # We use Python's **kwargs syntax to pass the whole dict as
          # keyword arguments
          anvil_table.add_row(**d)

@anvil.server.callable
def get_user_data():
    return app_tables.user_data.search()

@anvil.server.callable
def import_ExcelFile(file):
    global gMyFile
    with anvil.media.TempFile(file) as filename:
        gMyFile = v.BatchFile(filename)

    #gMyFile.data = gMyFile.data.reset_index().set_index('Label', drop=False) #turn Label column into normal column so that it's not lost
    add_df_to_anvil_table(gMyFile.data, app_tables.user_data)

    return("Your file has been loaded!")

@anvil.server.callable
def anvil_calculate_saturation_pressure(temperature, model, **kwargs):
    global gMyFile
    global gSatP_results_df
    gSatP_results_df = gMyFile.calculate_saturation_pressure(temperature=temperature, model=model)
    return("Successfully calculated saturation pressures")

@anvil.server.callable
def anvil_calculate_dissolved_volatiles(temperature, pressure, X_fluid, model, **kwargs):
    global gMyFile
    global gDissolved_results_df
    gDissolved_results_df = gMyFile.calculate_dissolved_volatiles(temperature=temperature, 
                                                                pressure=pressure,
                                                                X_fluid=X_fluid,
                                                                model=model)
    return("Successfully calculated dissolved volatile concentrations")

@anvil.server.callable
def anvil_calculate_equilibrium_fluid_comp(temperature, pressure, model, **kwargs):
    global gMyFile
    global gEqfluid_results_df
    gEqfluid_results_df = gMyFile.calculate_equilibrium_fluid_comp(temperature=temperature, 
                                                                pressure=pressure,
                                                                model=model)
    return("Successfully calculated equilibrium fluid compositions")

@anvil.server.callable
def download_satPs():
    global gSatP_results_df
    with pd.ExcelWriter('VESIcal_satPs.xlsx') as writer:
        gSatP_results_df.to_excel(writer)

    return anvil.media.from_file("./VESIcal_satPs.xlsx")

@anvil.server.callable
def download_dissolved():
    global gDissolved_results_df
    with pd.ExcelWriter('VESIcal_dissolved.xlsx') as writer:
        gDissolved_results_df.to_excel(writer)

    return anvil.media.from_file("./VESIcal_dissolved.xlsx")

@anvil.server.callable
def download_eqfluid():
    global gEqfluid_results_df
    with pd.ExcelWriter('VESIcal_eqfluid.xlsx') as writer:
        gEqfluid_results_df.to_excel(writer)

    return anvil.media.from_file("./VESIcal_eqfluid.xlsx")

@anvil.server.callable
def download_example_file():
    return anvil.media.from_file("./anvil/VESIcal_example_data.xlsx")

anvil.server.wait_forever()
