# ---------- DEFINE SOME CONSTANTS ------------- #
oxides = ['SiO2', 'TiO2', 'Al2O3', 'Fe2O3', 'Cr2O3', 'FeO', 'MnO', 'MgO', 'NiO', 'CoO', 'CaO',
          'Na2O', 'K2O', 'P2O5', 'H2O', 'CO2', 'F2O']
cations = ['Si', 'Ti', 'Al', 'Fe', 'Ca', 'Al', 'Na', 'K', 'Mn', 'Ti', 'P', 'Cr', 'Ni', 'Co',
           'Fe3', 'H', 'C', 'F']
magmasat_oxides = ['SiO2', 'TiO2', 'Al2O3', 'Fe2O3', 'Cr2O3', 'FeO', 'MnO', 'MgO', 'NiO', 'CoO',
                   'CaO', 'Na2O', 'K2O', 'P2O5', 'H2O', 'CO2']
anhydrous_oxides = ['SiO2', 'TiO2', 'Al2O3', 'Fe2O3', 'Cr2O3', 'FeO', 'MnO', 'MgO', 'NiO', 'CoO',
                    'CaO', 'Na2O', 'K2O', 'P2O5', 'F2O']
volatiles = ['H2O', 'CO2']
oxideMass = {'SiO2':  60.083,
             'MgO':   40.304,
             'FeO':   71.844,
             'CaO':   56.077,
             'Al2O3': 101.961,
             'Na2O':  61.979,
             'K2O':   94.195,
             'MnO':   70.937,
             'TiO2':  79.867,
             'P2O5':  141.943,
             'Cr2O3': 151.992,
             'NiO':   74.692,
             'CoO':   44.01,
             'Fe2O3': 159.687,
             'H2O':   18.02,
             'CO2':   44.01,
             'F2O':   37.997}

CationNum = {'SiO2': 1, 'MgO': 1, 'FeO': 1, 'CaO': 1, 'Al2O3': 2, 'Na2O': 2,
             'K2O': 2, 'MnO': 1, 'TiO2': 1, 'P2O5': 2, 'Cr2O3': 2,
             'NiO': 1, 'CoO': 1, 'Fe2O3': 2, 'H2O': 2, 'CO2': 1, 'F2O': 2}

OxygenNum = {'SiO2': 2, 'MgO': 1, 'FeO': 1, 'CaO': 1, 'Al2O3': 3, 'Na2O': 1,
             'K2O': 1, 'MnO': 1, 'TiO2': 2, 'P2O5': 5, 'Cr2O3': 3,
             'NiO': 1, 'CoO': 1, 'Fe2O3': 3, 'H2O': 1, 'CO2': 2, 'F2O': 1}

CationCharge = {'SiO2': 4, 'MgO': 2, 'FeO': 2, 'CaO': 2, 'Al2O3': 3, 'Na2O': 1,
                'K2O': 1, 'MnO': 2, 'TiO2': 4, 'P2O5': 5, 'Cr2O3': 3,
                'NiO': 2, 'CoO': 2, 'Fe2O3': 3, 'H2O': 1, 'CO2': 4, 'F2O': 1}

CationMass = {'SiO2': 28.085, 'MgO': 24.305, 'FeO': 55.845, 'CaO': 40.078, 'Al2O3': 26.982,
              'Na2O': 22.990, 'K2O': 39.098, 'MnO': 54.938, 'TiO2': 47.867, 'P2O5': 30.974,
              'Cr2O3': 51.996, 'NiO': 58.693, 'CoO': 28.01, 'Fe2O3': 55.845, 'H2O': 1.01,
              'CO2': 12.011, 'F2O': 18.998}

oxides_to_cations = {'SiO2': 'Si', 'MgO': 'Mg', 'FeO': 'Fe', 'CaO': 'Ca', 'Al2O3': 'Al',
                     'Na2O': 'Na', 'K2O': 'K', 'MnO': 'Mn', 'TiO2': 'Ti', 'P2O5': 'P',
                     'Cr2O3': 'Cr', 'NiO': 'Ni', 'CoO': 'Co', 'Fe2O3': 'Fe3', 'H2O': 'H',
                     'CO2': 'C', 'F2O': 'F'}
cations_to_oxides = {'Si': 'SiO2', 'Mg': 'MgO', 'Fe': 'FeO', 'Ca': 'CaO', 'Al': 'Al2O3',
                     'Na': 'Na2O', 'K': 'K2O', 'Mn': 'MnO', 'Ti': 'TiO2', 'P': 'P2O5',
                     'Cr': 'Cr2O3', 'Ni': 'NiO', 'Co': 'CoO', 'Fe3': 'Fe2O3', 'H': 'H2O',
                     'C': 'CO2', 'F': 'F2O'}


# ---------- DATA TRANSFORMATION FOR PANDAS DATAFRAMES --------- #
def fluid_molfrac_to_wt(data, H2O_colname='XH2O_fl_VESIcal', CO2_colname='XCO2_fl_VESIcal'):
    """
    Takes in a pandas dataframe object and converts only the fluid composition from mole fraction
    to wt%, leaving the melt composition in tact. The user must specify the names of the
    XH2O_fl and XCO2_fl columns.

    Parameters
    ----------
    data: pandas DataFrame
        Sample composition(s) containing columns for H2O and CO2 concentrations in the fluid.

    H2O_colname: str
        OPTIONAL. The default value is 'XH2O_fl', which is what is returned by BatchFile() core
        calculations. String containing the name of the column corresponding to the H2O
        concentration in the fluid, in mol fraction.

    CO2_colname: str
        OPTIONAL. The default value is 'XCO2_fl', which is what is returned by BatchFile() core
        calculations. String containing the name of the column corresponding to the CO2
        concentration in the fluid, in mol fraction.

    Returns
    -------
    pandas DataFrame
        Original data passed plus newly calculated values are returned.
    """
    convData = data.copy()

    MPO_H2O_list = []
    MPO_CO2_list = []
    for index, row in convData.iterrows():
        MPO_H2O_list.append(row[H2O_colname] * oxideMass["H2O"])
        MPO_CO2_list.append(row[CO2_colname] * oxideMass["CO2"])

    convData["MPO_H2O"] = MPO_H2O_list
    convData["MPO_CO2"] = MPO_CO2_list
    convData["H2O_fl_wt"] = 100 * convData["MPO_H2O"] / (convData["MPO_H2O"] + convData["MPO_CO2"])
    convData["CO2_fl_wt"] = 100 * convData["MPO_CO2"] / (convData["MPO_H2O"] + convData["MPO_CO2"])

    del convData["MPO_H2O"]
    del convData["MPO_CO2"]

    return convData


def fluid_wt_to_molfrac(data, H2O_colname='H2O_fl_wt', CO2_colname='CO2_fl_wt'):
    """
    Takes in a pandas dataframe object and converts only the fluid composition from wt% to mole
    fraction, leaving the melt composition in tact. The user must specify the names of the
    H2O_fl_wt and CO2_fl_wt columns.

    Parameters
    ----------
    data: pandas DataFrame
        DataFrame containing columns for H2O and CO2 concentrations in the fluid.

    H2O_colname: str
        OPTIONAL. The default value is 'H2O_fl_wt', which is what is returned by BatchFile() core
        calculations. String containing the name of the column corresponding to the H2O
        concentration in the fluid, in wt%.

    CO2_colname: str
        OPTIONAL. The default value is 'CO2_fl_wt', which is what is returned by BatchFile() core
        calculations. String containing the name of the column corresponding to the CO2
        concentration in the fluid, in wt%.

    Returns
    -------
    pandas DataFrame
        Original data passed plus newly calculated values are returned.
    """
    convData = data.copy()

    MPO_H2O_list = []
    MPO_CO2_list = []
    for index, row in convData.iterrows():
        MPO_H2O_list.append(row[H2O_colname] / oxideMass["H2O"])
        MPO_CO2_list.append(row[CO2_colname] / oxideMass["CO2"])

    convData["MPO_H2O"] = MPO_H2O_list
    convData["MPO_CO2"] = MPO_CO2_list
    convData["XH2O_fl"] = convData["MPO_H2O"] / (convData["MPO_H2O"] + convData["MPO_CO2"])
    convData["XCO2_fl"] = convData["MPO_CO2"] / (convData["MPO_H2O"] + convData["MPO_CO2"])

    del convData["MPO_H2O"]
    del convData["MPO_CO2"]

    return convData


class Error(Exception):
    """Base class for exceptions in this module."""
    pass


class InputError(Error):
    """Exception raised for errors in the input.

    Attributes:
        expression -- input expression in which the error occurred
        message -- explanation of the error
    """

    def __init__(self, message):
        self.message = message


class SaturationError(Error):
    """Exception raised for errors thrown when a sample does not reach saturation.

    Attributes:
        expression -- input expression in which the error occurred
        message -- explanation of the error
    """

    def __init__(self, message):
        self.message = message
