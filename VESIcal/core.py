"""
Single Source of Truth pattern for oxide-cation compositional factors.
All derived data structures are computed from a single master dictionary.
"""

# Master data structure - single source of truth
oxide_data = {
    'SiO2': {
        'mass': 60.083,
        'cation': 'Si',
        'cation_num': 1,
        'oxygen_num': 2,
        'cation_charge': 4,
        'cation_mass': 28.085,
        'groups': ['magmasat', 'anhydrous']
    },
    'TiO2': {
        'mass': 79.867,
        'cation': 'Ti',
        'cation_num': 1,
        'oxygen_num': 2,
        'cation_charge': 4,
        'cation_mass': 47.867,
        'groups': ['magmasat', 'anhydrous']
    },
    'Al2O3': {
        'mass': 101.961,
        'cation': 'Al',
        'cation_num': 2,
        'oxygen_num': 3,
        'cation_charge': 3,
        'cation_mass': 26.982,
        'groups': ['magmasat', 'anhydrous']
    },
    'Cr2O3': {
        'mass': 151.992,
        'cation': 'Cr',
        'cation_num': 2,
        'oxygen_num': 3,
        'cation_charge': 3,
        'cation_mass': 51.996,
        'groups': ['magmasat', 'anhydrous']
    },
    'FeO': {
        'mass': 71.844,
        'cation': 'Fe',
        'cation_num': 1,
        'oxygen_num': 1,
        'cation_charge': 2,
        'cation_mass': 55.845,
        'groups': ['magmasat', 'anhydrous']
    },
    'Fe2O3': {
        'mass': 159.687,
        'cation': 'Fe3',
        'cation_num': 2,
        'oxygen_num': 3,
        'cation_charge': 3,
        'cation_mass': 55.845,
        'groups': ['magmasat', 'anhydrous']
    },
    'MnO': {
        'mass': 70.937,
        'cation': 'Mn',
        'cation_num': 1,
        'oxygen_num': 1,
        'cation_charge': 2,
        'cation_mass': 54.938,
        'groups': ['magmasat', 'anhydrous']
    },
    'MgO': {
        'mass': 40.304,
        'cation': 'Mg',
        'cation_num': 1,
        'oxygen_num': 1,
        'cation_charge': 2,
        'cation_mass': 24.305,
        'groups': ['magmasat', 'anhydrous']
    },
    'CaO': {
        'mass': 56.077,
        'cation': 'Ca',
        'cation_num': 1,
        'oxygen_num': 1,
        'cation_charge': 2,
        'cation_mass': 40.078,
        'groups': ['magmasat', 'anhydrous']
    },
    'Na2O': {
        'mass': 61.979,
        'cation': 'Na',
        'cation_num': 2,
        'oxygen_num': 1,
        'cation_charge': 1,
        'cation_mass': 22.990,
        'groups': ['magmasat', 'anhydrous']
    },
    'K2O': {
        'mass': 94.195,
        'cation': 'K',
        'cation_num': 2,
        'oxygen_num': 1,
        'cation_charge': 1,
        'cation_mass': 39.098,
        'groups': ['magmasat', 'anhydrous']
    },
    'P2O5': {
        'mass': 141.943,
        'cation': 'P',
        'cation_num': 2,
        'oxygen_num': 5,
        'cation_charge': 5,
        'cation_mass': 30.974,
        'groups': ['magmasat', 'anhydrous']
    },
    'NiO': {
        'mass': 74.692,
        'cation': 'Ni',
        'cation_num': 1,
        'oxygen_num': 1,
        'cation_charge': 2,
        'cation_mass': 58.693,
        'groups': ['magmasat', 'anhydrous']
    },
    'CoO': {
        'mass': 44.01,
        'cation': 'Co',
        'cation_num': 1,
        'oxygen_num': 1,
        'cation_charge': 2,
        'cation_mass': 28.01,
        'groups': ['magmasat', 'anhydrous']
    },
    'H2O': {
        'mass': 18.02,
        'cation': 'H',
        'cation_num': 2,
        'oxygen_num': 1,
        'cation_charge': 1,
        'cation_mass': 1.01,
        'groups': ['magmasat', 'volatile']
    },
    'CO2': {
        'mass': 44.01,
        'cation': 'C',
        'cation_num': 1,
        'oxygen_num': 2,
        'cation_charge': 4,
        'cation_mass': 12.011,
        'groups': ['magmasat', 'volatile']
    },
    'F2O': {
        'mass': 37.997,
        'cation': 'F',
        'cation_num': 2,
        'oxygen_num': 1,
        'cation_charge': 1,
        'cation_mass': 18.998,
        'groups': ['anhydrous']
    }
}

# Derived data structures - all computed from oxide_data
# Lists
oxides = list(oxide_data.keys())
cations = list(set(data['cation'] for data in oxide_data.values()))
magmasat_oxides = [oxide for oxide, data in oxide_data.items() if 'magmasat' in data['groups']]
anhydrous_oxides = [oxide for oxide, data in oxide_data.items() if 'anhydrous' in data['groups']]
volatiles = [oxide for oxide, data in oxide_data.items() if 'volatile' in data['groups']]

# Property dictionaries
oxideMass = {oxide: data['mass'] for oxide, data in oxide_data.items()}
CationNum = {oxide: data['cation_num'] for oxide, data in oxide_data.items()}
OxygenNum = {oxide: data['oxygen_num'] for oxide, data in oxide_data.items()}
CationCharge = {oxide: data['cation_charge'] for oxide, data in oxide_data.items()}
CationMass = {oxide: data['cation_mass'] for oxide, data in oxide_data.items()}

# Conversion mappings
oxides_to_cations = {oxide: data['cation'] for oxide, data in oxide_data.items()}
cations_to_oxides = {data['cation']: oxide for oxide, data in oxide_data.items()}


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
