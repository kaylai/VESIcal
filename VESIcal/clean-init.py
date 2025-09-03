"""
VESIcal

A Python library for calculating and visualizing mixed volatile (H2O-CO2) solubility
in silicate melts, including saturation pressures, degassing paths, and related
thermodynamic and petrologic processes.

Documentation:
https://vesical.readthedocs.io

Peer-reviewed manuscripts:
Iacovino K., Matthews S., Wieser P.E., Moore G.M., and Begue F. (2021) VESIcal Part I: An
open-source thermodynamic model engine for mixed volatile (H2O-CO2) solubility in silicate melts,
Earth and Space Science, 8, e2020EA001584. https://doi.org/10.1029/2020EA001584.

Wieser P.E., Iacovino K., Matthews S., Moore G.M., Allison C.M. (2022) VESIcal Part II: A critical
approach to volatile solubility modelling using an open-source Python3 engine, Earth and Space
Science. https://doi.org/10.1029/2021EA001932

Main API
--------
    Sample
    BatchFile
    calculate_dissolved_volatiles
    calculate_saturation_pressure
    calculate_degassing_path
    calculate_liquid_density
    plot
    calib_plot
    save_to_file
"""

__version__ = "1.2.10"
__author__ = "Kayla Iacovino, Simon Matthews, and Penny Wieser"

# ----------------- WARNINGS CLEANUP ----------------- #
import warnings as w
w.filterwarnings("ignore", message="rubicon.objc.ctypes_patch has only been tested ")
w.filterwarnings("ignore", message="The handle")
w.filterwarnings("ignore", message="duanH2ODriver(b):")
w.filterwarnings("ignore", message="duanDriver-2:")
w.filterwarnings("ignore", message="Error for element")

# ----------------- IMPORTS ----------------- #
import VESIcal.core
from VESIcal.core import oxides, anhydrous_oxides, volatiles  # noqa F401
from VESIcal.core import fluid_molfrac_to_wt, fluid_wt_to_molfrac  # noqa F401
import VESIcal.activity_models
import VESIcal.batchfile
import VESIcal.batchmodel
import VESIcal.calculate_classes
import VESIcal.calibration_checks
import VESIcal.calibrations
import VESIcal.fugacity_models
import VESIcal.models
import VESIcal.sample_class
import VESIcal.vplot
import VESIcal.thermo

from VESIcal.sample_class import Sample
from VESIcal.batchmodel import BatchFile, BatchFile_from_DataFrame

from VESIcal.calculate_classes import (
    calculate_dissolved_volatiles,
    calculate_equilibrium_fluid_comp,
    calculate_isobars_and_isopleths,
    calculate_saturation_pressure,
    calculate_degassing_path,
)

from VESIcal.models import get_model_names

from VESIcal.thermo.thermo_calculate_classes import (
    calculate_liquid_density,
    calculate_liquid_viscosity,
)

from VESIcal import utils
from VESIcal.utils.save_util import save_results
from VESIcal.vplot import plot, calib_plot, show

import pandas as pd
import importlib as _importlib

__all__ = [
    # Core classes
    "Sample",
    "BatchFile",

    # Calculations
    "calculate_dissolved_volatiles",
    "calculate_equilibrium_fluid_comp",
    "calculate_isobars_and_isopleths",
    "calculate_saturation_pressure",
    "calculate_degassing_path",
    "calculate_liquid_density",
    "calculate_liquid_viscosity",

    # Models
    "models",
    
    # Utilities
    "BatchFile_from_DataFrame",
    "get_model_names",
    "get_oxides",
    "save_results",
    "utils",

    # Plotting
    "plot",
    "calib_plot",
    "show",
]

# ----------------- PUBLIC API RESOLUTION ----------------- #

def __dir__():
    return __all__


def __getattr__(name):
    if name in {
        "calculate_dissolved_volatiles",
        "calculate_equilibrium_fluid_comp",
        "calculate_isobars_and_isopleths",
        "calculate_saturation_pressure",
        "calculate_degassing_path",
    }:
        mod = _importlib.import_module("VESIcal.calculate_classes")
        return getattr(mod, name)

    elif name in {
        "calculate_liquid_density",
        "calculate_liquid_viscosity",
    }:
        mod = _importlib.import_module("VESIcal.thermo.thermo_calculate_classes")
        return getattr(mod, name)

    elif name in {
        "plot",
        "calib_plot",
        "show",
    }:
        mod = _importlib.import_module("VESIcal.vplot")
        return getattr(mod, name)

    elif name == "get_model_names":
        mod = _importlib.import_module("VESIcal.models")
        return getattr(mod, name)

    elif name == "save_results":
        mod = _importlib.import_module("VESIcal.utils.save_util")
        return getattr(mod, name)
    
    elif name == "utils":
        mod = _importlib.import_module("VESIcal.utils")
        return getattr(mod, name)

    elif name == "get_oxides":
        mod = _importlib.import_module("VESIcal.get_oxides")
        return getattr(mod, name)

    elif name == "Sample":
        mod = _importlib.import_module("VESIcal.sample_class")
        return getattr(mod, name)

    elif name == "BatchFile":
        mod = _importlib.import_module("VESIcal.batchmodel")
        return getattr(mod, name)

    elif name == "BatchFile_from_DataFrame":
        mod = _importlib.import_module("VESIcal.batchmodel")
        return getattr(mod, name)

    raise AttributeError(f"Module 'VESIcal' has no attribute '{name}'")

def get_oxides(sample):
    """
    Returns a sample composition with only compositional oxide data, removing
    any extranneous data. Useful when passing a self-defined sample (e.g. dict
    or pandas Series) to a some VESIcal function.

    Parameters
    ----------
    sample: pandas Series or dictionary
    A sample composition plus other sample information

    Returns
    -------
    Same type as passed sample (pandas Series or dictionary)
    Sample composition with extranneous information removed.
    """

    clean = {oxide: sample[oxide] for oxide in VESIcal.core.oxides}

    if isinstance(sample, dict):
        return clean
    if isinstance(sample, pd.core.series.Series):
        return pd.Series(clean)

"""
                        ,,,                                     .*****
                       ,***,*                                  ,* *****
                      ,***,,,,*                              ,* ,*******
                     .****,,,.,**                          ,,*.,*********
                     *****,,,,.****                     ,,,,*.,********* *
                     *,**,,,,,,,******                ,,,,,* ,,********,.*                                     .                # noqa: E501
                    .*,*,, //,,.,*******            ,,,,,**.//,*******,,.,                                   //                 # noqa: E501
                    ,**,, //// ,,******** ,,,,,,,, **,,,*** ((,,,,,,,,,,,*                                  //*####             # noqa: E501
                    .,* , ///((,,,******. ,,,,,,, *********(((,.#### ,,.,*                                   ///#######         # noqa: E501
                     ,,** .((((,******** ,,,*,,, *********.(((#####/.,,**    .,,                              ///#*#######      # noqa: E501
                     ,,**(((((.******* ,,**** **.*********.#(####(( .,**.******                                ///##/#######    # noqa: E501
                      ,,***(.( ****** ,*******#*********** #*###( ,,** ******,                              /////(###########   # noqa: E501
                    ,,, ,**,,,,,,. ,,,**,**.### ****.************.  .******,                        .,,,,,****, ##############  # noqa: E501
                    .,,,,,,.,,,,,,,,,* *,#### ************# *************,                    ,,,,********,***/#######(#######  # noqa: E501
                      ,,,,,, /(   (#,*##### *****###/ #### ************                 ,,,,,,**********,.**.#######(#########( # noqa: E501
                        ,,,,, ((((((# ### * * #..######## *****************.##        ,,,,,****************,####. (## #######   # noqa: E501
                      ,,,,,,,, *(((( *###.****,#  ##### *********** #### ###         ,,,,,*********************,######*#######  # noqa: E501
                     ./ ,,,,,,,   ,****########*****,       .**########### ,         ,,,*********************** (##./##.##/#*   # noqa: E501
                       /////*,  ,,  *****,*,******   ,,.,((( ########### ,****      .,,.*****************,######  //#/////,     # noqa: E501
                        //////  /(,. *****###.****  *((,*&((########## *****       .****,************ ###########(/////*  //    # noqa: E501
                          /////,*#((/,(((*########* *((#(*##########********** .,   *****,*******    .###########//////,        # noqa: E501
                             ///// ,* (((((####### %#  ##########,**********,************* **  **** /##### ///////*/ /          # noqa: E501
                                ./////((((((################# ###***************************** ,,,,,  ///////*  ////            # noqa: E501
                              ///// .////(((############ .(((###**#.******,********************     .    //////.                # noqa: E501
                                * //// ////(######## ///(#########*******,****.*****************                                # noqa: E501
                                 *//////.////### */////(#########.*****,,****** ****************                                # noqa: E501
                               /////////////////////##/##########* .* ,,*******,,************ **.                               # noqa: E501
                                   ,///////////////################(,,,********,*  ************ *,*                             # noqa: E501
                                        , ////////################,,,,*********  ..,******,** **.***                            # noqa: E501
                                        , //////##/############# ,,,,,******* *****.,,*.,,, *********,,                         # noqa: E501
                                       ,,,//////(#/############  ,,,,,************ ((( ,, .***********                          # noqa: E501
                                       ,,///  .  ///###########((,,,,,*********.* //   ,  ,** *  *,                             # noqa: E501
                                       ,,,,,,,,,,  //#########,(( ,,,,******** *       ,, ,,**   *,                             # noqa: E501
                                       ,,,,,,,,,.   *//###( #,/// ,,,,******* *        .,,, ,,***..                             # noqa: E501
                                       .,,,***,,      ///# *,//*,,.,,,****** ,           ,,,,,,.,,,,                            # noqa: E501
                                       ,,,*****         *  // ,,,. ,,,***** *              ,,,,,,,,*.                           # noqa: E501
                                       , ,*****              ,, ,, *,,******                 ,,,,,,**                           # noqa: E501
                                       ,,,.***               ,,,,, ** *.****                  , ,,,,,                           # noqa: E501
                                       ,,,****              ,,,,,  ********                    ,,,,,,                           # noqa: E501
                                       ,,****               ,,,,,  *******                      ,,,,,                           # noqa: E501
                                       ,,***.             .,,,,    ******                       ,,,,*                           # noqa: E501
                                       ,****                       *****,                       ,,**                            # noqa: E501
                                       .***                       ,*****                        , ,                             # noqa: E501
                                                                  *****
                                                                 ,****


                           ***** *      **        ***** **          *******            *****  *                       ***       # noqa: E501
                        ******  *    *****     ******  **** *     *       ***       ******  *                          ***      # noqa: E501
                       **   *  *       *****  **   *  * ****     *         **      **   *  *                            **      # noqa: E501
                      *    *  **       * **  *    *  *   **      **        *      *    *  *                             **      # noqa: E501
                          *  ***      *          *  *             ***                 *  *                              **      # noqa: E501
                         **   **      *         ** **            ** ***              ** **         ****       ****      **      # noqa: E501
                         **   **      *         ** **             *** ***            ** **        * ***  *   * ***  *   **      # noqa: E501
                         **   **     *          ** ******           *** ***        **** **       *   ****   *   ****    **      # noqa: E501
                         **   **     *          ** *****              *** ***     * *** **      **         **    **     **      # noqa: E501
                         **   **     *          ** **                   ** ***       ** **      **         **    **     **      # noqa: E501
                          **  **    *           *  **                    ** **  **   ** **      **         **    **     **      # noqa: E501
                           ** *     *              *                      * *  ***   *  *       **         **    **     **      # noqa: E501
                            ***     *          ****         *   ***        *    ***    *        ***     *  **    **     **      # noqa: E501
                             *******          *  ***********   *  *********      ******          *******    ***** **    *** *   # noqa: E501
                               ***           *     ******     *     *****          ***            *****      ***   **    ***    # noqa: E501
                                             *                *
                                              **               **

logo courtesy Twai (http://www.twitter.com/_twai)
"""

def WhatDoesTheFoxSay():
    print("https://www.youtube.com/watch?v=jofNR_WkoCE")
