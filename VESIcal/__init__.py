"""
VESIcal

A generalized python library for calculating and plotting various things related to mixed volatile (H2O-CO2) solubility in silicate melts.
"""

__version__ = "0.1.7"
__author__ = 'Kayla Iacovino, Simon Matthews, and Penny Wieser'

# -------------- TURN OFF WARNINGS ------------- #
import warnings as w
w.filterwarnings("ignore", message="rubicon.objc.ctypes_patch has only been tested ")
w.filterwarnings("ignore", message="The handle")

# ----------------- IMPORTS ----------------- #
import pandas as pd
import numpy as np
from scipy.optimize import root_scalar
from scipy.optimize import root
from abc import abstractmethod
import sys
import sympy
from copy import copy

# import matplotlib.font_manager as font_manager
# from cycler import cycler
# from scipy.optimize import minimize

from VESIcal.core import *
import VESIcal.activity_models
import VESIcal.batchfile
import VESIcal.batchmodel
import VESIcal.calculate_classes
import VESIcal.calibration_checks
import VESIcal.fugacity_models
import VESIcal.models
import VESIcal.sample_class
import VESIcal.vplot

# -------------- CALCULATION DEFINITIONS ----- #
class calculate_dissolved_volatiles(calculate_classes.calculate_dissolved_volatiles):
    pass

class calculate_equilibrium_fluid_comp(calculate_classes.calculate_equilibrium_fluid_comp):
    pass

class calculate_isobars_and_isopleths(calculate_classes.calculate_isobars_and_isopleths):
    pass

class calculate_saturation_pressure(calculate_classes.calculate_saturation_pressure):
    pass

class calculate_degassing_path(calculate_classes.calculate_degassing_path):
    pass

# -------------- ACCESS TO GET_MODEL_NAMES ----- #
def get_model_names(model='all'):
  return models.get_model_names(model=model)

# -------------- PLOTTING DEFINITIONS ----- #
def plot(**kwargs):
    return vplot.plot(**kwargs)

def show():
    return vplot.show()


# -------------- SAMPLE PROCESSING ---------- #
class Sample(sample_class.Sample):
    """
    Provides methods for working with rock compositions.
    """
    pass

# -------------- BATCH PROCESSING ------------ #
class BatchFile(batchmodel.BatchFile):
    """
    Provides methods for batch processing files of data.
    """
    pass

def BatchFile_from_DataFrame(dataframe, **kwargs):
    """
    Provides method for creating a BatchFile object from an existing pandas DataFrame.
    """
    return batchmodel.BatchFile_from_DataFrame(dataframe, **kwargs)

def test_BatchFile(filename=None):
    """
    A test routine that takes in an excel file and runs multiple calculations on that file.

    Parameters
    ----------
    filename: string
        OPTIONAL. Name with extension of file to be tested. If no file is passed, data will
        be generated automatially.
    """
    pd.set_option('display.max_rows', None, 'display.max_columns', None)

    print("\n================================\n= MAGMASATPLUS BATCH TESTING ROUTINE =\n================================")

    print("\n This routine will check that key methods run using typical values of variables or given user data. \
The routine does not check that the results are correct (though this may be obvious from the outputs),\
nor does it check every possible iteration of methods and input types. It will check that an update hasn't \
COMPLETELY broken the module.")

    #SET UP THE DATA
    if filename == None:
        fakedata = pd.DataFrame({'Label': ['Samp1', 'Samp2', 'Samp3'],
                                'SiO2': [47.95, 69.02, 55.4],
                                'TiO2': [1.67, 0.78, 1.01],
                                'Al2O3': [17.32, 15.2, 10.1],
                                'FeO': [10.24, 4.2, 8.9],
                                'Fe2O3': [0.1, 0.2, 0.5],
                                'MgO': [5.76, 0.3, 3.0],
                                'CaO': [10.93, 12.99, 10.9],
                                'Na2O': [3.45, 5.67, 4.22],
                                'K2O': [1.99, 3.2, 3.2],
                                'P2O5': [0.51, 0.2, 0.5],
                                'MnO': [0.1, 0.15, 0.1],
                                'CO2': [0.8, 0.2, 0.3],
                                'H2O': [4.0, 6.0, 2.0]})
        fakedata = fakedata.set_index('Label')
        myfile = BatchFile(filename=None, dataframe=fakedata)
    else:
        myfile = BatchFile(filename)

    test_temperature = 1000
    test_pressure = 2000
    test_X_fluid = 1

    #CALCULTE SATURATION PRESSURE
    print("\n Saturation Pressures:")
    print(" =====================")
    satPs = myfile.calculate_saturation_pressure(temperature=test_temperature)
    print(satPs)

    #CALCULATE DISSOLVED VOLATILES
    print("\n Dissolved Volatile Concentrations:")
    print(" ==================================")
    dissolved = myfile.calculate_dissolved_volatiles(temperature=test_temperature, pressure=test_pressure, X_fluid=test_X_fluid)
    print(dissolved)

    #CALCULATE EQUILIBRIUM FLUID COMPOSITIONS
    print("\n Equilibrium Fluid Compositions:")
    print(" ===============================")
    eqfluid = myfile.calculate_equilibrium_fluid_comp(temperature=test_temperature, pressure=test_pressure)
    print(eqfluid)

    print("\nTesting routine complete.\n")




"""
                        ,,,                                     .*****
                       ,***,*                                  ,* *****
                      ,***,,,,*                              ,* ,*******
                     .****,,,.,**                          ,,*.,*********
                     *****,,,,.****                     ,,,,*.,********* *
                     *,**,,,,,,,******                ,,,,,* ,,********,.*                                     .
                    .*,*,, //,,.,*******            ,,,,,**.//,*******,,.,                                   //
                    ,**,, //// ,,******** ,,,,,,,, **,,,*** ((,,,,,,,,,,,*                                  //*####
                    .,* , ///((,,,******. ,,,,,,, *********(((,.#### ,,.,*                                   ///#######
                     ,,** .((((,******** ,,,*,,, *********.(((#####/.,,**    .,,                              ///#*#######
                     ,,**(((((.******* ,,**** **.*********.#(####(( .,**.******                                ///##/#######
                      ,,***(.( ****** ,*******#*********** #*###( ,,** ******,                              /////(###########
                    ,,, ,**,,,,,,. ,,,**,**.### ****.************.  .******,                        .,,,,,****, ##############
                    .,,,,,,.,,,,,,,,,* *,#### ************# *************,                    ,,,,********,***/#######(#######
                      ,,,,,, /(   (#,*##### *****###/ #### ************                 ,,,,,,**********,.**.#######(#########(
                        ,,,,, ((((((# ### * * #..######## *****************.##        ,,,,,****************,####. (## #######
                      ,,,,,,,, *(((( *###.****,#  ##### *********** #### ###         ,,,,,*********************,######*#######
                     ./ ,,,,,,,   ,****########*****,       .**########### ,         ,,,*********************** (##./##.##/#*
                       /////*,  ,,  *****,*,******   ,,.,((( ########### ,****      .,,.*****************,######  //#/////,
                        //////  /(,. *****###.****  *((,*&((########## *****       .****,************ ###########(/////*  //
                          /////,*#((/,(((*########* *((#(*##########********** .,   *****,*******    .###########//////,
                             ///// ,* (((((####### %#  ##########,**********,************* **  **** /##### ///////*/ /
                                ./////((((((################# ###***************************** ,,,,,  ///////*  ////
                              ///// .////(((############ .(((###**#.******,********************     .    //////.
                                * //// ////(######## ///(#########*******,****.*****************
                                 *//////.////### */////(#########.*****,,****** ****************
                               /////////////////////##/##########* .* ,,*******,,************ **.
                                   ,///////////////################(,,,********,*  ************ *,*
                                        , ////////################,,,,*********  ..,******,** **.***
                                        , //////##/############# ,,,,,******* *****.,,*.,,, *********,,
                                       ,,,//////(#/############  ,,,,,************ ((( ,, .***********
                                       ,,///  .  ///###########((,,,,,*********.* //   ,  ,** *  *,
                                       ,,,,,,,,,,  //#########,(( ,,,,******** *       ,, ,,**   *,
                                       ,,,,,,,,,.   *//###( #,/// ,,,,******* *        .,,, ,,***..
                                       .,,,***,,      ///# *,//*,,.,,,****** ,           ,,,,,,.,,,,
                                       ,,,*****         *  // ,,,. ,,,***** *              ,,,,,,,,*.
                                       , ,*****              ,, ,, *,,******                 ,,,,,,**
                                       ,,,.***               ,,,,, ** *.****                  , ,,,,,
                                       ,,,****              ,,,,,  ********                    ,,,,,,
                                       ,,****               ,,,,,  *******                      ,,,,,
                                       ,,***.             .,,,,    ******                       ,,,,*
                                       ,****                       *****,                       ,,**
                                       .***                       ,*****                        , ,
                                                                  *****
                                                                 ,****


                           ***** *      **        ***** **          *******            *****  *                       ***
                        ******  *    *****     ******  **** *     *       ***       ******  *                          ***
                       **   *  *       *****  **   *  * ****     *         **      **   *  *                            **
                      *    *  **       * **  *    *  *   **      **        *      *    *  *                             **
                          *  ***      *          *  *             ***                 *  *                              **
                         **   **      *         ** **            ** ***              ** **         ****       ****      **
                         **   **      *         ** **             *** ***            ** **        * ***  *   * ***  *   **
                         **   **     *          ** ******           *** ***        **** **       *   ****   *   ****    **
                         **   **     *          ** *****              *** ***     * *** **      **         **    **     **
                         **   **     *          ** **                   ** ***       ** **      **         **    **     **
                          **  **    *           *  **                    ** **  **   ** **      **         **    **     **
                           ** *     *              *                      * *  ***   *  *       **         **    **     **
                            ***     *          ****         *   ***        *    ***    *        ***     *  **    **     **
                             *******          *  ***********   *  *********      ******          *******    ***** **    *** *
                               ***           *     ******     *     *****          ***            *****      ***   **    ***
                                             *                *
                                              **               **

logo courtesy Twai (http://www.twitter.com/_twai)
"""

def WhatDoesTheFoxSay():
    print("https://www.youtube.com/watch?v=jofNR_WkoCE")
