"""
VESIcal

A generalized python library for calculating and plotting various things related to mixed volatile (H2O-CO2) solubility in silicate melts.
"""

__version__ = "0.9.15"
__author__ = 'Kayla Iacovino, Simon Matthews, and Penny Wieser'

# -------------- TURN OFF WARNINGS ------------- #
import warnings as w
w.filterwarnings("ignore", message="rubicon.objc.ctypes_patch has only been tested ")
w.filterwarnings("ignore", message="The handle")

# ----------------- IMPORTS ----------------- #
import pandas as pd

# import matplotlib.font_manager as font_manager
# from cycler import cycler
# from scipy.optimize import minimize

import VESIcal.core
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

# -------------- CALCULATION DEFINITIONS ----- #
class calculate_dissolved_volatiles(VESIcal.calculate_classes.calculate_dissolved_volatiles):
    pass

class calculate_equilibrium_fluid_comp(VESIcal.calculate_classes.calculate_equilibrium_fluid_comp):
    pass

class calculate_isobars_and_isopleths(VESIcal.calculate_classes.calculate_isobars_and_isopleths):
    pass

class calculate_saturation_pressure(VESIcal.calculate_classes.calculate_saturation_pressure):
    pass

class calculate_degassing_path(VESIcal.calculate_classes.calculate_degassing_path):
    pass

# -------------- ACCESS TO GET_MODEL_NAMES ----- #
def get_model_names(model='all'):
  return VESIcal.models.get_model_names(model=model)

# -------------- PLOTTING DEFINITIONS ----- #
def plot(**kwargs):
    """
    Inherits from vplot.plot().
    """
    return VESIcal.vplot.plot(**kwargs)

def calib_plot(**kwargs):
    """
    Inherits from vplot.calib_plot()
    """
    return VESIcal.vplot.calib_plot(**kwargs)

def show():
    """
    Inherits from vplot.show()
    """
    return VESIcal.vplot.show()


# -------------- SAMPLE PROCESSING ---------- #
class Sample(VESIcal.sample_class.Sample):
    """
    Provides methods for working with rock compositions. Inherits from sample_class.Sample().
    """
    pass

def get_oxides(sample):
  """
  Returns a sample composition with only compositional oxide data, removing any extranneous data.
  Useful when passing a self-defined sample (e.g. dict or pandas Series) to a some VESIcal function.

  Parameters
  ----------
  sample: pandas Series or dictionary
    A sample composition plus other sample information

  Returns
  -------
  Same type as passed sample (pandas Series or dictionary)
    Sample composition with extranneous information removed.
  """

  clean = {oxide:  sample[oxide] for oxide in VESIcal.core.oxides}

  if isinstance(sample, dict):
    return clean
  if isinstance(sample, pd.core.series.Series):
    return pd.Series(clean)

# -------------- BATCH PROCESSING ------------ #
class BatchFile(VESIcal.batchmodel.BatchFile):
    """
    Provides methods for batch processing files of data. Inherits from batchfile.BatchFile().
    """
    pass

def BatchFile_from_DataFrame(dataframe, **kwargs):
    """
    Provides method for creating a BatchFile object from an existing pandas DataFrame. Inherits from batchfile.BatchFile().
    """
    return VESIcal.batchmodel.BatchFile_from_DataFrame(dataframe, **kwargs)



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
