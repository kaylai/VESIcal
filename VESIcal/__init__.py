"""
VESIcal

A generalized python library for calculating and plotting various things
related to mixed volatile (H2O-CO2) solubility in silicate melts.
"""

__version__ = "1.1.0"
__author__ = "Kayla Iacovino, Simon Matthews, and Penny Wieser"

# ----------------- IMPORTS ----------------- #
import warnings as w
import pandas as pd

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

# -------------- TURN OFF WARNINGS ------------- #
w.filterwarnings("ignore", message="rubicon.objc.ctypes_patch has only been tested ")
w.filterwarnings("ignore", message="The handle")


# -------------- CALCULATION DEFINITIONS ----- #
class calculate_dissolved_volatiles(
    VESIcal.calculate_classes.calculate_dissolved_volatiles
):
    pass


class calculate_equilibrium_fluid_comp(
    VESIcal.calculate_classes.calculate_equilibrium_fluid_comp
):
    pass


class calculate_isobars_and_isopleths(
    VESIcal.calculate_classes.calculate_isobars_and_isopleths
):
    pass


class calculate_saturation_pressure(
    VESIcal.calculate_classes.calculate_saturation_pressure
):
    pass


class calculate_degassing_path(VESIcal.calculate_classes.calculate_degassing_path):
    pass


class calculate_liquid_density(
    VESIcal.thermo.thermo_calculate_classes.calculate_liquid_density
):
    pass


class calculate_liquid_viscosity(
    VESIcal.thermo.thermo_calculate_classes.calculate_liquid_viscosity
):
    pass


# -------------- ACCESS TO GET_MODEL_NAMES ----- #
def get_model_names(model="all"):
    return VESIcal.models.get_model_names(model=model)


# -------------- PLOTTING DEFINITIONS ----- #
def plot(**kwargs):
    """
    Custom automatic plotting of model calculations in VESIcal.
    Isobars, isopleths, and degassing paths can be plotted. Labels can be
    specified for each. Any combination of isobars, isopleths, and degassing
    paths can be plotted.

    Parameters
    ----------
    isobars: pandas DataFrame or list
        OPTIONAL. DataFrame object containing isobar information as calculated
        by calculate_isobars_and_isopleths. Or a list of DataFrame objects.

    isopleths: pandas DataFrame or list
        OPTIONAL. DataFrame object containing isopleth information as
        calculated by calculate_isobars_and_isopleths. Or a list of DataFrame
        objects.

    degassing_paths: list
        OPTIONAL. List of DataFrames with degassing information as generated
        by calculate_degassing_path().

    custom_H2O: list
        OPTIONAL. List of groups of H2O values to plot as points. For example
        myfile.data['H2O'] is one group of H2O values. Must be passed with
        custom_CO2 and must be same length as custom_CO2.

    custom_CO2: list
        OPTIONAL. List of groups of CO2 values to plot as points.For example
        myfile.data['CO2'] is one group of CO2 values. Must be passed with
        custom_H2O and must be same length as custom_H2O.

    isobar_labels: list
        OPTIONAL. Labels for the plot legend. Default is None, in which case
        each plotted line will be given the generic legend name of
        "Isobars n", with n referring to the nth isobars passed. Isobar
        pressure is given in parentheses. The user can pass their own labels
        as a list of strings. If more than one set of isobars is passed, the
        labels should refer to each set of isobars, not each pressure.

    isopleth_labels: list
        OPTIONAL. Labels for the plot legend. Default is None, in which case
        each plotted isopleth will be given the generic legend name of
        "Isopleth n", with n referring to the nth isopleths passed. Isopleth
        XH2O values are given in parentheses. The user can pass their own
        labels as a list of strings. If more than one set of isopleths is
        passed, the labels should refer to each set of isopleths, not each
        XH2O value.

    degassing_path_labels: list
        OPTIONAL. Labels for the plot legend. Default is None, in which case
        each plotted line will be given the generic legend name of "Pathn",
        with n referring to the nth degassing path passed. The user can pass
        their own labels as a list of strings.

    custom_labels: list
        OPTIONAL. Labels for the plot legend. Default is None, in which case
        each group of custom points will be given the generic legend name of
        "Customn", with n referring to the nth degassing path passed. The user
        can pass their own labels as a list of strings.

    custom_colors: list
        OPTIONAL. Default value is "VESIcal", which uses VESIcal's color ramp.
        A list of color values readable by matplotlib can be passed here if
        custom symbol colors are desired. The length of this list must match
        that of custom_H2O and custom_CO2.

    custom_symbols: list
        OPTIONAL. Default value is None, in which case data are plotted as
        filled circles.. A list of symbol tyles readable by matplotlib can be
        passed here if custom symbol types are desired. The length of this
        list must match that of custom_H2O and custom_CO2.

    markersize: int
        OPTIONAL. Default value is 10. Same as markersize kwarg in matplotlib.
        Any numeric value passed here will set the marker size for
        (custom_H2O, custom_CO2) points.

    figsize: tuple
        OPTIONAL. Default value is (12,8). Sets the matplotlib.pyplot figsize
        value as (x_dimension, y_dimension)

    save_fig: False or str
        OPTIONAL. Default value is False, in which case the figure will not be
        saved. If a string is passed, the figure will be saved with the string
        as the filename. The string must include the file extension.

    extend_isobars_to_zero: bool
        OPTIONAL. If True (default), isobars will be extended to zero, even if
        there is a finite solubility at zero partial pressure.

    smooth_isobars: bool
        OPTIONAL. Default is False. If set to True, isobar data will be fit to
        a polynomial and plotted. If False, the raw input data will be plotted.

    smooth_isopleths: bool
        OPTIONAL. Default is False. If set to True, isopleth data will be fit
        to a polynomial and plotted. If False, the raw input data will be
        plotted.

    Returns
    -------
    fig, axes Matplotlib objects
        fig and axes matploblib objects defining a plot with x-axis as H2O wt%
        in the melt and y-axis as CO2 wt%in the melt. Isobars, or lines of
        constant pressure at which the sample magma composition is saturated,
        and isopleths, or lines of constant fluid composition at which the
        sample magma composition is saturated, are plotted if passed.
        Degassing paths, or the concentration of dissolved H2O and CO2 in a
        melt equilibrated along a path of decreasing pressure, is plotted if
        passed.
    """
    return VESIcal.vplot.plot(**kwargs)


def calib_plot(**kwargs):
    """
    Plots user data and calibration set of any or all models on any x-y plot
    or a total alkalis vs silica (TAS) diagram. TAS diagram boundaries
    provided by tasplot python module, copyright John A Stevenson.

    Parameters
    ----------
    user_data: BatchFile object, Sample object, pandas DataFrame, pandas Series,
        or dict.
        OPTIONAL. Default value is None, in which case only the model
        calibration set is plotted. User provided sample data describing the
        oxide composition of one or more samples. Multiple samples can be
        passed as an BatchFile object or pandas DataFrame. A single sample can
        be passed as a pandas Series.

    model: str or list
        OPTIONAL. Default value is 'all', in which case all model calibration
        datasets will be plotted. 'Mixed' can be used to plot all mixed fluid
        models. String of the name of the model calibration dataset to plot
        (e.g., 'Shishkina'). Multiple models can be plotted by passing them as
        strings within a list (e.g., ['Shishkina', 'Dixon']).

    plot_type: str
        OPTIONAL. Default value is 'TAS', which returns a total alkali vs
        silica (TAS) diagram. Any two oxides can be plotted as an x-y plot by
        setting plot_type='xy' and specifying x- and y-axis oxides, e.g.,
        x='SiO2', y='Al2O3'.

    zoom: str or list
        OPTIONAL. Default value is None in which case axes will be set to the
        default of 35<x<100 wt% and 0<y<25 wt% for TAS type plots and the best
        values to show the data for xy type plots. Can pass "user_data" to
        plot the figure where the x and y axes are scaled down to zoom in and
        only show the region surrounding the user_data. A list of tuples may
        be passed to manually specify x and y limits. Pass in data as
        [(x_min, x_max), (y_min, y_max)]. For example, the default limits here
        would be passed in as [(35,100), (0,25)].

    figsize: tuple
        OPTIONAL. Default value is (17,8). Sets the matplotlib.pyplot figsize
        value as (x_dimension, y_dimension).

    legend: bool
        OPTIONAL. Default value is True. Can be set to False in which case the
        legend will not be displayed.

    save_fig: False or str
        OPTIONAL. Default value is False, in which case the figure will not be
        saved. If a string is passed, the figure will be saved with the string
        as the filename. The string must include the file extension.

    Returns
    -------
    matplotlib object
    """
    return VESIcal.vplot.calib_plot(**kwargs)


def show():
    """
    Inherits from vplot.show()
    """
    return VESIcal.vplot.show()


# -------------- SAMPLE PROCESSING ---------- #
class Sample(VESIcal.sample_class.Sample):
    """The sample class stores compositional information for samples, and contains methods for
    normalization and other compositional calculations.

    The composition is stored as wtpt. If the composition is provided as wtpt, no
    normalization will be applied. If the composition is supplied as mols, the composition
    will be normalized to 100 wt%.

    Parameters
    ----------
    composition     dict or pandas.Series
        The composition of the sample in the format specified by the composition_type
        parameter. Default is oxides in wtpt.

    units     str
        Specifies the units and type of compositional information passed in the composition
        parameter. Choose from 'wtpt_oxides', 'mol_oxides', 'mol_cations'.

    default_normalization:     None or str
        The type of normalization to apply to the data by default. One of:
            - None (no normalization)
            - 'standard' (default): Normalizes an input composition to 100%.
            - 'fixedvolatiles': Normalizes major element oxides to 100 wt%, including
              volatiles. The volatile wt% will remain fixed, whilst the other major element
              oxides are reduced proportionally so that the total is 100 wt%.
            - 'additionalvolatiles': Normalises major element oxide wt% to 100%, assuming it
              is volatile-free. If H2O or CO2 are passed to the function, their un-normalized
              values will be retained in addition to the normalized non-volatile oxides,
              summing to >100%.

    default_units     str
        The type of composition to return by default, one of:
        - wtpt_oxides (default)
        - mol_oxides
        - mol_cations
        - mol_singleO
    """

    pass


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


# -------------- BATCH PROCESSING ------------ #
class BatchFile(VESIcal.batchmodel.BatchFile):
    """A batch file with sample names and oxide compositions

    Attributes
    ----------
        filename: str
            Path to the batch file, e.g., "my_file.xlsx". This always needs to
            be passed, even if the user is passing a pandas DataFrame rather
            than an batch file. If passing a DataFrame, filename should be set
            to None. File can be excel file (.xlsx) or .csv.

        sheet_name: str
            OPTIONAL. For Excel files. Default value is 0 which gets the first
            sheet in the batch spreadsheet file. This implements the pandas.
            read_excel() sheet_name parameter. But functionality to read in
            more than one sheet at a time (e.g.,
            pandas.read_excel(sheet_name=None)) is not yet imlpemented in
            VESIcal. From the pandas 1.0.4 documentation:

            Available cases:
            - Defaults to 0: 1st sheet as a DataFrame
            - 1: 2nd sheet as a DataFrame
            - "Sheet1": Load sheet with name “Sheet1”

        file_type: str
            OPTIONAL. Default is 'excel', which denotes that passed file has
            extension .xlsx. Other option is 'csv', which denotes that the
            passed file has extension .csv.

        units: str
            OPTIONAL. Default is 'wtpt_oxides'. String defining whether the
            oxide composition is given in wt percent ("wtpt_oxides", which is
            the default), mole oxides (mol_oxides) or mole cations
            (mol_cations).

        default_normalization:     None or str
            The type of normalization to apply to the data by default. One of:
            - None (no normalization)
            - 'standard' (default): Normalizes an input composition to 100%.
            - 'fixedvolatiles': Normalizes major element oxides to 100 wt%,
               including volatiles. The volatile wt% will remain fixed, whilst
               the other major element oxides are reduced proportionally so
               that the total is 100 wt%.
            - 'additionalvolatiles': Normalises major element oxide wt% to
              100%, assuming it is volatile-free. If H2O or CO2 are passed to
              the function, their un-normalized values will be retained in
              addition to the normalized non-volatile oxides, summing to >100%.

        default_units     str
            The type of composition to return by default, one of:
            - wtpt_oxides (default)
            - mol_oxides
            - mol_cations

        label: str
            OPTIONAL. Default is 'Label'. Name of the column within the passed
            file referring to sample names.

        dataframe: pandas DataFrame
            OPTIONAL. Default is None in which case this argument is ignored.
            This argument is used when the user wishes to turn a pandas
            DataFrame into an BatchFile object, for example when user data is
            already in python rather than being imported from a file. In this
            case set `dataframe` equal to the dataframe object being passed in.
            If using this option, pass None to filename.
    """

    pass


def BatchFile_from_DataFrame(dataframe, **kwargs):
    """
    Provides method for creating a BatchFile object from an existing pandas
    DataFrame. Inherits from batchfile.BatchFile().
    """
    return VESIcal.batchmodel.BatchFile_from_DataFrame(dataframe, **kwargs)


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
