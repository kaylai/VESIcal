from VESIcal import core
from VESIcal import calibrations
from VESIcal.tasplot import add_LeMaitre_fields

import pandas as pd
import numpy as np
import warnings as w
import matplotlib as mpl
import matplotlib.pyplot as plt


# ---------- DEFINE CUSTOM PLOTTING FORMATTING ------------ #
style = "seaborn-colorblind"
plt.style.use(style)
plt.rcParams["mathtext.default"] = "regular"
plt.rcParams["mathtext.fontset"] = "dejavusans"
mpl.rcParams['patch.linewidth'] = 1
mpl.rcParams['axes.linewidth'] = 1
plt.rcParams['axes.titlesize'] = 20
plt.rcParams['axes.labelsize'] = 18
plt.rcParams['xtick.labelsize'] = 14
plt.rcParams['ytick.labelsize'] = 14
plt.rcParams['legend.fontsize'] = 14
mpl.rcParams['lines.markersize'] = 10

# Define color cycler based on plot style set here
# get style formatting set by plt.style.use():
the_rc = plt.style.library[style]
# list of colors by hex code:
color_list = the_rc['axes.prop_cycle'].by_key()['color'] * 10
color_cyler = the_rc['axes.prop_cycle']  # get the cycler


# ----------- MAGMASAT PLOTTING FUNCTIONS ----------- #
def smooth_isobars_and_isopleths(isobars=None, isopleths=None):
    """
    Takes in a dataframe with calculated isobar and isopleth information
    (e.g., output from calculate_isobars_and_isopleths) and smooths the data
    for plotting.

    Parameters
    ----------
    isobars: pandas DataFrame
        OPTIONAL. DataFrame object containing isobar information as calculated
        by calculate_isobars_and_isopleths.

    isopleths: pandas DataFrame
        OPTIONAL. DataFrame object containing isopleth information as
        calculated by calculate_isobars_and_isopleths.

    Returns
    -------
    pandas DataFrame
        DataFrame with x and y values for all isobars and all isopleths.
        Useful if a user wishes to do custom plotting with isobar and isopleth
        data rather than using the built-in `plot_isobars_and_isopleths()`
        function.
    """
    np.seterr(divide='ignore', invalid='ignore')  # turn off numpy warning
    w.filterwarnings("ignore", message="Polyfit may be poorly conditioned")

    if isobars is not None:
        P_vals = isobars.Pressure.unique()
        isobars_lists = isobars.values.tolist()
        # add zero values to volatiles list
        isobars_lists.append([0.0, 0.0, 0.0, 0.0])

        isobars_pressure = []
        isobars_H2O_liq = []
        isobars_CO2_liq = []
        # do some data smoothing
        for pressure in P_vals:
            Pxs = [item[1] for item in isobars_lists if item[0] == pressure]
            Pys = [item[2] for item in isobars_lists if item[0] == pressure]

            try:
                # calcualte polynomial
                Pz = np.polyfit(Pxs, Pys, 3)
                Pf = np.poly1d(Pz)

                # calculate new x's and y's
                Px_new = np.linspace(Pxs[0], Pxs[-1], 50)
                Py_new = Pf(Px_new)

                # Save x's and y's
                Px_new_list = list(Px_new)
                isobars_H2O_liq += Px_new_list

                Py_new_list = list(Py_new)
                isobars_CO2_liq += Py_new_list

                pressure_vals_for_list = [pressure]*len(Px_new)
                isobars_pressure += pressure_vals_for_list

            except Exception:
                Px_list = list(Pxs)
                isobars_H2O_liq += Px_list

                Py_list = list(Pys)
                isobars_CO2_liq += Py_list

                pressure_vals_for_list = [pressure]*len(Pxs)
                isobars_pressure += pressure_vals_for_list

        isobar_df = pd.DataFrame({"Pressure": isobars_pressure,
                                  "H2O_liq": isobars_H2O_liq,
                                  "CO2_liq": isobars_CO2_liq})

    if isopleths is not None:
        XH2O_vals = isopleths.XH2O_fl.unique()
        isopleths_lists = isopleths.values.tolist()

        isopleths_XH2O_fl = []
        isopleths_H2O_liq = []
        isopleths_CO2_liq = []
        for Xfl in XH2O_vals:
            Xxs = [item[1] for item in isopleths_lists if item[0] == Xfl]
            Xys = [item[2] for item in isopleths_lists if item[0] == Xfl]

            try:
                # calculate polynomial
                Xz = np.polyfit(Xxs, Xys, 2)
                Xf = np.poly1d(Xz)

                # calculate new x's and y's
                Xx_new = np.linspace(Xxs[0], Xxs[-1], 50)
                Xy_new = Xf(Xx_new)

                # Save x's and y's
                Xx_new_list = list(Xx_new)
                isopleths_H2O_liq += Xx_new_list

                Xy_new_list = list(Xy_new)
                isopleths_CO2_liq += Xy_new_list

                XH2Ofl_vals_for_list = [Xfl]*len(Xx_new)
                isopleths_XH2O_fl += XH2Ofl_vals_for_list

            except Exception:
                Xx_list = list(Xxs)
                isopleths_H2O_liq += Xx_list

                Xy_list = list(Xys)
                isopleths_CO2_liq += Xy_list

                XH2Ofl_vals_for_list = [Xfl]*len(Xxs)
                isopleths_XH2O_fl += XH2Ofl_vals_for_list

        isopleth_df = pd.DataFrame({"XH2O_fl": isopleths_XH2O_fl,
                                    "H2O_liq": isopleths_H2O_liq,
                                    "CO2_liq": isopleths_CO2_liq})

    np.seterr(divide='warn', invalid='warn')  # turn numpy warning back on
    w.filterwarnings("always", message="Polyfit may be poorly conditioned")

    if isobars is not None:
        if isopleths is not None:
            return isobar_df, isopleth_df
        else:
            return isobar_df
    else:
        if isopleths is not None:
            return isopleth_df


def plot(isobars=None, isopleths=None, degassing_paths=None, custom_H2O=None,
         custom_CO2=None, isobar_labels=None, isopleth_labels=None,
         degassing_path_labels=None, custom_labels=None,
         custom_colors="VESIcal", custom_symbols=None, markersize=10,
         figsize=(12, 8), save_fig=False, extend_isobars_to_zero=True,
         smooth_isobars=False, smooth_isopleths=False, **kwargs):
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
    # Turn off warnings:
    np.seterr(divide='ignore', invalid='ignore')  # turn off numpy warning
    w.filterwarnings("ignore", message="Polyfit may be poorly conditioned")

    def check_inputs(custom_H2O, custom_CO2):
        if custom_H2O is not None:
            if custom_CO2 is None:
                raise core.InputError("If x data is passed, y data must also "
                                      "be passed.")
            else:
                if len(custom_H2O) == len(custom_CO2):
                    pass
                else:
                    raise core.InputError("x and y data must be same length")
        if custom_CO2 is not None:
            if custom_H2O is None:
                raise core.InputError("If y data is passed, x data must also "
                                      "be passed.")

    def check_colors(custom_colors):
        if custom_colors == "VESIcal":
            use_colors = color_list
        elif isinstance(custom_colors, list):
            use_colors = custom_colors
        else:
            raise core.InputError("Argument custom_colors must be type list. "
                                  "Just passing one item? Try putting square "
                                  "brackets, [], around it.")
        return use_colors

    def calc_extend_isobars_to_zero(Pxs, Pys):
        """
        Calculates new end-points for plotting isobars when
        extend_isobars_to_zero option is set to True.

        Parameters
        ----------
        Pxs, Pys: list
            List of x and y values corresponding to isobars.
        """
        if Pxs[0]*Pys[0] != 0.0:
            if Pxs[0] > Pys[0]:
                # create new array of length n+1 if n is the length of the
                # original array:
                Px_new = np.zeros(np.shape(Pxs)[0]+1)

                # set the first x value in the new array equal to 0:
                Px_new[0] = 0

                # fill the rest of the new array with the original array
                # values:
                Px_new[1:] = Pxs

                # overwrite the original array with the new one:
                Pxs = Px_new

                Py_new = np.zeros(np.shape(Pys)[0]+1)
                Py_new[0] = Pys[0]
                Py_new[1:] = Pys
                Pys = Py_new
            else:
                Px_new = np.zeros(np.shape(Pxs)[0]+1)
                Px_new[0] = Pxs[0]
                Px_new[1:] = Pxs
                Pxs = Px_new

                Py_new = np.zeros(np.shape(Pys)[0]+1)
                Py_new[0] = 0
                Py_new[1:] = Pys
                Pys = Py_new

        if Pxs[-1]*Pys[-1] != 0.0:
            if Pxs[-1] < Pys[-1]:
                Px_new = np.zeros(np.shape(Pxs)[0]+1)
                Px_new[-1] = 0
                Px_new[:-1] = Pxs
                Pxs = Px_new

                Py_new = np.zeros(np.shape(Pys)[0]+1)
                Py_new[-1] = Pys[-1]
                Py_new[:-1] = Pys
                Pys = Py_new
            else:
                Px_new = np.zeros(np.shape(Pxs)[0]+1)
                Px_new[-1] = Pxs[-1]
                Px_new[:-1] = Pxs
                Pxs = Px_new

                Py_new = np.zeros(np.shape(Pys)[0]+1)
                Py_new[-1] = 0
                Py_new[:-1] = Pys
                Pys = Py_new

        return Pxs, Pys

    # -------- HANDLE USER INPUT ERRORS, SET COLORS, SMOOTH LINES -------- ##
    check_inputs(custom_H2O=custom_H2O, custom_CO2=custom_CO2)
    use_colors = check_colors(custom_colors=custom_colors)

    if smooth_isobars:
        isobars = smooth_isobars_and_isopleths(isobars=isobars)
    if smooth_isopleths:
        isopleths = smooth_isobars_and_isopleths(isopleths=isopleths)

    # -------- CREATE FIGURE -------- ##
    fig, ax = plt.subplots(figsize=figsize)
    if 'custom_x' in kwargs:
        ax.set(xlabel=kwargs['xlabel'], ylabel=kwargs['ylabel'])
    else:
        ax.set(xlabel='H$_2$O wt%', ylabel='CO$_2$ wt%')

    labels = []

    # -------- PLOT ISOBARS -------- ##
    if isobars is not None:
        if isinstance(isobars, pd.DataFrame):
            isobars = [isobars]

        for i in range(len(isobars)):
            P_vals = isobars[i].Pressure.unique()
            isobars_lists = isobars[i].values.tolist()

            # add zero values to volatiles list
            isobars_lists.append([0.0, 0.0, 0.0, 0.0])

            P_iter = 0
            for pressure in P_vals:
                P_iter += 1
                Pxs = [item[1] for item in isobars_lists
                       if item[0] == pressure]
                Pys = [item[2] for item in isobars_lists
                       if item[0] == pressure]

                if extend_isobars_to_zero:
                    try:
                        Pxs, Pys = calc_extend_isobars_to_zero(Pxs, Pys)
                    except Exception:
                        pass
                else:
                    print(extend_isobars_to_zero)

                if len(isobars) > 1:
                    if P_iter == 1:
                        P_list = [int(i) for i in P_vals]
                        if isinstance(isobar_labels, list):
                            labels.append(str(isobar_labels[i]) + ' (' +
                                          ', '.join(map(str, P_list)) +
                                          " bars)")
                        else:
                            labels.append('Isobars ' + str(i+1) + ' (' +
                                          ', '.join(map(str, P_list)) +
                                          " bars)")
                    else:
                        labels.append('_nolegend_')

                if len(isobars) > 1:
                    ax.plot(Pxs, Pys, color=color_list[i])
                else:
                    ax.plot(Pxs, Pys)

            if len(isobars) == 1:
                labels += [str(P_val) + " bars" for P_val in P_vals]

    # -------- PLOT ISOPLETHS -------- ##
    if isopleths is not None:
        if isinstance(isopleths, pd.DataFrame):
            isopleths = [isopleths]

        for i in range(len(isopleths)):
            XH2O_vals = isopleths[i].XH2O_fl.unique()
            isopleths_lists = isopleths[i].values.tolist()

            H_iter = 0
            for Xfl in XH2O_vals:
                H_iter += 1
                Xxs = [item[1] for item in isopleths_lists if item[0] == Xfl]
                Xys = [item[2] for item in isopleths_lists if item[0] == Xfl]

                if len(isopleths) > 1:
                    if H_iter == 1:
                        H_list = [i for i in XH2O_vals]
                        if isinstance(isopleth_labels, list):
                            labels.append(str(isopleth_labels[i]) + ' (' +
                                          ', '.join(map(str, H_list)) +
                                          " XH2Ofluid)")
                        else:
                            labels.append('Isopleths ' + str(i+1) + ' (' +
                                          ', '.join(map(str, H_list)) +
                                          " XH2Ofluid)")
                    else:
                        labels.append('_nolegend_')
                    ax.plot(Xxs, Xys, ls='dashed', color=color_list[i])

                if len(isopleths) == 1:
                    H_list = [i for i in XH2O_vals]
                    if H_iter == 1:
                        labels.append('Isopleths (' +
                                      ', '.join(map(str, H_list)) +
                                      " XH2Ofluid)")
                    else:
                        labels.append('_nolegend_')
                    ax.plot(Xxs, Xys, ls='dashed', color='k')

    # -------- PLOT DEGASSING PATHS -------- ##
    if degassing_paths is not None:
        if isinstance(degassing_paths, pd.DataFrame):
            degassing_paths = [degassing_paths]

        degassing_colors = color_list.copy()
        iterno = 0
        plot_type = None
        for i in range(len(degassing_paths)):
            # handle whether labels are passed
            if degassing_path_labels is None:
                iterno += 1
                labels.append('Path%s' % iterno)
            else:
                labels.append(degassing_path_labels[iterno])
                iterno += 1

            # handle if only one volatile present
            if (degassing_paths[i]["H2O_liq"].max() > 0 and
                    degassing_paths[i]["CO2_liq"].max() > 0):
                ax.plot(degassing_paths[i]["H2O_liq"],
                        degassing_paths[i]["CO2_liq"],
                        ls='dotted', color=degassing_colors[i])
                ax.plot(degassing_paths[i]["H2O_liq"].max(),
                        degassing_paths[i]["CO2_liq"].max(),
                        'o', color=degassing_colors[i])
                # warn user if trying to plot mixed types on same figure
                if plot_type not in [None, "mixed"]:
                    w.warn("Warning: you are trying to plot two different plot types on the same"
                           " figure. Curves may not match axis labels.", stacklevel=2)
                plot_type = "mixed"
            elif (degassing_paths[i]["H2O_liq"].max() > 0 and
                  degassing_paths[i]["CO2_liq"].max() <= 0):
                ax.plot(degassing_paths[i]["Pressure_bars"],
                        degassing_paths[i]["H2O_liq"],
                        ls='dotted', color=degassing_colors[i])
                ax.plot(degassing_paths[i]["Pressure_bars"].max(),
                        degassing_paths[i]["H2O_liq"].max(),
                        'o', color=degassing_colors[i])
                ax.set(xlabel='Pressure (bars)', ylabel='H$_2$O wt%')
                # warn user if trying to plot mixed types on same figure
                if plot_type not in [None, "H2O"]:
                    w.warn("Warning: you are trying to plot two different plot types on the same"
                           " figure. Curves may not match axis labels.", stacklevel=2)
                plot_type = "H2O"
            elif (degassing_paths[i]["H2O_liq"].max() <= 0 and
                  degassing_paths[i]["CO2_liq"].max() > 0):
                ax.plot(degassing_paths[i]["Pressure_bars"],
                        degassing_paths[i]["CO2_liq"],
                        ls='dotted', color=degassing_colors[i])
                ax.plot(degassing_paths[i]["Pressure_bars"].max(),
                        degassing_paths[i]["CO2_liq"].max(),
                        'o', color=degassing_colors[i])
                ax.set(xlabel='Pressure (bars)', ylabel='CO$_2$ wt%')
                # warn user if trying to plot mixed types on same figure
                if plot_type not in [None, "CO2"]:
                    w.warn("Warning: you are trying to plot two different plot types on the same"
                           " figure. Curves may not match axis labels.", stacklevel=2)
                plot_type = "CO2"

            labels.append('_nolegend_')

    # -------- PLOT CUSTOM H2O-CO2 -------- ##
    if custom_H2O is not None and custom_CO2 is not None:
        if isinstance(custom_H2O, pd.DataFrame):
            custom_H2O = [custom_H2O]
        if isinstance(custom_CO2, pd.DataFrame):
            custom_CO2 = [custom_CO2]

        if custom_symbols is None:
            use_marker = ['o'] * len(custom_H2O)
        else:
            use_marker = custom_symbols

        iterno = 0
        for i in range(len(custom_H2O)):
            if custom_labels is None:
                iterno += 1
                labels.append('Custom%s' % iterno)
                ax.plot(custom_H2O[i], custom_CO2[i], use_marker[i],
                        color=use_colors[i], markersize=markersize)
            else:
                labels.append(custom_labels[iterno])
                ax.plot(custom_H2O[i], custom_CO2[i], use_marker[i],
                        color=use_colors[i], markersize=markersize)
                iterno += 1

    # -------- PLOT CUSTOM X-Y -------- ##
    if 'custom_x' in kwargs:
        custom_x = kwargs['custom_x']
        custom_y = kwargs['custom_y']

        if isinstance(custom_x, pd.core.series.Series):
            custom_x = [list(custom_x.values)]
        if isinstance(custom_y, pd.core.series.Series):
            custom_y = [list(custom_y.values)]

        if custom_symbols is None:
            use_marker = ['o'] * len(custom_x)
        else:
            use_marker = custom_symbols

        iterno = 0
        for i in range(len(custom_x)):
            if custom_labels is None:
                iterno += 1
                labels.append('Custom%s' % iterno)
                ax.plot(custom_x[i], custom_y[i], use_marker[i],
                        color=use_colors[i], markersize=markersize)
            else:
                labels.append(custom_labels[iterno])
                ax.plot(custom_x[i], custom_y[i], use_marker[i],
                        color=use_colors[i], markersize=markersize)
                iterno += 1

    # -------- PLOT LEGEND -------- ##
    ax.legend(labels, bbox_to_anchor=(1.01, 1), loc='upper left')

    if 'custom_x' not in kwargs:
        ax.set_xlim(left=0)
        ax.set_ylim(bottom=0)

    np.seterr(divide='warn', invalid='warn')  # turn numpy warning back on
    w.filterwarnings("always", message="Polyfit may be poorly conditioned")

    # -------- SAVE FIGURE IF DESIRED -------- ##
    if save_fig is not False:
        fig.savefig(save_fig)

    return fig, ax


def scatterplot(custom_x, custom_y, xlabel=None, ylabel=None, **kwargs):
    """
    Custom x-y plotting using VESIcal's built-in plot() function, built
    Matplotlib's plot and scatter functions.

    Parameters
    ----------
    custom_x: list
        List of groups of x-values to plot as points or lines

    custom_y: list
        List of groups of y-values to plot as points or lines

    xlabel: str
        OPTIONAL. What to display along the x-axis.

    ylabel: str
        OPTIONAL. What to display along the y-axis.

    kwargs:
        Can take in any key word agruments that can be passed to `plot()`.

    Returns
    -------
    fig, ax matplotlib objects
        X-y plot with custom x and y axis values and labels.
    """

    if isinstance(custom_x, list) and isinstance(custom_y, list):
        if len(custom_x) != len(custom_y):
            raise core.InputError("X and y lists must be same length")

    if xlabel is not None:
        if isinstance(xlabel, str):
            pass
        else:
            raise core.InputError("xlabel must be string")

    if ylabel is not None:
        if isinstance(ylabel, str):
            pass
        else:
            raise core.InputError("ylabel must be string")

    return plot(custom_x=custom_x, custom_y=custom_y, xlabel=xlabel,
                ylabel=ylabel, **kwargs)


# ------- Define custom plotting tools for checking calibrations ------- #

def calib_plot(user_data=None, model='all', plot_type='TAS', zoom=None,
               figsize=(17, 8), legend=True, save_fig=False, **kwargs):
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

    # Get x and y axis limits, if user passed them
    if zoom is None:
        user_xmin = 35
        user_xmax = 100
        user_ymin = 0
        user_ymax = 25
    elif zoom == 'user_data':
        if isinstance(user_data, pd.DataFrame):
            print("'user_data' type zoom for more than one sample is not ",
                  "implemented yet.")
            user_xmin = 35
            user_xmax = 100
            user_ymin = 0
            user_ymax = 25
        elif (isinstance(user_data, pd.core.series.Series) or
              isinstance(user_data, dict)):
            user_xmin = user_data['SiO2'] - 5
            user_xmax = user_data['SiO2'] + 5
            user_ymin = user_data['Na2O'] + user_data['K2O'] - 2
            if user_ymin < 0:
                user_ymin = 0
            user_ymax = user_data['Na2O'] + user_data['K2O'] + 2
    elif isinstance(zoom, list):
        user_xmin, user_xmax = zoom[0]
        user_ymin, user_ymax = zoom[1]
    else:
        raise core.InputError('Trying to pass zoom coords? Pass as ' +
                              '[(x, x), (y, y)]')

    # Create the figure
    fig, ax1 = plt.subplots(figsize=figsize)
    font = {'family': 'sans-serif',
            'color':  'black',
            'weight': 'normal',
            'size':   20,
            }

    # TAS figure
    if plot_type == 'TAS':
        # adjust x limits here if you want to focus on a specific part of
        # compostional space:
        ax1.set_xlim([user_xmin, user_xmax])

        # adjust y limits here
        ax1.set_ylim([user_ymin, user_ymax])
        plt.xlabel('SiO$_2$, wt%', fontdict=font, labelpad=15)
        plt.ylabel('Na$_2$O+K$_2$O, wt%', fontdict=font, labelpad=15)

        # add LeMaitre fields
        if zoom is None:
            add_LeMaitre_fields(ax1)

    elif plot_type == 'xy':
        if 'x' in kwargs and 'y' in kwargs:
            x = kwargs['x']
            y = kwargs['y']
            if zoom is not None:
                ax1.set_xlim([user_xmin, user_xmax])
                ax1.set_ylim([user_ymin, user_ymax])
            plt.xlabel(str(x)+", wt%", fontdict=font, labelpad=15)
            plt.ylabel(str(y)+", wt%", fontdict=font, labelpad=15)
        else:
            raise core.InputError("If plot_type is 'xy', then x and y "
                                  "values must be passed as strings. For "
                                  "example, x='SiO2', y='Al2O3'.")

    # Plot Calibration Data
    if model == 'all':
        model = ['MagmaSat',
                 'Shishkina',
                 'Dixon',
                 'IaconoMarziano',
                 'Liu',
                 'AllisonCarbon',
                 'MooreWater']
    if model == 'mixed':
        model = ['MagmaSat',
                 'Shishkina',
                 'Dixon',
                 'IaconoMarziano',
                 'Liu']

    if isinstance(model, str):
        model = [model]

    if isinstance(model, list):
        # set legends to false
        h2o_legend = False
        co2_h2oco2_legend = False

        # check which legends to turn to True
        for modelname in model:
            model_type = calibrations.return_calibration_type(modelname)
            if model_type['H2O']:
                h2o_legend = True
            if model_type['CO2'] or model_type['Mixed']:
                co2_h2oco2_legend = True

        if h2o_legend:
            plt.scatter([], [], marker='', label=r"$\bf{Pure \ H_2O:}$")

            for modelname in model:
                calibdata = calibrations.return_calibration(modelname)
                model_type = calibrations.return_calibration_type(modelname)
                if isinstance(calibdata, str):
                    w.warn(calibdata)
                else:
                    if model_type['H2O']:
                        if plot_type == 'TAS':
                            try:
                                plt.scatter(calibdata['H2O']['SiO2'],
                                            (calibdata['H2O']['Na2O'] +
                                                calibdata['H2O']['K2O']),
                                            marker='s', edgecolors='k',
                                            facecolors=calibdata['facecolor'],
                                            label=str(modelname))
                            except Exception:
                                plt.scatter(calibdata['H2O']['SiO2'],
                                            calibdata['H2O']['Na2O+K2O'],
                                            marker='s', edgecolors='k',
                                            facecolors=calibdata['facecolor'],
                                            label=str(modelname))
                        if plot_type == 'xy':
                            try:
                                plt.scatter(calibdata['H2O'][x],
                                            calibdata['H2O'][y],
                                            marker='s', edgecolors='k',
                                            facecolors=calibdata['facecolor'],
                                            label=str(modelname))
                            except Exception:
                                w.warn("The requested oxides were not found",
                                       "in the calibration dataset for " +
                                       str(modelname) + ".")

            if co2_h2oco2_legend:
                plt.scatter([], [], marker='', label=r"${\ }$")

        if co2_h2oco2_legend:
            plt.scatter([], [], marker='',
                        label=r"$\bf{\ CO_2 \ and \ H_2O\!-\!CO_2:}$")

        for modelname in model:
            calibdata = calibrations.return_calibration(modelname)
            model_type = calibrations.return_calibration_type(modelname)
            if isinstance(calibdata, str):
                w.warn(calibdata)
            else:
                if model_type['CO2'] and model_type['Mixed']:
                    frames = [calibdata['CO2'], calibdata['Mixed']]
                    co2_and_mixed = pd.concat(frames)
                    if plot_type == 'TAS':
                        try:
                            plt.scatter(co2_and_mixed['SiO2'],
                                        (co2_and_mixed['Na2O'] +
                                        co2_and_mixed['K2O']),
                                        marker='d', edgecolors='k',
                                        facecolors=calibdata['facecolor'],
                                        label=str(modelname))
                        except Exception:
                            plt.scatter(co2_and_mixed['SiO2'],
                                        co2_and_mixed['Na2O+K2O'],
                                        marker='d', edgecolors='k',
                                        facecolors=calibdata['facecolor'],
                                        label=str(modelname))
                    if plot_type == 'xy':
                        try:
                            plt.scatter(co2_and_mixed[x], co2_and_mixed[y],
                                        marker='d', edgecolors='k',
                                        facecolors=calibdata['facecolor'],
                                        label=str(modelname))
                        except Exception:
                            w.warn("The requested oxides were not found in ",
                                   "the calibration dataset for " +
                                   str(modelname) + ".")
                elif model_type['CO2'] or model_type['Mixed']:
                    if model_type['CO2']:
                        thistype = 'CO2'
                    if model_type['Mixed']:
                        thistype = 'Mixed'
                    if plot_type == 'TAS':
                        try:
                            plt.scatter(calibdata[thistype]['SiO2'],
                                        (calibdata[thistype]['Na2O'] +
                                        calibdata[thistype]['K2O']),
                                        marker='d', edgecolors='k',
                                        facecolors=calibdata['facecolor'],
                                        label=str(modelname))
                        except Exception:
                            plt.scatter(calibdata[thistype]['SiO2'],
                                        calibdata[thistype]['Na2O+K2O'],
                                        marker='d', edgecolors='k',
                                        facecolors=calibdata['facecolor'],
                                        label=str(modelname))
                    if plot_type == 'xy':
                        try:
                            plt.scatter(calibdata[thistype][x],
                                        calibdata[thistype][y],
                                        marker='d', edgecolors='k',
                                        facecolors=calibdata['facecolor'],
                                        label=str(modelname))
                        except Exception:
                            w.warn("The requested oxides were not found in ",
                                   "the calibration dataset for "
                                   + str(modelname) + ".")
    else:
        raise core.InputError("model must be of type str or list")

    # Plot user data
    if user_data is None:
        pass
    else:
        if ((user_data.__class__.__module__, user_data.__class__.__name__) ==
           ('VESIcal', 'BatchFile')):
            user_data = user_data.get_data()
            # batchfile and VESIcal (__init__) are not imported to avoid
            # circular imports
            # use above notation to interrogate datatype
        if ((user_data.__class__.__module__, user_data.__class__.__name__) ==
           ('VESIcal', 'Sample')):
            user_data = user_data.get_composition()
            # batchfile and VESIcal (__init__) are not imported to avoid
            # circular imports
            # use above notation to interrogate datatype
        if plot_type == 'TAS':
            _sample = user_data.copy()
            try:
                _sample["TotalAlkalis"] = _sample["Na2O"] + _sample["K2O"]
            except Exception:
                core.InputError("Na2O and K2O data must be in user_data")
            plt.scatter(_sample['SiO2'], _sample['TotalAlkalis'],
                        s=150, edgecolors='w', facecolors='red', marker='P',
                        label='User Data')
        if plot_type == 'xy':
            _sample = user_data.copy()
            plt.scatter(_sample[x], _sample[y],
                        s=150, edgecolors='w', facecolors='red', marker='P',
                        label='User Data')

    if legend:
        plt.legend(bbox_to_anchor=(1.04, 1), loc="upper left")
    fig.tight_layout()
    if isinstance(save_fig, str):
        fig.savefig(save_fig)

    return fig, ax1


def show():
    """
    Local implementation of pyplot.show(). For displaying created plots.
    """
    plt.show()
