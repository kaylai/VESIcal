from abc import abstractmethod
import numpy as np
import pandas as pd
import warnings as w
from scipy.optimize import root_scalar
from scipy.optimize import root
from copy import deepcopy

from VESIcal import activity_models
from VESIcal import core
from VESIcal import fugacity_models


class Model(object):
    """The model object implements a volatile solubility model. It is composed
    of the methods needed to evaluate
    :func:`VESIcal.calculate_dissolved_volatiles`,
    :func:`VESIcal.calculate_equilibrium_fluid_comp`, and
    :func:`calculate_saturation_pressure`. The fugacity and activity models for
    the volatiles species must be specified, defaulting to ideal.
    """

    def __init__(self):
        self.set_volatile_species(None)
        self.set_fugacity_model(fugacity_models.fugacity_idealgas())
        self.set_activity_model(activity_models.activity_idealsolution())
        self.set_calibration_ranges([])
        self.set_solubility_dependence(False)

    def set_volatile_species(self, volatile_species):
        if type(volatile_species) == str:
            volatile_species = [volatile_species]
        elif type(volatile_species) != list:
            raise core.InputError("volatile_species must be a str or list.")
        self.volatile_species = volatile_species

    def set_fugacity_model(self, fugacity_model):
        self.fugacity_model = fugacity_model

    def set_activity_model(self, activity_model):
        self.activity_model = activity_model

    def set_calibration_ranges(self, calibration_ranges):
        self.calibration_ranges = calibration_ranges

    def set_solubility_dependence(self, solubility_dependence):
        self.solubility_dependence = solubility_dependence

    def get_calibration_values(self, variable_names):
        """ Returns the values stored as the calibration range for the given
        variable(s). However, for checks where there is a single value- i.e.
        cr_GreaterThan or crf_LessThan, the logical operation will remain a
        mystery until someone figures out an elegant way of communicating it.

        Parameters
        ----------
        variable_names     str or list
            The name(s) of the variables you want the calibration ranges for.

        Returns
        -------
        list
            A list of values or tuples for the calibration ranges in the order
            given.
        """

        # Check if the variable name is passed as a string, and if so put it in
        # a list:
        if type(variable_names) == str:
            variable_names = [variable_names]

        calibration_values = []

        for var in variable_names:
            found_var = False
            for cr in self.calibration_ranges:
                if found_var is False:
                    if cr.parameter_name == var:
                        found_var = True
                        calibration_values.append(cr.value)
            if found_var is False:
                calibration_values.append(np.nan)

        return calibration_values

    @abstractmethod
    def calculate_dissolved_volatiles(self, **kwargs):
        pass

    @abstractmethod
    def calculate_equilibrium_fluid_comp(self, **kwargs):
        pass

    @abstractmethod
    def calculate_saturation_pressure(self, **kwargs):
        pass

    # @abstractmethod
    # def preprocess_sample(self,**kwargs):
    #     pass

    def check_calibration_range(self, parameters, report_nonexistance=True):
        """ Checks whether the given parameters are within the ranges defined
        by the CalibrationRange objects for the model and its fugacity and
        activity models. An empty string will be returned if all parameters are
        within the calibration range. If a parameter is not within the
        calibration range, a description of the problem will be returned in the
        string.

        Parameters
        ----------
        parameters     dict
            Dictionary keys are the names of the parameters to be checked,
            e.g., pressure temperature, SiO2, etc. Values are the values of
            each parameter. A complete set need not be given.

        Returns
        -------
        str
            String description of any parameters falling outside of the
            calibration range.
        """
        s = ''
        for cr in self.calibration_ranges:
            if cr.check(parameters) is False:
                s += cr.string(parameters, report_nonexistance)
        for cr in self.fugacity_model.calibration_ranges:
            if cr.check(parameters) is False:
                s += cr.string(parameters, report_nonexistance)
        for cr in self.activity_model.calibration_ranges:
            if cr.check(parameters) is False:
                s += cr.string(parameters, report_nonexistance)
        return s

    def get_calibration_range(self):
        """ Returns a string describing the calibration ranges defined by the
        CalibrationRange objects for each model, and its associated fugacity
        and activity models.

        Returns
        -------
        str
            String description of the calibration range objects."""
        s = ''
        for cr in self.calibration_ranges:
            s += cr.string(None)
        for cr in self.fugacity_model.calibration_ranges:
            s += cr.string(None)
        for cr in self.activity_model.calibration_ranges:
            s += cr.string(None)
        return s


# ------------ MIXED FLUID MODELS ------------------------------- #
class MixedFluid(Model):
    """
    Implements the generic framework for mixed fluid solubility. Any set of
    pure fluid solubility models may be specified.
    """
    def __init__(self, models):
        """
        Initializes the mixed fluid model.

        Parameters
        ----------
        models     dictionary
            Dictionary with names of volatile species as keys, and the model
            objects as values.
        """
        self.models = tuple(model for model in models.values())
        self.set_volatile_species(list(models.keys()))

    def get_calibration_values(self, variable_names):
        """ Placeholder method to prevent an error when this generic method is
        called for a MixedFluid model.

        Returns
        -------
        np.nan
        """
        return np.nan

    def calculate_dissolved_volatiles(self, pressure, X_fluid,
                                      returndict=False, **kwargs):
        """
        Calculates the dissolved volatile concentrations in wt%, using each
        model's calculate_dissolved_volatiles method. At present the volatile
        concentrations are not propagated through.

        Parameters
        ----------
        pressure     float
            The total pressure in bars.
        X_fluid     float, numpy.ndarry, dict, pandas Series
            The mole fraction of each species in the fluid. If the mixed fluid
            model contains only two species (e.g. CO2 and H2O), the value of
            the first species in self.volatile_species may be passed on its own
            as a float.
        returndict         bool
            If True, the results will be returned in a dict, otherwise they
            will be returned as a tuple.

        Returns
        -------
        tuple
            Dissolved volatile concentrations of each species in the model, in
            the order set by self.volatile_species.
        """
        if ((type(X_fluid) == float or type(X_fluid) == int) and
           len(self.volatile_species) == 2):
            X_fluid = (X_fluid, 1-X_fluid)
        elif len(X_fluid) != len(self.volatile_species):
            raise core.InputError("X_fluid must have the same length as the "
                                  "number of volatile species in the "
                                  "MixedFluids Model class, or it may have "
                                  "length 1 if two species are present in "
                                  "the MixedFluids Model class.")

        if np.sum(X_fluid) != 1.0:
            raise core.InputError("X_fluid must sum to 1.0")
        if any(val < 0 for val in X_fluid) or any(val > 1 for val in X_fluid):
            raise core.InputError("Each mole fraction in X_fluid must have a "
                                  "value between 0 and 1.")

        if type(X_fluid) == dict or type(X_fluid) == pd.core.series.Series:
            X_fluid = tuple(X_fluid[species] for species in
                            self.volatile_species)

        # If the models don't depend on the concentration of volatiles,
        # themselves:
        if all(model.solubility_dependence is False for model in self.models):
            result = tuple(model.calculate_dissolved_volatiles(
                        pressure=pressure, X_fluid=Xi, **kwargs) for model, Xi
                        in zip(self.models, X_fluid))
        # If one of the models depends on the other volatile concentration
        elif (len(self.models) == 2 and
              self.models[0].solubility_dependence is False and
              'sample' in kwargs):
            result0 = self.models[0].calculate_dissolved_volatiles(
                              pressure=pressure, X_fluid=X_fluid[0], **kwargs)
            samplecopy = kwargs['sample'].change_composition(
                             {self.volatile_species[0]: result0},
                             inplace=False)
            kwargs['sample'] = samplecopy
            result1 = self.models[1].calculate_dissolved_volatiles(
                            pressure=pressure, X_fluid=X_fluid[1], **kwargs)
            result = (result0, result1)
        elif(len(self.models) == 2 and
             self.models[1].solubility_dependence is False and
             'sample' in kwargs):
            result1 = self.models[1].calculate_dissolved_volatiles(
                               pressure=pressure, X_fluid=X_fluid[1], **kwargs)
            samplecopy = kwargs['sample'].change_composition(
                            {self.volatile_species[1]: result1}, inplace=False)
            kwargs['sample'] = samplecopy
            result0 = self.models[0].calculate_dissolved_volatiles(
                               pressure=pressure, X_fluid=X_fluid[0], **kwargs)
            result = (result0, result1)
        else:
            raise core.InputError("The solubility dependence of the models "
                                  "is not currently supported by the "
                                  "MixedFluid model.")

        if returndict:
            resultsdict = {}
            for i, v in zip(range(len(self.volatile_species)),
                            self.volatile_species):
                resultsdict.update({v+'_liq': result[i]})
            return resultsdict
        else:
            return result

    def calculate_equilibrium_fluid_comp(self, pressure, sample,
                                         return_dict=True, **kwargs):
        """ Calculates the composition of the fluid in equilibrium with the
        dissolved volatile concentrations passed. If a fluid phase is
        undersaturated at the chosen pressure (0,0) will be returned. Note,
        this currently assumes the given H2O and CO2 concentrations are the
        system total, not the total dissolved. If one of the volatile species
        has a zero or negative concentration, the pure fluid model for the
        other volatile species will be used.

        Parameters
        ----------
        pressure     float
            The total pressure in bars.
        sample    Sample class
            Magma major element composition.
        return_dict     bool
            Set the return type, if true a dict will be returned, if False two
            floats will be returned. Default is True.

        Returns
        -------
        dict or floats
            Mole fractions of the volatile species in the fluid, in the order
            given by self.volatile_species if floats.
        """
        if len(self.volatile_species) != 2:
            raise core.InputError("Currently equilibrium fluid compositions "
                                  "can only be calculated when two volatile "
                                  "species are present.")

        dissolved_at_0bar = [self.models[0].calculate_dissolved_volatiles(
                                    sample=sample, pressure=0.0, **kwargs),
                             self.models[1].calculate_dissolved_volatiles(
                                    sample=sample, pressure=0.0, **kwargs)]

        if (sample.get_composition(self.volatile_species[0]) <= 0.0 or
           (sample.get_composition(self.volatile_species[0]) <=
           dissolved_at_0bar[0])):
            Xv0 = 0.0
            Xv1 = self.models[1].calculate_equilibrium_fluid_comp(
                                    pressure=pressure, sample=sample, **kwargs)
        elif (sample.get_composition(self.volatile_species[1]) <= 0.0 or
              (sample.get_composition(self.volatile_species[1]) <=
              dissolved_at_0bar[1])):
            Xv1 = 0.0
            Xv0 = self.models[0].calculate_equilibrium_fluid_comp(
                                    pressure=pressure, sample=sample, **kwargs)
        else:
            satP = self.calculate_saturation_pressure(sample, **kwargs)

            if satP < pressure:
                if return_dict:
                    return {self.volatile_species[0]: 0,
                            self.volatile_species[1]: 0}
                else:
                    return (0, 0)

            molfracs = sample.get_composition(units='mol_oxides')
            (Xt0, Xt1) = (molfracs[self.volatile_species[0]],
                          molfracs[self.volatile_species[1]])

            try:
                Xv0 = root_scalar(self.root_for_fluid_comp,
                                  bracket=[1e-15, 1-1e-15],
                                  args=(pressure, Xt0, Xt1, sample,
                                        kwargs)).root
                Xv1 = 1 - Xv0
            except Exception:
                try:
                    Xv0 = root_scalar(self.root_for_fluid_comp, x0=0.5, x1=0.1,
                                      args=(pressure, Xt0, Xt1, sample,
                                            kwargs)).root
                    Xv1 = 1 - Xv0
                except Exception:
                    raise core.SaturationError("Equilibrium fluid not found. "
                                               "Likely an issue with the "
                                               "numerical solver.")

        if return_dict:
            return {self.volatile_species[0]: Xv0,
                    self.volatile_species[1]: Xv1}
        else:
            return Xv0, Xv1

    def calculate_saturation_pressure(self, sample, **kwargs):
        """
        Calculates the pressure at which a fluid will be saturated, given the
        dissolved volatile concentrations. If one of the volatile species has a
        zero or negative concentration the pure fluid model for the other
        species will be used. If one of the volatile species has a
        concentration lower than the concentration dissolved at 0 bar, the pure
        fluid model for the other species will be used.

        Parameters
        ----------
        sample     Sample class
            Magma major element composition (including volatiles).

        Returns
        -------
        float
            The saturation pressure in bars.
        """
        dissolved_at_0bar = [self.models[0].calculate_dissolved_volatiles(
                                        sample=sample, pressure=0.0, **kwargs),
                             self.models[1].calculate_dissolved_volatiles(
                                        sample=sample, pressure=0.0, **kwargs)]

        if (sample.get_composition(self.volatile_species[0]) <= 0.0 or
            (sample.get_composition(self.volatile_species[0]) <=
             dissolved_at_0bar[0])):
            satP = self.models[1].calculate_saturation_pressure(sample=sample,
                                                                **kwargs)
        elif (sample.get_composition(self.volatile_species[1]) <= 0.0 or
              (sample.get_composition(self.volatile_species[1]) <=
               dissolved_at_0bar[1])):
            satP = self.models[0].calculate_saturation_pressure(sample=sample,
                                                                **kwargs)
        else:
            volatile_concs = np.array(tuple(sample.get_composition(species) for
                                            species in self.volatile_species))

            x0 = 0
            for model in self.models:
                xx0 = model.calculate_saturation_pressure(sample=sample,
                                                          **kwargs)
                if np.isfinite(xx0):
                    x0 += xx0

            try:
                satP = root(self.root_saturation_pressure, x0=[x0, 0.5],
                            args=(volatile_concs, sample, kwargs)).x[0]
            except Exception:
                w.warn("Saturation pressure not found.", RuntimeWarning,
                       stacklevel=2)
                satP = np.nan

        return satP

    def calculate_isobars_and_isopleths(self, pressure_list,
                                        isopleth_list=[0, 1], points=51,
                                        return_dfs=True, extend_to_zero=True,
                                        **kwargs):
        """
        Calculates isobars and isopleths. Isobars can be calculated for any
        number of pressures. Variables required by each of the pure fluid
        models must be passed, e.g. sample, temperature, etc.

        Parameters
        ----------
        pressure_list     list or float
            List of all pressure values at which to calculate isobars, in bars.
        isopleth_list     list
            Default value is None, in which case only isobars will be
            calculated. List of all fluid compositions in mole fraction (of the
            first species in self.volatile_species) at which to calcualte
            isopleths. Values can range from 0 to 1.
        points     int
            The number of points in each isobar and isopleth. Default value is
            101.
        return_dfs     bool
            If True, the results will be returned as two pandas DataFrames, as
            produced by the MagmaSat method. If False the results will be
            returned as lists of numpy arrays.

        Returns
        -------
        pandas DataFrame object(s) or list(s)
            If isopleth_list is not None, two objects will be returned, one
            with the isobars and the second withthe isopleths. If return_dfs is
            True, two pandas DataFrames will be returned with column names
            'Pressure' or 'XH2O_fl', 'H2O_liq', and 'CO2_liq'. If return_dfs is
            False, two lists of numpy arrays will be returned. Each array is an
            individual isobar or isopleth, in the order passed via
            pressure_list or isopleth_list. The arrays are the concentrations
            of H2O and CO2 in the liquid, in the order of the species in
            self.volatile_species.

        """
        if (len(self.volatile_species) != 2 or
                'H2O' not in self.volatile_species or
                'CO2' not in self.volatile_species):
            raise core.InputError("calculate_isobars_and_isopleths may only "
                                  "be used with a H2O-CO2 fluid.")

        H2O_id = self.volatile_species.index('H2O')
        CO2_id = self.volatile_species.index('CO2')

        if isinstance(pressure_list, list):
            pass
        elif (isinstance(pressure_list, int) or
              isinstance(pressure_list, float)):
            pressure_list = [pressure_list]
        else:
            raise core.InputError("pressure_list must be a single float "
                                  "(1000.0), int (1000), or list of those "
                                  "[1000, 2000.0, 3000].")

        has_isopleths = True
        if isopleth_list is None:
            has_isopleths = False

        isobars_df = pd.DataFrame(columns=['Pressure', 'H2O_liq', 'CO2_liq'])
        isobars = []
        for pressure in pressure_list:
            dissolved = np.zeros([2, points])
            Xv0 = np.linspace(0.0, 1.0, points)
            for i in range(points):
                dissolved[:, i] = self.calculate_dissolved_volatiles(
                       pressure=pressure, X_fluid=(Xv0[i], 1-Xv0[i]), **kwargs)
                isobars_df = isobars_df.append({
                                               'Pressure': pressure,
                                               'H2O_liq': dissolved[H2O_id, i],
                                               'CO2_liq': dissolved[CO2_id, i]
                                               }, ignore_index=True)
            isobars.append(dissolved)

        if has_isopleths:
            isopleths_df = pd.DataFrame(columns=['XH2O_fl', 'H2O_liq',
                                                 'CO2_liq'])
            isopleths = []
            for isopleth in isopleth_list:
                dissolved = np.zeros([2, points])
                pmin = np.nanmin(pressure_list)
                pmax = np.nanmax(pressure_list)
                if pmin == pmax:
                    pmin = 0.0
                pressure = np.linspace(pmin, pmax, points)
                for i in range(points):
                    dissolved[:, i] = self.calculate_dissolved_volatiles(
                          pressure=pressure[i], X_fluid=(isopleth, 1-isopleth),
                          **kwargs)
                    isopleths_df = isopleths_df.append({
                       'XH2O_fl': [isopleth, 1-isopleth][H2O_id],
                       'H2O_liq': dissolved[H2O_id, i],
                       'CO2_liq': dissolved[CO2_id, i]}, ignore_index=True)
                isopleths.append(dissolved)

        if return_dfs:
            if has_isopleths:
                return (isobars_df, isopleths_df)
            else:
                return isobars_df
        else:
            if has_isopleths:
                return (isobars, isopleths)
            else:
                return isobars

    def calculate_degassing_path(self, sample, pressure='saturation',
                                 fractionate_vapor=0.0, final_pressure=100.0,
                                 steps=101, return_dfs=True,
                                 round_to_zero=True, **kwargs):
        """
        Calculates the dissolved volatiles in a progressively degassing sample.

        Parameters
        ----------
        sample     Sample class
            Magma major element composition (including volatiles).
        pressure     string, float, int, list, or numpy array
            Defaults to 'saturation', the calculation will begin at the
            saturation pressure. If a number is passed as either a float or
            int, this will be the starting pressure. If a list of numpy array
            is passed, the pressure values in the list or array will define the
            degassing path, i.e. final_pressure and steps variables will be
            ignored. Units are bars.
        fractionate_vapor     float
            What proportion of vapor should be removed at each step. If 0.0
            (default), the degassing path will correspond to closed-system
            degassing. If 1.0, the degassing path will correspond to
            open-system degassing.
        final_pressure         float
            The final pressure on the degassing path, in bars. Ignored if a
            list or numpy array is passed as the pressure variable. Default is
            1 bar.
        steps     int
            The number of steps in the degassing path. Ignored if a list or
            numpy array are passed as the pressure variable.
        return_dfs     bool
            If True, the results will be returned in a pandas DataFrame, if
            False, two numpy arrays will be returned.
        round_to_zero   bool
            If True, the first entry of FluidProportion_wt will be rounded to
            zero, rather than being a value within numerical error of zero.
            Default is True.

        Returns
        -------
        pandas DataFrame or numpy arrays
            If return_dfs is True (default), a DataFrame with columns
            'Pressure_bars', 'H2O_liq', 'CO2_liq', 'H2O_fl', 'CO2_fl', and
            'FluidProportion_wt', is returned. Dissolved volatiles are in wt%,
            the proportions of volatiles in the fluid are in mole fraction.
            Otherwise a numpy array containing the dissolved volatile
            concentrations, and a numpy array containing the mole fractions of
            volatiles in the fluid is returned. The columns are in the order of
            the volatiles in self.volatile_species.
        """

        # Create a copy of the sample so that initial volatile concentrations
        # are not overwritten.
        sample = deepcopy(sample)

        # Its imperative that normalization doesn't change the volatile
        # concentrations throughout the calculation.
        if sample.default_normalization not in ['fixedvolatiles', 'none']:
            sample.set_default_normalization('fixedvolatiles')
            w.warn('Sample normalization changed to fixedvolatiles.')

        wtptoxides = sample.get_composition(units='wtpt_oxides')
        wtm0s, wtm1s = (wtptoxides[self.volatile_species[0]],
                        wtptoxides[self.volatile_species[1]])

        if pressure == 'saturation':
            p0 = self.calculate_saturation_pressure(sample, **kwargs)
            pressures = np.linspace(p0, final_pressure, steps)
        elif type(pressure) == float or type(pressure) == int:
            pressures = np.linspace(pressure, final_pressure, steps)
        elif type(pressure) == list or type(pressure) == np.ndarray:
            pressures = pressure

        Xv = np.zeros([2, len(pressures)])
        wtm = np.zeros([2, len(pressures)])

        for i in range(len(pressures)):
            try:
                wtptoxides = sample.get_composition(units='wtpt_oxides')
                X_fluid = self.calculate_equilibrium_fluid_comp(
                    pressure=pressures[i], sample=sample, return_dict=False,
                    **kwargs)
                Xv[:, i] = X_fluid
                if X_fluid == (0, 0):
                    wtm[:, i] = (wtptoxides[self.volatile_species[0]],
                                 wtptoxides[self.volatile_species[1]])
                else:
                    if X_fluid[0] == 0:
                        wtm[0, i] = wtptoxides[self.volatile_species[0]]
                        wtm[1, i] = self.calculate_dissolved_volatiles(
                            pressure=pressures[i], sample=sample,
                            X_fluid=X_fluid, **kwargs)[1]
                    elif X_fluid[1] == 0:
                        wtm[1, i] = wtptoxides[self.volatile_species[1]]
                        wtm[0, i] = self.calculate_dissolved_volatiles(
                            pressure=pressures[i], sample=sample,
                            X_fluid=X_fluid, **kwargs)[0]
                    else:
                        wtm[:, i] = self.calculate_dissolved_volatiles(
                            pressure=pressures[i], sample=sample,
                            X_fluid=X_fluid, **kwargs)

                    sample.change_composition({
                        self.volatile_species[0]: (wtm[0, i] +
                                                   (1-fractionate_vapor) *
                                                   (wtm0s-wtm[0, i])),
                        self.volatile_species[1]: (wtm[1, i] +
                                                   (1-fractionate_vapor) *
                                                   (wtm1s-wtm[1, i]))})

            except Exception:
                Xv[:, i] = [np.nan]*np.shape(Xv)[0]
                wtm[:, i] = wtm[:, i-1]

        if return_dfs:
            exsolved_degassing_df = pd.DataFrame()
            exsolved_degassing_df['Pressure_bars'] = pressures
            exsolved_degassing_df['H2O_liq'] = (
                                    wtm[self.volatile_species.index('H2O'), :])
            exsolved_degassing_df['CO2_liq'] = (
                                    wtm[self.volatile_species.index('CO2'), :])
            exsolved_degassing_df['H2O_fl'] = (
                                     Xv[self.volatile_species.index('H2O'), :])
            exsolved_degassing_df['CO2_fl'] = (
                                     Xv[self.volatile_species.index('CO2'), :])
            exsolved_degassing_df['FluidProportion_wt'] = (
                (wtm0s+wtm1s) - exsolved_degassing_df['H2O_liq'] -
                exsolved_degassing_df['CO2_liq'])

            if (round_to_zero is True and np.round(
                  exsolved_degassing_df.loc[0, 'FluidProportion_wt'], 2) == 0):
                exsolved_degassing_df.loc[0, 'FluidProportion_wt'] = 0.0

            return exsolved_degassing_df

        else:
            return (wtm, Xv)

    def root_saturation_pressure(self, x, volatile_concs, sample, kwargs):
        """ Function called by scipy.root when finding the saturation pressure
        using calculate_saturation_pressure.

        Parameters
        ----------
        x     numpy array
            The guessed value for the root. x[0] is the pressure (in bars) and
            x[1] is the mole fraction of the first volatile in
            self.volatile_species.
        volatile_concs     numpy array
            The dissolved volatile concentrations, in the same order as
            self.volatile_species.
        sample:     Sample class
            Magma major element composition (including volatiles).
        kwargs     dictionary
            Dictionary of keyword arguments, which may be required by the
            pure-fluid models.

        Returns
        -------
        numpy array
            The difference in the dissolved volatile concentrations, and those
            predicted with the pressure and fluid composition specified by x.
        """
        if x[1] < 0:
            x[1] = 0
        elif x[1] > 1:
            x[1] = 1
        if x[0] <= 0:
            x[0] = 1e-15
        misfit = (np.array(self.calculate_dissolved_volatiles(
            pressure=x[0], X_fluid=(x[1], 1-x[1]), sample=sample, **kwargs))
            - volatile_concs)
        return misfit

    def root_for_fluid_comp(self, Xv0, pressure, Xt0, Xt1, sample, kwargs):
        """ Function called by scipy.root_scalar when calculating the
        composition of equilibrium fluid in the
        calculate_equilibrium_fluid_comp method.

        Parameters
        ----------
        Xv0     float
            The guessed mole fraction of the first volatile species in
            self.volatile_species.
        pressure     float
            The total pressure in bars.
        Xt0     float
            The total mole fraction of the first volatile species in
            self.volatile_species.
        Xt1        float
            The total mole fraction of the second volatile species in
            self.volatile_species.
        sample     Sample class
            Magma major element composition.
        kwargs     dictionary
            A dictionary of keyword arguments that may be required by the pure
            fluid models.

        Returns
        -------
        float
            The differene in the LHS and RHS of the mass balance equation. Eq X
            in manuscript.
        """

        wtt0 = sample.get_composition(self.volatile_species[0])
        wtt1 = sample.get_composition(self.volatile_species[1])

        wtm0, wtm1 = self.calculate_dissolved_volatiles(
            pressure=pressure, X_fluid=(Xv0, 1-Xv0), sample=sample, **kwargs)

        Xm0 = Xt0 / wtt0 * wtm0
        Xm1 = Xt1 / wtt1 * wtm1

        if self.volatile_species[0] == 'CO2' and Xv0 != Xm0:
            f = (Xt0 - Xm0) / (Xv0 - Xm0)
            return (1 - f) * Xm1 + f * (1 - Xv0) - Xt1
        else:
            f = (Xt1 - Xm1) / ((1 - Xv0) - Xm1)
            return (1 - f) * Xm0 + f * Xv0 - Xt0

    def check_calibration_range(self, parameters, report_nonexistance=True):
        """ Checks whether the given parameters are within the ranges defined
        by the CalibrationRange objects for each model and its fugacity and
        activity models. An empty string will be returned if all parameters are
        within the calibration range. If a parameter is not within the
        calibration range, a description of the problem will be returned in the
        string.

        Parameters
        ----------
        parameters     dict
            Dictionary keys are the names of the parameters to be checked,
            e.g., pressure temperature, SiO2, etc. Values are the values of
            each parameter. A complete set need not be given.

        Returns
        -------
        str
            String description of any parameters falling outside of the
            calibration range.
        """
        s = ''
        for model in self.models:
            for cr in model.calibration_ranges:
                if cr.check(parameters) is False:
                    s += cr.string(parameters, report_nonexistance)
            for cr in model.fugacity_model.calibration_ranges:
                if cr.check(parameters) is False:
                    s += cr.string(parameters, report_nonexistance)
            for cr in model.activity_model.calibration_ranges:
                if cr.check(parameters) is False:
                    s += cr.string(parameters, report_nonexistance)
        return s

    def get_calibration_range(self):
        """ Returns a string describing the calibration ranges defined by the
        CalibrationRange objects for each model, and its associated fugacity
        and activity models.

        Returns
        -------
        str
            String description of the calibration range objects."""
        s = ''
        for model in self.models:
            for cr in model.calibration_ranges:
                s += cr.string(None)
            for cr in model.fugacity_model.calibration_ranges:
                s += cr.string(None)
            for cr in model.activity_model.calibration_ranges:
                s += cr.string(None)
        return s
