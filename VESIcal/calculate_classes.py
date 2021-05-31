from abc import abstractmethod
import warnings as w
import numpy as np

from VESIcal import core
from VESIcal import models
from VESIcal.models import magmasat

from copy import deepcopy


class Calculate(object):
    """ The Calculate object is a template for implementing user-friendly
    methods for running calculations using the volatile solubility models.
    All Calculate methods have a common workflow- sample is read in,
    preprocessed, the calculation is performed, the calibration range is
    checked, and the results stored.
    """
    def __init__(self, sample, model='MagmaSat', silence_warnings=False,
                 **kwargs):
        """
        Initializes the calculation.

        Parameters
        ----------
        sample:    Sample class
            The rock composition as a Sample object.
        model:     string or Model class
            Which model to use for the calculation. If passed a string, it
            will look up the name in the default_models dictionary. Default is
            MagmaSat.
        silence_warnings:     bool
            Silence warnings about calibration ranges. Default is False.
        preprocess_sample:     bool
            Before running the calculation, run the sample through the
            preprocessing routine. As of Feb 2021 this functionality should be
            redundant.
        """
        self.model_name = model
        if model == 'MagmaSat':
            self.model = magmasat.MagmaSat()
        elif type(model) == str:
            if model in models.default_models.keys():
                self.model = models.default_models[model]
            else:
                raise core.InputError("The model name given is not recognised."
                                      " Run the method get_model_names() to "
                                      "find allowed names.")
        else:
            self.model = model

        self.sample = sample

        self.result = self.calculate(sample=self.sample, **kwargs)
        self.calib_check = self.check_calibration_range(sample=self.sample,
                                                        **kwargs)

        if self.calib_check is not None and silence_warnings is False:
            if self.calib_check != '':
                w.warn(self.calib_check, RuntimeWarning)

    @abstractmethod
    def calculate(self):
        """ """

    @abstractmethod
    def check_calibration_range(self):
        """ """


class calculate_dissolved_volatiles(Calculate):
    """ Calculates the dissolved volatile concentration using a chosen model
    (default is MagmaSat). Using this interface will preprocess the sample,
    run the calculation, and then check the calibration ranges. All parameters
    required by the chosen model must be passed.

    Parameters
    ----------
    sample:    Sample class
        The rock composition as a Sample object.
    pressure:   float
        Total pressure in bars.
    model:  string or Model object
        Model to be used. If using one of the default models, this can be
        the string corresponding to the model in the default_models dict.
    silence_warnings     bool
        If set to True, no warnings will be raised automatically when
        calibration checks fail.
    preprocess_sample     bool
        If True (default), the sample will be preprocessed according to the
        preprocessing operations within the models. If you obtain unexpected
        results, try setting to False.

    Returns
    -------
    Calculate object
        Calculate object, access results by fetching the result property.
        Dissolved volatile concentrations (in wt%), in order (CO2, H2O, if
        using a mixed fluid default model).
    """
    def return_default_units(self, sample, calc_result, **kwargs):
        """ Checkes the default units set for the sample_class.Sample
        object and returns the result of a calculation in those units.

        Parameters
        ----------
        sample: Sample class
            The rock composition as a Sample object.

        calc_result: dict or float
            Result of a calculate_dissolved_volatiles() calculation on a
            sample.
        """
        default_units = sample.default_units

        # get the composition of
        bulk_comp = deepcopy(sample)

        # check if calculation result is H2O-only, CO2-only, or mixed
        # H2O-CO2
        if isinstance(self.model_name, str):
            if self.model_name in models.get_model_names(model='mixed'):
                # set dissolved H2O and CO2 values.
                # this action assumes they are input as wt% but updates the
                # composition in its default units.
                bulk_comp.change_composition({'H2O': calc_result['H2O_liq'],
                                             'CO2': calc_result['CO2_liq']})
                return {'H2O_liq': bulk_comp.get_composition(
                                                      species='H2O',
                                                      units=default_units),
                        'CO2_liq': bulk_comp.get_composition(
                                                      species='CO2',
                                                      units=default_units)}
            elif 'Water' in self.model_name:
                bulk_comp.change_composition({'H2O': calc_result})
                return bulk_comp.get_composition(species='H2O',
                                                 units=default_units)
            elif 'Carbon' in self.model_name:
                bulk_comp.change_composition({'CO2': calc_result})
                return bulk_comp.get_composition(species='CO2',
                                                 units=default_units)
            elif self.model_name == 'MagmaSat':
                bulk_comp.change_composition({'H2O': calc_result['H2O_liq'],
                                             'CO2': calc_result['CO2_liq']})
                # check if verbose method has been chosen
                if 'verbose' in kwargs and kwargs['verbose']:
                    return {'H2O_liq': bulk_comp.get_composition(
                                                      species='H2O',
                                                      units=default_units),
                            'CO2_liq': bulk_comp.get_composition(
                                                      species='CO2',
                                                      units=default_units),
                            'XH2O_fl': calc_result['XH2O_fl'],
                            'XCO2_fl': calc_result['XCO2_fl'],
                            'FluidProportion_wt': calc_result[
                                                     'FluidProportion_wt']}
                else:
                    return {'H2O_liq': bulk_comp.get_composition(
                                                      species='H2O',
                                                      units=default_units),
                            'CO2_liq': bulk_comp.get_composition(
                                                      species='CO2',
                                                      units=default_units)}
        else:
            # if self.model_name is not a string, most likely this is a
            # user-created model class, and so we cannot interrogate it for
            # model type (H2O, CO2, or mixed).
            # TODO: create model type parameter so that we can ID model type
            # this way instead of via strings in a list!
            return calc_result

    def calculate(self, sample, pressure, **kwargs):
        dissolved = self.model.calculate_dissolved_volatiles(pressure=pressure,
                                                             sample=sample,
                                                             returndict=True,
                                                             **kwargs)
        dissolved_default_units = self.return_default_units(sample, dissolved,
                                                            **kwargs)
        return dissolved_default_units

    def check_calibration_range(self, sample, pressure, **kwargs):
        parameters = kwargs
        parameters['sample'] = sample
        parameters.update(dict(sample.get_composition(units='wtpt_oxides',
                                                      normalization='none')))
        parameters['pressure'] = pressure
        if len(self.model.volatile_species) == 1:
            volspec = self.model.volatile_species[0]
            volconc = self.result
            parameters.update({volspec: volconc})
        else:
            parameters.update(self.result)

        calib_check = self.model.check_calibration_range(parameters)
        return calib_check


class calculate_equilibrium_fluid_comp(Calculate):
    """ Calculates the equilibrium fluid composition using a chosen model
    (default is MagmaSat). Using this interface will preprocess the sample,
    run the calculation, and then check the calibration ranges. All parameters
    required by the chosen model must be passed.

    Parameters
    ----------
    sample:    Sample class
        The rock composition as a Sample object.
    pressure:   float or None
        Total pressure in bars. If None, the saturation pressure will be used.
    model:  string or Model object
        Model to be used. If using one of the default models, this can be
        the string corresponding to the model in the default_models dict.
    silence_warnings     bool
        If set to True, no warnings will be raised automatically when
        calibration checks fail.
    preprocess_sample     bool
        If True (default), the sample will be preprocessed according to the
        preprocessing operations within the models. If you obtain unexpected
        results, try setting to False.

    Returns
    -------
    Calculate object
        Calculate object, access result by fetching the result property. Mole
        fractions of each volatile species, in order (CO2, then H2O, if using
        a mixed-fluid default model).
    """
    def calculate(self, sample, pressure=None, **kwargs):
        if pressure is None:
            pressure = float(self.model.calculate_saturation_pressure(
                                                                 sample=sample,
                                                                 verbose=False,
                                                                 **kwargs))

        fluid_comp = self.model.calculate_equilibrium_fluid_comp(
                                                             pressure=pressure,
                                                             sample=sample,
                                                             **kwargs)
        return fluid_comp

    def check_calibration_range(self, sample, pressure=None, **kwargs):
        if pressure is None:
            pressure = self.model.calculate_saturation_pressure(sample=sample,
                                                                **kwargs)
        parameters = kwargs
        parameters.update(dict(sample.get_composition(units='wtpt_oxides',
                                                      normalization='none')))
        parameters['sample'] = sample
        parameters['pressure'] = pressure
        if len(self.model.volatile_species) == 1:
            volspec = self.model.volatile_species
            volconc = {volspec[0]: self.result}
            parameters.update(volconc)
        elif type(self.model.volatile_species) == list:
            parameters.update(self.result)

        calib_check = self.model.check_calibration_range(parameters)
        return calib_check


class calculate_isobars_and_isopleths(Calculate):
    """ Calculates isobars and isopleths using a chosen model (default is
    MagmaSat). Using this interface will preprocess the sample, run the
    calculation, and then check the calibration ranges. All parameters
    required by the chosen model must be passed.

    Parameters
    ----------
    sample:    Sample class
        The rock composition as a Sample object.
    pressure_list:   list
        List of all pressure values at which to calculate isobars, in bars.
    isopleth_list:   list
        OPTIONAL: Default value is None, in which case only isobars will be
        calculated. List of all fluid compositions in mole fraction (of the
        first species in self.volatile_species) at which to calcualte
        isopleths. Values can range from 0 to 1.
    points:     int
        The number of points in each isobar and isopleth. Default value is 101.
    model:  string or Model object
        Model to be used. If using one of the default models, this can be
        the string corresponding to the model in the default_models dict.
    silence_warnings     bool
        If set to True, no warnings will be raised automatically when
        calibration checks fail.
    preprocess_sample     bool
        If True (default), the sample will be preprocessed according to the
        preprocessing operations within the models. If you obtain unexpected
        results, try setting to False.

    Returns
    -------
    Calculate object
        Calculate object, access results by fetching the result property.
        If isopleth_list is not None, two objects will be returned, one with
        the isobars and the second with the isopleths. If return_dfs is True,
        two pandas DataFrames will be returned with column names
        'Pressure' or 'XH2O_fl', 'H2O_liq', and 'CO2_liq'. If return_dfs is
        False, two lists of numpy arrays will be returned. Each array is an
        individual isobar or isopleth, in the order passed via pressure_list
        or isopleth_list. The arrays are the concentrations of H2O and CO2 in
        the liquid, in the order of the species in self.volatile_species.
    """
    def calculate(self, sample, pressure_list, isopleth_list=[0, 1],
                  points=101, **kwargs):
        check = getattr(self.model, "calculate_isobars_and_isopleths", None)
        if callable(check):
            # samplenorm = sample.copy()
            # samplenorm = normalize_AdditionalVolatiles(samplenorm)
            isobars, isopleths = self.model.calculate_isobars_and_isopleths(
                               sample=self.sample, pressure_list=pressure_list,
                               isopleth_list=isopleth_list, points=points,
                               **kwargs)
            return isobars, isopleths
        else:
            raise core.InputError("This model does not have a "
                                  "calculate_isobars_and_isopleths method "
                                  "built in, most likely because it is a pure "
                                  "fluid model.")

    def check_calibration_range(self, sample, pressure_list, **kwargs):
        parameters = kwargs
        parameters.update(dict(sample.get_composition(units='wtpt_oxides')))
        parameters['sample'] = sample
        s = ''
        s += self.model.check_calibration_range(parameters)
        parameters = {}
        if isinstance(pressure_list, list):
            pass
        else:
            pressure_list = [pressure_list]
        for pressure in pressure_list:
            parameters['pressure'] = pressure
            s += self.model.check_calibration_range(parameters,
                                                    report_nonexistance=False)
        return s


class calculate_saturation_pressure(Calculate):
    """
    Calculates the pressure at which a fluid will be saturated, given the
    dissolved volatile concentrations. Using this interface will preprocess
    the sample, run the calculation, and then check the calibration ranges.
    All parameters required by the chosen model must be passed.

    Parameters
    ----------
    sample:    Sample class
        The rock composition as a Sample object.
    model:  string or Model object
        Model to be used. If using one of the default models, this can be
        the string corresponding to the model in the default_models dict.
    silence_warnings     bool
        If set to True, no warnings will be raised automatically when
        calibration checks fail.
    preprocess_sample     bool
        If True (default), the sample will be preprocessed according to the
        preprocessing operations within the models. If you obtain unexpected
        results, try setting to False.

    Returns
    -------
    Calculate object
        Calculate object, access results by fetching the result property.
        The saturation pressure in bars as a float.
    """
    def calculate(self, sample, **kwargs):
        satP = self.model.calculate_saturation_pressure(sample=sample,
                                                        **kwargs)
        return satP

    def check_calibration_range(self, sample, **kwargs):
        parameters = kwargs
        if isinstance(self.result, dict):  # handles cases where verbose=True
            parameters['pressure'] = next(iter(self.result.values()))
        else:
            parameters['pressure'] = self.result
        parameters.update(dict(sample.get_composition(units='wtpt_oxides')))
        parameters['sample'] = sample
        s = self.model.check_calibration_range(parameters)
        return s


class calculate_degassing_path(Calculate):
    """
    Calculates the dissolved volatiles in a progressively degassing sample.

    Parameters
    ----------
    sample:    Sample class
        The rock composition as a Sample object.
    pressure     string, float, int, list, or numpy array
        Defaults to 'saturation', the calculation will begin at the saturation
        pressure. If a number is passed as either a float or int, this will be
        the starting pressure. If a list of numpy array is passed, the
        pressure values in the list or array will define the degassing path,
        i.e. final_pressure and steps variables will be ignored. Units are
        bars.
    fractionate_vapor     float
        What proportion of vapor should be removed at each step. If 0.0
        (default), the degassing path will correspond to closed-system
        degassing. If 1.0, the degassing path will correspond to open-system
        degassing.
    final_pressure         float
        The final pressure on the degassing path, in bars. Ignored if a list
        or numpy array is passed as the pressure variable. Default is 1 bar.
    steps     int
        The number of steps in the degassing path. Ignored if a list or
        numpy array are passed as the pressure variable.
    model:  string or Model object
        Model to be used. If using one of the default models, this can be
        the string corresponding to the model in the default_models dict.
    silence_warnings     bool
        If set to True, no warnings will be raised automatically when
        calibration checks fail.
    preprocess_sample     bool
        If True (default), the sample will be preprocessed according to the
        preprocessing operations within the models. If you obtain unexpected
        results, try setting to False.

    Returns
    -------
    Calculate object
        Calculate object, access results by fetching the result property.
        A DataFrame with columns 'Pressure', 'H2O_liq', 'CO2_liq',
        'H2O_fl', 'CO2_fl', and 'FluidProportion_wt', is returned. Dissolved
        volatiles are in wt%, the proportions of volatiles in the fluid are
        in mole fraction.
    """
    def calculate(self, sample, pressure='saturation', fractionate_vapor=0.0,
                  final_pressure=100.0, **kwargs):
        check = getattr(self.model, "calculate_degassing_path", None)
        if callable(check):
            data = self.model.calculate_degassing_path(
                                           sample=sample, pressure=pressure,
                                           fractionate_vapor=fractionate_vapor,
                                           final_pressure=final_pressure,
                                           **kwargs)
            return data
        else:
            raise core.InputError("This model does not have a "
                                  "calculate_isobars_and_isopleths method "
                                  "built in, most likely because it is a pure "
                                  "fluid model.")

    def check_calibration_range(self, sample, **kwargs):
        parameters = kwargs
        parameters.update(dict(sample.get_composition()))
        parameters['sample'] = sample
        s = self.model.check_calibration_range(parameters)
        parameters = {}
        parameters['pressure'] = np.nanmax(self.result.Pressure_bars)
        s += self.model.check_calibration_range(parameters,
                                                report_nonexistance=False)
        return s
