from VESIcal import activity_models
from VESIcal import calibration_checks
from VESIcal import core
from VESIcal import fugacity_models
from VESIcal import model_classes
from VESIcal import sample_class

import numpy as np
import warnings as w
from scipy.optimize import root_scalar


class carbon(model_classes.Model):
    """
    Implementation of the Allison et al. (2019) CO2 solubility model. Which type of fit, and
    which composition must be selected when the Model is initialized. The fit may be either
    thermodynamic or power-law. The composition may be chosen from sunset, sfvf, erebus, vesuvius,
    etna, or stromboli. Default is the power-law fit to sunset.
    """

    def __init__(self, model_loc='sunset', model_fit='thermodynamic'):
        """
        Initialize the model.

        Parameters
        ----------
        model_fit     str
            Either 'power' for the power-law fits, or 'thermodynamic' for the
            thermodynamic fits.
        model_loc     str
            One of 'sunset', 'sfvf', 'erebus', 'vesuvius', 'etna', 'stromboli'.
        """

        self.set_volatile_species(['CO2'])
        self.set_fugacity_model(fugacity_models.fugacity_HB_co2())
        self.set_activity_model(activity_models.activity_idealsolution())
        self.set_calibration_ranges([
            calibration_checks.CalibrationRange(
                'temperature', 1200, calibration_checks.crf_EqualTo, 'oC',
                'Allison et al. (2019) carbon', fail_msg=crmsg_Temp,
                pass_msg=calibration_checks.crmsg_EqualTo_pass,
                description_msg=calibration_checks.crmsg_EqualTo_description),
            calibration_checks.CalibrationRange(
                'temperature', [1000, 1400], calibration_checks.crf_Between, 'oC',
                'Allison et al. (2019) carbon', fail_msg=crmsg_Between_Temp,
                pass_msg=calibration_checks.crmsg_EqualTo_pass,
                description_msg=calibration_checks.crmsg_EqualTo_description),
            calibration_checks.CalibrationRange(
                'H2O', 0.5, calibration_checks.crf_LessThan, 'wt%', 'Allison et al. (2019) carbon',
                fail_msg=crmsg_H2O, pass_msg=calibration_checks.crmsg_LessThan_pass)])
        self.set_solubility_dependence(False)
        self.model_loc = model_loc
        self.model_fit = model_fit

    def calculate_dissolved_volatiles(self, pressure, temperature=1200, sample=None, X_fluid=1.0,
                                      **kwargs):
        """
        Calclates the dissolved CO2 concentration using (Eqns) 2-7 or 10-11 from Allison et al.
        (2019).

        Parameters
        ----------
        pressure     float
            Pressure in bars.
        temperature     float
            Temperature in C.
        sample      NoneType or Sample class
            Magma major element composition. Not required for this model, therefore None may be
            passed.
        X_fluid     float
            The mole fraction of CO2 in the fluid. Default is 1.0.

        Returns
        -------
        float
            Dissolved CO2 concentration in wt%.
        """
        # temperature = 1200 #temp in degrees C
        temperature = temperature + 273.15  # translate T from C to K

        if pressure < 0.0:
            raise core.InputError("Pressure must be positive.")
        if X_fluid < 0 or X_fluid > 1:
            raise core.InputError("X_fluid must have a value between 0 and 1.")

        if self.model_fit not in ['power', 'thermodynamic']:
            raise core.InputError("model_fit must be one of 'power', or 'thermodynamic'.")
        if self.model_loc not in ['sunset', 'sfvf', 'erebus', 'vesuvius', 'etna', 'stromboli']:
            raise core.InputError("model_loc must be one of 'sunset', 'sfvf', 'erebus', ",
                                  "'vesuvius', 'etna', or 'stromboli'.")

        if pressure == 0:
            return 0

        if self.model_fit == 'thermodynamic':
            P0 = 1000  # bar
            params = dict({'sunset':    [16.4, -14.67],
                           'sfvf':      [15.02, -14.87],
                           'erebus':    [15.83, -14.65],
                           'vesuvius':  [24.42, -14.04],
                           'etna':      [21.59, -14.28],
                           'stromboli': [14.93, -14.68]})

            DV = params[self.model_loc][0]
            lnK0 = params[self.model_loc][1]

            lnK = lnK0 - (pressure-P0)*DV/(10*8.3141*temperature)
            fCO2 = self.fugacity_model.fugacity(pressure=pressure, temperature=temperature-273.15,
                                                X_fluid=X_fluid, **kwargs)
            Kf = np.exp(lnK)*fCO2
            XCO3 = Kf/(1-Kf)

            FWone = 36.594
            wtCO2 = (44.01*XCO3)/((44.01*XCO3)+(1-XCO3)*FWone)*100
            return wtCO2

        if self.model_fit == 'power':
            params = dict({'stromboli': [1.05, 0.883],
                           'etna':      [2.831, 0.797],
                           'vesuvius':  [4.796, 0.754],
                           'sfvf':      [3.273, 0.74],
                           'sunset':    [4.32, 0.728],
                           'erebus':    [5.145, 0.713]})

            fCO2 = self.fugacity_model.fugacity(pressure=pressure, temperature=temperature-273.15,
                                                X_fluid=X_fluid, **kwargs)

            return params[self.model_loc][0]*fCO2**params[self.model_loc][1]/1e4

    def calculate_equilibrium_fluid_comp(self, pressure, sample, temperature=1200, **kwargs):
        """ Returns 1.0 if a pure CO2 fluid is saturated. Returns 0.0 if a pure CO2 fluid is
        undersaturated.

        Parameters
        ----------
        pressure     float
            The total pressure of the system in bars.
        temperature     float
            The temperature of the system in C.
        sample:        Sample class
            Magma major element composition (including H2O).

        Returns
        -------
        float
            1.0 if CO2-fluid saturated, 0.0 otherwise.
        """

        satP = self.calculate_saturation_pressure(temperature=temperature, sample=sample,
                                                  X_fluid=1.0, **kwargs)
        if pressure < satP:
            return 1.0
        else:
            return 0.0

    def calculate_saturation_pressure(self, sample, temperature=1200, X_fluid=1.0, **kwargs):
        """
        Calculates the pressure at which a pure CO2 fluid is saturated, for the given sample
        composition and CO2 concentration. Calls the scipy.root_scalar routine, which makes
        repeated called to the calculate_dissolved_volatiles method.

        Parameters
        ----------
        temperature     float
            The temperature of the system in C.
        sample:        Sample class
            Magma major element composition (including CO2).
        X_fluid     float
            The mole fraction of H2O in the fluid. Default is 1.0.

        Returns
        -------
        float
            Calculated saturation pressure in bars.
        """

        if X_fluid < 0 or X_fluid > 1:
            raise core.InputError("X_fluid must have a value between 0 and 1.")
        if isinstance(sample, sample_class.Sample) is False:
            raise core.InputError("Sample must be an instance of the Sample class.")
        if sample.check_oxide('CO2') is False:
            raise core.InputError("sample must contain CO2.")
        if sample.get_composition('CO2') < 0.0:
            raise core.InputError("Dissolved CO2 concentration must be greater than 0 wt%.")

        try:
            satP = root_scalar(self.root_saturation_pressure,
                               args=(temperature, sample, X_fluid, kwargs),
                               x0=1000.0, x1=2000.0).root
        except Exception:
            w.warn("Saturation pressure not found.", RuntimeWarning, stacklevel=2)
            satP = np.nan
        return satP

    def root_saturation_pressure(self, pressure, temperature, sample, X_fluid, kwargs):
        """ Function called by scipy.root_scalar when finding the saturation pressure using
        calculate_saturation_pressure.

        Parameters
        ----------
        pressure     float
            Pressure guess in bars
        temperature     float
            The temperature of the system in C.
        sample:     Sample class
            Magma major element composition, including CO2.
        kwargs         dictionary
            Additional keyword arguments supplied to calculate_saturation_pressure. Might be
            required for the fugacity or activity models.

        Returns
        -------
        float
            The differece between the dissolved CO2 at the pressure guessed, and the CO2
            concentration passed in the sample variable.
        """
        return (sample.get_composition('CO2') -
                self.calculate_dissolved_volatiles(pressure=pressure, temperature=temperature,
                                                   sample=sample, X_fluid=X_fluid, **kwargs))


# ALLISON COMPOSITIONAL LIMITS -DEFINED AS MIN AND MAX OF CALIBRATION DATASET (-5% AND +5%
# RESPECTIVELY)

def crf_generic(calibval=None, sample=sample_class.Sample({}), bounds={}):
    comp = sample.get_composition(units='wtpt_oxides')
    testresults = []
    for ox in sfvfCompRange:
        testresults.append(comp[ox] >= bounds[ox][0])
        testresults.append(comp[ox] <= bounds[ox][1])
    return all(testresults)


sfvfCompRange = {'SiO2':  [50.0, 55.97],
                 'TiO2':  [1.08, 1.25],
                 'Al2O3': [15.76, 18.15],
                 'FeO':   [7.07, 8.06],
                 'MgO':   [5.89, 7.09],
                 'CaO':   [8.80, 9.91],
                 'Na2O':  [3.05, 3.53],
                 'K2O':   [1.25, 1.49]
                 }


def crf_sfvf(calibval=None, sample=sample_class.Sample({})):
    return crf_generic(calibval, sample, sfvfCompRange)


sunsetCompRange = {'SiO2':  [45.72, 50.62],
                   'TiO2':  [1.75, 1.96],
                   'Al2O3': [15.62, 17.45],
                   'FeO':   [9.13, 10.42],
                   'MgO':   [8.13, 9.19],
                   'CaO':   [9.56, 10.59],
                   'Na2O':  [3.29, 3.64],
                   'K2O':   [0.77, 0.86]
                   }


def crf_sunset(calibval=None, sample={}):
    return crf_generic(calibval, sample, sunsetCompRange)


erebusCompRange = {'SiO2':  [45.96, 51.04],
                   'TiO2':  [2.67, 3.00],
                   'Al2O3': [18.31, 20.63],
                   'FeO':   [7.46, 9.37],
                   'MgO':   [3.03, 3.43],
                   'CaO':   [6.58, 7.42],
                   'Na2O':  [5.8, 6.49],
                   'K2O':   [2.75, 3.13]
                   }


def crf_erebus(calibval=None, sample={}):
    return crf_generic(calibval, sample, erebusCompRange)


vesuviusCompRange = {'SiO2':  [45.65, 53.29],
                     'TiO2':  [0.93, 1.13],
                     'Al2O3': [13.75, 16.28],
                     'FeO':   [4.98, 7.48],
                     'MgO':   [6.41, 7.76],
                     'CaO':   [11.12, 14.16],
                     'Na2O':  [1.74, 2.07],
                     'K2O':   [5.48, 6.35]
                     }


def crf_vesuvius(calibval=None, sample={}):
    return crf_generic(calibval, sample, vesuviusCompRange)


etnaCompRange = {'SiO2':  [46.03, 52.97],
                 'TiO2':  [1.61, 1.89],
                 'Al2O3': [15.87, 18.24],
                 'FeO':   [6.75, 10.21],
                 'MgO':   [5.9, 7.0],
                 'CaO':   [9.38, 11.99],
                 'Na2O':  [3.4, 3.94],
                 'K2O':   [1.69, 2.25]
                 }


def crf_etna(calibval=None, sample={}):
    return crf_generic(calibval, sample, etnaCompRange)


stromboliCompRange = {'SiO2':  [47.23, 55.15],
                      'TiO2':  [0.74, 0.94],
                      'Al2O3': [14.87, 17.73],
                      'FeO':   [5.08, 7.51],
                      'MgO':   [7.52, 9.26],
                      'CaO':   [11.99, 13.46],
                      'Na2O':  [2.29, 2.67],
                      'K2O':   [1.79, 2.17]
                      }


def crf_stromboli(calibval=None, sample={}):
    return crf_generic(calibval, sample, stromboliCompRange)


crmsg_Comp_pass = ("The sample appears to be similar in composition to the compositional dataset "
                   "for the selected Carbon model of Allison et al. (2019).")
crmsg_Comp_fail = (" These calibration limits were selected based on the minimum and maximum "
                   "values of these oxides (+-5%) in the calibration dataset. As the Allison et "
                   "al. model incorperates no term for compositional dependence, users must take "
                   "extreme care when extrapolating this model to compositions which differ "
                   "significantly from the calibration dataset. These warnings are simply a "
                   "guide;  we suggest that users carefully compare their major element data to"
                   " the calibration dataset to check for suitability ")
crmsg_Comp_description = ("The Allison et al. (2019) Carbon model is defined for 6 different "
                          "alkali compositions.")
crmsg_Temp = ("All calculations for {model_name} are performed at 1200 C "
              "(inputted Temp={param_val:.1f} {units}). Allison et al. (2019) suggest the "
              "results are likely applicable between 1000-1400°C). ")
crmsg_PressureGreater = ("{param_name} ({param_val:.1f} {units}) is less than the lowest P "
                         "experiment in the calibration dataset ({calib_val:.1f} {units}). ")
crmsg_PressureLess = ("{param_name} ({param_val:.1f} {units}) exceeds the upper limit of "
                      "{calib_val:.1f} {units} suggested by Allison et al. (2019) - Their "
                      "spreadsheet would return 7000 bar for this input. ")
crmsg_H2O = ("{param_name} ({param_val:.1f} {units}) is > {param_val:.1f} {units}: this model "
             "does not account for the effect of H$_2$O on volatile solubility. VESIcal allows "
             "you to combine Allison Carbon with a variety of H$_2$O models. ")
crmsg_Between_Temp = ("{param_name} ({param_val:.1f} {units} is outside the recomended ",
                      "temperature range for {model_name} (1000-1400°C). ")


# Create objects for each model location in order to set their calibration ranges
sunset = carbon(model_loc='sunset')
sfvf = carbon(model_loc='sfvf')
erebus = carbon(model_loc='erebus')
vesuvius = carbon(model_loc='vesuvius')
etna = carbon(model_loc='etna')
stromboli = carbon(model_loc='stromboli')

_crs_to_update = sunset.calibration_ranges

cr_oxide_list = []
for ox in sunsetCompRange.keys():
    cr_oxide_list.append(
        calibration_checks.CalibrationRange(
                         ox, sunsetCompRange[ox], calibration_checks.crf_Between, 'wt%',
                         'Allison et al. (2019) sunset carbon',
                         fail_msg=calibration_checks.crmsg_BC_fail,
                         pass_msg=calibration_checks.crmsg_BC_pass,
                         description_msg=calibration_checks.crmsg_Between_description))

sunset.set_calibration_ranges(
    _crs_to_update + [
        calibration_checks.CalibrationRange(
            'pressure', 7000, calibration_checks.crf_LessThan, 'bar',
            'Allison et al. (2019) sunset carbon',
            fail_msg=crmsg_PressureLess,
            pass_msg=calibration_checks.crmsg_LessThan_pass),
        calibration_checks.CalibrationRange(
            'pressure', 4071, calibration_checks.crf_GreaterThan, 'bar',
            'Allison et al. (2019) sunset carbon',
            fail_msg=crmsg_PressureGreater,
            pass_msg=calibration_checks.crmsg_GreaterThan_pass),
        calibration_checks.CalibrationRange(
            'sample', None, crf_sunset, None, None,
            fail_msg=crmsg_Comp_fail,
            pass_msg=crmsg_Comp_pass)] + cr_oxide_list)

# SFVF
_crs_to_update = sfvf.calibration_ranges
cr_oxide_list = []
for ox in sfvfCompRange.keys():
    cr_oxide_list.append(
        calibration_checks.CalibrationRange(
            ox, sfvfCompRange[ox], calibration_checks.crf_Between, 'wt%',
            'Allison et al. (2019) sfvf carbon',
            fail_msg=calibration_checks.crmsg_BC_fail,
            pass_msg=calibration_checks.crmsg_BC_pass,
            description_msg=calibration_checks.crmsg_Between_description))

sfvf.set_calibration_ranges(_crs_to_update + [
    calibration_checks.CalibrationRange(
        'pressure', 7000, calibration_checks.crf_LessThan, 'bar',
        'Allison et al. (2019) sfvf carbon',
        fail_msg=crmsg_PressureLess,
        pass_msg=calibration_checks.crmsg_LessThan_pass),
    calibration_checks.CalibrationRange(
        'pressure', 4133, calibration_checks.crf_GreaterThan, 'bar',
        'Allison et al. (2019) sfvf carbon',
        fail_msg=crmsg_PressureGreater,
        pass_msg=calibration_checks.crmsg_GreaterThan_pass),
    calibration_checks.CalibrationRange(
        'sample', None, crf_sfvf, None, None,
        fail_msg=crmsg_Comp_fail,
        pass_msg=crmsg_Comp_pass)] + cr_oxide_list)


# Erebus
_crs_to_update = erebus.calibration_ranges
cr_oxide_list = []
for ox in erebusCompRange.keys():
    cr_oxide_list.append(
        calibration_checks.CalibrationRange(
            ox, erebusCompRange[ox], calibration_checks.crf_Between, 'wt%',
            'Allison et al. (2019) erebus carbon',
            fail_msg=calibration_checks.crmsg_BC_fail,
            pass_msg=calibration_checks.crmsg_BC_pass,
            description_msg=calibration_checks.crmsg_Between_description))

erebus.set_calibration_ranges(_crs_to_update + [
    calibration_checks.CalibrationRange(
        'pressure', 7000, calibration_checks.crf_LessThan, 'bar',
        'Allison et al. (2019) erebus carbon',
        fail_msg=crmsg_PressureLess,
        pass_msg=calibration_checks.crmsg_LessThan_pass),
    calibration_checks.CalibrationRange(
        'pressure', 4078, calibration_checks.crf_GreaterThan, 'bar',
        'Allison et al. (2019) erebus carbon',
        fail_msg=crmsg_PressureGreater,
        pass_msg=calibration_checks.crmsg_GreaterThan_pass),
    calibration_checks.CalibrationRange(
        'sample', None, crf_erebus, None, None,
        fail_msg=crmsg_Comp_fail,
        pass_msg=crmsg_Comp_pass)] + cr_oxide_list)


_crs_to_update = etna.calibration_ranges
cr_oxide_list = []
for ox in etnaCompRange.keys():
    cr_oxide_list.append(
        calibration_checks.CalibrationRange(
            ox, etnaCompRange[ox], calibration_checks.crf_Between, 'wt%',
            'Allison et al. (2019) etna carbon',
            fail_msg=calibration_checks.crmsg_BC_fail,
            pass_msg=calibration_checks.crmsg_BC_pass,
            description_msg=calibration_checks.crmsg_Between_description))


# Etna
etna.set_calibration_ranges(_crs_to_update + [
    calibration_checks.CalibrationRange(
        'pressure', 7000, calibration_checks.crf_LessThan, 'bar',
        'Allison et al. (2019) etna carbon',
        fail_msg=crmsg_PressureLess,
        pass_msg=calibration_checks.crmsg_LessThan_pass),
    calibration_checks.CalibrationRange(
        'pressure', 485, calibration_checks.crf_GreaterThan, 'bar',
        'Allison et al. (2019) etna carbon',
        fail_msg=crmsg_PressureGreater,
        pass_msg=calibration_checks.crmsg_GreaterThan_pass),
    calibration_checks.CalibrationRange(
        'sample', None, crf_etna, None, None,
        fail_msg=crmsg_Comp_fail,
        pass_msg=crmsg_Comp_pass)] + cr_oxide_list)


# Vesuvius
_crs_to_update = vesuvius.calibration_ranges
cr_oxide_list = []
for ox in vesuviusCompRange.keys():
    cr_oxide_list.append(
        calibration_checks.CalibrationRange(
            ox, vesuviusCompRange[ox], calibration_checks.crf_Between, 'wt%',
            'Allison et al. (2019) vesuvius carbon',
            fail_msg=calibration_checks.crmsg_BC_fail,
            pass_msg=calibration_checks.crmsg_BC_pass,
            description_msg=calibration_checks.crmsg_Between_description))

vesuvius.set_calibration_ranges(_crs_to_update + [
    calibration_checks.CalibrationRange(
        'pressure', 7000, calibration_checks.crf_LessThan, 'bar',
        'Allison et al. (2019) vesuvius carbon',
        fail_msg=crmsg_PressureLess,
        pass_msg=calibration_checks.crmsg_LessThan_pass),
    calibration_checks.CalibrationRange(
        'pressure', 269, calibration_checks.crf_GreaterThan, 'bar',
        'Allison et al. (2019) vesuvius carbon',
        fail_msg=crmsg_PressureGreater,
        pass_msg=calibration_checks.crmsg_GreaterThan_pass),
    calibration_checks.CalibrationRange(
        'sample', None, crf_vesuvius, None, None,
        fail_msg=crmsg_Comp_fail,
        pass_msg=crmsg_Comp_pass)] + cr_oxide_list)


# Stromboli
_crs_to_update = stromboli.calibration_ranges
cr_oxide_list = []
for ox in stromboliCompRange.keys():
    cr_oxide_list.append(
        calibration_checks.CalibrationRange(
            ox, stromboliCompRange[ox], calibration_checks.crf_Between, 'wt%',
            'Allison et al. (2019) stromboli carbon',
            fail_msg=calibration_checks.crmsg_BC_fail,
            pass_msg=calibration_checks.crmsg_BC_pass,
            description_msg=calibration_checks.crmsg_Between_description))

stromboli.set_calibration_ranges(_crs_to_update + [
    calibration_checks.CalibrationRange(
        'pressure', 7000, calibration_checks.crf_LessThan, 'bar',
        'Allison et al. (2019) stromboli carbon',
        fail_msg=crmsg_PressureLess,
        pass_msg=calibration_checks.crmsg_LessThan_pass),
    calibration_checks.CalibrationRange(
        'pressure', 524, calibration_checks.crf_GreaterThan, 'bar',
        'Allison et al. (2019) stromboli carbon',
        fail_msg=crmsg_PressureGreater,
        pass_msg=calibration_checks.crmsg_GreaterThan_pass),
    calibration_checks.CalibrationRange(
        'sample', None, crf_stromboli, None, None,
        fail_msg=crmsg_Comp_fail,
        pass_msg=crmsg_Comp_pass)] + cr_oxide_list)
