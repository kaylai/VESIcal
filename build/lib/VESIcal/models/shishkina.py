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
    """ Implementation of the Shishkina et al. (2014) carbon solubility model, as a Model class.
    """
    def __init__(self):
        self.set_volatile_species(['CO2'])
        self.set_fugacity_model(fugacity_models.fugacity_idealgas())
        self.set_activity_model(activity_models.activity_idealsolution())
        self.set_solubility_dependence(False)
        self.set_calibration_ranges([
            calibration_checks.CalibrationRange(
                'pressure', [500.0, 5000.0], calibration_checks.crf_Between, 'bar',
                'Shishkina et al. carbon',
                fail_msg=calibration_checks.crmsg_Between_fail,
                pass_msg=calibration_checks.crmsg_Between_pass,
                description_msg=calibration_checks.crmsg_Between_description),
            calibration_checks.CalibrationRange(
                'temperature', [1200.0, 1300.0], calibration_checks.crf_Between, 'oC',
                'Shishkina et al. carbon',
                fail_msg=calibration_checks.crmsg_Between_fail,
                pass_msg=calibration_checks.crmsg_Between_pass,
                description_msg=calibration_checks.crmsg_Between_description),
            calibration_checks.CalibrationRange(
                'SiO2', [40, 57], calibration_checks.crf_Between, 'wt%',
                'Shishkina et al. carbon',
                fail_msg=crmsg_BC_fail_C1,
                pass_msg=calibration_checks.crmsg_Between_pass,
                description_msg=calibration_checks.crmsg_Between_description)])

    def PiStar(self, sample):
        """Shishkina et al. (2014) Eq (11)

        Calculates the Pi* parameter for use in calculating CO2 solubility.

        Parameters
        ----------
        sample:        Sample class
            The magma composition stored in a Sample class.

        Returns
        -------
        float
            The value of the Pi* compositional parameter.
        """

        _mols = sample.get_composition(units='mol_cations')

        if all(cation in _mols for cation in ['Ca', 'K', 'Na', 'Mg', 'Fe', 'Si', 'Al']) is False:
            raise core.InputError("To calculate PiStar, values for CaO, K2O, Na2O, MgO, FeO,"
                                  " SiO2, and Al2O3 must be provided in sample.")

        # Calculate assuming all Fe in Fe+2, as done during calibration
        if 'Fe3' in _mols:
            _mols['Fe'] += _mols['Fe3']

        _pi = ((_mols['Ca'] + 0.8*_mols['K'] + 0.7*_mols['Na'] + 0.4*_mols['Mg'] +
               0.4*_mols['Fe'])/(_mols['Si']+_mols['Al']))

        return _pi

    def calculate_dissolved_volatiles(self, pressure, sample, X_fluid=1, **kwargs):
        """ Calculates the dissolved CO2 concentration in wt%, using equation (13) of Shishkina et
        al. (2014).

        Parameters
        ----------
        pressure:    float
            (Total) pressure in bars.
        sample:        Sample class
            Magma composition.
        X_fluid:    float
            The mol-fraction of the fluid that is CO2. Default is 1, i.e. a pure CO2 fluid.

        Returns
        -------
        float
            The dissolved CO2 concentration in wt%.
        """

        if X_fluid < 0 or X_fluid > 1:
            raise core.InputError("X_fluid must have a value between 0 and 1.")
        if pressure < 0:
            raise core.InputError("pressure must be a positive value.")

        PiStar = self.PiStar(sample)
        fugacity = self.fugacity_model.fugacity(pressure=pressure, X_fluid=X_fluid, **kwargs)

        A = 1.150
        B = 6.71
        C = -1.345

        if fugacity == 0:
            return 0
        else:
            return np.exp(A*np.log(fugacity/10)+B*PiStar+C)/1e4

    def calculate_equilibrium_fluid_comp(self, pressure, sample, **kwargs):
        """ Returns 1.0 if a pure CO2 fluid is saturated. Returns 0.0 if a pure CO2 fluid is
        undersaturated.

        Parameters
        ----------
        pressure     float
            The total pressure of the system in bars.
        sample         Sample class
            Magma major element composition.

        Returns
        -------
        float
            1.0 if CO2-fluid saturated, 0.0 otherwise.
        """
        if self.calculate_saturation_pressure(sample=sample, **kwargs) < pressure:
            return 0.0
        else:
            return 1.0

    def calculate_saturation_pressure(self, sample, **kwargs):
        """ Calculates the pressure at which a pure CO2 fluid is saturated, for the given
        sample composition and CO2 concentration. Calls the scipy.root_scalar routine, which makes
        repeated calls to the calculate_dissolved_volatiles method.

        Parameters
        ----------
        sample         Sample class
            Magma major element composition.

        Returns
        -------
        float
            Saturation pressure in bar
        """

        if sample.check_oxide('CO2') is False:
            raise core.InputError("sample must contain CO2.")
        if sample.get_composition('CO2') < 0:
            raise core.InputError("CO2 concentration must be greater than 0 wt%.")

        try:
            satP = root_scalar(self.root_saturation_pressure, bracket=[1e-15, 1e5],
                               args=(sample, kwargs)).root
        except Exception:
            w.warn("Saturation pressure not found.", RuntimeWarning, stacklevel=2)
            satP = np.nan
        return satP

    def root_saturation_pressure(self, pressure, sample, kwargs):
        """ Function called by scipy.root_scalar when finding the saturation pressure using
        calculate_saturation_pressure.

        Parameters
        ----------
        pressure     float
            Pressure guess in bars
        sample         Sample class
            Magma major element composition.
        kwargs         dictionary
            Additional keyword arguments supplied to calculate_saturation_pressure. Might be
            required for the fugacity or activity models.

        Returns
        -------
        float
            The differece between the dissolved CO2 at the pressure guessed, and the CO2
            concentration passed in the sample variable.
        """
        return (self.calculate_dissolved_volatiles(pressure=pressure, sample=sample, **kwargs) -
                sample.get_composition('CO2'))


class water(model_classes.Model):
    """ Implementation of the Shishkina et al. (2014) H2O solubility model as a Model class.
    """
    def __init__(self):
        self.set_volatile_species(['H2O'])
        self.set_fugacity_model(fugacity_models.fugacity_idealgas())
        self.set_activity_model(activity_models.activity_idealsolution())
        self.set_solubility_dependence(False)
        self.set_calibration_ranges([
            calibration_checks.CalibrationRange(
                'pressure', [500.0, 5000.0], calibration_checks.crf_Between, 'bar',
                'Shishkina et al. water',
                fail_msg=calibration_checks.crmsg_Between_fail,
                pass_msg=calibration_checks.crmsg_Between_pass,
                description_msg=calibration_checks.crmsg_Between_description),
            calibration_checks.CalibrationRange(
                'temperature', [1050, 1400], calibration_checks.crf_Between, 'oC',
                'Shishkina et al. water',
                fail_msg=crmsg_BC_fail_T1,
                pass_msg=calibration_checks.crmsg_Between_pass,
                description_msg=calibration_checks.crmsg_Between_description),
            calibration_checks.CalibrationRange(
                'SiO2', 65, calibration_checks.crf_LessThan, 'wt%',
                'Shishkina et al. water',
                fail_msg=crmsg_LessThan_fail_WaterSi,
                pass_msg=calibration_checks.crmsg_LessThan_pass,
                description_msg=calibration_checks.crmsg_LessThan_description),
            calibration_checks.CalibrationRange(
                'SiO2', 40, calibration_checks.crf_GreaterThan, 'wt%',
                'Shishkina et al. water',
                fail_msg=crmsg_GreaterThan_fail_WaterSi,
                pass_msg=calibration_checks.crmsg_GreaterThan_pass,
                description_msg=calibration_checks.crmsg_GreaterThan_description)])

    def calculate_dissolved_volatiles(self, pressure, sample, X_fluid=1.0, **kwargs):
        """Calculates the dissolved H2O concentration using Eqn (9) of Shishkina et al. (2014).

        Parameters
        ----------
        pressure     float
            Total pressure in bars
        sample         Sample class
            Magma major element composition.
        X_fluid     float
            The mol fraction of H2O in the fluid

        Returns
        -------
        float
            The H2O concentration in wt%
        """
        if isinstance(sample, sample_class.Sample) is False:
            raise core.InputError("Sample must be an instance of the Sample class.")

        if sample.check_oxide('Na2O') is False or sample.check_oxide('K2O') is False:
            raise core.InputError("Na2O and K2O must be present in sample.")

        if pressure < 0:
            raise core.InputError("Pressure must be positive.")

        _mols = sample.get_composition(units='mol_cations')
        _mol_volatiles = 0
        if 'H' in _mols:
            _mol_volatiles += _mols['H']
        if 'C' in _mols:
            _mol_volatiles += _mols['C']

        total_alkalis = (_mols['Na'] + _mols['K'])/(1-_mol_volatiles)

        fugacity = self.fugacity_model.fugacity(pressure, X_fluid=X_fluid, **kwargs)

        a = 3.36e-7 * (fugacity/10)**3 - 2.33e-4*(fugacity/10)**2 + 0.0711*(fugacity/10) - 1.1309
        b = -1.2e-5*(fugacity/10)**2 + 0.0196*(fugacity/10)+1.1297

        return a*total_alkalis + b

    def calculate_equilibrium_fluid_comp(self, pressure, sample, **kwargs):
        """ Returns 1.0 if a pure H2O fluid is saturated.
        Returns 0.0 if a pure H2O fluid is undersaturated.

        Parameters
        ----------
        pressure     float
            The total pressure of the system in bars.
        sample       Sample class
            Magma major element composition (including H2O).

        Returns
        -------
        float
            1.0 if H2O-fluid saturated, 0.0 otherwise.
        """
        if self.calculate_saturation_pressure(sample=sample, **kwargs) < pressure:
            return 0.0
        else:
            return 1.0

    def calculate_saturation_pressure(self, sample, **kwargs):
        """ Calculates the pressure at which a pure H2O fluid is saturated, for the given
        sample composition and H2O concentration. Calls the scipy.root_scalar routine, which makes
        repeated calls to the calculate_dissolved_volatiles method.

        Parameters
        ----------
        sample         Sample class
            Magma major element composition (including H2O).

        Returns
        -------
        float
            Saturation pressure in bar
        """
        if sample.check_oxide('H2O') is False:
            raise core.InputError("sample must contain H2O")
        if sample.get_composition('H2O') < 0:
            raise core.InputError("H2O concentration must be greater than 0 wt%.")

        if sample.get_composition('H2O') < self.calculate_dissolved_volatiles(sample=sample,
                                                                              pressure=0,
                                                                              **kwargs):
            return np.nan

        try:
            satP = root_scalar(self.root_saturation_pressure, bracket=[1e-15, 1e5],
                               args=(sample, kwargs)).root
        except Exception:
            w.warn("Saturation pressure not found.", RuntimeWarning, stacklevel=2)
            satP = np.nan
        return satP

    def root_saturation_pressure(self, pressure, sample, kwargs):
        """ Function called by scipy.root_scalar when finding the saturation pressure using
        calculate_saturation_pressure.

        Parameters
        ----------
        pressure     float
            Pressure guess in bars
        sample         Sample class
            Magma major element composition (including H2O).
        kwargs         dictionary
            Additional keyword arguments supplied to calculate_saturation_pressure. Might be
            required for the fugacity or activity models.

        Returns
        -------
        float
            The differece between the dissolved H2O at the pressure guessed, and the H2O
            concentration passed in the sample variable.
        """
        return (self.calculate_dissolved_volatiles(pressure=pressure, sample=sample, **kwargs) -
                sample.get_composition('H2O'))


crmsg_BC_fail_T1 = ("{param_name} ({param_val:.1f} {units}) is outside the calibration range of "
                    "{model_name} ({calib_val0:.1f}-{calib_val1:.1f} {units}). Note, the authors "
                    "recomend that this model is optimally calibrated between 1150-1250C. ")
crmsg_BC_fail_C1 = ("{param_name} ({param_val:.1f} {units}) is outside the calibration range of "
                    "{model_name} (calculated from the max and minimum concentrations in the "
                    "calibration dataset +-5%; {calib_val0:.1f}-{calib_val1:.1f} {units}).")
# This check warns users if SiO2<40 in the water model
crmsg_GreaterThan_fail_WaterSi = ("{param_name}{units} exceeds the lower calibration limit of "
                                  "{model_name} ({calib_val:.1f} {units}, based on the minimum "
                                  "concentration minus 5% in the calibration dataset). ")
# This check warns users if SiO2>65 wt% , the upper limit suggested by Shishkina et al. (2014)
crmsg_LessThan_fail_WaterSi = ("{param_name} ({param_val:.1f} {units}) exceeds the upper limit of"
                               " {calib_val:.1f} {units} suggested by Shishkina et al. for their "
                               "H2O model. ")


mixed = model_classes.MixedFluid({'H2O': water(), 'CO2': carbon()})

mixed.models[0].calibration_ranges.append(
    calibration_checks.CalibrationRange(
        'X_fluid', 0.0, calibration_checks.crf_MixedFluidWarning, '', 'Shishkina et al. (2014)',
        fail_msg=calibration_checks.crmsg_MixedFluidWarning_fail,
        pass_msg=calibration_checks.crmsg_MixedFluidWarning_pass,
        description_msg=calibration_checks.crmsg_MixedFluidWarning_description))
