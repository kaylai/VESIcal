from VESIcal import activity_models
from VESIcal import calibration_checks
from VESIcal import core
from VESIcal import fugacity_models
from VESIcal import model_classes
from VESIcal import sample_class

import numpy as np
import warnings as w
from scipy.optimize import root_scalar


class water(model_classes.Model):
    """
    Implementation of the Moore et al. (1998) H2O solubility model for magmas up to 3,000 bars.
    """

    def __init__(self):
        """
        Initialize the model.
        """
        self.set_volatile_species(['H2O'])
        self.set_fugacity_model(fugacity_models.fugacity_HB_h2o())
        self.set_activity_model(activity_models.activity_idealsolution())
        self.set_solubility_dependence(False)
        self.set_calibration_ranges([
            calibration_checks.CalibrationRange(
                'pressure', 3000.0, calibration_checks.crf_LessThan, 'bar',
                'Moore et al. (1998) water',
                fail_msg=calibration_checks.crmsg_LessThan_fail,
                pass_msg=calibration_checks.crmsg_LessThan_pass,
                description_msg=calibration_checks.crmsg_LessThan_description),
            calibration_checks.CalibrationRange(
                'temperature', [700.0, 1200], calibration_checks.crf_Between, 'oC',
                'Moore et al. (1998) water',
                fail_msg=calibration_checks.crmsg_Between_fail,
                pass_msg=calibration_checks.crmsg_Between_pass,
                description_msg=calibration_checks.crmsg_Between_description)])

    def calculate_dissolved_volatiles(self, sample, pressure, temperature, X_fluid=1.0, **kwargs):
        """
        Parameters
        ----------
        sample:     Sample class
            Magma major element chemistry.

        pressure float
            Pressure in bars.

        temperature float
            Temperature in degrees C.

        X_fluid float
            OPTIONAL. Default is 1.0. Mole fraction of H2O in the H2O-CO2 fluid.

        Returns
        -------
        float
            Calculated dissolved H2O concentration in wt%.
        """

        _sample = sample.change_composition({'H2O': 0.0, 'CO2': 0.0}, inplace=False)

        fH2O = self.fugacity_model.fugacity(pressure=pressure, temperature=temperature,
                                            X_fluid=X_fluid, **kwargs)
        aParam = 2565.0
        bParam_Al2O3 = -1.997
        bParam_FeOt = -0.9275
        bParam_Na2O = 2.736
        cParam = 1.171
        dParam = -14.21

        temperatureK = temperature + 273.15

        sample_molfrac = _sample.get_composition(units='mol_oxides')
        FeOtot = sample_molfrac['FeO'] + sample_molfrac['Fe2O3']*0.8998

        b_x_sum = ((bParam_Al2O3 * sample_molfrac['Al2O3']) + (bParam_FeOt * FeOtot) +
                   (bParam_Na2O * sample_molfrac['Na2O']))
        two_ln_XH2Omelt = ((aParam / temperatureK) + b_x_sum * (pressure/temperatureK) +
                           cParam * np.log(fH2O) + dParam)
        ln_XH2Omelt = two_ln_XH2Omelt / 2.0
        XH2Omelt = np.exp(ln_XH2Omelt)
        sample_molfrac['H2O'] = XH2Omelt

        # Normalize mol fractions to sum to 1, while preserving XH2O
        sample_molfrac = dict(sample_molfrac)
        for key, value in sample_molfrac.items():
            if key != 'H2O':
                sample_molfrac.update({key: value/((1/(1-sample_molfrac['H2O'])))})

        _sample.change_composition(sample_molfrac, units='mol_oxides')

        return _sample.get_composition('H2O')

    def calculate_equilibrium_fluid_comp(self, sample, pressure, temperature, **kwargs):
        """
        Parameters
        ----------
        sample:     Sample class
            Magma major element composition.

        pressure float
            Pressure in bars.

        temperature float
            Temperature in degrees C.

        Returns
        -------
        float
            Calculated equilibrium fluid concentration in XH2Ofluid mole fraction.
        """

        sample_anhy = sample.change_composition({'H2O': 0.0, 'CO2': 0.0}, inplace=False)

        aParam = 2565.0
        bParam_Al2O3 = -1.997
        bParam_FeOt = -0.9275
        bParam_Na2O = 2.736
        cParam = 1.171
        dParam = -14.21

        temperatureK = temperature + 273.15

        sample_molfrac_anhy = sample_anhy.get_composition(units='mol_oxides')
        sample_molfrac_hy = sample.get_composition(units='mol_oxides')
        FeOtot = sample_molfrac_anhy['FeO'] + sample_molfrac_anhy['Fe2O3']*0.8998

        b_x_sum = ((bParam_Al2O3 * sample_molfrac_anhy['Al2O3']) + (bParam_FeOt * FeOtot) +
                   (bParam_Na2O * sample_molfrac_anhy['Na2O']))
        ln_fH2O = ((2 * np.log(sample_molfrac_hy['H2O']) - (aParam/temperatureK) -
                   b_x_sum * (pressure/temperatureK) - dParam) / cParam)
        fH2O = np.exp(ln_fH2O)
        XH2O_fl = fH2O / pressure

        return XH2O_fl

    def calculate_saturation_pressure(self, temperature, sample, X_fluid=1.0, **kwargs):
        """
        Calculates the pressure at which a an H2O-bearing fluid is saturated. Calls the
        scipy.root_scalar routine, which makes repeated called to the
        calculate_dissolved_volatiles method.

        Parameters
        ----------
        sample:     Sample class
            Magma major element composition.

        temperature float
            Temperature in degrees C.

        X_fluid float
            OPTIONAL. Default is 1.0. Mole fraction of H2O in the H2O-CO2 fluid.

        Returns
        -------
        float
            Calculated saturation pressure in bars.
        """

        temperatureK = temperature + 273.15
        if temperatureK <= 0.0:
            raise core.InputError("Temperature must be greater than 0K.")
        if X_fluid < 0 or X_fluid > 1:
            raise core.InputError("X_fluid must have a value between 0 and 1.")
        if isinstance(sample, sample_class.Sample) is False:
            raise core.InputError("Sample must be an instance of the Sample class.")
        if sample.check_oxide('H2O') is False:
            raise core.InputError("sample must contain H2O.")
        if sample.get_composition('H2O') < 0.0:
            raise core.InputError("Dissolved H2O concentration must be greater than 0 wt%.")

        try:
            satP = root_scalar(self.root_saturation_pressure,
                               args=(temperature, sample, X_fluid, kwargs),
                               x0=100.0, x1=2000.0).root
        except Exception:
            w.warn("Saturation pressure not found.", RuntimeWarning, stacklevel=2)
            satP = np.nan
        return np.real(satP)

    def root_saturation_pressure(self, pressure, temperature, sample, X_fluid, kwargs):
        """ Function called by scipy.root_scalar when finding the saturation pressure using
        calculate_saturation_pressure.

        Parameters
        ----------
        pressure     float
            Pressure guess in bars
        temperature     float
            The temperature of the system in C.
        sample         Sample class
            Magma major element composition.
        kwargs         dictionary
            Additional keyword arguments supplied to calculate_saturation_pressure. Might be
            required for the fugacity or activity models.

        Returns
        -------
        float
            The differece between the dissolved H2O at the pressure guessed, and the H2O
            concentration passed in the sample variable.
        """
        return (self.calculate_dissolved_volatiles(pressure=pressure, temperature=temperature,
                                                   sample=sample, X_fluid=X_fluid, **kwargs) -
                sample.get_composition('H2O'))
