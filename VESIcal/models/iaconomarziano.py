from VESIcal import activity_models
from VESIcal import calibration_checks
from VESIcal import core
from VESIcal import fugacity_models
from VESIcal import model_classes
from VESIcal import sample_class

import numpy as np
from scipy.optimize import root_scalar
import warnings as w


class water(model_classes.Model):
    """
    Implementation of the Iacono-Marziano et al. (2012) water solubility model, as a Model class.
    Three calibrations are provided- the one incorporating the H2O content as a parameter
    (hydrous), and the one that does not (anhydrous), in addition to the coefficient values given
    incorrect and do not align with the web versions of the model. Which model should be used is
    specified when the methods are called. The default choice is the hydrous model.
    """

    def __init__(self):
        """
        Initialise the model.

        """
        self.set_volatile_species(['H2O'])
        self.set_fugacity_model(fugacity_models.fugacity_idealgas())
        self.set_activity_model(activity_models.activity_idealsolution())
        self.set_calibration_ranges([
            calibration_checks.CalibrationRange(
                'temperature', [1000, 1400], calibration_checks.crf_Between, 'oC',
                'IaconoMarzianoWater',
                fail_msg=crmsg_BC_T,
                pass_msg=calibration_checks.crmsg_Between_pass,
                description_msg=calibration_checks.crmsg_Between_description),
            calibration_checks.CalibrationRange(
                'pressure', [100, 10000], calibration_checks.crf_Between, 'bars',
                'IaconoMarzianoWater',
                fail_msg=crmsg_BC_P,
                pass_msg=calibration_checks.crmsg_Between_pass,
                description_msg=calibration_checks.crmsg_Between_description)])
        # Not dependent on CO2 conc, H2O dependence dealt with within model.
        self.set_solubility_dependence(False)

        # The oxide masses used in the IM webapp.
        self.IM_oxideMasses = {'Al2O3': 101.96,
                               'CaO':    56.08,
                               'FeO':    71.85,
                               'K2O':    94.2,
                               'MgO':    40.32,
                               'Na2O':   61.98,
                               'SiO2':   60.09,
                               'TiO2':   79.9,
                               'H2O':    18.01}

    def calculate_dissolved_volatiles(self, pressure, temperature, sample, X_fluid=1.0,
                                      coeffs='webapp', **kwargs):
        """
        Calculates the dissolved H2O concentration, using Eq (13) of Iacono-Marziano et al. (2012).
        If using the hydrous parameterization, it will use the scipy.root_scalar routine to find
        the root of the root_dissolved_volatiles method.

        Parameters
        ----------
        pressure    float
            Total pressure in bars.
        temperature     float
            Temperature in C
        sample     pandas Series or dict
            Major element oxides in wt%.
        X_fluid      float
            Mole fraction of H2O in the fluid. Default is 1.0.
        coeffs  str
            Which set of coefficients should be used in the calculations:
            - 'webapp' (default) for the hydrous NBO/O parameterisation coefficients used in
              the Iacono-Marziano webapp.
            - 'manuscript' for the hydrous NBO/O parameterisation coefficients given in the
              Iacono-Marziano et al. (2012) manuscript, but were incorrect (pers.comm.).
            - 'anhydrous' for the anhydrous NBO/O parameterisation coefficients given in the
              Iacono-Marziano et al. (2012) manuscript.

        Returns
        -------
        float
            Dissolved H2O concentration in wt%.
        """

        if coeffs not in ['webapp', 'manuscript', 'anhydrous']:
            raise core.InputError("The coeffs argument must be one of 'webapp', 'manuscript', "
                                  "or 'anhydrous'")

        temperature = temperature + 273.15  # translate T from C to K

        if isinstance(sample, sample_class.Sample) is False:
            raise core.InputError("Sample must be an instance of the Sample class.")
        if pressure < 0:
            raise core.InputError("Pressure must be positive.")
        if X_fluid < 0 or X_fluid > 1:
            raise core.InputError("X_fluid must have a value between 0 and 1.")

        if pressure == 0:
            return 0

        if coeffs == 'webapp' or coeffs == 'manuscript':
            if X_fluid == 0:
                return 0
            H2O = root_scalar(
                self.root_dissolved_volatiles,
                args=(pressure, temperature, sample, X_fluid, coeffs, kwargs),
                x0=1.0, x1=2.0).root
            return H2O
        else:
            a = 0.54
            b = 1.24
            B = -2.95
            C = 0.02

            fugacity = self.fugacity_model.fugacity(pressure=pressure, X_fluid=X_fluid,
                                                    temperature=temperature-273.15, **kwargs)
            if fugacity == 0:
                return 0
            NBO_O = self.NBO_O(sample=sample, coeffs=coeffs)

            H2O = np.exp(a*np.log(fugacity) + b*NBO_O + B + C*pressure/temperature)

            return H2O

    def calculate_equilibrium_fluid_comp(self, pressure, temperature, sample, **kwargs):
        """ Returns 1.0 if a pure H2O fluid is saturated. Returns 0.0 if a pure H2O fluid is
        undersaturated.

        Parameters
        ----------
        pressure     float
            The total pressure of the system in bars.
        temperature     float
            The temperature of the system in C.
        sample         pandas Series or dict
            Major element oxides in wt% (including H2O).

        Returns
        -------
        float
            1.0 if H2O-fluid saturated, 0.0 otherwise.
        """

        if pressure > self.calculate_saturation_pressure(temperature=temperature,
                                                         sample=sample, **kwargs):
            return 0.0
        else:
            return 1.0

    def calculate_saturation_pressure(self, temperature, sample, **kwargs):
        """
        Calculates the pressure at which a pure H2O fluid is saturated, for the given sample
        composition and H2O concentration. Calls the scipy.root_scalar routine, which makes
        repeated called to the calculate_dissolved_volatiles method.

        Parameters
        ----------
        temperature     float
            The temperature of the system in C.
        sample         pandas Series or dict
            Major element oxides in wt% (including H2O).
        X_fluid     float
            The mole fraction of H2O in the fluid. Default is 1.0.

        Returns
        -------
        float
            Calculated saturation pressure in bars.
        """

        if isinstance(sample, sample_class.Sample) is False:
            raise core.InputError("Sample must be an instance of the Sample class.")
        if sample.check_oxide('H2O') is False:
            raise core.InputError("sample must contain H2O.")
        if sample.get_composition('H2O') < 0.0:
            raise core.InputError("Dissolved H2O must be greater than 0 wt%.")

        try:
            # Checks whether the upper bound for the numerical solver returns a positive H2O
            # concentration, and if it doesn't, it will progressively decrease the bound until
            # it does.
            upperbound = 1e5
            while self.calculate_dissolved_volatiles(upperbound, temperature,
                                                     sample, **kwargs) < 0:
                upperbound = upperbound*0.9

            satP = root_scalar(self.root_saturation_pressure, args=(temperature, sample, kwargs),
                               bracket=[1e-15, upperbound]).root
        except Exception:
            w.warn("Saturation pressure not found.", RuntimeWarning, stacklevel=2)
            satP = np.nan
        return satP

    def root_saturation_pressure(self, pressure, temperature, sample, kwargs):
        """ Function called by scipy.root_scalar when finding the saturation pressure using
        calculate_saturation_pressure.

        Parameters
        ----------
        pressure     float
            Pressure guess in bars
        temperature     float
            The temperature of the system in C.
        sample         pandas Series or dict
            Major elements in wt% (normalized to 100%), including H2O.
        kwargs         dictionary
            Additional keyword arguments supplied to calculate_saturation_pressure. Might be
            required for the fugacity or activity models.

        Returns
        -------
        float
            The differece between the dissolved H2O at the pressure guessed, and the H2O
            concentration passed in the sample variable.
        """

        return (sample.get_composition('H2O') -
                self.calculate_dissolved_volatiles(pressure=pressure, temperature=temperature,
                                                   sample=sample, **kwargs))

    def root_dissolved_volatiles(self, h2o, pressure, temperature, sample, X_fluid, coeffs,
                                 kwargs):
        """ Function called by calculate_dissolved_volatiles method when the hydrous
        parameterization is being used.

        Parameters
        ----------
        h2o     float
            Guess for the H2O concentration in wt%.
        pressure     float
            Total pressure in bars.
        temperature     float
            Temperature in K.
        sample         pandas Series or dict
            Major element oxides in wt%.
        X_fluid     float
            Mole fraction of H2O in the fluid.
        coeffs  str
            One of 'webapp','manuscript','anhydrous'.
        kwargs     dictionary
            Keyword arguments

        Returns
        -------
        float
            Difference between H2O guessed and the H2O calculated.
        """

        if coeffs == 'manuscript':
            a = 0.53
            b = 2.35
            B = -3.37
            C = -0.02
        else:
            a = 0.52096846
            b = 2.11575907
            B = -3.24443335
            C = -0.02238884

        sample_copy = sample.change_composition({'H2O': h2o}, inplace=False)

        NBO_O = self.NBO_O(sample=sample_copy, coeffs=coeffs)
        fugacity = self.fugacity_model.fugacity(pressure=pressure, X_fluid=X_fluid,
                                                temperature=temperature, **kwargs)

        return h2o - np.exp(a*np.log(fugacity) + b*NBO_O + B + C*pressure/(temperature+273.15))

    def NBO_O(self, sample, coeffs='webapp'):
        """
        Calculates NBO/O according to Appendix A.1. of Iacono-Marziano et al. (2012). NBO/O
        is calculated on either a hydrous or anhyrous basis, as set when initialising the
        Model class.

        Parameters
        ----------
        sample     pandas Series or dict
            Major element oxides in wt% (including H2O if using the hydrous parameterization).

        coeffs  str
            One of:
            - 'webapp' or 'manuscript' to include H2O in NBO/O
            - 'anhydrous' to exclude H2O from NBO/O

        Returns
        -------
        float
            NBO/O.
        """
        if isinstance(sample, sample_class.Sample) is False:
            raise core.InputError("Sample must be an instance of the Sample class.")
        if all(sample.check_oxide(ox) for ox in ['K2O', 'Na2O', 'CaO', 'MgO',
                                                 'FeO', 'Al2O3', 'SiO2', 'TiO2']) is False:
            raise core.InputError("Sample must contain K2O, Na2O, CaO, MgO, FeO, Al2O3, SiO2, "
                                  "and TiO2.")

        X = sample.get_composition(units='mol_oxides', oxide_masses=self.IM_oxideMasses)

        if 'Fe2O3' in X:
            Fe2O3 = X['Fe2O3']
        else:
            Fe2O3 = 0

        NBO = 2*(X['K2O'] + X['Na2O'] + X['CaO'] + X['MgO'] + X['FeO'] + 2*Fe2O3 - X['Al2O3'])
        Ox = (2*X['SiO2'] + 2*X['TiO2'] + 3*X['Al2O3'] + X['MgO'] + X['FeO'] + 2*Fe2O3 +
              X['CaO'] + X['Na2O'] + X['K2O'])

        if coeffs == 'webapp' or coeffs == 'manuscript':
            if 'H2O' not in X:
                raise core.InputError("sample must contain H2O if using the hydrous"
                                      " parameterization.")
            NBO = NBO + 2*X['H2O']
            Ox = Ox + X['H2O']

        return NBO/Ox


class carbon(model_classes.Model):
    """
    Implementation of the Iacono-Marziano et al. (2012) carbon solubility model, as a Model class.
    """

    def __init__(self):
        """
        Initialise the model.

        """
        self.set_volatile_species(['CO2'])
        self.set_fugacity_model(fugacity_models.fugacity_idealgas())
        self.set_activity_model(activity_models.activity_idealsolution())
        self.set_calibration_ranges([
            calibration_checks.CalibrationRange(
                'temperature', [1000, 1400], calibration_checks.crf_Between, 'oC',
                'IaconoMarzianoCarbon',
                fail_msg=crmsg_BC_T,
                pass_msg=calibration_checks.crmsg_Between_pass,
                description_msg=calibration_checks.crmsg_Between_description),
            calibration_checks.CalibrationRange(
                'pressure', [100, 10000], calibration_checks.crf_Between, 'bars',
                'IaconoMarzianoCarbon',
                fail_msg=crmsg_BC_P,
                pass_msg=calibration_checks.crmsg_Between_pass,
                description_msg=calibration_checks.crmsg_Between_description)])
        self.set_solubility_dependence(True)

        # The oxide masses used in the IM webapp.
        self.IM_oxideMasses = {'Al2O3': 101.96,
                               'CaO':    56.08,
                               'FeO':    71.85,
                               'K2O':    94.2,
                               'MgO':    40.32,
                               'Na2O':   61.98,
                               'SiO2':   60.09,
                               'TiO2':   79.9,
                               'H2O':    18.01}

    def calculate_dissolved_volatiles(self, pressure, temperature, sample, X_fluid=1,
                                      coeffs='webapp', **kwargs):
        """
        Calculates the dissolved CO2 concentration, using Eq (12) of Iacono-Marziano et al. (2012).
        If using the hydrous parameterization, it will use the scipy.root_scalar routine to find
        the root of the root_dissolved_volatiles method.

        Parameters
        ----------
        pressure    float
            Total pressure in bars.
        temperature     float
            Temperature in C
        sample     Sample class
            Magma major element composition.
        X_fluid      float
            Mole fraction of H2O in the fluid. Default is 1.0.
        coeffs  str
            Which set of coefficients should be used for H2O calculations:
            - 'webapp' (default) for the hydrous NBO/O parameterisation coefficients used in the
              Iacono-Marziano webapp.
            - 'manuscript' for the hydrous NBO/O parameterisation coefficients given in the
              Iacono-Marziano et al. (2012) manuscript.
            - 'anhydrous' for the anhydrous NBO/O parameterisation coefficients given in the
              Iacono-Marziano et al. (2012) manuscript.

        Returns
        -------
        float
            Dissolved H2O concentration in wt%.
        """

        if coeffs not in ['webapp', 'manuscript', 'anhydrous']:
            raise core.InputError("The coeffs argument must be one of 'webapp', 'manuscript', "
                                  "or 'anhydrous'")

        temperature = temperature + 273.15  # translate T from C to K

        if isinstance(sample, sample_class.Sample) is False:
            raise core.InputError("Sample must be an instance of the Sample class.")
        if pressure < 0:
            raise core.InputError("Pressure must be positive.")
        if temperature <= 0:
            raise core.InputError("Temperature must be greater than 0K.")
        if X_fluid < 0 or X_fluid > 1:
            raise core.InputError("X_fluid must have a value between 0 and 1.")

        if pressure == 0:
            return 0

        if coeffs == 'webapp' or coeffs == 'manuscript':
            im_h2o_model = water()
            h2o = im_h2o_model.calculate_dissolved_volatiles(pressure=pressure,
                                                             temperature=temperature-273.15,
                                                             sample=sample, X_fluid=1-X_fluid,
                                                             coeffs=coeffs, **kwargs)

            sample_h2o = sample.change_composition({'H2O': h2o}, inplace=False)

            d = np.array([-16.4, 4.4, -17.1, 22.8])
            a = 1.0
            b = 17.3
            B = -6.0
            C = 0.12

            NBO_O = self.NBO_O(sample=sample_h2o, coeffs=coeffs)

        else:
            im_h2o_model = water()
            h2o = im_h2o_model.calculate_dissolved_volatiles(pressure=pressure,
                                                             temperature=temperature-273.15,
                                                             sample=sample, X_fluid=1-X_fluid,
                                                             coeffs=coeffs, **kwargs)

            sample_h2o = sample.change_composition({'H2O': h2o}, inplace=False)

            d = np.array([2.3, 3.8, -16.3, 20.1])
            a = 1.0
            b = 15.8
            B = -5.3
            C = 0.14

            NBO_O = self.NBO_O(sample=sample, coeffs=coeffs)

        fugacity = self.fugacity_model.fugacity(pressure=pressure, X_fluid=X_fluid,
                                                temperature=temperature-273.15, **kwargs)

        if fugacity == 0:
            return 0

        molarProps = sample_h2o.get_composition(units='mol_oxides',
                                                oxide_masses=self.IM_oxideMasses)

        if all(ox in molarProps for ox in ['Al2O3', 'CaO', 'K2O', 'Na2O',
                                           'FeO', 'MgO', 'Na2O', 'K2O']) is False:
            raise core.InputError("sample must contain Al2O3, CaO, K2O, Na2O, FeO, MgO, Na2O, "
                                  "and K2O.")
        if 'Fe2O3' in molarProps:
            Fe2O3 = molarProps['Fe2O3']
        else:
            Fe2O3 = 0

        x = list()
        if 'H2O' in molarProps:
            x.append(molarProps['H2O'])
        else:
            x.append(0.0)
        x.append(molarProps['Al2O3']/(molarProps['CaO']+molarProps['K2O']+molarProps['Na2O']))
        x.append((molarProps['FeO']+Fe2O3*2+molarProps['MgO']))
        x.append((molarProps['Na2O']+molarProps['K2O']))
        x = np.array(x)

        CO3 = np.exp(np.sum(x*d) + a*np.log(fugacity) + b*NBO_O + B + C*pressure/temperature)
        CO2 = CO3/1e4

        return CO2

    def calculate_equilibrium_fluid_comp(self, pressure, temperature, sample, **kwargs):
        """ Returns 1.0 if a pure CO2 fluid is saturated. Returns 0.0 if a pure CO2 fluid is
        undersaturated.

        Parameters
        ----------
        pressure     float
            The total pressure of the system in bars.
        temperature     float
            The temperature of the system in C.
        sample         Sample class
            Magma major element composition (including H2O).

        Returns
        -------
        float
            1.0 if CO2-fluid saturated, 0.0 otherwise.
        """

        if pressure > self.calculate_saturation_pressure(temperature=temperature, sample=sample,
                                                         **kwargs):
            return 0.0
        else:
            return 1.0

    def calculate_saturation_pressure(self, temperature, sample, **kwargs):
        """
        Calculates the pressure at which a pure CO2 fluid is saturated, for the given sample
        composition and CO2 concentration. Calls the scipy.root_scalar routine, which makes
        repeated called to the calculate_dissolved_volatiles method.

        Parameters
        ----------
        temperature     float
            The temperature of the system in C.
        sample         Sample class
            Magma major element composition (including CO2).

        Returns
        -------
        float
            Calculated saturation pressure in bars.
        """

        if temperature <= 0:
            raise core.InputError("Temperature must be greater than 0K.")
        if isinstance(sample, sample_class.Sample) is False:
            raise core.InputError("Sample must be an instance of the Sample class.")
        if sample.check_oxide('CO2') is False:
            raise core.InputError("sample must contain CO2")
        if sample.get_composition('CO2') < 0:
            raise core.InputError("Dissolved CO2 must be greater than 0 wt%.")

        try:
            satP = root_scalar(self.root_saturation_pressure, args=(temperature, sample, kwargs),
                               bracket=[1e-15, 1e5]).root
        except Exception:
            w.warn("Saturation pressure not found.", RuntimeWarning, stacklevel=2)
            satP = np.nan
        return satP

    def root_saturation_pressure(self, pressure, temperature, sample, kwargs):
        """ Function called by scipy.root_scalar when finding the saturation pressure using
        calculate_saturation_pressure.

        Parameters
        ----------
        pressure     float
            Pressure guess in bars
        temperature     float
            The temperature of the system in C.
        sample       Sample class
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
                                                   sample=sample, **kwargs))

    def NBO_O(self, sample, coeffs='webapp'):
        """
        Calculates NBO/O according to Appendix A.1. of Iacono-Marziano et al. (2012). NBO/O is
        calculated on either a hydrous or anhyrous basis, as set when initialising the Model class.

        Parameters
        ----------
        sample     pandas Series or dict
            Major element oxides in wt% (including H2O if using the hydrous parameterization).

        coeffs  str
            One of:
            - 'webapp' or 'manuscript' to include H2O in NBO/O
            - 'anhydrous' to exclude H2O from NBO/O

        Returns
        -------
        float
            NBO/O.
        """
        if isinstance(sample, sample_class.Sample) is False:
            raise core.InputError("Sample must be an instance of the Sample class.")
        if all(sample.check_oxide(ox) for ox in ['K2O', 'Na2O', 'CaO', 'MgO', 'FeO',
                                                 'Al2O3', 'SiO2', 'TiO2']) is False:
            raise core.InputError("sample must contain K2O, Na2O, CaO, MgO, FeO, Al2O3, SiO2, "
                                  "and TiO2.")

        X = sample.get_composition(units='mol_oxides', oxide_masses=self.IM_oxideMasses)

        if 'Fe2O3' in X:
            Fe2O3 = X['Fe2O3']
        else:
            Fe2O3 = 0

        NBO = 2*(X['K2O']+X['Na2O']+X['CaO']+X['MgO']+X['FeO']+2*Fe2O3-X['Al2O3'])
        Ox = (2*X['SiO2'] + 2*X['TiO2'] + 3*X['Al2O3'] + X['MgO'] + X['FeO'] + 2*Fe2O3 +
              X['CaO'] + X['Na2O'] + X['K2O'])

        if coeffs == 'webapp' or coeffs == 'manuscript':
            if 'H2O' not in X:
                raise core.InputError("sample must contain H2O if using the hydrous "
                                      "parameterization.")
            NBO = NBO + 2*X['H2O']
            Ox = Ox + X['H2O']

        return NBO/Ox


crmsg_BC_T = ("{param_name} ({param_val:.1f} {units}) is outside the broad range suggested by "
              "Iacono-Marziano ({calib_val0:.1f}-{calib_val1:.1f} {units}, although they note "
              "that this model is best calibrated at 1200-1300C). ")
# Warning for Iacono-Marziano pressure
crmsg_BC_P = ("{param_name} ({param_val:.1f} {units}) is outside the broad range suggested by "
              "Iacono-Marziano ({calib_val0:.1f}-{calib_val1:.1f} {units}, although most "
              "calibration experiments were conducted at <5000 bars).  ")


mixed = model_classes.MixedFluid({'H2O': water(), 'CO2': carbon()})
