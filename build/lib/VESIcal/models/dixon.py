from VESIcal import activity_models
from VESIcal import calibration_checks
from VESIcal import core
from VESIcal import fugacity_models
from VESIcal import model_classes
from VESIcal import sample_class

import numpy as np
from scipy.optimize import root_scalar
import warnings as w


class carbon(model_classes.Model):
    """
    Implementation of the Dixon (1997) carbon solubility model, as a Model class.
    """

    def __init__(self):
        self.set_volatile_species(['CO2'])
        self.set_fugacity_model(fugacity_models.fugacity_MRK_co2())
        self.set_activity_model(activity_models.activity_idealsolution())
        self.set_calibration_ranges([])
        self.set_solubility_dependence(False)

    def calculate_dissolved_volatiles(self, pressure, sample, X_fluid=1.0, **kwargs):
        """Calculates the dissolved CO2 concentration using Eqn (3) of Dixon (1997).

        Parameters
        ----------
        pressure  float
            Total pressure in bars.
        sample      Sample class
            Magma major element composition.
        X_fluid      float
            The mol fraction of CO2 in the fluid.

        Returns
        -------
        float
            The CO2 concentration in wt%.
        """

        if X_fluid < 0 or X_fluid > 1:
            raise core.InputError("X_fluid must have a value between 0 and 1.")
        if pressure < 0:
            raise core.InputError("Pressure must be positive.")
        if sample.check_oxide('SiO2') is False:
            raise core.InputError("sample must contain SiO2.")

        if pressure == 0:
            return 0

        XCO3 = self.molfrac_molecular(pressure=pressure, sample=sample, X_fluid=X_fluid, **kwargs)
        return (4400 * XCO3) / (36.594 - 44*XCO3)  # Following Dixon 1997 setting Mr as constant

    def calculate_equilibrium_fluid_comp(self, pressure, sample, **kwargs):
        """ Returns 1.0 if a pure H2O fluid is saturated. Returns 0.0 if a pure H2O fluid is
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

    def calculate_saturation_pressure(self, sample, X_fluid=1.0, **kwargs):
        """
        Calculates the pressure at which a pure CO2 fluid is saturated, for the given sample
        composition and CO2 concentration. Calls the scipy.root_scalar routine, which makes
        repeated called to the calculate_dissolved_volatiles method.

        Parameters
        ----------
        sample         Sample class
            Magma major element composition (including CO2).
        X_fluid     float
            The mole fraction of CO2 in the fluid. Default is 1.0.

        Returns
        -------
        float
            Calculated saturation pressure in bars.
        """
        if isinstance(sample, sample_class.Sample) is False:
            raise core.InputError("Sample must be an instance of the Sample class.")
        if sample.check_oxide('CO2') is False:
            raise core.InputError("sample must contain CO2.")
        if sample.get_composition('CO2') < 0:
            raise core.InputError("Dissolved CO2 concentration must be greater than 0 wt%.")

        try:
            satP = root_scalar(
                self.root_saturation_pressure, x0=100.0, x1=1000.0, args=(sample, kwargs)).root
        except Exception:
            w.warn("Saturation pressure not found.", RuntimeWarning, stacklevel=2)
            satP = np.nan
        return np.real(satP)

    def molfrac_molecular(self, pressure, sample, X_fluid=1.0, **kwargs):
        """Calculates the mole fraction of CO3(-2) dissolved when in equilibrium with a pure CO2
        fluid at 1200C, using Eqn (1) of Dixon (1997).

        Parameters
        ----------
        pressure      float
            Total pressure in bars.
        sample        Sample class
            Magma major element composition.
        X_fluid     float
            Mole fraction of CO2 in the fluid.

        Returns
        -------
        float
            Mole fraction of CO3(2-) dissolved."""

        DeltaVr = 23  # Changed to match dixon spreadsheet.14 (cm3 mole-1)
        P0 = 1
        R = 83.15
        T0 = 1473.15

        fugacity = self.fugacity_model.fugacity(pressure=pressure, X_fluid=X_fluid, **kwargs)

        XCO3Std = self.XCO3_Std(sample)

        return XCO3Std * fugacity * np.exp(-DeltaVr * (pressure-P0)/(R*T0))

    def XCO3_Std(self, sample):
        """ Calculates the mole fraction of CO3(2-) dissolved when in equilibrium with pure CO2
        vapour at 1200C and 1 bar, using Eq (8) of Dixon (1997).

        Parameters
        ----------
        sample    Sample class
            Magma major element chemistry.

        Returns
        -------
        float
            Mole fraction of CO3(2-) dissolved at 1 bar and 1200C.
        """
        if sample.get_composition('SiO2') > 48.9:
            return 3.817e-7
        else:
            return 8.697e-6 - 1.697e-7*sample.get_composition('SiO2')

    def root_saturation_pressure(self, pressure, sample, kwargs):
        """ The function called by scipy.root_scalar when finding the saturation pressure using
        calculate_saturation_pressure.

        Parameters
        ----------
        pressure     float
            Total pressure in bars.
        sample         Sample class
            Magma major element composition.

        Returns
        -------
        float
            The difference between the dissolved CO2 the pressure guessed, and the CO2
            concentration passed in the sample variable.
        """
        return (self.calculate_dissolved_volatiles(pressure=pressure, sample=sample, **kwargs) -
                sample.get_composition('CO2'))


class water(model_classes.Model):
    """
    Implementation of the Dixon (1997) water solubility model, as a Model class.
    """

    def __init__(self):
        self.set_volatile_species(['H2O'])
        self.set_fugacity_model(fugacity_models.fugacity_MRK_h2o())
        self.set_activity_model(activity_models.activity_idealsolution())
        self.set_calibration_ranges([
            calibration_checks.CalibrationRange(
                'pressure', 1000, calibration_checks.crf_LessThan, 'bar',
                'Dixon (1997, Pi-SiO2 simpl.) Water',
                fail_msg=crmsg_1000bar_fail,
                pass_msg=calibration_checks.crmsg_LessThan_pass,
                description_msg=calibration_checks.crmsg_LessThan_description),
            calibration_checks.CalibrationRange(
                'pressure', 2000, calibration_checks.crf_LessThan, 'bar',
                'Dixon (1997, Pi-SiO2 simpl.) Water',
                fail_msg=crmsg_2000bar_fail,
                pass_msg=calibration_checks.crmsg_LessThan_pass,
                description_msg=calibration_checks.crmsg_LessThan_description),
            calibration_checks.CalibrationRange(
                'SiO2', 49, calibration_checks.crf_LessThan, 'wt%',
                'Dixon (1997, Pi-SiO2 simpl.) Water',
                fail_msg=crmsg_49_fail,
                pass_msg=calibration_checks.crmsg_LessThan_pass,
                description_msg=calibration_checks.crmsg_LessThan_description),
            calibration_checks.CalibrationRange(
                'SiO2', 40, calibration_checks.crf_GreaterThan, 'wt%',
                'Dixon (1997, Pi-SiO2 simpl.) Water',
                fail_msg=crmsg_40_fail,
                pass_msg=calibration_checks.crmsg_GreaterThan_pass,
                description_msg=calibration_checks.crmsg_GreaterThan_description),
            calibration_checks.CalibrationRange(
                'temperature', [1000, 1400], calibration_checks.crf_Between, 'oC',
                'Dixon (1997, Pi-SiO2 simpl.) Water',
                fail_msg=crmsg_BC_T,
                pass_msg=calibration_checks.crmsg_Between_pass,
                description_msg=calibration_checks.crmsg_Between_description)])
        self.set_solubility_dependence(False)

    def calculate_dissolved_volatiles(self, pressure, sample, X_fluid=1.0, **kwargs):
        """Calculates the dissolved H2O concentration using Eqns (5) and (6) of Dixon (1997).

        Parameters
        ----------
        pressure  float
            Total pressure in bars.
        sample      Sample class
            Magma major element composition.
        X_fluid      float
            The mol fraction of H2O in the fluid.

        Returns
        -------
        float
            The H2O concentration in wt%.
        """
        if isinstance(sample, sample_class.Sample) is False:
            raise core.InputError("Sample must be an instance of the Sample class.")
        if sample.check_oxide('SiO2') is False:
            raise core.InputError("sample must contain SiO2.")
        if pressure < 0:
            raise core.InputError("Pressure must be positive")
        if X_fluid < 0 or X_fluid > 1:
            raise core.InputError("X_fluid must have a value between 0 and 1.")

        if pressure == 0:
            return 0

        XH2O = self.molfrac_molecular(pressure=pressure, sample=sample, X_fluid=X_fluid, **kwargs)
        XOH = self.XOH(pressure=pressure, sample=sample, X_fluid=X_fluid, **kwargs)

        XB = XH2O + 0.5*XOH
        return 1801.5*XB/(36.594-18.579*XB)  # Following Dixon spreadsheet

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

    def calculate_saturation_pressure(self, sample, X_fluid=1.0, **kwargs):
        """
        Calculates the pressure at which a pure H2O fluid is saturated, for the given sample
        composition and H2O concentration. Calls the scipy.root_scalar routine, which makes
        repeated called to the calculate_dissolved_volatiles method.

        Parameters
        ----------
        sample      Sample class
            Magma major element composition (including H2O).
        X_fluid     float
            The mole fraction of H2O in the fluid. Default is 1.0.

        Returns
        -------
        float
            Calculated saturation pressure in bars.
        """
        if sample.check_oxide('H2O') is False:
            raise core.InputError("sample must contain H2O")
        if sample.get_composition('H2O') < 0:
            raise core.InputError("H2O concentration must be greater than 0 wt%.")
        try:
            satP = root_scalar(
                self.root_saturation_pressure, x0=100.0, x1=1000.0, args=(sample, kwargs)).root
        except Exception:
            w.warn("Saturation pressure not found.", RuntimeWarning, stacklevel=2)
            satP = np.nan
        return np.real(satP)

    def molfrac_molecular(self, pressure, sample, X_fluid=1.0, **kwargs):
        """Calculates the mole fraction of molecular H2O dissolved when in equilibrium with
        a pure H2O fluid at 1200C, using Eqn (2) of Dixon (1997).

        Parameters
        ----------
        pressure      float
            Total pressure in bars.
        sample      Sample class
            Magma major element composition.
        X_fluid     float
            Mole fraction of H2O in the fluid.

        Returns
        -------
        float
            Mole fraction of molecular H2O dissolved.
        """

        VH2O = 12  # cm3 mole-1
        P0 = 1
        R = 83.15
        T0 = 1473.15

        XH2OStd = self.XH2O_Std(sample)

        fugacity = self.fugacity_model.fugacity(pressure=pressure, X_fluid=X_fluid, **kwargs)

        return XH2OStd * fugacity * np.exp(-VH2O * (pressure-P0)/(R*T0))

    def XH2O_Std(self, sample):
        """ Calculates the mole fraction of molecular H2O dissolved when in equilibrium with pure
        H2O vapour at 1200C and 1 bar, using Eq (9) of Dixon (1997).

        Parameters
        ----------
        sample    Sample class
            Magma major element composition.

        Returns
        -------
        float
            Mole fraction of molecular water dissolved at 1 bar and 1200C.
        """
        if sample.get_composition('SiO2') > 48.9:
            return 3.28e-5
        else:
            return -3.04e-5 + 1.29e-6*sample.get_composition('SiO2')

    def XOH(self, pressure, sample, X_fluid=1.0, **kwargs):
        """
        Calculates the mole fraction of hydroxyl groups dissolved by solving Eq (4) of Dixon
        (1997). Calls scipy.root_scalar to find the root of the XOH_root method.

        Parameters
        ----------
        pressure     float
            Total pressure in bars.
        sample         pandas Series or dict
            Major element oxides in wt%.
        X_fluid     float
            Mole fraction of H2O in the fluid.

        Returns
        -------
        float
            Mole fraction of hydroxyl groups dissolved.
        """

        XH2O = self.molfrac_molecular(pressure=pressure, sample=sample, X_fluid=X_fluid, **kwargs)
        if XH2O < 1e-14:
            return 0
        return np.exp(root_scalar(self.XOH_root, x0=np.log(0.5), x1=np.log(0.1), args=(XH2O)).root)

    def XOH_root(self, XOH, XH2O):
        """
        Method called by scipy.root_scalar when finding the saturation pressure using the
        calculate_saturation_pressure method. Implements Eq (4) of Dixon (1997).

        Parameters
        ----------
        XOH         float
            Guess for the mole fraction of hydroxyl groups dissolved in melt.
        XH2O    float
            Mole fraction of molecular water dissolved in melt.

        Returns
        -------
        float
            The difference between the RHS and LHS of Eq (4) of Dixon (1997) for the
            guessed value of XOH.
        """

        A = 0.403
        B = 15.333
        C = 10.894

        XOH = np.exp(XOH)

        term = (XOH)**2.0/(XH2O*(1.0-XOH-XH2O))

        lhs = - np.log(term)

        rhs = A + B*XOH + C*XH2O

        return rhs - lhs

    def root_saturation_pressure(self, pressure, sample, kwargs):
        """ Function called by scipy.root_scalar when finding the saturation pressure using
        calculate_saturation_pressure.

        Parameters
        ----------
        pressure     float
            Pressure guess in bars
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
        return (self.calculate_dissolved_volatiles(pressure=pressure, sample=sample, **kwargs) -
                sample.get_composition('H2O'))


crmsg_40_fail = ("{param_name} ({param_val:.1f} {units})<40 wt%, which is the lower calibration "
                 "limit of the Dixon (1997, Pi-SiO2 simpl.) Model. VESIcal has performed "
                 "calculations assuming SiO2=40wt%. ")
crmsg_1000bar_fail = ("{param_name} exceeds 1000 bar, which Iacono-Marziano et al. (2012) suggest"
                      " as an upper calibration limit of the Dixon (1997, Pi-SiO2 simpl.) Model, ")
crmsg_2000bar_fail = ("as well as the upper calibration limit of 2000 bar suggested by Lesne et"
                      " al. (2011), ")
crmsg_5000bar_fail = ("and the upper calibration limit of 5000 bar suggested by Newman and"
                      " Lowenstern, (2002). ")
crmsg_49_fail = ("{param_name} ({param_val:.1f} {units})>49 wt%, which is the calibration limit "
                 "of the Dixon (1997, Pi-SiO2 simpl.) Model. VESIcal has performed calculations "
                 "assuming SiO2=49wt% for this sample. ")
# Warning for Dixon temperature
crmsg_BC_T = ("{param_name} ({param_val:.1f} {units}) lies more than 200°C away from the "
              "temperature the Dixon (1997, Pi-SiO2 simpl.) model was calibrated for (the range "
              "suggested by Newman and Lowenstern, 2002; VolatileCalc). ")


mixed = model_classes.MixedFluid({'H2O': water(), 'CO2': carbon()})
