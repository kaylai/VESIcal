from VESIcal import activity_models
from VESIcal import calibration_checks
from VESIcal import core
from VESIcal import fugacity_models
from VESIcal import model_classes
from VESIcal import sample_class

import numpy as np
import warnings as w
import sympy
from scipy.optimize import root_scalar


class water(model_classes.Model):
    """
    Implementation of the Liu et al. (2005) H2O solubility model for metaluminous high-silica
    rhyolitic melts.
    """

    def __init__(self):
        """
        Initialize the model.
        """
        self.set_volatile_species(['H2O'])
        self.set_fugacity_model(fugacity_models.fugacity_idealgas())
        self.set_activity_model(activity_models.activity_idealsolution())
        self.set_solubility_dependence(False)

        # Generate calibration range objects for each oxide
        cr_oxide_list = []
        for ox in watercomprange.keys():
            cr_oxide_list.append(
                calibration_checks.CalibrationRange(
                    ox, watercomprange[ox], calibration_checks.crf_Between, 'wt%',
                    'Liu et al. (2005) water',
                    fail_msg=calibration_checks.crmsg_BC_fail,
                    pass_msg=calibration_checks.crmsg_BC_pass,
                    description_msg=calibration_checks.crmsg_Between_description))

        self.set_calibration_ranges([
            calibration_checks.CalibrationRange(
                'pressure', [0, 5000.0], calibration_checks.crf_Between, 'bar',
                'Liu et al. (2005) water',
                fail_msg=calibration_checks.crmsg_Between_fail,
                pass_msg=calibration_checks.crmsg_Between_pass,
                description_msg=calibration_checks.crmsg_Between_description),
            calibration_checks.CalibrationRange(
                'temperature', [700.0, 1200], calibration_checks.crf_Between, 'oC',
                'Liu et al. (2005) water',
                fail_msg=calibration_checks.crmsg_Between_fail,
                pass_msg=calibration_checks.crmsg_Between_pass,
                description_msg=calibration_checks.crmsg_Between_description),
            calibration_checks.CalibrationRange(
                'sample', None, crf_WaterComp, None, None,
                fail_msg=crmsg_WaterComp_fail,
                pass_msg=crmsg_Comp_pass,
                description_msg=crmsg_Comp_description)] + cr_oxide_list)

    def calculate_dissolved_volatiles(self, sample, pressure, temperature, X_fluid=1.0, **kwargs):
        """
        Parameters
        ----------
        sample:     Sample class
            Magma major element composition.

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
        pressureMPa = pressure / 10.0
        Pw = pressureMPa * X_fluid
        PCO2 = pressureMPa * (1 - X_fluid)

        temperatureK = temperature + 273.15

        H2Ot = ((354.94*Pw**(0.5) + 9.623*Pw - 1.5223*Pw**(1.5)) / temperatureK +
                0.0012439*Pw**(1.5) + PCO2*(-1.084*10**(-4)*Pw**(0.5) - 1.362*10**(-5)*Pw))

        return H2Ot

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
        temperatureK = temperature + 273.15
        pressureMPa = pressure / 10.0

        H2Ot = sample.get_composition("H2O")

        # calculate saturation pressure and assert that input P <= SatP
        satP = self.calculate_saturation_pressure(temperature, sample)
        is_saturated = satP - pressure
        if is_saturated >= 0:
            pass
        else:
            w.warn("{:.1f} bars is above the saturation pressure ({:.1f} bars) for this sample. "
                   "Results from this calculation may be nonsensical.".format(pressure, satP))

        # Use sympy to solve solubility equation for XH2Ofluid
        XH2Ofluid = sympy.symbols('XH2Ofluid')  # XH2Ofluid is the variable to solve for

        equation = ((354.94*(XH2Ofluid*pressureMPa)**(0.5) + 9.623*(XH2Ofluid*pressureMPa)
                    - 1.5223*(XH2Ofluid*pressureMPa)**(1.5)) / temperatureK
                    + 0.0012439*(XH2Ofluid*pressureMPa)**(1.5)
                    + pressureMPa*(1-XH2Ofluid)*(-1.084*10**(-4)*(XH2Ofluid*pressureMPa)**(0.5)
                    - 1.362*10**(-5)*(XH2Ofluid*pressureMPa)) - H2Ot)

        XH2Ofluid = sympy.solve(equation, XH2Ofluid)[0]
        if XH2Ofluid > 1:
            XH2Ofluid = 1
        if XH2Ofluid < 0:
            XH2Ofluid = 0

        return XH2Ofluid

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
                               x0=1.0, x1=2.0).root
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
        sample:    Sample class
            Magma major element composition, including H2O.
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


class carbon(model_classes.Model):
    """
    Implementation of the Liu et al. (2005) H2O-CO2 solubility model for metaluminous high-silica
    rhyolitic melts.
    """

    def __init__(self):
        """
        Initialize the model.
        """
        self.set_volatile_species(['CO2'])
        self.set_fugacity_model(fugacity_models.fugacity_idealgas())
        self.set_activity_model(activity_models.activity_idealsolution())
        self.set_solubility_dependence(False)

        # Generate calibration range objects for each oxide
        cr_oxide_list = []
        for ox in carboncomprange.keys():
            cr_oxide_list.append(
                calibration_checks.CalibrationRange(
                    ox, carboncomprange[ox], calibration_checks.crf_Between, 'wt%',
                    'Liu et al. (2005) carbon',
                    fail_msg=calibration_checks.crmsg_BC_fail,
                    pass_msg=calibration_checks.crmsg_BC_pass,
                    description_msg=calibration_checks.crmsg_Between_description))

        self.set_calibration_ranges([
            calibration_checks.CalibrationRange(
                'pressure', [0, 5000.0], calibration_checks.crf_Between, 'bar',
                'Liu et al. (2005) Carbon',
                fail_msg=calibration_checks.crmsg_Between_fail,
                pass_msg=calibration_checks.crmsg_Between_pass,
                description_msg=calibration_checks.crmsg_Between_description),
            calibration_checks.CalibrationRange(
                'temperature', [700.0, 1200], calibration_checks.crf_Between, 'oC',
                'Liu et al. (2005) Carbon',
                fail_msg=calibration_checks.crmsg_Between_fail,
                pass_msg=calibration_checks.crmsg_Between_pass,
                description_msg=calibration_checks.crmsg_Between_description),
            calibration_checks.CalibrationRange(
                'sample', None, crf_CarbonComp, None, None,
                fail_msg=crmsg_CarbonComp_fail,
                pass_msg=crmsg_Comp_pass,
                description_msg=crmsg_Comp_description)] + cr_oxide_list)

    def calculate_dissolved_volatiles(self, sample, pressure, temperature, X_fluid=1, **kwargs):
        """
        Parameters
        ----------
        sample:        Sample class
            Magma major element composition.

        pressure float
            Pressure in bars.

        temperature float
            Temperature in degrees C.

        X_fluid float
            OPTIONAL. Default is 1. Mole fraction of CO2 in the H2O-CO2 fluid.

        Returns
        -------
        float
            Calculated dissolved CO2 concentration in wt%.
        """
        pressureMPa = pressure / 10.0
        Pw = pressureMPa * (1 - X_fluid)
        PCO2 = pressureMPa * X_fluid

        temperatureK = temperature + 273.15

        CO2melt_ppm = (PCO2*(5668 - 55.99*Pw)/temperatureK
                       + PCO2*(0.4133*Pw**(0.5) + 2.041*10**(-3)*Pw**(1.5)))

        CO2melt = CO2melt_ppm / 10000

        return CO2melt

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
            Calculated equilibrium fluid concentration in XCO2fluid mole fraction.
        """
        temperatureK = temperature + 273.15
        pressureMPa = pressure / 10.0

        if temperatureK <= 0.0:
            raise core.InputError("Temperature must be greater than 0K.")
        if isinstance(sample, sample_class.Sample) is False:
            raise core.InputError("Sample must be an instance of the Sample class.")
        if sample.check_oxide('CO2') is False:
            raise core.InputError("sample must contain CO2.")
        if sample.get_composition('CO2') < 0.0:
            raise core.InputError("Dissolved CO2 concentration must be greater than 0 wt%.")

        CO2melt_wt = sample.get_composition("CO2")
        CO2melt_ppm = CO2melt_wt * 10000

        # calculate saturation pressure and assert that input P <= SatP
        satP = self.calculate_saturation_pressure(temperature, sample)
        is_saturated = satP - pressure
        if is_saturated >= 0:
            pass
        else:
            w.warn(str(pressure) + " bars is above the saturation pressure (" + str(satP) +
                   " bars) for this sample. Results from this calculation may be nonsensical.")

        # Use sympy to solve solubility equation for XH2Ofluid
        XCO2fluid = sympy.symbols('XCO2fluid')  # XCO2fluid is the variable to solve for

        equation = ((XCO2fluid*pressureMPa*(5668 - 55.99*(pressureMPa*(1-XCO2fluid)))/temperatureK
                    + (XCO2fluid*pressureMPa)*(0.4133*(pressureMPa*(1-XCO2fluid))**(0.5)
                    + 2.041*10**(-3)*(pressureMPa*(1-XCO2fluid))**(1.5))) - CO2melt_ppm)

        XCO2fluid = sympy.solve(equation, XCO2fluid, real=True)[0]

        if type(XCO2fluid) != float:
            w.warn("Could not find equilibrium fluid composition.")
            return 0
        else:
            if XCO2fluid > 1:
                XCO2fluid = 1
            if XCO2fluid < 0:
                XCO2fluid = 0

            return XCO2fluid

    def calculate_saturation_pressure(self, temperature, sample, X_fluid=1.0, **kwargs):
        """
        Calculates the pressure at which a an CO2-bearing fluid is saturated. Calls the
        scipy.root_scalar routine, which makes repeated called to the
        calculate_dissolved_volatiles method.

        Parameters
        ----------
        sample:        Sample class
            Magma major element composition.

        temperature float
            Temperature in degrees C.

        X_fluid float
            OPTIONAL. Default is 0. Mole fraction of CO2 in the H2O-CO2 fluid.

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
        if sample.check_oxide('CO2') is False:
            raise core.InputError("sample must contain CO2.")
        if sample.get_composition('CO2') < 0.0:
            raise core.InputError("Dissolved CO2 concentration must be greater than 0 wt%.")

        try:
            satP = root_scalar(self.root_saturation_pressure,
                               args=(temperature, sample, X_fluid, kwargs),
                               x0=10.0, x1=2000.0).root
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
        sample:        Sample class
            Magma major element composition, including H2O.
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
                sample.get_composition('CO2'))


# Defining compositional ranges for Liu - based on the Max value of the calibration dataset +-5%
# of that value
watercomprange = {'SiO2':  [71, 82],
                  'TiO2':  [0, 0.21],
                  'FeO':   [0, 1.5],
                  'Al2O3': [11.5, 14.2],
                  'MgO':   [0, 0.18],
                  'CaO':   [0, 1.2],
                  'Na2O':  [3.2, 4.9],
                  'K2O':   [3.4, 6.0]
                  }

carboncomprange = {'SiO2':  [73, 82],
                   'TiO2':  [0.07, 0.12],
                   'FeO':   [0.36, 1.1],
                   'Al2O3': [11.9, 13.7],
                   'MgO':   [0.05, 0.08],
                   'CaO':   [0.23, 0.6],
                   'Na2O':  [3.9, 4.4],
                   'K2O':   [4.0, 5.0]
                   }


def crf_WaterComp(calibval=None, sample=sample_class.Sample({})):
    comp = sample.get_composition(units='wtpt_oxides')

    test_results = []
    for ox in watercomprange.keys():
        test_results.append(comp[ox] >= watercomprange[ox][0] and
                            comp[ox] <= watercomprange[ox][1])

    return all(test_results)


def crf_CarbonComp(calibval=None, sample=sample_class.Sample({})):
    comp = sample.get_composition(units='wtpt_oxides')

    test_results = []
    for ox in carboncomprange.keys():
        test_results.append(comp[ox] >= carboncomprange[ox][0] and
                            comp[ox] <= carboncomprange[ox][1])

    return all(test_results)


def crf_MixedComp(calibval=None, sample=sample_class.Sample({})):
    comp = sample.get_composition(units='wtpt_oxides')

    test_results = []
    for ox in carboncomprange.keys():
        test_results.append(comp[ox] >= np.min([watercomprange[ox][0], carboncomprange[ox][0]]))
        test_results.append(comp[ox] <= np.max([watercomprange[ox][1], carboncomprange[ox][1]]))

    return all(test_results)


crmsg_Comp_pass = ("The sample appears to be similar in composition to the rhyolites and "
                   "haplogranites used to calibrate the Liu et al. model.")
crmsg_WaterComp_fail = (" These calibration limits were selected based on the minimum and "
                        "maximum values of these oxides (+-5%) in the Water calibration dataset. "
                        "As the Liu et al. model incorperates no term for compositional "
                        "dependence, users must take extreme care when extrapolating this model "
                        "to compositions which differ significantly from the haplogranites and "
                        "rhyolites in the calibration dataset. These warnings are simply a guide;"
                        " we suggest that users carefully compare their major element data to the "
                        "calibration dataset to check for suitability ")
crmsg_CarbonComp_fail = (" These calibration limits were selected based on the minimum and "
                         "maximum values of these oxides (+-5%) in the Carbon calibration dataset."
                         " As the Liu et al. model incorperates no term for compositional "
                         "dependence, users must take extreme care when extrapolating this model "
                         "to compositions which differ significantly from the haplogranites and "
                         "rhyolites in the calibration dataset. These warnings are simply a guide;"
                         " we suggest that users carefully compare their major element data to the"
                         " calibration dataset to check for suitability ")
crmsg_Comp_fail = (" These calibration limits were selected based on the minimum and maximum "
                   "values of these oxides (+-5%) in the combined Water and Carbon calibration"
                   " dataset. As the Liu et al. model incorperates no term for compositional "
                   "dependence, users must take extreme care when extrapolating this model to "
                   "compositions which differ significantly from the haplogranites and rhyolites "
                   "in the calibration dataset. These warnings are simply a guide; we suggest that"
                   " users carefully compare their major element data to the calibration dataset "
                   "to check for suitability ")
crmsg_Comp_description = "The Liu et al. model is suitable for haplogranites and rhyolites."


mixed = model_classes.MixedFluid({'H2O': water(), 'CO2': carbon()})

# Prevent the same error being returned for both H2O and CO2 when using the mixed model.
mixed.models[1].set_calibration_ranges([])
_crs_to_update = mixed.models[0].calibration_ranges
for _cr in _crs_to_update:
    _cr.model_name = 'Liu et al. (2005)'

cr_oxide_list = []
for ox in carboncomprange.keys():
    lo = np.min([watercomprange[ox][0], carboncomprange[ox][0]])
    hi = np.max([watercomprange[ox][1], carboncomprange[ox][1]])
    cr_oxide_list.append(
        calibration_checks.CalibrationRange(
            ox, [lo, hi], calibration_checks.crf_Between, 'wt%', 'Liu et al. (2005)',
            fail_msg=calibration_checks.crmsg_BC_fail,
            pass_msg=calibration_checks.crmsg_BC_pass,
            description_msg=calibration_checks.crmsg_Between_description))

mixed.models[0].set_calibration_ranges([
    calibration_checks.CalibrationRange(
        'pressure', [0, 5000.0], calibration_checks.crf_Between, 'bar', 'Liu et al. (2005)',
        fail_msg=calibration_checks.crmsg_Between_fail,
        pass_msg=calibration_checks.crmsg_Between_pass,
        description_msg=calibration_checks.crmsg_Between_description),
    calibration_checks.CalibrationRange(
        'temperature', [700.0, 1200], calibration_checks.crf_Between, 'oC', 'Liu et al. (2005)',
        fail_msg=calibration_checks.crmsg_Between_fail,
        pass_msg=calibration_checks.crmsg_Between_pass,
        description_msg=calibration_checks.crmsg_Between_description),
    calibration_checks.CalibrationRange(
        'sample', None, crf_MixedComp, None, None,
        fail_msg=crmsg_Comp_fail,
        pass_msg=crmsg_Comp_pass,
        description_msg=crmsg_Comp_description)] + cr_oxide_list)
