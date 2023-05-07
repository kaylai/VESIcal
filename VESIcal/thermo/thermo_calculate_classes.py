from abc import abstractmethod
import warnings as w

from VESIcal.thermo import densityx, giordano


class ThermoCalculate(object):
    """
    The ThermoCalculate object is designed to emulate the VESIcal Calculate
    object. It provides a common workflow for all additional thermodynamic
    calculations one might wish to perform on a sample. It is separate from
    the VESIcal Calculate class as that is specifically for volatile
    solubility models
    """
    def __init__(self, sample, silence_warnings=False, **kwargs):
        """
        Initializes the calculation.

        Parameters
        ----------
        sample:    Sample class
            The rock composition as a Sample object.
        silence_warnings:     bool
            Silence warnings about calibration ranges. Default is False.
        """
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


class calculate_liquid_density(ThermoCalculate):
    """
    Calculates the density of the liquid using the DensityX model. Using
    this interface will preprocess the sample, run the calculation, and then
    check the calibration ranges. All parameters required by the chosen model
    must be passed.

    Parameters
    ----------
    sample:    Sample class
        The rock composition as a Sample object.
    pressure:   float
        Total pressure in bars.
    temperature:  float
        Temperature in degrees C.
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
        Density in g/L.
    """

    def calculate(self, sample, pressure, temperature, **kwargs):
        density_g_per_L = densityx.calculate_liquid_density(sample=sample,
                                                            pressure=pressure,
                                                            temperature=temperature,
                                                            **kwargs)
        return density_g_per_L

    def check_calibration_range(self, sample, pressure, temperature, **kwargs):
        pass


class calculate_liquid_viscosity(ThermoCalculate):
    """
    Calculates the density of the liquid using the Giordano et al. (2008)
    model. Using this interface will preprocess the sample, run the calculation,
    and then check the calibration ranges. All parameters required by the chosen
    model must be passed.

    Parameters
    ----------
    sample:    Sample class
        The rock composition as a Sample object.
    temperature:  float
        Temperature in degrees C.
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
        Log viscosity in Pa*s.
    """

    def calculate(self, sample, temperature, **kwargs):
        viscosity_log_nu_Pas = giordano.calculate_liquid_viscosity(sample=sample,
                                                                   temperature=temperature,
                                                                   **kwargs)
        return viscosity_log_nu_Pas

    def check_calibration_range(self, sample, temperature, **kwargs):
        pass
