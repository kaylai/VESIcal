from VESIcal import calibration_checks
from VESIcal import core
from VESIcal import model_classes
from VESIcal import sample_class
from VESIcal import vplot
from VESIcal import batchfile  # needed for status_bar functions

from copy import deepcopy
import numpy as np
import pandas as pd
import warnings as w
import sys
from contextlib import redirect_stdout
import io

from thermoengine import equilibrate

# Variable to send H2O driver warnings into the void
_f = io.StringIO()

w.filterwarnings("ignore", message="rubicon.objc.ctypes_patch has only been tested ")

# -------------- MELTS preamble --------------- #
# instantiate thermoengine equilibrate MELTS instance
melts = equilibrate.MELTSmodel("1.2.0")

# Suppress phases not required in the melts simulation
phases = melts.get_phase_names()
for phase in phases:
    melts.set_phase_inclusion_status({phase: False})
melts.set_phase_inclusion_status({"Fluid": True, "Liquid": True})
# --------------------------------------------- #


class MagmaSat(model_classes.Model):
    """
    An object to instantiate a thermoengine equilibrate class
    """

    def __init__(self):
        self.melts_version = (
            "1.2.0"  # just here so users can see which version is being used
        )

        self.set_volatile_species(["H2O", "CO2"])
        self.set_calibration_ranges(
            [
                calibration_checks.CalibrationRange(
                    "pressure",
                    [0.0, 20000.0],
                    calibration_checks.crf_Between,
                    "bar",
                    "MagmaSat",
                    fail_msg=calibration_checks.crmsg_Between_fail,
                    pass_msg=calibration_checks.crmsg_Between_pass,
                    description_msg=calibration_checks.crmsg_Between_description,
                ),
                calibration_checks.CalibrationRange(
                    "temperature",
                    [800, 1400],
                    calibration_checks.crf_Between,
                    "oC",
                    "MagmaSat",
                    fail_msg=calibration_checks.crmsg_Between_fail,
                    pass_msg=calibration_checks.crmsg_Between_pass,
                    description_msg=calibration_checks.crmsg_Between_description,
                ),
            ]
        )
        self.model_type = "MagmaSat"

    def preprocess_sample(self, sample):
        """
        Returns sample with 0.0 values for any oxides not passed.

        Parameters
        ----------
        sample:     Sample class
            Magma major element composition.

        Returns
        -------
        Sample class object

        """
        _sample = deepcopy(sample)
        _sample = _sample.get_composition(
            units="wtpt_oxides",
            normalization=_sample.default_normalization,
            asSampleClass=True,
        )

        # set any necessary magmasat oxides to 0 if no value given
        for oxide in core.magmasat_oxides:
            if oxide in _sample.get_composition():
                pass
            else:
                _sample.change_composition({oxide: 0.0})

        # remove any non-magmasat oxides
        for oxide in _sample.get_composition().index:
            if oxide in core.magmasat_oxides:
                pass
            else:
                _sample.delete_oxide(oxide)

        self.bulk_comp_orig = _sample.get_composition(units="wtpt_oxides")

        return _sample

    def check_calibration_range(self, parameters, **kwargs):
        """Checks whether supplied parameters and calculated results are within the calibration
        range of the model, defined by the CalibrationRange objects. An empty string will be
        returned if all parameters are within the calibration range. If a parameter is not within
        the calibration range, a description of the problem will be returned in the string.

        Parameters
        ----------
        parameters     dict
            Dictionary keys are the names of the parameters to be checked, e.g., pressure
            temperature, SiO2, etc. Values are the values of each parameter. A complete set
            need not be given.

        Returns
        -------
        str
            String description of any parameters falling outside of the calibration range.
        """
        s = ""
        for cr in self.calibration_ranges:
            if cr.check(parameters) is False:
                s += cr.string(parameters, report_nonexistance=False)
            if "notsaturated" in kwargs:
                s += "Sample not saturated at these conditions."
        return s

    def get_calibration_range(self):
        """Returns a string describing the calibration ranges defined by the CalibrationRange
        objects for the model.

        Returns
        -------
        str
            String description of the calibration range objects."""
        s = ""
        for cr in self.calibration_ranges:
            s += cr.string(None)
        return s

    def get_fluid_mass(self, sample, temperature, pressure, H2O, CO2):
        """An internally used function to calculate fluid mass.

        Parameters
        ----------
        sample:     Sample class
            Magma major element composition.

        temperature: float
            Temperature in degrees C.

        pressure: float
            Pressure in bars

        H2O: float
            wt% H2O in the system

        CO2: float
            wt% CO2 in the system

        Returns
        -------
        float
            mass of the fluid in grams
        """
        pressureMPa = pressure / 10.0

        bulk_comp_dict = sample.get_composition(units="wtpt_oxides")
        bulk_comp_dict = {oxide: bulk_comp_dict[oxide] for oxide in
                          core.magmasat_oxides}
        bulk_comp_dict["H2O"] = H2O
        bulk_comp_dict["CO2"] = CO2
        melts.set_bulk_composition(bulk_comp_dict)

        with redirect_stdout(_f):
            output = melts.equilibrate_tp(temperature, pressureMPa, initialize=True)
        (status, temperature, pressureMPa, xmlout) = output[0]
        fluid_mass = melts.get_mass_of_phase(xmlout, phase_name="Fluid")

        return fluid_mass

    def get_XH2O_fluid(self, sample, temperature, pressure, H2O, CO2):
        """An internally used function to calculate fluid composition.

        Parameters
        ----------
        sample:     Sample class
            Magma major element composition.

        temperature: float
            Temperature in degrees C.

        pressure: float
            Pressure in bars

        H2O: float
            wt% H2O in the system

        CO2: float
            wt% CO2 in the system

        Returns
        -------
        float
            Mole fraction of H2O in the H2O-CO2 fluid

        """
        _sample = self.preprocess_sample(sample)

        pressureMPa = pressure / 10.0

        _sample.change_composition({"H2O": H2O, "CO2": CO2})
        melts.set_bulk_composition(
            _sample.get_composition(units="wtpt_oxides", normalization="none")
        )

        with redirect_stdout(_f):
            output = melts.equilibrate_tp(temperature, pressureMPa, initialize=True)
        (status, temperature, pressureMPa, xmlout) = output[0]
        fluid_comp = melts.get_composition_of_phase(
            xmlout, phase_name="Fluid", mode="component"
        )
        # NOTE mode='component' returns endmember component keys with values in mol fraction.

        if "Water" in fluid_comp:
            H2O_fl = fluid_comp["Water"]
        else:
            H2O_fl = 0.0

        return H2O_fl

    def calculate_dissolved_volatiles(
        self,
        sample,
        temperature,
        pressure,
        X_fluid=1,
        H2O_guess=0.0,
        verbose=False,
        **kwargs
    ):
        """
        Calculates the amount of H2O and CO2 dissolved in a magma at saturation at the given P/T
        conditions and fluid composition. Fluid composition will be matched to within 0.0001 mole
        fraction.

        Parameters
        ----------
        sample:     Sample class
            Magma major element composition.

        temperature: float or int
            Temperature, in degrees C.

        presure: float or int
            Pressure, in bars.

        X_fluid: float or int
            The default value is 1. The mole fraction of H2O in the H2O-CO2 fluid. X_fluid=1 is a
            pure H2O fluid. X_fluid=0 is a pure CO2 fluid.

        verbose: bool
            OPTIONAL: Default is False. If set to True, returns H2O and CO2 concentration in the
            melt, H2O and CO2 concentration in the fluid, mass of the fluid in grams, and
            proportion of fluid in the system in wt%.

        Returns
        -------
        dict
            A dictionary of dissolved volatile concentrations in wt% with keys H2O and CO2.
        """
        _sample = self.preprocess_sample(sample)

        if isinstance(X_fluid, int) or isinstance(X_fluid, float):
            pass
        else:
            raise core.InputError("X_fluid must be type int or float")

        if isinstance(H2O_guess, int) or isinstance(H2O_guess, float):
            pass
        else:
            raise core.InputError("H2O_guess must be type int or float")

        pressureMPa = pressure / 10.0

        if X_fluid != 0 and X_fluid != 1:
            if X_fluid < 0.001 or X_fluid > 0.999:
                raise core.InputError(
                    "X_fluid is calculated to a precision of 0.0001 mole "
                    "fraction. Value for X_fluid must be between 0.0001 and "
                    "0.9999."
                )

        H2O_val = H2O_guess
        CO2_val = 0.0
        fluid_mass = 0.0
        while fluid_mass <= 0:
            if X_fluid == 0:
                CO2_val += 0.1
            elif X_fluid >= 0.5:
                H2O_val += 0.2
                # NOTE this is setting XH2Owt of the system (not of the fluid) to X_fluid
                CO2_val = (H2O_val / X_fluid) - H2O_val
                # TODO this is what needs to be higher for higher XH2O. Slows down computation
                # by a second or two
            else:
                H2O_val += 0.1
                # NOTE this is setting XH2Owt of the system (not of the fluid) to X_fluid
                CO2_val = (H2O_val / X_fluid) - H2O_val
                # TODO this is what needs to be higher for higher XH2O. Slows down computation
                # by a second or two

            fluid_mass = self.get_fluid_mass(
                _sample, temperature, pressure, H2O_val, CO2_val
            )

        _sample.change_composition(
            {"H2O": H2O_val, "CO2": CO2_val}, units="wtpt_oxides"
        )
        melts.set_bulk_composition(
            _sample.get_composition(units="wtpt_oxides", normalization="none")
        )

        with redirect_stdout(_f):
            output = melts.equilibrate_tp(temperature, pressureMPa, initialize=True)
        (status, temperature, pressureMPa, xmlout) = output[0]
        liquid_comp = melts.get_composition_of_phase(
            xmlout, phase_name="Liquid", mode="oxide_wt"
        )
        fluid_comp = melts.get_composition_of_phase(
            xmlout, phase_name="Fluid", mode="component"
        )

        if "Water" in fluid_comp:
            H2O_fl = fluid_comp["Water"]
        else:
            H2O_fl = 0.0

        XH2O_fluid = H2O_fl

        # ------ Coarse Check ------ #
        while XH2O_fluid < X_fluid - 0.1:  # too low coarse check
            H2O_val += 0.2
            XH2O_fluid = self.get_XH2O_fluid(
                sample, temperature, pressure, H2O_val, CO2_val
            )

        while XH2O_fluid > X_fluid + 0.1:  # too high coarse check
            CO2_val += 0.1
            XH2O_fluid = self.get_XH2O_fluid(
                sample, temperature, pressure, H2O_val, CO2_val
            )

        # ------ Refinement 1 ------ #
        while XH2O_fluid < X_fluid - 0.01:  # too low refinement 1
            H2O_val += 0.05
            XH2O_fluid = self.get_XH2O_fluid(
                sample, temperature, pressure, H2O_val, CO2_val
            )

        while XH2O_fluid > X_fluid + 0.01:  # too high refinement 1
            CO2_val += 0.01
            XH2O_fluid = self.get_XH2O_fluid(
                sample, temperature, pressure, H2O_val, CO2_val
            )

        # ------ Refinement 2 ------ #
        while XH2O_fluid < X_fluid - 0.001:  # too low refinement 2
            H2O_val += 0.005
            XH2O_fluid = self.get_XH2O_fluid(
                sample, temperature, pressure, H2O_val, CO2_val
            )

        while XH2O_fluid > X_fluid + 0.001:  # too high refinement 2
            CO2_val += 0.001
            XH2O_fluid = self.get_XH2O_fluid(
                sample, temperature, pressure, H2O_val, CO2_val
            )

        # ------ Final refinement ------ #
        while XH2O_fluid < X_fluid - 0.0001:  # too low final refinement
            H2O_val += 0.001
            XH2O_fluid = self.get_XH2O_fluid(
                sample, temperature, pressure, H2O_val, CO2_val
            )

        while XH2O_fluid > X_fluid + 0.0001:  # too high final refinement
            CO2_val += 0.0001
            XH2O_fluid = self.get_XH2O_fluid(
                sample, temperature, pressure, H2O_val, CO2_val
            )

        # ------ Get calculated values ------ #
        _sample.change_composition(
            {"H2O": H2O_val, "CO2": CO2_val}, units="wtpt_oxides"
        )
        melts.set_bulk_composition(
            _sample.get_composition(units="wtpt_oxides", normalization="none")
        )

        with redirect_stdout(_f):
            output = melts.equilibrate_tp(temperature, pressureMPa, initialize=True)
        (status, temperature, pressureMPa, xmlout) = output[0]
        fluid_mass = melts.get_mass_of_phase(xmlout, phase_name="Fluid")
        system_mass = melts.get_mass_of_phase(xmlout, phase_name="System")
        liquid_comp = melts.get_composition_of_phase(
            xmlout, phase_name="Liquid", mode="oxide_wt"
        )
        fluid_comp = melts.get_composition_of_phase(
            xmlout, phase_name="Fluid", mode="component"
        )

        if "H2O" in liquid_comp:
            H2O_liq = liquid_comp["H2O"]
        else:
            H2O_liq = 0

        if "CO2" in liquid_comp:
            CO2_liq = liquid_comp["CO2"]
        else:
            CO2_liq = 0

        if "Water" in fluid_comp:
            H2O_fl = fluid_comp["Water"]
        else:
            H2O_fl = 0.0
        if "Carbon Dioxide" in fluid_comp:
            CO2_fl = fluid_comp["Carbon Dioxide"]
        else:
            CO2_fl = 0.0

        XH2O_fluid = H2O_fl

        if verbose:
            return {
                "temperature": temperature,
                "pressure": pressure,
                "H2O_liq": H2O_liq,
                "CO2_liq": CO2_liq,
                "XH2O_fl": H2O_fl,
                "XCO2_fl": CO2_fl,
                "FluidProportion_wt": 100 * fluid_mass / system_mass,
            }

        if verbose is False:
            return {"CO2_liq": CO2_liq, "H2O_liq": H2O_liq}

    def calculate_equilibrium_fluid_comp(
        self, sample, temperature, pressure, verbose=False, **kwargs
    ):
        """
        Returns H2O and CO2 concentrations in wt% in a fluid in equilibrium with the given sample
        at the given P/T condition.

        Parameters
        ----------
        sample:     Sample class
            Magma major element composition.

        temperature: float or int
            Temperature, in degrees C.

        presure: float or int
            Pressure, in bars.

        verbose: bool
            OPTIONAL: Default is False. If set to True, returns H2O and CO2 concentration in the
            fluid, mass of the fluid in grams, and proportion of fluid in the system in wt%.

        Returns
        -------
        dict
            A dictionary of fluid composition in wt% with keys 'H2O' and 'CO2' is returned.
        """
        _sample = self.preprocess_sample(sample)
        bulk_comp_dict = _sample.get_composition(
            units="wtpt_oxides", normalization="none"
        )

        if isinstance(temperature, float) or isinstance(temperature, int):
            pass
        else:
            raise core.InputError("temp must be type float or int")

        if isinstance(pressure, float) or isinstance(pressure, int):
            pass
        else:
            raise core.InputError("presure must be type float or int")

        # Check if only single volatile species is passed. If so, can skip calculations.
        if bulk_comp_dict["H2O"] == 0:
            if bulk_comp_dict["CO2"] == 0:
                if verbose is False:
                    return {"CO2": 0.0, "H2O": 0.0}
                if verbose:
                    return {
                        "CO2": 0.0,
                        "H2O": 0.0,
                        "FluidMass_grams": 0.0,
                        "FluidProportion_wt": 0.0,
                    }
            else:
                if verbose is False:
                    return {"CO2": 1.0, "H2O": 0.0}
        else:
            if bulk_comp_dict["CO2"] == 0:
                if verbose is False:
                    return {"CO2": 0.0, "H2O": 1.0}

        pressureMPa = pressure / 10.0

        melts.set_bulk_composition(bulk_comp_dict)

        with redirect_stdout(_f):
            output = melts.equilibrate_tp(temperature, pressureMPa, initialize=True)
        (status, temperature, pressureMPa, xmlout) = output[0]
        fluid_mass = melts.get_mass_of_phase(xmlout, phase_name="Fluid")
        flsystem_wtper = (
            100
            * fluid_mass
            / (fluid_mass + melts.get_mass_of_phase(xmlout, phase_name="Liquid"))
        )

        if fluid_mass > 0.0:
            fluid_comp = melts.get_composition_of_phase(
                xmlout, phase_name="Fluid", mode="component"
            )
            fluid_comp_H2O = fluid_comp["Water"]
            fluid_comp_CO2 = fluid_comp["Carbon Dioxide"]
        else:
            fluid_comp_H2O = 0
            fluid_comp_CO2 = 0

        melts.set_bulk_composition(self.bulk_comp_orig)  # reset

        if verbose is False:
            return {"CO2": fluid_comp_CO2, "H2O": fluid_comp_H2O}

        if verbose:
            return {
                "CO2": fluid_comp_CO2,
                "H2O": fluid_comp_H2O,
                "FluidMass_grams": fluid_mass,
                "FluidProportion_wt": flsystem_wtper,
            }

    def calculate_saturation_pressure(
        self, sample, temperature, verbose=False, **kwargs
    ):
        """
        Calculates the saturation pressure of a sample composition.

        Parameters
        ----------
        sample:     Sample class
            Magma major element composition.

        temperature: flaot or int
            Temperature of the sample in degrees C.

        verbose: bool
            OPTIONAL: Default is False. If set to False, only the saturation pressure is returned.
            If set to True, the saturation pressure, mass of fluid in grams, proportion of fluid
            in wt%, and H2O and CO2 concentrations in the fluid in mole fraction are all returned
            in a dict.

        Returns
        -------
        float or dict
            If verbose is set to False: Saturation pressure in bars.
            If verbose is set to True: dict of all calculated values.
        """
        _sample = self.preprocess_sample(sample)

        melts.set_bulk_composition(
            _sample.get_composition(units="wtpt_oxides", normalization="fixedvolatiles")
        )
        # Coarse search
        # NOTE that pressure is in MPa for MagmaSat calculations but reported in bars.
        pressureMPa = 2000

        # Check if saturated at 2000 MPa (rare, for deep samples)
        with redirect_stdout(_f):
            output = melts.equilibrate_tp(temperature, pressureMPa, initialize=True)
        (status, temperature, pressureMPa, xmlout) = output[0]
        fluid_mass = melts.get_mass_of_phase(xmlout, phase_name="Fluid")

        if fluid_mass <= 0:  # if not sat'd at 2000 MPa
            while fluid_mass <= 0:
                pressureMPa -= 100
                if pressureMPa <= 0:
                    break

                with redirect_stdout(_f):
                    output = melts.equilibrate_tp(
                        temperature, pressureMPa, initialize=True
                    )
                (status, temperature, pressureMPa, xmlout) = output[0]
                fluid_mass = melts.get_mass_of_phase(xmlout, phase_name="Fluid")

            fluid_mass = 0
            pressureMPa += 100

        elif fluid_mass > 0:  # if sat'd at 2000 MPa, add pressure
            while fluid_mass > 0:
                pressureMPa += 100

                melts.set_bulk_composition(
                    _sample.get_composition(
                        units="wtpt_oxides", normalization="fixedvolatiles"
                    )
                )
                with redirect_stdout(_f):
                    output = melts.equilibrate_tp(
                        temperature, pressureMPa, initialize=True
                    )
                (status, temperature, pressureMPa, xmlout) = output[0]
                fluid_mass = melts.get_mass_of_phase(xmlout, phase_name="Fluid")

            fluid_mass = 1.0
            pressureMPa -= 100

        # Refined search 1
        melts.set_bulk_composition(
            _sample.get_composition(units="wtpt_oxides", normalization="fixedvolatiles")
        )

        if fluid_mass <= 0:  # proceed down pressure search
            while fluid_mass <= 0:
                pressureMPa -= 10
                if pressureMPa <= 0:
                    break

                with redirect_stdout(_f):
                    output = melts.equilibrate_tp(
                        temperature, pressureMPa, initialize=True
                    )
                (status, temperature, pressureMPa, xmlout) = output[0]
                fluid_mass = melts.get_mass_of_phase(xmlout, phase_name="Fluid")

            fluid_mass = 0
            pressureMPa += 10

        elif fluid_mass > 0:  # proceed upward pressure search
            while fluid_mass > 0:
                pressureMPa += 10

                melts.set_bulk_composition(
                    _sample.get_composition(
                        units="wtpt_oxides", normalization="fixedvolatiles"
                    )
                )
                with redirect_stdout(_f):
                    output = melts.equilibrate_tp(
                        temperature, pressureMPa, initialize=True
                    )
                (status, temperature, pressureMPa, xmlout) = output[0]
                fluid_mass = melts.get_mass_of_phase(xmlout, phase_name="Fluid")

            fluid_mass = 1.0
            pressureMPa -= 10

        # Refined search 2
        melts.set_bulk_composition(
            _sample.get_composition(units="wtpt_oxides", normalization="fixedvolatiles")
        )

        if fluid_mass <= 0:  # proceed down pressure search
            while fluid_mass <= 0:
                pressureMPa -= 1
                if pressureMPa <= 0:
                    break

                with redirect_stdout(_f):
                    output = melts.equilibrate_tp(
                        temperature, pressureMPa, initialize=True
                    )
                (status, temperature, pressureMPa, xmlout) = output[0]
                fluid_mass = melts.get_mass_of_phase(xmlout, phase_name="Fluid")

        elif fluid_mass > 0:  # proceed upward pressure search
            while fluid_mass > 0:
                pressureMPa += 1

                melts.set_bulk_composition(
                    _sample.get_composition(
                        units="wtpt_oxides", normalization="fixedvolatiles"
                    )
                )
                with redirect_stdout(_f):
                    output = melts.equilibrate_tp(
                        temperature, pressureMPa, initialize=True
                    )
                (status, temperature, pressureMPa, xmlout) = output[0]
                fluid_mass = melts.get_mass_of_phase(xmlout, phase_name="Fluid")

        if pressureMPa != np.nan:
            satP = pressureMPa * 10  # convert pressure to bars
            flmass = fluid_mass
            flsystem_wtper = (
                100
                * fluid_mass
                / (fluid_mass + melts.get_mass_of_phase(xmlout, phase_name="Liquid"))
            )
            flcomp = melts.get_composition_of_phase(
                xmlout, phase_name="Fluid", mode="component"
            )
            try:
                flH2O = flcomp["Water"]
            except Exception:
                flH2O = 0.0
            try:
                flCO2 = flcomp["Carbon Dioxide"]
            except Exception:
                flCO2 = 0.0
        else:
            flmass = np.nan
            flsystem_wtper = np.nan
            flH2O = np.nan
            flCO2 = np.nan
            warnmessage = "Calculation failed."

        melts.set_bulk_composition(
            self.bulk_comp_orig
        )  # this needs to be reset always!

        if verbose is False:
            try:
                w.warn(warnmessage)
            except Exception:
                pass
            return satP

        elif verbose:
            try:
                w.warn(warnmessage)
            except Exception:
                pass
            return {
                "SaturationP_bars": satP,
                "FluidMass_grams": flmass,
                "FluidProportion_wt": flsystem_wtper,
                "XH2O_fl": flH2O,
                "XCO2_fl": flCO2,
            }

    def calculate_isobars_and_isopleths(
        self,
        sample,
        temperature,
        pressure_list,
        isopleth_list=None,
        smooth_isobars=True,
        smooth_isopleths=True,
        print_status=True,
        **kwargs
    ):
        """
        Calculates isobars and isopleths at a constant temperature for a given sample. Isobars can
        be calculated for any number of pressures. Isobars are calculated using 5 XH2O values
        (0, 0.25, 0.5, 0.75, 1).

        Parameters
        ----------
        sample:     Sample class
            Magma major element composition.

        temperature: float
            Temperature in degrees C.

        pressure_list: list or float
            List of all pressure values at which to calculate isobars, in bars. If only one value
            is passed it can be as float instead of list.

        isopleth_list: list or float
            OPTIONAL: Default value is None in which case only isobars will be calculated.
            List of all fluid compositions in mole fraction H2O (XH2Ofluid) at which to calcualte
            isopleths. Values can range from 0-1. If only one value is passed it can be as float
            instead of list.

        smooth_isobars: bool
            OPTIONAL. Default is True. If set to True, polynomials will be fit to the computed
            isobar points.

        smooth_isopleths: bool
            OPTIONAL. Default is True. If set to True, polynomials will be fit to the computed
            isopleth points.

        print_status: bool
            OPTIONAL: Default is True. If set to True, progress of the calculations will be
            printed to the terminal.

        Returns
        -------
        pandas DataFrame objects
            Two pandas DataFrames are returned; the first has isobar data, and the second has
            isopleth data. Columns in the isobar dataframe are 'Pressure', 'H2Omelt', and
            'CO2melt', correpsonding to pressure in bars and dissolved H2O and CO2 in the liquid
            in wt%. Columns in the isopleth dataframe are 'Pressure', 'H2Ofl', and 'CO2fl',
            corresponding to pressure in bars and H2O and CO2 concentration in the H2O-CO2 fluid,
            in wt%.
        """

        if isinstance(pressure_list, list):
            P_vals = pressure_list
        elif isinstance(pressure_list, int) or isinstance(pressure_list, float):
            P_vals = [pressure_list]
        else:
            raise core.InputError("pressure_list must be a single float (1000.0), int (1000), or "
                                  "list of those [1000, 2000.0, 3000].")

        if isopleth_list is None:
            pass
        elif isinstance(isopleth_list, list):
            iso_vals = isopleth_list
        else:
            iso_vals = [isopleth_list]

        _sample = self.preprocess_sample(sample)

        required_iso_vals = [0, 0.25, 0.5, 0.75, 1]
        all_iso_vals = iso_vals + required_iso_vals
        all_iso_vals = list(dict.fromkeys(all_iso_vals))  # remove duplicates
        all_iso_vals.sort()  # sort from smallest to largest

        isobar_data = []
        isopleth_data = []
        for X in iso_vals:
            isopleth_data.append([X, 0.0, 0.0])
        # Calculate equilibrium phase assemblage for all P/T conditions, check if saturated in
        # fluid...
        for i in P_vals:
            guess = 0.0
            if print_status:
                print("Calculating isobar at " + str(i) + " bars")
            X_iter = 0
            for X in all_iso_vals:
                X_iter += 1
                if print_status:
                    if isopleth_list is not None and X in iso_vals:
                        sys.stdout.write(
                            "\r Calculating isopleth at XH2Ofluid = "
                            + str(X)
                            + "               "
                        )
                    if X not in iso_vals:
                        sys.stdout.write(
                            "\r Calculating isobar control point at XH2Ofluid = "
                            + str(X)
                            + "               "
                        )
                    if X_iter == len(all_iso_vals):
                        sys.stdout.write(
                            "\r done.                                               "
                            "                                                       "
                            "                     \n"
                        )
                saturated_vols = self.calculate_dissolved_volatiles(
                    sample=_sample,
                    temperature=temperature,
                    pressure=i,
                    H2O_guess=guess,
                    X_fluid=X,
                )

                if X in required_iso_vals:
                    isobar_data.append([i, saturated_vols["H2O_liq"], saturated_vols["CO2_liq"]])
                if X in iso_vals:
                    isopleth_data.append([X, saturated_vols["H2O_liq"], saturated_vols["CO2_liq"]])

                guess = saturated_vols["H2O_liq"]

        if print_status:
            print("Done!")

        isobars_df = pd.DataFrame(isobar_data, columns=["Pressure", "H2O_liq", "CO2_liq"])
        isopleths_df = pd.DataFrame(isopleth_data, columns=["XH2O_fl", "H2O_liq", "CO2_liq"])

        melts.set_bulk_composition(self.bulk_comp_orig)  # reset

        if smooth_isobars:
            isobars_smoothed = vplot.smooth_isobars_and_isopleths(isobars=isobars_df)
            res_isobars = isobars_smoothed.copy()
        else:
            res_isobars = isobars_df.copy()

        if smooth_isopleths:
            isopleths_smoothed = vplot.smooth_isobars_and_isopleths(isopleths=isopleths_df)
            res_isopleths = isopleths_smoothed.copy()
        else:
            res_isopleths = isopleths_df.copy()

        return res_isobars, res_isopleths

    def calculate_degassing_path(self, sample, temperature, pressure="saturation",
                                 fractionate_vapor=0.0, init_vapor=0.0, steps=50, **kwargs):
        """
        Calculates degassing path for one sample

        Parameters
        ----------
        sample:     Sample class
            Magma major element composition.

            Legacy info follows, might still be useful?
            If pulling from an uploaded file
            with data for many samples, first call get_sample_composition() to get the sample
            desired. Then pass the result into this function.

        temperature: float
            Temperature at which to calculate degassing paths, in degrees C.

        pressure: float
            OPTIONAL. The perssure at which to begin the degassing calculations. Default value is
            'saturation', which runs the calculation with the initial pressure at the saturation
            pressure. If a pressure greater than the saturation pressure is input, the calculation
            will start at saturation, since this is the first pressure at which any degassing will
            occur.

        fractionate_vapor: float
            OPTIONAL. Proportion of vapor removed at each pressure step.
            Default value is 0.0 (completely closed-system degassing). Specifies the type of
            calculation performed, either closed system (0.0) or open system (1.0) degassing. If
            any value between <1.0 is chosen, user can also specify the 'init_vapor' argument
            (see below). A value in between 0 and 1 will remove that proportion of vapor at each
            step. For example, for a value of 0.2, the calculation will remove 20% of the vapor
            and retain 80% of the vapor at each pressure step.

        init_vapor: float
            OPTIONAL. Default value is 0.0. Specifies the amount of vapor (in wt%) coexisting
            with the melt before degassing.

        steps: int
            OPTIONAL. Default value is 50. Specifies the number of steps in pressure space at
            which dissolved volatile concentrations are calculated.

        Returns
        -------
        pandas DataFrame object

        """
        sys.stdout.write("Finding saturation point... ")  # print start of calculation to terminal
        _sample = self.preprocess_sample(sample)

        # Normalize sample composition
        _normed_comp = _sample.get_composition(normalization="standard")
        _sample.change_composition(_normed_comp)
        _sample_dict = _sample.get_composition()

        # ------ RESET MELTS ------ #
        # MELTS needs to be reloaded here. If an unfeasible composition gets set inside of MELTS,
        # which can happen when running open-system degassing path calcs, the following calls to
        # MELTS will fail. This prevents that from happening.
        melts = equilibrate.MELTSmodel("1.2.0")

        # Suppress phases not required in the melts simulation
        phases = melts.get_phase_names()
        for phase in phases:
            melts.set_phase_inclusion_status({phase: False})
        melts.set_phase_inclusion_status({"Fluid": True, "Liquid": True})

        melts.set_bulk_composition(_sample_dict)
        # ------------------------- #

        # Get saturation pressure
        data = self.calculate_saturation_pressure(sample=_sample, temperature=temperature,
                                                  verbose=True)

        if pressure == "saturation" or pressure >= data["SaturationP_bars"]:
            SatP_MPa = data["SaturationP_bars"] / 10.0
        else:
            SatP_MPa = pressure / 10.0

        # convert number of steps to step size
        MPa_step = SatP_MPa / steps
        if MPa_step < 1:
            MPa_step = 1

        P_array = np.arange(1.0, SatP_MPa, MPa_step)
        P_array = -np.sort(-P_array)
        # add last few MPa steps
        P_array = np.append(P_array, 0.5)
        P_array = np.append(P_array, 0.1)
        fl_wtper = data["FluidProportion_wt"]

        while fl_wtper <= init_vapor:
            with redirect_stdout(_f):
                output = melts.equilibrate_tp(temperature, SatP_MPa, initialize=True)
            (status, temperature, p, xmlout) = output[0]
            fl_mass = melts.get_mass_of_phase(xmlout, phase_name="Fluid")
            liq_mass = melts.get_mass_of_phase(xmlout, phase_name="Liquid")
            fl_comp = melts.get_composition_of_phase(xmlout, phase_name="Fluid")
            fl_wtper = 100 * fl_mass / (fl_mass + liq_mass)
            try:
                _sample_dict["H2O"] += fl_comp["H2O"] * 0.0005
            except Exception:
                _sample_dict["H2O"] = _sample_dict["H2O"] * 1.1
            try:
                _sample_dict["CO2"] += fl_comp["CO2"] * 0.0005
            except Exception:
                _sample_dict["CO2"] = _sample_dict["CO2"] * 1.1
            _sample = sample_class.Sample(_sample_dict)
            _sample_dict = _sample.get_composition(normalization="standard", units="wtpt_oxides")
            melts.set_bulk_composition(_sample_dict)  # reset MELTS

        pressure = []
        H2Oliq = []
        CO2liq = []
        H2Ofl = []
        CO2fl = []
        fluid_wtper = []
        iterno = 0
        sys.stdout.write("\r")  # carriage return to remove previous printed text
        for i in P_array:
            # Handle status_bar
            iterno += 1
            percent = iterno / len(P_array)
            batchfile.status_bar.status_bar(percent, btext="Calculating degassing path...")

            fl_mass = 0.0
            melts.set_bulk_composition(_sample_dict)
            with redirect_stdout(_f):
                output = melts.equilibrate_tp(temperature, i, initialize=True)
            (status, temperature, p, xmlout) = output[0]
            liq_comp = melts.get_composition_of_phase(xmlout, phase_name="Liquid")
            fl_comp = melts.get_composition_of_phase(xmlout, phase_name="Fluid", mode="component")
            liq_mass = melts.get_mass_of_phase(xmlout, phase_name="Liquid")
            fl_mass = melts.get_mass_of_phase(xmlout, phase_name="Fluid")
            fl_wtper = 100 * fl_mass / (fl_mass + liq_mass)

            if fl_mass > 0:
                pressure.append(p * 10.0)
                try:
                    H2Oliq.append(liq_comp["H2O"])
                except Exception:
                    H2Oliq.append(0)
                try:
                    CO2liq.append(liq_comp["CO2"])
                except Exception:
                    CO2liq.append(0)
                try:
                    H2Ofl.append(fl_comp["Water"])
                except Exception:
                    H2Ofl.append(0)
                try:
                    CO2fl.append(fl_comp["Carbon Dioxide"])
                except Exception:
                    CO2fl.append(0)
                fluid_wtper.append(fl_wtper)

                try:
                    _sample_dict["H2O"] = (liq_comp["H2O"] +
                                           (_sample_dict["H2O"] - liq_comp["H2O"]) *
                                           (1.0 - fractionate_vapor))
                except Exception:
                    _sample_dict["H2O"] = 0
                try:
                    _sample_dict["CO2"] = (liq_comp["CO2"] +
                                           (_sample_dict["CO2"] - liq_comp["CO2"]) *
                                           (1.0 - fractionate_vapor))
                except Exception:
                    _sample_dict["CO2"] = 0
            _sample = sample_class.Sample(_sample_dict)
            _sample_dict = _sample.get_composition(normalization="standard", units="wtpt_oxides")

        melts.set_bulk_composition(self.bulk_comp_orig)  # this needs to be reset always!
        open_degassing_df = pd.DataFrame(list(zip(pressure, H2Oliq, CO2liq, H2Ofl, CO2fl,
                                                  fluid_wtper)),
                                         columns=["Pressure_bars",
                                                  "H2O_liq",
                                                  "CO2_liq",
                                                  "XH2O_fl",
                                                  "XCO2_fl",
                                                  "FluidProportion_wt"])

        open_degassing_df = open_degassing_df[open_degassing_df.CO2_liq >= 0.0]
        open_degassing_df = open_degassing_df[open_degassing_df.H2O_liq >= 0.0]

        return open_degassing_df
