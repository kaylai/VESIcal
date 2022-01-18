from VESIcal import core
from VESIcal import models
from VESIcal import calculate_classes
from VESIcal import batchfile

import numpy as np
import warnings as w
import sys

from thermoengine import equilibrate

w.filterwarnings("ignore", message="rubicon.objc.ctypes_patch has only been "
                                   "tested ")

# -------------- MELTS preamble --------------- #
# instantiate thermoengine equilibrate MELTS instance
melts = equilibrate.MELTSmodel('1.2.0')

# Suppress phases not required in the melts simulation
phases = melts.get_phase_names()
for phase in phases:
    melts.set_phase_inclusion_status({phase: False})
melts.set_phase_inclusion_status({'Fluid': True, 'Liquid': True})
# --------------------------------------------- #


# -------------- BATCH PROCESSING ----------- #
class BatchFile(batchfile.BatchFile):
    """Performs model functions on a batchfile.BatchFile object
    """
    pass

    def get_XH2O_fluid(self, sample, temperature, pressure, H2O, CO2):
        """An internally used function to calculate fluid composition.

        Parameters
        ----------
        sample: dictionary
            Sample composition in wt% oxides

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
        pressureMPa = pressure / 10.0

        bulk_comp = {oxide:  sample[oxide] for oxide in core.oxides}
        bulk_comp["H2O"] = H2O
        bulk_comp["CO2"] = CO2
        melts.set_bulk_composition(bulk_comp)

        output = melts.equilibrate_tp(temperature, pressureMPa,
                                      initialize=True)
        (status, temperature, pressureMPa, xmlout) = output[0]
        fluid_comp = melts.get_composition_of_phase(xmlout, phase_name='Fluid',
                                                    mode='component')
        # NOTE mode='component' returns endmember component keys with values
        # in mol fraction.

        if "Water" in fluid_comp:
            H2O_fl = fluid_comp["Water"]
        else:
            H2O_fl = 0.0

        return H2O_fl

    def calculate_dissolved_volatiles(self, temperature, pressure, X_fluid=1,
                                      print_status=True, model='MagmaSat',
                                      record_errors=False, **kwargs):
        """
        Calculates the amount of H2O and CO2 dissolved in a magma at the given
        P/T conditions and fluid composition. Fluid composition will be
        matched to within 0.0001 mole fraction.

        Parameters
        ----------
        temperature: float, int, or str
            Temperature, in degrees C. Can be passed as float, in which case
            the passed value is used as the temperature for all samples.
            Alternatively, temperature information for each individual
            sample may already be present in the BatchFile object. If so, pass
            the str value corresponding to the column title in the BatchFile
            object.

        pressure: float, int, or str
            Pressure, in bars. Can be passed as float or int, in which case the
            passed value is used as the pressure for all samples.
            Alternatively, pressure information for each individual sample may
            already be present in the BatchFile object. If so, pass the str
            value corresponding to the column title in the BatchFile object.

        X_fluid: float, int, or str
            OPTIONAL: Default value is 1. The mole fraction of H2O in the
            H2O-CO2 fluid. X_fluid=1 is a pure H2O fluid. X_fluid=0 is a pure
            CO2 fluid. Can be passed as a float or int, in which case the
            passed value is used as the X_fluid for all samples.
            Alternatively, X_fluid information for each individual sample may
            already be present in the BatchFile object. If so, pass the str
            value corresponding to the column title in the BatchFile object.

        print_status: bool
            OPTIONAL: The default value is True, in which case the progress of
            the calculation will be printed to the terminal. If set to False,
            nothing will be printed. MagmaSat calculations tend to be slow,
            and so a value of True is recommended for most use cases.

        model: string
            OPTIONAL: Default is 'MagmaSat'. Any other model name can be
            passed here.

        record_errors: bool
            OPTIONAL: If True, any errors arising during the calculation will
            be recorded as a column.

        Returns
        -------
        pandas DataFrame
            Original data passed plus newly calculated values are returned.
        """
        dissolved_data = self.get_data().copy()

        if isinstance(temperature, str):
            file_has_temp = True
            temp_name = temperature
        elif isinstance(temperature, float) or isinstance(temperature, int):
            file_has_temp = False
        else:
            raise core.InputError("temp must be type str or float or int")

        if isinstance(pressure, str):
            file_has_press = True
            press_name = pressure
        elif isinstance(pressure, float) or isinstance(pressure, int):
            file_has_press = False
        else:
            raise core.InputError("pressure must be type str or float or int")

        if isinstance(X_fluid, str):
            file_has_X = True
            X_name = X_fluid
        elif isinstance(X_fluid, float) or isinstance(X_fluid, int):
            file_has_X = False
            if X_fluid != 0 and X_fluid != 1:
                if X_fluid < 0.001 or X_fluid > 0.999:
                    raise core.InputError("X_fluid is calculated to a "
                                          "precision of 0.0001 mole fraction. "
                                          "Value for X_fluid must be between "
                                          "0.0001 and 0.9999.")
        else:
            raise core.InputError("X_fluid must be type str or float or int")

        # Check if the model passed as the attribute "model_type"
        # Currently only implemented for MagmaSat type models
        if hasattr(model, 'model_type') is True:
            model = model.model_type

        H2Ovals = []
        CO2vals = []
        warnings = []
        errors = []
        if model in models.get_model_names(model='mixed'):
            for index, row in dissolved_data.iterrows():
                try:
                    if file_has_temp:
                        temperature = row[temp_name]

                    if file_has_press:
                        pressure = row[press_name]

                    if file_has_X:
                        X_fluid = row[X_name]

                    # Get sample comp as Sample class with defaults
                    bulk_comp = self.get_sample_composition(
                                      index,
                                      normalization=self.default_normalization,
                                      units='wtpt_oxides', asSampleClass=True)
                    bulk_comp.set_default_units(self.default_units)
                    bulk_comp.set_default_normalization(
                                                    self.default_normalization)
                    calc = calculate_classes.calculate_dissolved_volatiles(
                                           sample=bulk_comp, pressure=pressure,
                                           temperature=temperature,
                                           X_fluid=(X_fluid, 1-X_fluid),
                                           model=model, silence_warnings=True,
                                           **kwargs)
                    H2Ovals.append(calc.result['H2O_liq'])
                    CO2vals.append(calc.result['CO2_liq'])
                    warnings.append(calc.calib_check)
                    errors.append('')
                except Exception:
                    H2Ovals.append(np.nan)
                    CO2vals.append(np.nan)
                    warnings.append('Calculation Failed.')
                    errors.append(sys.exc_info()[0])
            dissolved_data["H2O_liq_VESIcal"] = H2Ovals
            dissolved_data["CO2_liq_VESIcal"] = CO2vals

            if file_has_temp is False:
                dissolved_data["Temperature_C_VESIcal"] = temperature
            if file_has_press is False:
                dissolved_data["Pressure_bars_VESIcal"] = pressure
            if file_has_X is False:
                dissolved_data["X_fluid_input_VESIcal"] = X_fluid
            dissolved_data["Model"] = model
            dissolved_data["Warnings"] = warnings
            if record_errors:
                dissolved_data["Errors"] = errors

            return dissolved_data

        elif model == 'MagmaSat':
            XH2Ovals = []
            XCO2vals = []
            FluidProportionvals = []
            iterno = 0
            for index, row in dissolved_data.iterrows():
                iterno += 1
                if print_status:
                    percent = iterno/len(dissolved_data.index)
                    batchfile.status_bar.status_bar(percent, index)

                if file_has_temp:
                    temperature = row[temp_name]
                if temperature <= 0:
                    H2Ovals.append(np.nan)
                    CO2vals.append(np.nan)
                    XH2Ovals.append(np.nan)
                    XCO2vals.append(np.nan)
                    FluidProportionvals.append(np.nan)
                    warnings.append("Sample skipped. Bad temperature.")
                    errors.append(sys.exc_info()[0])
                    w.warn("Temperature for sample " + str(index) +
                           " is <=0. Skipping sample.", stacklevel=2)

                if file_has_press:
                    pressure = row[press_name]
                if temperature > 0 and pressure <= 0:
                    H2Ovals.append(np.nan)
                    CO2vals.append(np.nan)
                    XH2Ovals.append(np.nan)
                    XCO2vals.append(np.nan)
                    FluidProportionvals.append(np.nan)
                    warnings.append("Sample skipped. Bad pressure.")
                    errors.append(sys.exc_info()[0])
                    w.warn("Pressure for sample " + str(index) +
                           " is <=0. Skipping sample.", stacklevel=2)

                if file_has_X:
                    X_fluid = row[X_name]
                if temperature > 0 and pressure > 0 and X_fluid < 0:
                    H2Ovals.append(np.nan)
                    CO2vals.append(np.nan)
                    XH2Ovals.append(np.nan)
                    XCO2vals.append(np.nan)
                    FluidProportionvals.append(np.nan)
                    warnings.append("Sample skipped. Bad X_fluid.")
                    errors.append(sys.exc_info()[0])
                    w.warn("X_fluid for sample " + str(index) +
                           " is <0. Skipping sample.", stacklevel=2)

                if temperature > 0 and pressure > 0 and X_fluid > 1:
                    H2Ovals.append(np.nan)
                    CO2vals.append(np.nan)
                    XH2Ovals.append(np.nan)
                    XCO2vals.append(np.nan)
                    FluidProportionvals.append(np.nan)
                    warnings.append("Sample skipped. Bad X_fluid.")
                    errors.append(sys.exc_info()[0])
                    w.warn("X_fluid for sample " + str(index) +
                           " is >1. Skipping sample.", stacklevel=2)

                if (temperature > 0 and pressure > 0 and
                   X_fluid >= 0 and X_fluid <= 1):
                    try:
                        # Get sample comp as Sample class with defaults
                        bulk_comp = self.get_sample_composition(
                                      index,
                                      normalization=self.default_normalization,
                                      units='wtpt_oxides',
                                      asSampleClass=True)
                        bulk_comp.set_default_units(self.default_units)
                        bulk_comp.set_default_normalization(
                                                    self.default_normalization)
                        calc = calculate_classes.calculate_dissolved_volatiles(
                                           sample=bulk_comp, pressure=pressure,
                                           temperature=temperature,
                                           X_fluid=X_fluid, model=model,
                                           silence_warnings=True, verbose=True)
                        H2Ovals.append(calc.result['H2O_liq'])
                        CO2vals.append(calc.result['CO2_liq'])
                        XH2Ovals.append(calc.result['XH2O_fl'])
                        XCO2vals.append(calc.result['XCO2_fl'])
                        FluidProportionvals.append(
                                             calc.result['FluidProportion_wt'])
                        warnings.append(calc.calib_check)
                        errors.append('')
                    except Exception:
                        H2Ovals.append(np.nan)
                        CO2vals.append(np.nan)
                        XH2Ovals.append(np.nan)
                        XCO2vals.append(np.nan)
                        FluidProportionvals.append(np.nan)
                        warnings.append('Calculation Failed.')
                        errors.append(sys.exc_info()[0])
            dissolved_data["H2O_liq_VESIcal"] = H2Ovals
            dissolved_data["CO2_liq_VESIcal"] = CO2vals

            if file_has_temp is False:
                dissolved_data["Temperature_C_VESIcal"] = temperature
            if file_has_press is False:
                dissolved_data["Pressure_bars_VESIcal"] = pressure
            if file_has_X is False:
                dissolved_data["X_fluid_input_VESIcal"] = X_fluid
            dissolved_data["Model"] = model
            dissolved_data["Warnings"] = warnings
            if record_errors:
                dissolved_data["Errors"] = errors

            return dissolved_data

        else:
            XH2Ovals = []
            XCO2vals = []
            FluidProportionvals = []
            for index, row in dissolved_data.iterrows():
                if file_has_temp:
                    temperature = row[temp_name]
                if file_has_press:
                    pressure = row[press_name]
                if file_has_X:
                    X_fluid = row[X_name]
                # Get sample comp as Sample class with defaults
                bulk_comp = self.get_sample_composition(
                                      index,
                                      normalization=self.default_normalization,
                                      units='wtpt_oxides',
                                      asSampleClass=True)
                bulk_comp.set_default_units(self.default_units)
                bulk_comp.set_default_normalization(self.default_normalization)
                if 'Water' in model:
                    try:
                        calc = calculate_classes.calculate_dissolved_volatiles(
                                           sample=bulk_comp, pressure=pressure,
                                           temperature=temperature,
                                           X_fluid=X_fluid, model=model,
                                           silence_warnings=True)
                        H2Ovals.append(calc.result)
                        warnings.append(calc.calib_check)
                    except Exception:
                        H2Ovals.append(0)
                        warnings.append('Calculation Failed #001')
                if 'Carbon' in model:
                    try:
                        calc = calculate_classes.calculate_dissolved_volatiles(
                                           sample=bulk_comp, pressure=pressure,
                                           temperature=temperature,
                                           X_fluid=X_fluid, model=model,
                                           silence_warnings=True)
                        CO2vals.append(calc.result)
                        warnings.append(calc.calib_check)
                    except Exception:
                        CO2vals.append(0)
                        warnings.append('Calculation Failed #002')

            if 'Water' in model:
                dissolved_data["H2O_liq_VESIcal"] = H2Ovals
            if 'Carbon' in model:
                dissolved_data["CO2_liq_VESIcal"] = CO2vals
            if file_has_temp is False:
                dissolved_data["Temperature_C_VESIcal"] = temperature
            if file_has_press is False:
                dissolved_data["Pressure_bars_VESIcal"] = pressure
            if file_has_X is False:
                dissolved_data["X_fluid_input_VESIcal"] = X_fluid
            dissolved_data["Model"] = model
            dissolved_data["Warnings"] = warnings

            return dissolved_data

    def calculate_equilibrium_fluid_comp(self, temperature, pressure=None,
                                         print_status=False, model='MagmaSat',
                                         **kwargs):
        """
        Returns H2O and CO2 concentrations in wt% or mole fraction in a fluid
        in equilibrium with the given sample(s) at the given P/T condition.

        Parameters
        ----------
        sample: BatchFile object
            Compositional information on samples in oxides.

        temperature: float, int, or str
            Temperature, in degrees C. Can be passed as float, in which case
            the passed value is used as the temperature for all samples.
            Alternatively, temperature information for each individual
            sample may already be present in the BatchFile object. If so, pass
            the str value corresponding to the column title in the  BatchFile
            object.

        presure: float, int, or str
            Pressure, in bars. Can be passed as float or int, in which case
            the passed value is used as the pressure for all samples.
            Alternatively, pressure information for each individual sample may
            already be present in the BatchFile object. If so, pass the str
            value corresponding to the column title in the BatchFile object.

        model: string
            OPTIONAL: Default is 'MagmaSat'. Any other model name can be
            passed here.

        Returns
        -------
        pandas DataFrame
            Original data passed plus newly calculated values are returned.
        """
        fluid_data = self.get_data().copy()

        # Check if the model passed as the attribute "model_type"
        # Currently only implemented for MagmaSat type models
        if hasattr(model, 'model_type') is True:
            model = model.model_type

        if isinstance(temperature, str):
            file_has_temp = True
            temp_name = temperature
        elif isinstance(temperature, float) or isinstance(temperature, int):
            file_has_temp = False
        else:
            raise core.InputError("temp must be type str or float or int")

        if isinstance(pressure, str):
            file_has_press = True
            press_name = pressure
        elif (isinstance(pressure, float) or isinstance(pressure, int) or
              pressure is None):
            file_has_press = False
        else:
            raise core.InputError("pressure must be type str or float or int")

        H2Ovals = []
        CO2vals = []
        warnings = []
        if kwargs.get('verbose') is True:
            FluidMass_grams_vals = []
            FluidProportion_wt_vals = []
        if (model in models.get_model_names(model='mixed') or
           model == "MooreWater"):
            for index, row in fluid_data.iterrows():
                try:
                    if file_has_temp:
                        temperature = row[temp_name]
                    if file_has_press:
                        pressure = row[press_name]
                    # Get sample comp as Sample class with defaults
                    bulk_comp = self.get_sample_composition(
                                      index,
                                      normalization=self.default_normalization,
                                      units='wtpt_oxides', asSampleClass=True)
                    bulk_comp.set_default_units(self.default_units)
                    bulk_comp.set_default_normalization(
                                                    self.default_normalization)

                    calc = calculate_classes.calculate_equilibrium_fluid_comp(
                                           sample=bulk_comp, pressure=pressure,
                                           temperature=temperature,
                                           model=model, silence_warnings=True,
                                           **kwargs)

                    H2Ovals.append(calc.result['H2O'])
                    CO2vals.append(calc.result['CO2'])
                    if calc.result['H2O'] == 0 and calc.result['CO2'] == 0:
                        warnings.append(calc.calib_check + "Sample not " +
                                        "saturated at these conditions")
                    else:
                        warnings.append(calc.calib_check)
                except Exception:
                    H2Ovals.append(np.nan)
                    CO2vals.append(np.nan)
                    warnings.append("Calculation Failed.")
            fluid_data["XH2O_fl_VESIcal"] = H2Ovals
            fluid_data["XCO2_fl_VESIcal"] = CO2vals
            if file_has_temp is False:
                fluid_data["Temperature_C_VESIcal"] = temperature
            if file_has_press is False:
                fluid_data["Pressure_bars_VESIcal"] = pressure
            fluid_data["Model"] = model
            fluid_data["Warnings"] = warnings

            return fluid_data
        elif model == 'MagmaSat':
            iterno = 0
            for index, row in fluid_data.iterrows():
                iterno += 1
                if print_status:
                    percent = iterno/len(fluid_data.index)
                    batchfile.status_bar.status_bar(percent, index)

                if file_has_temp:
                    temperature = row[temp_name]
                if temperature <= 0:
                    H2Ovals.append(np.nan)
                    CO2vals.append(np.nan)
                    warnings.append("Calculation skipped. Bad temperature.")
                    w.warn("Temperature for sample " + str(index) +
                           " is <=0. Skipping sample.",
                           stacklevel=2)

                if file_has_press:
                    pressure = row[press_name]
                if temperature > 0 and pressure <= 0:
                    H2Ovals.append(np.nan)
                    CO2vals.append(np.nan)
                    warnings.append("Calculation skipped. Bad pressure.")
                    w.warn("Pressure for sample " + str(index) +
                           " is <=0. Skipping sample.", stacklevel=2)

                if temperature > 0 and pressure > 0:
                    try:
                        # Get sample comp as Sample class with defaults
                        bulk_comp = self.get_sample_composition(
                                      index,
                                      normalization=self.default_normalization,
                                      units='wtpt_oxides', asSampleClass=True)
                        bulk_comp.set_default_units(self.default_units)
                        bulk_comp.set_default_normalization(
                                                    self.default_normalization)

                        calc = (
                            calculate_classes.calculate_equilibrium_fluid_comp(
                                           sample=bulk_comp, pressure=pressure,
                                           temperature=temperature,
                                           model=model, silence_warnings=True,
                                           **kwargs))

                        H2Ovals.append(calc.result['H2O'])
                        CO2vals.append(calc.result['CO2'])
                        if kwargs.get('verbose') is True:
                            FluidMass_grams_vals.append(calc.result['FluidMass_grams'])
                            FluidProportion_wt_vals.append(calc.result['FluidProportion_wt'])
                        if calc.result['H2O'] == 0 and calc.result['CO2'] == 0:
                            warnings.append(calc.calib_check + "Sample not " +
                                            "saturated at these conditions")
                        else:
                            warnings.append(calc.calib_check)
                    except Exception:
                        H2Ovals.append(np.nan)
                        CO2vals.append(np.nan)
                        warnings.append("Calculation Failed.")
                        if kwargs.get('verbose') is True:
                            FluidMass_grams_vals.append(np.nan)
                            FluidProportion_wt_vals.append(np.nan)
            fluid_data["XH2O_fl_VESIcal"] = H2Ovals
            fluid_data["XCO2_fl_VESIcal"] = CO2vals
            if kwargs.get('verbose') is True:
                fluid_data["FluidMass_grams"] = FluidMass_grams_vals
                fluid_data["FluidProportion_wt"] = FluidProportion_wt_vals
            if file_has_temp is False:
                fluid_data["Temperature_C_VESIcal"] = temperature
            if file_has_press is False:
                fluid_data["Pressure_bars_VESIcal"] = pressure
            fluid_data["Model"] = model
            fluid_data["Warnings"] = warnings

            return fluid_data

        else:
            saturated = []
            for index, row in fluid_data.iterrows():
                try:
                    if file_has_temp:
                        temperature = row[temp_name]
                    if file_has_press:
                        pressure = row[press_name]
                    # Get sample comp as Sample class with defaults
                    bulk_comp = self.get_sample_composition(
                                      index,
                                      normalization=self.default_normalization,
                                      units='wtpt_oxides', asSampleClass=True)
                    bulk_comp.set_default_units(self.default_units)
                    bulk_comp.set_default_normalization(
                                                    self.default_normalization)

                    calc = calculate_classes.calculate_equilibrium_fluid_comp(
                                           sample=bulk_comp, pressure=pressure,
                                           temperature=temperature,
                                           model=model, silence_warnings=True)
                    saturated.append(calc.result)
                    warnings.append(calc.calib_check)
                except Exception:
                    saturated.append(np.nan)
                    warnings.append("Calculation Failed.")
            fluid_data["Saturated_VESIcal"] = saturated
            if file_has_temp is False:
                fluid_data["Temperature_C_VESIcal"] = temperature
            if file_has_press is False:
                fluid_data["Pressure_bars_VESIcal"] = pressure
            fluid_data["Model"] = model
            fluid_data["Warnings"] = warnings

            return fluid_data

    def calculate_saturation_pressure(self, temperature, print_status=None,
                                      model='MagmaSat', **kwargs):
        """
        Calculates the saturation pressure of multiple sample compositions in
        the BatchFile.

        Parameters
        ----------
        temperature: float, int, or str
            Temperature at which to calculate saturation pressures, in
            degrees C. Can be passed as float or int, in which case the
            passed value is used as the temperature for all samples.
            Alternatively, temperature information for each individual sample
            may already be present in the passed BatchFile object. If so, pass
            the str value corresponding to the column title in the passed
            BatchFile object.

        print_status: bool
            OPTIONAL: The default value for MagmaSat is True and the default
            for all other models is False. If set to True, the progress of the
            calculation will be printed to the terminal. If set to False,
            nothing will be printed. MagmaSat calculations tend to be slow,
            and so a value of True is recommended more most use cases.

        model: string
            OPTIONAL: Default is 'MagmaSat'. Any other model name can be
            passed here.

        Returns
        -------
        pandas DataFrame object
            Values returned are saturation pressure in bars, the mass of
            fluid present, and the composition of the fluid present.
        """
        satp_data = self.get_data().copy()

        # Check if the model passed has the attribute "model_type"
        # Currently only implemented for MagmaSat type models
        if hasattr(model, 'model_type') is True:
            model = model.model_type

        # set default print_status to True for MagmaSat, False for other
        # models if user doesn't pass any option
        if print_status is None:
            if model == 'MagmaSat':
                print_status = True
            else:
                print_status = False

        if isinstance(temperature, str):
            file_has_temp = True
            temp_name = temperature
        elif isinstance(temperature, float) or isinstance(temperature, int):
            file_has_temp = False
        else:
            raise core.InputError("temperature must be type str or float or "
                                  "int")

        if model != 'MagmaSat':
            satP = []
            warnings = []
            piStar = []
            iterno = 0
            for index, row in satp_data.iterrows():
                iterno += 1
                if print_status:
                    percent = iterno/len(satp_data.index)
                    batchfile.status_bar.status_bar(percent, index)
                if file_has_temp:
                    temperature = row[temp_name]
                # Get sample comp as Sample class with defaults
                bulk_comp = self.get_sample_composition(
                                      index,
                                      normalization=self.default_normalization,
                                      units='wtpt_oxides', asSampleClass=True)
                bulk_comp.set_default_units(self.default_units)
                bulk_comp.set_default_normalization(self.default_normalization)

                calc = calculate_classes.calculate_saturation_pressure(
                                     sample=bulk_comp, temperature=temperature,
                                     model=model, silence_warnings=True,
                                     **kwargs)
                satP.append(calc.result)
                warnings.append(calc.calib_check)
                if model == 'ShishkinaIdealMixing':
                    piStar.append(
                           models.default_models['ShishkinaIdealMixing'
                                                 ].models[1].PiStar(bulk_comp))

            satp_data["SaturationP_bars_VESIcal"] = satP
            if file_has_temp is False:
                satp_data["Temperature_C_VESIcal"] = temperature
            satp_data["Model"] = model
            satp_data["Warnings"] = warnings
            if len(piStar) > 0:
                satp_data['PiStar_VESIcal'] = piStar

            return satp_data

        elif model == 'MagmaSat':
            satP = []
            flmass = []
            flH2O = []
            flCO2 = []
            flsystem_wtper = []
            warnings = []
            iterno = 0
            for index, row in satp_data.iterrows():
                iterno += 1
                if print_status:
                    percent = iterno/len(satp_data.index)
                    batchfile.status_bar.status_bar(percent, index)

                if file_has_temp:
                    temperature = row[temp_name]
                if temperature <= 0:
                    satP.append(np.nan)
                    flmass.append(np.nan)
                    flsystem_wtper.append(np.nan)
                    flH2O.append(np.nan)
                    flCO2.append(np.nan)
                    warnings.append("Calculation skipped. Bad temperature.")
                    w.warn("Temperature for sample " + str(index) +
                           " is <=0. Skipping sample.", stacklevel=2)

                if temperature > 0:
                    try:
                        # Get sample comp as Sample class with defaults
                        bulk_comp = self.get_sample_composition(
                                      index,
                                      normalization=self.default_normalization,
                                      units='wtpt_oxides', asSampleClass=True)
                        bulk_comp.set_default_units(self.default_units)
                        bulk_comp.set_default_normalization(
                                                    self.default_normalization)

                        calc = calculate_classes.calculate_saturation_pressure(
                                     sample=bulk_comp, temperature=temperature,
                                     model=model, verbose=True,
                                     silence_warnings=True)
                        satP.append(calc.result["SaturationP_bars"])
                        flmass.append(calc.result["FluidMass_grams"])
                        flsystem_wtper.append(
                                             calc.result["FluidProportion_wt"])
                        flH2O.append(calc.result["XH2O_fl"])
                        flCO2.append(calc.result["XCO2_fl"])
                        warnings.append(calc.calib_check)
                    except Exception:
                        satP.append(np.nan)
                        flmass.append(np.nan)
                        flsystem_wtper.append(np.nan)
                        flH2O.append(np.nan)
                        flCO2.append(np.nan)
                        warnings.append("Calculation Failed")

            satp_data["SaturationP_bars_VESIcal"] = satP
            if file_has_temp is False:
                satp_data["Temperature_C_VESIcal"] = temperature
            satp_data["XH2O_fl_VESIcal"] = flH2O
            satp_data["XCO2_fl_VESIcal"] = flCO2
            satp_data["FluidMass_grams_VESIcal"] = flmass
            satp_data["FluidSystem_wt_VESIcal"] = flsystem_wtper
            satp_data["Model"] = model
            satp_data["Warnings"] = warnings

            return satp_data


def BatchFile_from_DataFrame(dataframe, units='wtpt_oxides', label=None):
    """
    Transforms any pandas DataFrame object into a VESIcal BatchFile object.

    Parameters
    ----------
    Same as batchfile.BatchFile_from_DataFrame()

    Returns
    -------
    VESIcal.BatchFile object
    """
    return BatchFile(filename=None, dataframe=dataframe, units=units,
                     label=label)
