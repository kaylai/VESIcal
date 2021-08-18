############################################
Writing new calculations for existing models
############################################
.. contents::

Download the Jupyter notebook `here <https://github.com/kaylai/VESIcal/blob/master/docs/jupyter_notebooks/adv_newcalcs.ipynb>`_.

Any h2llo model within VESIcal can be duplicated and then modified. New calculation methods that rely on the pre-existing structure within VESIcal, can be easily added. An example is shown here for a new method, :py:meth:`calculate_dissolved_CO2()`. This method will be used to calculate the dissolved CO2 concentration, in wt%, of a mixed fluid (H2O-CO2) saturated magma with a known dissolved H2O concentration, pressure, temperature, and bulk composition. This magma would thus be undersaturated in pure H2O but saturated in mixed H2O-CO2 fluid.

Creating and modifying model objects
====================================

First, we duplicate the MagmaSat model object.

.. code-block:: python

	# duplicate the MagmaSat model into a new object called mymodel
	mymodel = v.models.magmasat.MagmaSat()

We now have a new model object, called mymodel, that behaves exactly as the MagmaSat model would. You can test this by performing any VESIcal calculation and setting the argument 'model' to be 'mymodel'. Let's first create a sample to test this with:

.. code-block:: python

	# create a sample to test new method
	testsample = v.Sample({'SiO2':    47.95,
	                       'TiO2':    1.67,
	                       'Al2O3':   17.32,
	                       'FeO':     10.24,
	                       'Fe2O3':   0.1,
	                       'MgO':     5.76,
	                       'CaO':     10.93,
	                       'Na2O':    3.45,
	                       'K2O':     1.99,
	                       'P2O5':    0.51,
	                       'MnO':     0.1,
	                       'H2O':     2.0
	                    })

Now, let's try calculating the saturation pressure of testsample at 1000 degrees C using the mymodel object we've just created.

.. code-block:: python

	satP = v.calculate_saturation_pressure(sample=testsample, temperature=1000, model=mymodel).result
	print("The saturation pressure is " + str(satP) + " bars.")

.. code-block:: python

	The saturation pressure is 430 bars.

Creating new calculations using an existing model
-------------------------------------------------
Now, let's create a new calculation (called a method) within the framework of mymodel, which is a duplicate of MagmaSat. In the case of MagmaSat, the MELTS model also needs to be defined separately (with other models, this step is unnecessary).

.. code-block:: python

	# import MELTS preamble (copy-paste from magmasat.py)
	from thermoengine import equilibrate

	# -------------- MELTS preamble --------------- #
	# instantiate thermoengine equilibrate MELTS instance
	melts = equilibrate.MELTSmodel('1.2.0')

	# Suppress phases not required in the melts simulation
	phases = melts.get_phase_names()
	for phase in phases:
	    melts.set_phase_inclusion_status({phase: False})
	melts.set_phase_inclusion_status({'Fluid': True, 'Liquid': True})
	# --------------------------------------------- #

Now that MELTS is defined, we can get to work creating a new calculation. Below, we have taken the code from MagmaSat's :py:meth:`calculate_dissolved_volatiles` method and adjusted it to meet our needs. We begin with a magma that we know to be pure-H2O undersaturated but mixed H2O-CO2 saturated, with known dissolved H2O concentration in wt%, pressure in bars, and temperature in degrees C. The goal is to write a method to calculate the concentration of dissolved CO2 necessary to acheive mixed H2O-CO2 saturation at the given conditions.

.. code-block:: python

	# write a new method to add to our mymodel object
	# this method will hold H2O constant and calculate dissolved CO2 at given P, T
	# this is a modification of the existing calculate_dissolved_volatiles() method

	def calculate_dissolved_CO2(self, sample, temperature, pressure,
	                            H2O_liq, verbose=False, **kwargs):
	        """
	        Calculates the amount of CO2 dissolved in a magma at saturation at the given P/T
	        conditions and given dissolved H2O. 

	        Parameters
	        ----------
	        sample:     Sample class
	            Magma major element composition.

	        temperature: float or int
	            Temperature, in degrees C.

	        presure: float or int
	            Pressure, in bars.

	        H2O_liq: float or int
	            Dissolved H2O concentration, in wt%

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

	        if isinstance(H2O_liq, int) or isinstance(H2O_liq, float):
	            pass
	        else:
	            raise core.InputError("H2O_liq must be type int or float")

	        pressureMPa = pressure / 10.0

	        # coarse search
	        H2O_bulk = H2O_liq
	        CO2_bulk = 0.0
	        fluid_mass = 0.0
	        while fluid_mass <= 0:
	            CO2_bulk += 0.01
	            fluid_mass = self.get_fluid_mass(_sample, temperature, pressure, H2O_bulk, CO2_bulk)
	        
	        # calculated dissolved H2O, then increment up
	        H2O_diss = 0
	        while H2O_diss < H2O_liq:
	            _sample.change_composition({'H2O': H2O_bulk, 'CO2': CO2_bulk}, units='wtpt_oxides')
	            melts.set_bulk_composition(_sample.get_composition(units='wtpt_oxides',
	                                                               normalization='none'))

	            output = melts.equilibrate_tp(temperature, pressureMPa, initialize=True)
	            (status, temperature, pressureMPa, xmlout) = output[0]
	            liquid_comp = melts.get_composition_of_phase(xmlout, phase_name='Liquid', mode='oxide_wt')

	            if "H2O" in liquid_comp:
	                H2O_diss = liquid_comp["H2O"]
	            else:
	                H2O_diss = 0
	            # changing this value changes how close to the original
	            # known H2O value the resulting H2O_liquid wt% will be
	            H2O_bulk += 0.001
	        
	        H2O_val = H2O_bulk
	        CO2_val = CO2_bulk

	        # ------ Get calculated values ------ #
	        _sample.change_composition({'H2O': H2O_val, 'CO2': CO2_val}, units='wtpt_oxides')
	        melts.set_bulk_composition(_sample.get_composition(units='wtpt_oxides',
	                                                           normalization='none'))

	        output = melts.equilibrate_tp(temperature, pressureMPa, initialize=True)
	        (status, temperature, pressureMPa, xmlout) = output[0]
	        fluid_mass = melts.get_mass_of_phase(xmlout, phase_name='Fluid')
	        system_mass = melts.get_mass_of_phase(xmlout, phase_name='System')
	        liquid_comp = melts.get_composition_of_phase(xmlout, phase_name='Liquid', mode='oxide_wt')
	        fluid_comp = melts.get_composition_of_phase(xmlout, phase_name='Fluid', mode='component')

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
	            return {"temperature": temperature, "pressure": pressure,
	                    "H2O_liq": H2O_liq, "CO2_liq": CO2_liq,
	                    "XH2O_fl": H2O_fl, "XCO2_fl": CO2_fl,
	                    "FluidProportion_wt": 100*fluid_mass/system_mass}

	        if verbose is False:
	            return {"CO2_liq": CO2_liq, "H2O_liq": H2O_liq}	

In order to bind our newly created method to mymodel (in other words, in order to allow mymodel to access and execute the code we have just written), we use python's universal .get method.

.. code-block:: python

	# add our new method to mymodel
	mymodel.calculate_dissolved_CO2 = calculate_dissolved_CO2.__get__(mymodel)

Now, let's test our new method.

.. code-block:: python

	mymodel.calculate_dissolved_CO2(testsample, pressure=5000.0,
					temperature=1000.0,
					H2O_liq=testsample.get_composition()['H2O'],
					verbose=True)

.. code-block:: python

	{'temperature': 1000.0,
	 'pressure': 5000.0,
	 'H2O_liq': 2.00137179156045,
	 'CO2_liq': 0.520077091019401,
	 'XH2O_fl': 0.146938900738546,
	 'XCO2_fl': 0.853061099261454,
	 'FluidProportion_wt': 0.006638675057019597}

Some notes about our new calculation
------------------------------------
Notice that the final output dissolved H2O concentration matches our given H2O concentration of 2.0 wt% to within ~0.001 wt%. This final output value can be made to match much more closely to the given H2O concentration by adjusting one line of code. See comment "changing this value changes how close to the original known H2O value the resulting H2O_liquid wt% will be".


