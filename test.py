from VESIcal import *
from VESIcal import activity_models
from VESIcal import calculate_classes
from VESIcal import calibration_checks
from VESIcal import core
from VESIcal import fugacity_models
from VESIcal import models
from VESIcal import sample_class

# THIS ROUTINE SHOULD BE REMOVED ONCE THE 'PROPER' TESTING ROUTINES HAVE BEEN IMPLEMENTED.
def test():
    """
    This is a set of tests for the module, firstly to ensure all of the functions that should run, do run.
    Secondly, the output can be used to check the module reproduces the results in each of the manuscripts
    describing the models.
    """

    sample_comp = {'SiO2':47.95,
                 'TiO2':1.67,
                 'Al2O3':17.32,
                 'FeO':10.24,
                 'Fe2O3':0.1,
                 'MgO':5.76,
                 'CaO':10.93,
                 'Na2O':3.45,
                 'K2O':1.99,
                 'P2O5':0.51,
                 'MnO':0.1,
                 'CO2':0.08,
                 'H2O':4.0}

    test_sample = sample_class.Sample(sample_comp)

    test_pressure = 2000.0
    test_temperature = 1473.15
    test_pressure_list = [1000.0,2000.0,5000.0]
    test_isopleth_list = [0.0,0.5,1.0]


    print("\n================================\n= MAGMASATPLUS TESTING ROUTINE =\n================================")

    print("\n This routine will check that key methods run using typical values of variables. \
The routine does not check that the results are correct (though this may be obvious from the outputs),\
nor does it check every possible iteration of methods and input types. It will check that an update hasn't \
COMPLETELY broken the module.")

    for i, model_name in zip(range(len(models.default_models)),list(models.default_models.keys())):
        print("\nTesting model {:d} of {:d}: {:s}".format(i+1,len(models.default_models),model_name))

        ### calculate_dissolved_volatiles
        model = models.default_models[model_name]
        print("Model contains "+" ".join(model.volatile_species))

        print("Testing calculate_dissolved_volatiles method...")

        if len(model.volatile_species) == 1:
            dissolved = model.calculate_dissolved_volatiles(pressure=test_pressure,temperature=test_temperature,
                                                            sample=test_sample)
        else:
            X_fluid = 1.0/len(model.volatile_species)
            X_fluid = tuple(X_fluid for x in range(len(model.volatile_species)))
            print("Setting X_fluid to "+str(X_fluid))
            dissolved = model.calculate_dissolved_volatiles(pressure=test_pressure,temperature=test_temperature,
                                                            sample=test_sample,X_fluid=X_fluid)

        if len(model.volatile_species) == 1:
            print("  {:s} solubility at {:.0f} bars and {:.0f} K is {:.3f} wt%".format(model.volatile_species[0],
                                                                                     test_pressure,test_temperature,
                                                                                     dissolved))
        else:
            for i,volatile in zip(range(len(model.volatile_species)),model.volatile_species):
                print("  {:s} solubility at {:.0f} bars and {:.0f} K is {:.3f} wt%".format(volatile,
                                                                                         test_pressure,test_temperature,
                                                                                         dissolved[i]))

        print("Testing calculate_dissolved_volatiles class interface...")
        if len(model.volatile_species) == 1:
            result = calculate_dissolved_volatiles(sample=test_sample,pressure=test_pressure,
                                                   temperature=test_temperature,model=model_name)
        else:
            result = calculate_dissolved_volatiles(sample=test_sample,pressure=test_pressure,
                                                   temperature=test_temperature,model=model_name,
                                                   X_fluid=X_fluid)

        if len(model.volatile_species) == 1:
            print("  {:s} solubility at {:.0f} bars and {:.0f} K is {:.3f} wt%".format(model.volatile_species[0],
                                                                                     test_pressure,test_temperature,
                                                                                     result.result))
        else:
            for volatile in model.volatile_species:
                print("  {:s} solubility at {:.0f} bars and {:.0f} K is {:.3f} wt%".format(volatile,
                                                                                         test_pressure,test_temperature,
                                                                                         result.result[volatile+'_liq']))



        ### calculate_saturation_pressure
        print("Testing calculate_saturation_pressure method...")
        satP = model.calculate_saturation_pressure(sample=test_sample,temperature=test_temperature)
        if len(model.volatile_species) == 1:
            print("  A concentration of {:.2f} wt% {:s} is saturated at {:.0f} bars at {:.0f} K".format(test_sample.get_composition(model.volatile_species[0]),
                                                                                                       model.volatile_species[0],satP,
                                                                                                       test_temperature))
        else:
            concstr = ""
            for volatile in model.volatile_species:
                concstr += "{:.2f}".format(test_sample.get_composition(volatile))
                concstr += " wt% "
                concstr += volatile
                concstr += ", "
            print("  Concentrations of "+concstr[:-1]+" are saturated at {:.0f} bars at {:.0f} K".format(satP,test_temperature))

        print("Testing calculate_saturation_pressure class interface...")
        satP = calculate_saturation_pressure(model=model_name,sample=test_sample,temperature=test_temperature).result
        if len(model.volatile_species) == 1:
            print("  A concentration of {:.2f} wt% {:s} is saturated at {:.0f} bars at {:.0f} K".format(test_sample.get_composition(model.volatile_species[0]),
                                                                                                       model.volatile_species[0],satP,
                                                                                                       test_temperature))
        else:
            concstr = ""
            for volatile in model.volatile_species:
                concstr += "{:.2f}".format(test_sample.get_composition(volatile))
                concstr += " wt% "
                concstr += volatile
                concstr += ", "
            print("  Concentrations of "+concstr[:-1]+" are saturated at {:.0f} bars at {:.0f} K".format(satP,test_temperature))


        ### calculate_equilibrium_fluid_comp
        print("Testing calculate_equilibrium_fluid_comp method...")
        fluid = model.calculate_equilibrium_fluid_comp(sample=test_sample,temperature=test_temperature,pressure=test_pressure)

        if len(model.volatile_species) == 1:
            print("  A mole fraction of {:.2f} of {:s} is present in the fluid.".format(fluid,model.volatile_species[0]))
        else:
            fluidstr = ""
            for i,volatile in zip(range(len(model.volatile_species)),model.volatile_species):
                fluidstr += "{:.2f}".format(fluid[model.volatile_species[i]])
                fluidstr += " "
                fluidstr += volatile
                fluidstr += ", "
            print("  Mole fractions of "+fluidstr +"are present in the fluid.")
            if np.sum(list(fluid.values())) != 0.0 and np.sum(list(fluid.values())) != 1.0:
                print("  WARNING: MOLE FRACTIONS DO NOT SUM TO 1.0")

        print("Testing calculate_equilibrium_fluid_comp class interface...")
        fluid = model.calculate_equilibrium_fluid_comp(model=model_name,sample=test_sample,
                                                       temperature=test_temperature,pressure=test_pressure)
        if len(model.volatile_species) == 1:
            print("  A mole fraction of {:.2f} of {:s} is present in the fluid.".format(fluid,model.volatile_species[0]))
        else:
            fluidstr = ""
            for i,volatile in zip(range(len(model.volatile_species)),model.volatile_species):
                fluidstr += "{:.2f}".format(fluid[model.volatile_species[i]])
                fluidstr += " "
                fluidstr += volatile
                fluidstr += ", "
            print("  Mole fractions of "+fluidstr +"are present in the fluid.")
            if np.sum(list(fluid.values())) != 0.0 and np.sum(list(fluid.values())) != 1.0:
                print("  WARNING: MOLE FRACTIONS DO NOT SUM TO 1.0")

        ### calculate_isobars_and_isopleths
        if len(model.volatile_species) > 1:
            print("Testing calculate_isobars_and_isopleths method...")
            isobars, isopleths = model.calculate_isobars_and_isopleths(pressure_list=test_pressure_list,
                                                                       isopleth_list=test_isopleth_list,
                                                                       sample=test_sample,
                                                                       temperature=test_temperature)
            print("Isobars:")
            print(isobars)
            print("\nIsopleths:")
            print(isopleths)

            print("Testing calculate_isobars_and_isopleths class interface...")
            isobars, isopleths = calculate_isobars_and_isopleths(model=model_name,pressure_list=test_pressure_list,
                                                                 isopleth_list=test_isopleth_list,
                                                                 sample=test_sample,
                                                                 temperature=test_temperature).result
            print("Isobars:")
            print(isobars)
            print("\nIsopleths:")
            print(isopleths)

        ### calculate_degassing_path
        test_sample = Sample(sample_comp)
        if len(model.volatile_species) > 1:
            print("Testing calculate_degassing_path method...")
            degassing = model.calculate_degassing_path(sample=test_sample,temperature=test_temperature)
            print("  Degassing path:")
            print(degassing)

            print("Testing calculate_degassing_path class interface...")
            degassing = calculate_degassing_path(model=model_name,sample=test_sample,
                                                        temperature=test_temperature).result
            print("  Degassing path:")
            print(degassing)


    print("\nTesting routine complete.\n")




test()
