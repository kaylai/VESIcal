import unittest
import VESIcal as v

class TestDissolvedVolatiles(unittest.TestCase):
    def setUp(self):
        # Sample with units as wtpt_oxides
        self.majors_wtpt = {'SiO2':    47.95,
                         'TiO2':    1.67,
                         'Al2O3':   17.32,
                         'FeO':     10.24,
                         'Fe2O3':   0.1,
                         'MgO':     5.76,
                         'CaO':     10.93,
                         'Na2O':    3.45,
                         'K2O':     1.99,
                         'P2O5':    0.51,
                         'MnO':     0.1}

        # create Sample object and set default units to wtpt_oxides
        self.sample_wtpt = v.Sample(self.majors_wtpt)
        self.sample_wtpt.set_default_units("wtpt_oxides")

        # create sample and set default units to mol_oxides
        self.sample_molox = v.Sample(self.majors_wtpt)
        self.sample_molox.set_default_units("mol_oxides")

        # BatchFile with test sample as defined above in wtpt_oxides
        try:
            self.batch_wtpt = v.BatchFile('BatchTest.xlsx', units='wtpt_oxides')
        except:
            self.batch_wtpt = v.BatchFile('tests/BatchTest.xlsx', units='wtpt_oxides')
        self.batch_wtpt.set_default_units("wtpt_oxides")

        # BatchFile with test sample as defined above in mol_oxides
        try:
            self.batch_molox = v.BatchFile('BatchTest.xlsx')
        except:
            self.batch_molox = v.BatchFile('tests/BatchTest.xlsx')
        self.batch_molox.set_default_units("mol_oxides")

        # solubilities calculated with VESIcal in wt% and translated externally to mol fraction
        self.shishkinaMixed_wtpt =      {'H2O_liq': 2.24019,
                                        'CO2_liq': 0.028104173}
        self.shishkinaMixed_molox =     {'H2O_liq': 0.07415139,
                                        'CO2_liq': 0.000379788}

        self.dixonMixed_wtpt =          {'H2O_liq': 2.06748958,
                                        'CO2_liq': 0.035434893}
        self.dixonMixed_molox =         {'H2O_liq': 0.0688214,
                                        'CO2_liq': 0.000481216}

        self.iaconomarzianoMixed_wtpt = {'H2O_liq': 2.0853469,
                                        'CO2_liq': 0.0500067}
        self.iaconomarzianoMixed_molox ={'H2O_liq': 0.0693591,
                                        'CO2_liq': 0.000680900}

        self.liuMixed_wtpt =            {'H2O_liq': 2.293923234,
                                        'CO2_liq': 0.029485739}
        self.liuMixed_molox =           {'H2O_liq': 0.075793677,
                                        'CO2_liq': 0.000398905}

        self.magmasat_wtpt =            {'H2O_liq': 1.982045999,
                                        'CO2_liq': 0.044939669}
        self.magmasat_molox =           {'H2O_liq': 0.06615768,
                                        'CO2_liq': 0.000614167}

        self.allisonCarbon_wtpt             = 0.06101661
        self.allisonCarbon_molox            = 0.00089276
        self.allisonCarbon_sunset_wtpt      = 0.03249080
        self.allisonCarbon_sunset_molox     = 0.00047558
        self.allisonCarbon_sfvf_wtpt        = 0.02660018
        self.allisonCarbon_sfvf_molox       = 0.00038939
        self.allisonCarbon_erebus_wtpt      = 0.03314730
        self.allisonCarbon_erebus_molox     = 0.00048519
        self.allisonCarbon_vesuvius_wtpt    = 0.06101661
        self.allisonCarbon_vesuvius_molox   = 0.00089276
        self.allisonCarbon_etna_wtpt        = 0.04799322
        self.allisonCarbon_etna_molox       = 0.00070234
        self.allisonCarbon_stromboli_wtpt   = 0.03216744
        self.allisonCarbon_stromboli_molox  = 0.00047085
        self.mooreWater_wtpt                = 2.24343860
        self.mooreWater_molox               = 0.07427734

        self.mixed_dict =   {"ShishkinaIdealMixing": {'wtpt_oxides': self.shishkinaMixed_wtpt, 'mol_oxides': self.shishkinaMixed_molox},
                             "Dixon"               : {'wtpt_oxides': self.dixonMixed_wtpt, 'mol_oxides': self.dixonMixed_molox},
                             "IaconoMarziano"      : {'wtpt_oxides': self.iaconomarzianoMixed_wtpt, 'mol_oxides': self.iaconomarzianoMixed_molox},
                             "Liu"                 : {'wtpt_oxides': self.liuMixed_wtpt, 'mol_oxides': self.liuMixed_molox},
                             "MagmaSat"            : {'wtpt_oxides': self.magmasat_wtpt, 'mol_oxides': self.magmasat_molox}
                             }

        self.nonmixed_dict = {"AllisonCarbon"           : {'wtpt_oxides': self.allisonCarbon_wtpt, 'mol_oxides': self.allisonCarbon_molox},
                              "AllisonCarbon_sunset"    : {'wtpt_oxides': self.allisonCarbon_sunset_wtpt, 'mol_oxides': self.allisonCarbon_sunset_molox},
                              'AllisonCarbon_sfvf'      : {'wtpt_oxides': self.allisonCarbon_sfvf_wtpt, 'mol_oxides': self.allisonCarbon_sfvf_molox},
                              'AllisonCarbon_erebus'    : {'wtpt_oxides': self.allisonCarbon_erebus_wtpt, 'mol_oxides': self.allisonCarbon_erebus_molox},
                              'AllisonCarbon_vesuvius'  : {'wtpt_oxides': self.allisonCarbon_vesuvius_wtpt, 'mol_oxides': self.allisonCarbon_vesuvius_molox},
                              'AllisonCarbon_etna'      : {'wtpt_oxides': self.allisonCarbon_etna_wtpt, 'mol_oxides': self.allisonCarbon_etna_molox},
                              'AllisonCarbon_stromboli' : {'wtpt_oxides': self.allisonCarbon_stromboli_wtpt, 'mol_oxides': self.allisonCarbon_stromboli_molox},
                              'MooreWater'              : {'wtpt_oxides': self.mooreWater_wtpt, 'mol_oxides': self.mooreWater_molox},
                            }

    def test_calculate_single_wtpt_mixed(self):
        for model in self.mixed_dict.keys():
            calcd_result = v.calculate_dissolved_volatiles(self.sample_wtpt, pressure=1000, temperature=1000, X_fluid=0.5, model=model, verbose=False).result
            known_result = self.mixed_dict[model]['wtpt_oxides']
            for k in calcd_result.keys():
                self.assertAlmostEqual(calcd_result[k], known_result[k], places=4)

    # def test_calculate_batch_wtpt_mixed(self):
    #     for model in self.mixed_dict.keys():
    #         calcd_result = self.batch_wtpt.calculate_dissolved_volatiles(self.sample_wtpt, pressure=1000, temperature=1000, X_fluid=0.5, model=model, verbose=False).result
    #         known_result = self.mixed_dict[model]['wtpt_oxides']
    #         for k in calcd_result.keys():
    #             self.assertAlmostEqual(calcd_result[k], known_result[k], places=4)

    def test_calculate_single__wtpt_nonmixed(self):
        for model in self.nonmixed_dict.keys():
            calcd_result = v.calculate_dissolved_volatiles(self.sample_wtpt, pressure=1000, temperature=1000, X_fluid=0.5, model=model, verbose=False).result
            known_result = self.nonmixed_dict[model]['wtpt_oxides']
            self.assertAlmostEqual(calcd_result, known_result, places=4)

    def test_calculate_single__molox_mixed(self):
        for model in self.mixed_dict.keys():
            calcd_result = v.calculate_dissolved_volatiles(self.sample_molox, pressure=1000, temperature=1000, X_fluid=0.5, model=model, verbose=False).result
            known_result = self.mixed_dict[model]['mol_oxides']
            for k in calcd_result.keys():
                self.assertAlmostEqual(calcd_result[k], known_result[k], places=4)

    def test_calculate_single__molox_nonmixed(self):
        for model in self.nonmixed_dict.keys():
            calcd_result = v.calculate_dissolved_volatiles(self.sample_molox, pressure=1000, temperature=1000, X_fluid=0.5, model=model, verbose=False).result
            known_result = self.nonmixed_dict[model]['mol_oxides']
            self.assertAlmostEqual(calcd_result, known_result, places=4)

class TestSaturationPressure(unittest.TestCase):
    def setUp(self):
        # Sample with units as wtpt_oxides
        self.majors_wtpt = {'SiO2':    47.95,
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
                         'H2O':     2.0,
                         'CO2':     0.1}

        # Set conditions of calculation
        self.temperature = 1000

        # create Sample object and set default units to wtpt_oxides
        self.sample_wtpt = v.Sample(self.majors_wtpt)
        self.sample_wtpt.set_default_units("wtpt_oxides")

        # create sample and set default units to mol_oxides
        self.sample_molox = v.Sample(self.majors_wtpt)
        self.sample_molox.set_default_units("mol_oxides")

        # BatchFile with test sample as defined above in wtpt_oxides
        try:
            self.batch_wtpt = v.BatchFile('BatchTest.xlsx', units='wtpt_oxides')
        except:
            self.batch_wtpt = v.BatchFile('tests/BatchTest.xlsx', units='wtpt_oxides')
        self.batch_wtpt.set_default_units("wtpt_oxides")

        # BatchFile with test sample as defined above in mol_oxides
        try:
            self.batch_molox = v.BatchFile('BatchTest.xlsx')
        except:
            self.batch_molox = v.BatchFile('tests/BatchTest.xlsx')
        self.batch_molox.set_default_units("mol_oxides")

        # saturation pressures calculated with VESIcal
        self.shishkinaMixed =      1903.1192979127836
        self.dixonMixed =          1847.1637265676327
        self.iaconomarzianoMixed = 1437.268470833670
        self.liuMixed =            2157.191497153953
        self.magmasat =            1630  # Generated with MagmaSat app

        self.shishkinaCarbon           = 1507.647195272281
        self.dixonCarbon               = 1375.61109469857
        self.iaconomarzianoCarbon      = 1215.0221712464863
        self.allisonCarbon             = 821.9428189121331
        self.allisonCarbon_sunset      = 1468.1852834512158
        self.allisonCarbon_sfvf        = 1728.744343059253
        self.allisonCarbon_erebus      = 1439.8801634549727
        self.allisonCarbon_vesuvius    = 821.9428189121331
        self.allisonCarbon_etna        = 1039.5129481979354
        self.allisonCarbon_stromboli   = 1472.5760950857439
        self.liuCarbon                 = 2246.2067748765

        self.shishkinaWater            = 395.4721026237296
        self.dixonWater                = 431.1140172567279
        self.iaconomarzianoWater       = 459.7132843026062
        self.mooreWater                = 366.7939178950552
        self.liuWater                  = 374.4357814145504

        self.mixed_dict =   {"ShishkinaIdealMixing": self.shishkinaMixed,
                             "Dixon"               : self.dixonMixed,
                             "IaconoMarziano"      : self.iaconomarzianoMixed,
                             "Liu"                 : self.liuMixed,
                             "MagmaSat"            : self.magmasat
                             }

        self.carbon_dict = {"ShishkinaCarbon"        : self.shishkinaCarbon,
                           "DixonCarbon"             : self.dixonCarbon,
                           "IaconoMarzianoCarbon"    : self.iaconomarzianoCarbon,
                           "AllisonCarbon"           : self.allisonCarbon,
                           "AllisonCarbon_sunset"    : self.allisonCarbon_sunset,
                           'AllisonCarbon_sfvf'      : self.allisonCarbon_sfvf,
                           'AllisonCarbon_erebus'    : self.allisonCarbon_erebus,
                           'AllisonCarbon_vesuvius'  : self.allisonCarbon_vesuvius,
                           'AllisonCarbon_etna'      : self.allisonCarbon_etna,
                           'AllisonCarbon_stromboli' : self.allisonCarbon_stromboli,
                           'LiuCarbon'               : self.liuCarbon
                           }

        self.water_dict = {"ShishkinaWater"       : self.shishkinaWater,
                          "DixonWater"            : self.dixonWater,
                          "IaconoMarzianoWater"   : self.iaconomarzianoWater,
                          "MooreWater"            : self.mooreWater,
                          "LiuWater"              : self.liuWater
                          }

    def test_calculate_single_wtpt_mixed(self):
        for model in self.mixed_dict.keys():
            calcd_result = v.calculate_saturation_pressure(self.sample_wtpt, temperature=self.temperature, model=model, verbose=False).result
            known_result = self.mixed_dict[model]
            self.assertAlmostEqual(calcd_result, known_result, places=4)

    def test_calculate_batch_wtpt_mixed(self):
        for model in self.mixed_dict.keys():
            batch_result = self.batch_wtpt.calculate_saturation_pressure(temperature=self.temperature, model=model, verbose=True)
            calcd_result = batch_result['SaturationP_bars_VESIcal'].loc['test_samp']
            known_result = self.mixed_dict[model]
            #print("Pi* wtpt = " + str(batch_result['PiStar_VESIcal'].loc['test_samp']))
            self.assertAlmostEqual(calcd_result, known_result, places=4)

    def test_calculate_single_wtpt_carbon(self):
        for model in self.carbon_dict.keys():
            calcd_result = v.calculate_saturation_pressure(self.sample_wtpt, temperature=self.temperature, model=model, verbose=False).result
            known_result = self.carbon_dict[model]
            self.assertAlmostEqual(calcd_result, known_result, places=4)

    def test_calculate_batch_wtpt_carbon(self):
        for model in self.carbon_dict.keys():
            batch_result = self.batch_wtpt.calculate_saturation_pressure(temperature=self.temperature, model=model, verbose=False)
            calcd_result = batch_result['SaturationP_bars_VESIcal'].loc['test_samp']
            known_result = self.carbon_dict[model]
            self.assertAlmostEqual(calcd_result, known_result, places=4)

    def test_calculate_single_wtpt_water(self):
        for model in self.water_dict.keys():
            calcd_result = v.calculate_saturation_pressure(self.sample_wtpt, temperature=self.temperature, model=model, verbose=False).result
            known_result = self.water_dict[model]
            self.assertAlmostEqual(calcd_result, known_result, places=4)

    def test_calculate_batch_wtpt_water(self):
        for model in self.water_dict.keys():
            batch_result = self.batch_wtpt.calculate_saturation_pressure(temperature=self.temperature, model=model, verbose=False)
            calcd_result = batch_result['SaturationP_bars_VESIcal'].loc['test_samp']
            known_result = self.water_dict[model]
            self.assertAlmostEqual(calcd_result, known_result, places=4)

    def test_calculation_single_molox_mixed(self):
        for model in self.mixed_dict.keys():
            calcd_result = v.calculate_saturation_pressure(self.sample_molox, temperature=self.temperature, model=model, verbose=False).result
            known_result = self.mixed_dict[model]
            self.assertAlmostEqual(calcd_result, known_result, places=4)

    def test_calculation_batch_molox_mixed(self):
        for model in self.mixed_dict.keys():
            batch_result = self.batch_molox.calculate_saturation_pressure(temperature=self.temperature, model=model, verbose=True)
            calcd_result = batch_result['SaturationP_bars_VESIcal'].loc['test_samp']
            known_result = self.mixed_dict[model]
            #print("Pi* molox = " + str(batch_result['PiStar_VESIcal'].loc['test_samp']))
            self.assertAlmostEqual(calcd_result, known_result, places=4)

    def test_calculate_single_molox_carbon(self):
        for model in self.carbon_dict.keys():
            calcd_result = v.calculate_saturation_pressure(self.sample_molox, temperature=self.temperature, model=model, verbose=False).result
            known_result = self.carbon_dict[model]
            self.assertAlmostEqual(calcd_result, known_result, places=4)

    def test_calculate_batch_molox_carbon(self):
        for model in self.carbon_dict.keys():
            batch_result = self.batch_molox.calculate_saturation_pressure(temperature=self.temperature, model=model, verbose=False)
            calcd_result = batch_result['SaturationP_bars_VESIcal'].loc['test_samp']
            known_result = self.carbon_dict[model]
            self.assertAlmostEqual(calcd_result, known_result, places=4)

    def test_calculate_single_molox_water(self):
        for model in self.water_dict.keys():
            calcd_result = v.calculate_saturation_pressure(self.sample_molox, temperature=self.temperature, model=model, verbose=False).result
            known_result = self.water_dict[model]
            self.assertAlmostEqual(calcd_result, known_result, places=4)

    def test_calculate_batch_molox_water(self):
        for model in self.water_dict.keys():
            batch_result = self.batch_molox.calculate_saturation_pressure(temperature=self.temperature, model=model, verbose=False)
            calcd_result = batch_result['SaturationP_bars_VESIcal'].loc['test_samp']
            known_result = self.water_dict[model]
            self.assertAlmostEqual(calcd_result, known_result, places=4)

class TestEquilibriumFluidComp(unittest.TestCase):
    def setUp(self):
        # Sample with units as wtpt_oxides
        self.majors_wtpt = {'SiO2':    47.95,
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
                         'H2O':     2.0,
                         'CO2':     0.1}

        # Set conditions of calculation
        self.temperature = 1000
        self.pressure = 500

        # create Sample object and set default units to wtpt_oxides
        self.sample_wtpt = v.Sample(self.majors_wtpt)
        self.sample_wtpt.set_default_units("wtpt_oxides")

        # create sample and set default units to mol_oxides
        self.sample_molox = v.Sample(self.majors_wtpt)
        self.sample_molox.set_default_units("mol_oxides")

        # equilibrium fluid comps calculated with VESIcal
        self.shishkinaMixed =      {'H2O': 0.7158408854774682, 'CO2': 0.28415911452253184}
        self.dixonMixed =          {'H2O': 0.7750859842655139, 'CO2': 0.2249140157344861}
        self.iaconomarzianoMixed = {'H2O': 0.8103088488568084, 'CO2': 0.18969115114319157}
        self.liuMixed =            {'H2O': 0.7066707740811572, 'CO2': 0.2933292259188428}
        self.magmasat =            {'CO2': 0.219354223233457, 'H2O': 0.780645776766543}

        self.shishkinaCarbon           = 1.0
        self.dixonCarbon               = 1.0
        self.iaconomarzianoCarbon      = 1.0
        self.allisonCarbon             = 1.0
        self.allisonCarbon_sunset      = 1.0
        self.allisonCarbon_sfvf        = 1.0
        self.allisonCarbon_erebus      = 1.0
        self.allisonCarbon_vesuvius    = 1.0
        self.allisonCarbon_etna        = 1.0
        self.allisonCarbon_stromboli   = 1.0
        self.liuCarbon                 = 0

        self.shishkinaWater            = 0.0
        self.dixonWater                = 0.0
        self.iaconomarzianoWater       = 0.0
        self.mooreWater                = 0.7002733835379757
        self.liuWater                  = 0.758191771357082

        self.mixed_dict =   {"ShishkinaIdealMixing": self.shishkinaMixed,
                             "Dixon"               : self.dixonMixed,
                             "IaconoMarziano"      : self.iaconomarzianoMixed,
                             "Liu"                 : self.liuMixed,
                             "MagmaSat"            : self.magmasat
                             }

        self.carbon_dict = {"ShishkinaCarbon"        : self.shishkinaCarbon,
                           "DixonCarbon"             : self.dixonCarbon,
                           "IaconoMarzianoCarbon"    : self.iaconomarzianoCarbon,
                           "AllisonCarbon"           : self.allisonCarbon,
                           "AllisonCarbon_sunset"    : self.allisonCarbon_sunset,
                           'AllisonCarbon_sfvf'      : self.allisonCarbon_sfvf,
                           'AllisonCarbon_erebus'    : self.allisonCarbon_erebus,
                           'AllisonCarbon_vesuvius'  : self.allisonCarbon_vesuvius,
                           'AllisonCarbon_etna'      : self.allisonCarbon_etna,
                           'AllisonCarbon_stromboli' : self.allisonCarbon_stromboli,
                           'LiuCarbon'               : self.liuCarbon
                           }

        self.water_dict = {"ShishkinaWater"       : self.shishkinaWater,
                          "DixonWater"            : self.dixonWater,
                          "IaconoMarzianoWater"   : self.iaconomarzianoWater,
                          "MooreWater"            : self.mooreWater,
                          "LiuWater"              : self.liuWater
                          }

    def test_calculate_wtpt_mixed(self):
        for model in self.mixed_dict.keys():
            calcd_result = v.calculate_equilibrium_fluid_comp(self.sample_wtpt,
                                                              temperature=self.temperature,
                                                              pressure=self.pressure,
                                                              model=model,
                                                              verbose=False).result
            known_result = self.mixed_dict[model]
            for key in known_result.keys():
                self.assertAlmostEqual(calcd_result[key], known_result[key], places=4)

    def test_calculate_wtpt_carbon(self):
        for model in self.carbon_dict.keys():
            calcd_result = v.calculate_equilibrium_fluid_comp(self.sample_wtpt,
                                                              temperature=self.temperature,
                                                              pressure=self.pressure,
                                                              model=model,
                                                              verbose=False).result
            known_result = self.carbon_dict[model]
            self.assertAlmostEqual(calcd_result, known_result, places=4)

    def test_calculate_wtpt_water(self):
        for model in self.water_dict.keys():
            calcd_result = v.calculate_equilibrium_fluid_comp(self.sample_wtpt,
                                                              temperature=self.temperature,
                                                              pressure=self.pressure,
                                                              model=model,
                                                              verbose=False).result
            known_result = self.water_dict[model]
            self.assertAlmostEqual(calcd_result, known_result, places=4)

    def test_calculation_molox_mixed(self):
        for model in self.mixed_dict.keys():
            calcd_result = v.calculate_equilibrium_fluid_comp(self.sample_molox,
                                                              temperature=self.temperature,
                                                              pressure=self.pressure,
                                                              model=model,
                                                              verbose=False).result
            known_result = self.mixed_dict[model]
            for key in known_result.keys():
                self.assertAlmostEqual(calcd_result[key], known_result[key], places=4)

    def test_calculate_molox_carbon(self):
        for model in self.carbon_dict.keys():
            calcd_result = v.calculate_equilibrium_fluid_comp(self.sample_molox,
                                                              temperature=self.temperature,
                                                              pressure=self.pressure,
                                                              model=model,
                                                              verbose=False).result
            known_result = self.carbon_dict[model]
            self.assertAlmostEqual(calcd_result, known_result, places=4)

    def test_calculate_molox_water(self):
        for model in self.water_dict.keys():
            calcd_result = v.calculate_equilibrium_fluid_comp(self.sample_molox,
                                                              temperature=self.temperature,
                                                              pressure=self.pressure,
                                                              model=model,
                                                              verbose=False).result
            known_result = self.water_dict[model]
            self.assertAlmostEqual(calcd_result, known_result, places=4)

class TestLiquidViscosity(unittest.TestCase):
    def setUp(self):
        # Sample with units as wtpt_oxides
        # Default sample in Giordano excel spreadsheet
        self.majors_wtpt = {'SiO2':    62.40,
                            'TiO2':    0.55,
                            'Al2O3':   20.01,
                            'FeOT':    0.03,
                            'MnO':     0.02,
                            'MgO':     3.22,
                            'CaO':     9.08,
                            'Na2O':    3.52,
                            'K2O':     0.93,
                            'P2O5':    0.12,
                            'H2O':     2.0,
                            'F2O':     0.5}

        self.majors_Fe_Fe2O3_wtpt = {'SiO2':   62.40,
                                    'TiO2':    0.55,
                                    'Al2O3':   20.01,
                                    'FeO':     5.0,
                                    'Fe2O3':   0.5, 
                                    'MnO':     0.02,
                                    'MgO':     3.22,
                                    'CaO':     9.08,
                                    'Na2O':    3.52,
                                    'K2O':     0.93,
                                    'P2O5':    0.12,
                                    'H2O':     2.0,
                                    'F2O':     0.5}

        # Set conditions of calculation
        self.temperature = 1200

        # create Sample object and set default units to wtpt_oxides
        self.sample_wtpt = v.Sample(self.majors_wtpt)
        self.sample_wtpt.set_default_units("wtpt_oxides")
        self.sample_Fe_Fe2O3_wtpt = v.Sample(self.majors_Fe_Fe2O3_wtpt)
        self.sample_wtpt.set_default_units("wtpt_oxides")

        # create sample and set default units to mol_oxides
        self.sample_molox = v.Sample(self.majors_wtpt)
        self.sample_molox.set_default_units("mol_oxides")

        # BatchFile with test sample as defined above in wtpt_oxides
        try:
            self.batch_wtpt = v.BatchFile('BatchTest.xlsx', 
                                          sheet_name='giordano_test',
                                          units='wtpt_oxides')
        except:
            self.batch_wtpt = v.BatchFile('tests/BatchTest.xlsx',
                                            sheet_name='giordano_test', 
                                            units='wtpt_oxides')
        self.batch_wtpt.set_default_units("wtpt_oxides")

        # BatchFile with test sample as defined above in mol_oxides
        try:
            self.batch_molox = v.BatchFile('BatchTest.xlsx',
                                            sheet_name='giordano_test')
        except:
            self.batch_molox = v.BatchFile('tests/BatchTest.xlsx',
                                            sheet_name='giordano_test')
        self.batch_molox.set_default_units("mol_oxides")
        
        # viscosities calculated with Giordano excel spreadsheet
        self.giordano_default_viscosity = 2.1733
        self.batch_test_samp_viscosity = 0.7721
        self.batch_test_giordano_spreadsheet_default_comp = 2.1733

        # Sample majors_Fe_Fe2O3_wtpt normalized w/Giordano spreadsheet
        # This is in terms of mol% oxides
        self.giordano_normalized = {'SiO2':    59.003,
                                    'TiO2':    0.391,
                                    'Al2O3':   11.150,
                                    'FeO':     4.310,
                                    'MnO':     0.016,
                                    'MgO':     4.539,
                                    'CaO':     9.199,
                                    'Na2O':    3.227,
                                    'K2O':     0.561,
                                    'P2O5':    0.048,
                                    'H2O':     6.809,
                                    'F2O':     0.748}

    def test_calculate_single_wtpt(self):
        calcd_result = v.calculate_liquid_viscosity(self.sample_wtpt,
                                            temperature=self.temperature).result
        known_result = self.giordano_default_viscosity
        self.assertAlmostEqual(calcd_result, known_result, places=3)

    def test_calculate_batch_wtpt_1(self):
        batch_result = self.batch_wtpt.calculate_liquid_viscosity(
                                                temperature=self.temperature)
        calcd_result_1 = batch_result['Viscosity_liq_VESIcal'].loc['test_samp']
        known_result_1 = self.batch_test_samp_viscosity
        self.assertAlmostEqual(calcd_result_1, known_result_1, places=3)

    def test_calculate_batch_wtpt_2(self):
        batch_result = self.batch_wtpt.calculate_liquid_viscosity(
                                                temperature=self.temperature)
        calcd_result_2 = batch_result['Viscosity_liq_VESIcal'].loc[
                                            'giordano_spreadsheet_default_comp']
        known_result_2 = self.batch_test_giordano_spreadsheet_default_comp
        self.assertAlmostEqual(calcd_result_2, known_result_2, places=3)


    def test_calculate_single_molox(self):
        calcd_result = v.calculate_liquid_viscosity(self.sample_molox,
                                            temperature=self.temperature).result
        known_result = self.giordano_default_viscosity
        self.assertAlmostEqual(calcd_result, known_result, places=3)

    def test_calculate_batch_molox_1(self):
        batch_result = self.batch_molox.calculate_liquid_viscosity(
                                                temperature=self.temperature)
        calcd_result_1 = batch_result['Viscosity_liq_VESIcal'].loc['test_samp']
        known_result_1 = self.batch_test_samp_viscosity
        self.assertAlmostEqual(calcd_result_1, known_result_1, places=3)

    def test_calculate_batch_molox_2(self):
        batch_result = self.batch_molox.calculate_liquid_viscosity(
                                                temperature=self.temperature)
        calcd_result_2 = batch_result['Viscosity_liq_VESIcal'].loc[
                                            'giordano_spreadsheet_default_comp']
        known_result_2 = self.batch_test_giordano_spreadsheet_default_comp
        self.assertAlmostEqual(calcd_result_2, known_result_2, places=3)

    def test_normalize_giordano(self):
        known_result = self.giordano_normalized
        calcd_result = {}
        calculation = v.thermo.giordano._normalize_Giordano(
                                                      self.sample_Fe_Fe2O3_wtpt)
        for key, val in calculation.items():
            calcd_result[key] = round(val,3)
        self.assertDictEqual(calcd_result, known_result)

class TestLiquidDensity(unittest.TestCase):
    def setUp(self):
        # Sample with units as wtpt_oxides
        self.majors_wtpt = {'SiO2':    47.95,
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
                         'H2O':     2.0,
                         'CO2':     0.1}

        # Set conditions of calculation
        self.temperature = 1000
        self.pressure = 1000

        # create Sample object and set default units to wtpt_oxides
        self.sample_wtpt = v.Sample(self.majors_wtpt)
        self.sample_wtpt.set_default_units("wtpt_oxides")

        # create sample and set default units to mol_oxides
        self.sample_molox = v.Sample(self.majors_wtpt)
        self.sample_molox.set_default_units("mol_oxides")

        # BatchFile with test sample as defined above in wtpt_oxides
        try:
            self.batch_wtpt = v.BatchFile('BatchTest.xlsx', units='wtpt_oxides')
        except:
            self.batch_wtpt = v.BatchFile('tests/BatchTest.xlsx', 
                                            units='wtpt_oxides')
        self.batch_wtpt.set_default_units("wtpt_oxides")

        # BatchFile with test sample as defined above in mol_oxides
        try:
            self.batch_molox = v.BatchFile('BatchTest.xlsx')
        except:
            self.batch_molox = v.BatchFile('tests/BatchTest.xlsx')
        self.batch_molox.set_default_units("mol_oxides")
        
        # densities calculated with DensityX
        self.densityx = 2620.832

    def test_calculate_single_wtpt(self):
        calcd_result = v.calculate_liquid_density(self.sample_wtpt,
                                                  temperature=self.temperature, 
                                                  pressure=self.pressure).result
        known_result = self.densityx
        self.assertAlmostEqual(calcd_result, known_result, places=4)

    def test_calculate_batch_wtpt(self):
        batch_result = self.batch_wtpt.calculate_liquid_density(
                                            temperature=self.temperature,
                                            pressure=self.pressure)
        calcd_result = batch_result['Density_liq_VESIcal'].loc['test_samp']
        known_result = self.densityx
        self.assertAlmostEqual(calcd_result, known_result, places=4)

    def test_calculate_single_molox(self):
        calcd_result = v.calculate_liquid_density(self.sample_molox,
                                                  temperature=self.temperature, 
                                                  pressure=self.pressure).result
        known_result = self.densityx
        self.assertAlmostEqual(calcd_result, known_result, places=4)

    def test_calculate_batch_molox(self):
        batch_result = self.batch_molox.calculate_liquid_density(
                                            temperature=self.temperature,
                                            pressure=self.pressure)
        calcd_result = batch_result['Density_liq_VESIcal'].loc['test_samp']
        known_result = self.densityx
        self.assertAlmostEqual(calcd_result, known_result, places=4)

if __name__ == '__main__':
    unittest.main()
