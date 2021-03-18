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

        # solubilities calculated with VESIcal in wt% and translated externally to mol fraction
        self.shishkinaMixed_wtpt =      {'H2O_liq': 2.24019,
                                        'CO2_liq': 0.028022292} 
        self.shishkinaMixed_molox =     {'H2O_liq': 0.07415139,
                                        'CO2_liq': 0.000379788}

        self.dixonMixed_wtpt =          {'H2O_liq': 2.061908389,
                                        'CO2_liq': 0.035300022}
        self.dixonMixed_molox =         {'H2O_liq': 0.068648519,
                                        'CO2_liq': 0.000481216}

        self.iaconomarzianoMixed_wtpt = {'H2O_liq': 2.085306353,
                                        'CO2_liq': 0.049996964}
        self.iaconomarzianoMixed_molox ={'H2O_liq': 0.069359596,
                                        'CO2_liq': 0.000680901}

        self.liuMixed_wtpt =            {'H2O_liq': 2.293923234,
                                        'CO2_liq': 0.029485739}
        self.liuMixed_molox =           {'H2O_liq': 0.075793677,
                                        'CO2_liq': 0.000398905}

        self.magmasat_wtpt =            {'H2O_liq': 1.982045999,
                                        'CO2_liq': 0.044939669}
        self.magmasat_molox =           {'H2O_liq': 0.066156818,
                                        'CO2_liq': 0.000614178}

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
        self.mooreWater_wtpt                = 2.24337989
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

    def test_calculate_wtpt_mixed(self):
        for model in self.mixed_dict.keys():
            calcd_result = v.calculate_dissolved_volatiles(self.sample_wtpt, pressure=1000, temperature=1000, X_fluid=0.5, model=model, verbose=False).result
            known_result = self.mixed_dict[model]['wtpt_oxides']
            for k in calcd_result.keys():
                self.assertAlmostEqual(calcd_result[k], known_result[k], places=4)

    def test_calculate_wtpt_nonmixed(self):
        for model in self.nonmixed_dict.keys():
            calcd_result = v.calculate_dissolved_volatiles(self.sample_wtpt, pressure=1000, temperature=1000, X_fluid=0.5, model=model, verbose=False).result
            known_result = self.nonmixed_dict[model]['wtpt_oxides']
            self.assertAlmostEqual(calcd_result, known_result, places=4)

    def test_calculation_molox_mixed(self):
        for model in self.mixed_dict.keys():
            calcd_result = v.calculate_dissolved_volatiles(self.sample_molox, pressure=1000, temperature=1000, X_fluid=0.5, model=model, verbose=False).result
            known_result = self.mixed_dict[model]['mol_oxides']
            for k in calcd_result.keys():
                self.assertAlmostEqual(calcd_result[k], known_result[k], places=4)

    def test_calculation_molox_nonmixed(self):
        for model in self.nonmixed_dict.keys():
            calcd_result = v.calculate_dissolved_volatiles(self.sample_molox, pressure=1000, temperature=1000, X_fluid=0.5, model=model, verbose=False).result
            known_result = self.nonmixed_dict[model]['mol_oxides']
            self.assertAlmostEqual(calcd_result, known_result, places=4)

if __name__ == '__main__':
    unittest.main()




