import unittest
import VESIcal as v
import numpy as np
import pandas as pd

class TestCreateSample(unittest.TestCase):
    def setUp(self):
        self.majors = pd.Series({'SiO2':    47.95,
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
                    })

        self.majors_normed = pd.Series({'SiO2':    47.94,
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
                    })

        # Mol fractions calculated externally
        self.majors_moloxides = pd.Series({ 'SiO2':  51.43,
                                            'TiO2':  1.35,
                                            'Al2O3': 10.95,
                                            'FeO':   9.19,
                                            'Fe2O3': 0.04,
                                            'MnO':   0.09,
                                            'MgO':   9.21,
                                            'CaO':   12.56,
                                            'Na2O':  3.59,
                                            'K2O':   1.36,
                                            'P2O5': 0.23})
        self.majors_moloxides = self.majors_moloxides/100

        self.majors_molcations = pd.Series({'Si':   0.4427,
                                            'Ti':   0.0116,
                                            'Al':   0.1885,
                                            'Fe':   0.0791,
                                            'Fe3':  0.0007,
                                            'Mn':   0.0008,
                                            'Mg':   0.0793,
                                            'Ca':   0.1081,
                                            'Na':   0.0618,
                                            'K':    0.0234,
                                            'Cr':   0.0040})

        self.majorsv = pd.Series({'SiO2':   47.95,
                                  'TiO2':   1.67,
                                  'Al2O3':  17.32,
                                  'FeO':    10.24,
                                  'Fe2O3':  0.1,
                                  'MgO':    5.76,
                                  'CaO':    10.93,
                                  'Na2O':   3.45,
                                  'K2O':    1.99,
                                  'P2O5':   0.51,
                                  'MnO':    0.1,
                                  'CO2':    0.08,
                                  'H2O':    4.0
                        })

        self.majorsv_normed = pd.Series({'SiO2':   46.06,
                                         'TiO2':   1.60,
                                         'Al2O3':  16.64,
                                         'FeO':    9.84,
                                         'Fe2O3':  0.1,
                                         'MgO':    5.53,
                                         'CaO':    10.50,
                                         'Na2O':   3.31,
                                         'K2O':    1.91,
                                         'P2O5':   0.49,
                                         'MnO':    0.1,
                                         'CO2':    0.08,
                                         'H2O':    3.84
                        })

        self.majorsv_moloxides = pd.Series({'SiO2':  44.95,
                                            'TiO2':  1.178,
                                            'Al2O3': 9.568,
                                            'FeO':   8.028,
                                            'Fe2O3': 0.035,
                                            'MnO':   0.079,
                                            'MgO':   8.050,
                                            'CaO':   10.978,
                                            'Na2O':  3.135,
                                            'K2O':   1.190,
                                            'P2O5':  0.202,
                                            'H2O':   12.503,
                                            'CO2':   0.102})
        self.majorsv_moloxides = self.majorsv_moloxides/100

        self.majorsv_molcations = pd.Series({'Si':   0.35496,
                                             'Ti':   0.00930,
                                             'Al':   0.15111,
                                             'Fe':   0.06340,
                                             'Fe3':  0.00056,
                                             'Mn':   0.00063,
                                             'Mg':   0.06357,
                                             'Ca':   0.08669,
                                             'Na':   0.04952,
                                             'K':    0.01879,
                                             'P':    0.00320,
                                             'H':    0.19747,
                                             'C':    0.00081})

        self.sample = v.Sample(self.majorsv)

    def test_createSample(self):
        for ox in self.majors.index:
            self.assertEqual(self.sample._composition[ox],self.majors[ox])

    def test_setdefault_noargs(self):
        self.assertEqual(self.sample.default_normalization,'none')
        self.assertEqual(self.sample.default_units,'wtpt_oxides')

    def test_setdefaults_none_wtpt(self):
        sample = v.Sample(self.majorsv, default_normalization='none',default_units='wtpt_oxides')
        self.assertEqual(sample.default_normalization,'none')
        self.assertEqual(sample.default_units,'wtpt_oxides')

    def test_setdefaults_standard_molox(self):
        sample = v.Sample(self.majorsv, default_normalization='standard',default_units='mol_oxides')
        self.assertEqual(sample.default_normalization,'standard')
        self.assertEqual(sample.default_units,'mol_oxides')

    def test_setdefaults_fixvol_molcat(self):
        sample = v.Sample(self.majorsv, default_normalization='fixedvolatiles',default_units='mol_cations')
        self.assertEqual(sample.default_normalization,'fixedvolatiles')
        self.assertEqual(sample.default_units,'mol_cations')

    def test_setdefaults_addvol_molSinO(self):
        sample = v.Sample(self.majorsv, default_normalization='additionalvolatiles',default_units='mol_singleO')
        self.assertEqual(sample.default_normalization,'additionalvolatiles')
        self.assertEqual(sample.default_units,'mol_singleO')

    def test_setdefaults_garbageNorm(self):
        with self.assertRaises(v.core.InputError):
            v.Sample(composition=self.majorsv,default_normalization='garbage')

    def test_setdefaults_garbageType(self):
        with self.assertRaises(v.core.InputError):
            v.Sample(composition=self.majorsv,default_units='garbage')

    def test_type_garbage(self):
        with self.assertRaises(v.core.InputError):
            v.Sample(composition=self.majorsv,units='garbage')

    def test_type_wtptoxides(self):
        sample = v.Sample(self.majors,units='wtpt_oxides')
        for ox in self.majors.index:
            self.assertEqual(self.sample._composition[ox],self.majors[ox])

    def test_type_moloxides(self):
        sample = v.Sample(self.majorsv_moloxides,units='mol_oxides')
        for ox in self.majorsv.index:
            self.assertEqual(np.round(sample._composition[ox],2),np.round(self.majorsv_normed[ox],2))

    def test_type_molcations(self):
        sample = v.Sample(self.majorsv_molcations,units='mol_cations')
        for ox in self.majorsv.index:
            self.assertEqual(np.round(sample._composition[ox],2),np.round(self.majorsv_normed[ox],2))


class TestGetComposition(unittest.TestCase):
    def setUp(self):
        self.majors = pd.Series({'SiO2':    47.95,
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
                    })

        self.majors_normed = pd.Series({'SiO2':    47.94,
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
                    })

        # Mol fractions calculated externally
        self.majors_moloxides = pd.Series({ 'SiO2':  51.434,
                                            'TiO2':  1.348,
                                            'Al2O3': 10.948,
                                            'FeO':   9.186,
                                            'Fe2O3': 0.004,
                                            'MnO':   0.091,
                                            'MgO':   9.211,
                                            'CaO':   12.562,
                                            'Na2O':  3.588,
                                            'K2O':   1.362,
                                            'P2O5': 0.232})
        self.majors_moloxides = self.majors_moloxides/100

        self.majors_molcations = pd.Series({'Si':   0.44275,
                                            'Ti':   0.01160,
                                            'Al':   0.18848,
                                            'Fe':   0.07908,
                                            'Fe3':  0.00069,
                                            'Mn':   0.00078,
                                            'Mg':   0.07929,
                                            'Ca':   0.10813,
                                            'Na':   0.06176,
                                            'K':    0.02344,
                                            'P':    0.00399})

        self.majors_molSingleO = pd.Series({ 'Si':   0.29276,
                                             'Ti':   0.00767,
                                             'Al':   0.12463,
                                             'Fe':   0.05229,
                                             'Fe3':  0.00046,
                                             'Mn':   0.00052,
                                             'Mg':   0.05243,
                                             'Ca':   0.07150,
                                             'Na':   0.04084,
                                             'K':    0.01550,
                                             'P':    0.00264,
                                            })

        self.majors_fw = 36.6928

        self.majorsv = pd.Series({'SiO2':   47.95,
                                  'TiO2':   1.67,
                                  'Al2O3':  17.32,
                                  'FeO':    10.24,
                                  'Fe2O3':  0.1,
                                  'MgO':    5.76,
                                  'CaO':    10.93,
                                  'Na2O':   3.45,
                                  'K2O':    1.99,
                                  'P2O5':   0.51,
                                  'MnO':    0.1,
                                  'CO2':    0.08,
                                  'H2O':    4.0
                        })

        self.majorsv_normed = pd.Series({'SiO2':   46.06,
                                         'TiO2':   1.60,
                                         'Al2O3':  16.64,
                                         'FeO':    9.84,
                                         'Fe2O3':  0.1,
                                         'MgO':    5.53,
                                         'CaO':    10.50,
                                         'Na2O':   3.31,
                                         'K2O':    1.91,
                                         'P2O5':   0.49,
                                         'MnO':    0.1,
                                         'CO2':    0.08,
                                         'H2O':    3.84
                        })

        self.majorsv_moloxides = pd.Series({'SiO2':  44.95,
                                            'TiO2':  1.178,
                                            'Al2O3': 9.568,
                                            'FeO':   8.028,
                                            'Fe2O3': 0.035,
                                            'MnO':   0.079,
                                            'MgO':   8.050,
                                            'CaO':   10.978,
                                            'Na2O':  3.135,
                                            'K2O':   1.190,
                                            'P2O5':  0.202,
                                            'H2O':   12.503,
                                            'CO2':   0.102})
        self.majorsv_moloxides = self.majorsv_moloxides/100

        self.majorsv_molcations = pd.Series({'Si':   0.35496,
                                             'Ti':   0.00930,
                                             'Al':   0.15111,
                                             'Fe':   0.06340,
                                             'Fe3':  0.00056,
                                             'Mn':   0.00063,
                                             'Mg':   0.06357,
                                             'Ca':   0.08669,
                                             'Na':   0.04952,
                                             'K':    0.01879,
                                             'P':    0.00320,
                                             'H':    0.19747,
                                             'C':    0.00081})

        self.majorsv_molSingleO = pd.Series({'Si':   0.27038,
                                             'Ti':   0.00708,
                                             'Al':   0.11510,
                                             'Fe':   0.04829,
                                             'Fe3':  0.00042,
                                             'Mn':   0.00048,
                                             'Mg':   0.04842,
                                             'Ca':   0.06604,
                                             'Na':   0.03772,
                                             'K':    0.01432,
                                             'P':    0.00243,
                                             'H':    0.15042,
                                             'C':    0.00062})

        self.majorsv_fw = 35.2704

        self.sample = v.Sample(self.majorsv)

    def test_default(self):
        composition = self.sample.get_composition()
        for ox in composition.index:
            self.assertEqual(composition[ox],self.majorsv[ox])

    def test_wtptoxides_none(self):
        composition = self.sample.get_composition(normalization='none')
        for ox in composition.index:
            self.assertEqual(composition[ox],self.majorsv[ox])

    def test_wtptoxides_none_exclV(self):
        composition = self.sample.get_composition(normalization='none',exclude_volatiles=True)
        for ox in composition.index:
            self.assertEqual(composition[ox],self.majors[ox])

    def test_wtptoxides_std(self):
        composition = self.sample.get_composition(normalization='standard')
        for ox in composition.index:
            self.assertEqual(np.round(composition[ox],2),np.round(self.majorsv_normed[ox],2))

    def test_wtptoxides_std_exclV(self):
        composition = self.sample.get_composition(normalization='standard',exclude_volatiles=True)
        for ox in composition.index:
            self.assertEqual(np.round(composition[ox],2),np.round(self.majors_normed[ox],2))

    def test_wtptoxides_fixedVolatiles(self):
        composition = self.sample.get_composition(normalization='fixedvolatiles')
        for ox in composition.index:
            if ox not in ['CO2','H2O']:
                self.assertEqual(np.round(composition[ox],2),np.round(self.majors_normed[ox]*(100-self.majorsv['CO2']-self.majorsv['H2O'])/100,2))
            else:
                self.assertEqual(np.round(composition[ox],2),np.round(self.majorsv[ox],2))

    def test_wtptoxides_fixedVolatiles_exclV(self):
        composition = self.sample.get_composition(normalization='fixedvolatiles',exclude_volatiles=True)
        for ox in composition.index:
            self.assertEqual(np.round(composition[ox],2),np.round(self.majors_normed[ox],2))

    def test_wtptoxides_additionalVolatiles(self):
        composition = self.sample.get_composition(normalization='additionalvolatiles')
        for ox in composition.index:
            if ox not in ['CO2','H2O']:
                self.assertEqual(np.round(composition[ox],2),np.round(self.majors_normed[ox],2))
            else:
                self.assertEqual(np.round(composition[ox],2),np.round(self.majorsv[ox],2))

    def test_wtptoxides_additionalVolatiles_exclV(self):
        composition = self.sample.get_composition(normalization='additionalvolatiles',exclude_volatiles=True)
        for ox in composition.index:
            self.assertEqual(np.round(composition[ox],2),np.round(self.majors_normed[ox],2))

    def test_moloxides(self):
        composition = self.sample.get_composition(units='mol_oxides')
        for ox in composition.index:
            self.assertEqual(np.round(composition[ox],3),np.round(self.majorsv_moloxides[ox],3))

    def test_moloxides_exclV(self):
        composition = self.sample.get_composition(units='mol_oxides',exclude_volatiles=True)
        for ox in composition.index:
            self.assertEqual(np.round(composition[ox],3),np.round(self.majors_moloxides[ox],3))

    def test_molcations(self):
        composition = self.sample.get_composition(units='mol_cations')
        for el in composition.index:
            self.assertEqual(np.round(composition[el],3),np.round(self.majorsv_molcations[el],3))

    def test_molcations_exclV(self):
        composition = self.sample.get_composition(units='mol_cations',exclude_volatiles=True)
        for el in composition.index:
            self.assertEqual(np.round(composition[el],3),np.round(self.majors_molcations[el],3))

    def test_molsingleO(self):
        composition = self.sample.get_composition(units='mol_singleO')
        for el in composition.index:
            self.assertEqual(np.round(composition[el],3),np.round(self.majorsv_molSingleO[el],3))

    def test_molsingleO_exclV(self):
        composition = self.sample.get_composition(units='mol_singleO',exclude_volatiles=True)
        for el in composition.index:
            self.assertEqual(np.round(composition[el],3),np.round(self.majors_molSingleO[el],3))

    def test_formulawt(self):
        fw = self.sample.get_formulaweight()
        self.assertEqual(np.round(self.majorsv_fw,2),np.round(fw,2))

    def test_formulawt_exclV(self):
        fw = self.sample.get_formulaweight(exclude_volatiles=True)
        self.assertEqual(np.round(self.majors_fw,2),np.round(fw,2))
