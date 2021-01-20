import unittest
import VESIcal as v
import numpy as np
import pandas as pd

class TestCreateSample(unittest.TestCase):
    def setUp(self):
        self.majors = { 'SiO2':47.95,
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
                    }
        self.majors = pd.Series(self.majors)

        self.majorsv = {'SiO2':47.95,
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
                        'H2O':4.0
                        }
        self.majorsv =pd.Series(self.majorsv)

        self.majors_comparison = self.majors.copy()

        self.sample = v.sample(self.majorsv)

    def test_createSample(self):
        for ox in self.majors.index:
            self.assertEqual(self.sample._composition[ox],self.majors_comparison[ox])

    def test_setdefault_noargs(self):
        self.assertEqual(self.sample.default_normalization,'none')
        self.assertEqual(self.sample.default_type,'wtpt_oxides')

    def test_setdefaults_none_wtpt(self):
        sample = v.sample(self.majorsv, default_normalization='none',default_type='wtpt_oxides')
        self.assertEqual(sample.default_normalization,'none')
        self.assertEqual(sample.default_type,'wtpt_oxides')

    def test_setdefaults_standard_molox(self):
        sample = v.sample(self.majorsv, default_normalization='standard',default_type='mol_oxides')
        self.assertEqual(sample.default_normalization,'standard')
        self.assertEqual(sample.default_type,'mol_oxides')

    def test_setdefaults_fixvol_molcat(self):
        sample = v.sample(self.majorsv, default_normalization='fixedvolatiles',default_type='mol_cations')
        self.assertEqual(sample.default_normalization,'fixedvolatiles')
        self.assertEqual(sample.default_type,'mol_cations')

    def test_setdefaults_addvol_molSinO(self):
        sample = v.sample(self.majorsv, default_normalization='additionalvolatiles',default_type='mol_singleO')
        self.assertEqual(sample.default_normalization,'additionalvolatiles')
        self.assertEqual(sample.default_type,'mol_singleO')

    def test_setdefaults_garbageNorm(self):
        with self.assertRaises(v.InputError):
            v.sample(composition=self.majorsv,default_normalization='garbage')
