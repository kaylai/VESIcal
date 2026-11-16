from copy import deepcopy
import unittest

import pandas as pd

import VESIcal as v
from VESIcal import core


# CIAAW abridged atomic weights, at the library's three-decimal precision:
# Co = 58.933, O = 15.999, Mg = 24.305, Si = 28.085.
# https://www.ciaaw.org/abridged-atomic-weights.htm
# These literal wt% mixtures contain equal moles of the two named oxides.
CO_MG_WT = {"CoO": 65.02481863306606, "MgO": 34.975181366933946}
CO_SI_WT = {"CoO": 55.4990186275599, "SiO2": 44.5009813724401}
MG_SI_WT = {"MgO": 40.148624821939094, "SiO2": 59.851375178060906}
OVERRIDE_WT = {"CoO": 52.19773703062362, "MgO": 47.80226296937638}


class TestCobaltMolarMasses(unittest.TestCase):
    def make_sample(self, composition, units="wtpt_oxides"):
        before = deepcopy(composition)
        try:
            return v.Sample(composition, units=units)
        finally:
            if isinstance(composition, pd.Series):
                pd.testing.assert_series_equal(composition, before)
            else:
                self.assertEqual(composition, before)
                self.assertEqual(list(composition), list(before))

    def assert_series_result(self, sample, expected, **kwargs):
        before = sample.get_composition()
        try:
            actual = sample.get_composition(**kwargs)
        finally:
            pd.testing.assert_series_equal(sample.get_composition(), before)
        self.assertIs(type(actual), pd.Series)
        pd.testing.assert_series_equal(
            actual, pd.Series(expected, dtype="float64"),
            check_exact=False, rtol=1e-12, atol=1e-12,
        )

    def assert_formula_weight(self, sample, expected):
        before = sample.get_composition()
        try:
            actual = sample.get_formulaweight()
        finally:
            pd.testing.assert_series_equal(sample.get_composition(), before)
        self.assertAlmostEqual(actual, expected, places=12)

    def test_weight_percent_to_mole_oxides(self):
        sample = self.make_sample(CO_MG_WT)
        self.assert_series_result(sample, {"CoO": 0.5, "MgO": 0.5}, units="mol_oxides")

    def test_weight_percent_to_mole_cations(self):
        composition = pd.Series(CO_MG_WT, dtype="float64", name="Co-Mg mixture")
        sample = self.make_sample(composition)
        self.assert_series_result(sample, {"Co": 0.5, "Mg": 0.5}, units="mol_cations")

    def test_weight_percent_to_single_oxygen(self):
        sample = self.make_sample(CO_SI_WT)
        # One mole each of CoO and SiO2 contains three moles of oxygen atoms.
        self.assert_series_result(sample, {"Co": 1 / 3, "Si": 1 / 3}, units="mol_singleO")

    def test_mole_oxide_input_to_weight_percent(self):
        composition = pd.Series({"CoO": 0.5, "MgO": 0.5}, dtype="float64", name="molar oxides")
        sample = self.make_sample(composition, units="mol_oxides")
        self.assert_series_result(sample, CO_MG_WT)

    def test_mole_cation_input_to_weight_percent(self):
        sample = self.make_sample({"Co": 0.5, "Mg": 0.5}, units="mol_cations")
        self.assert_series_result(sample, CO_MG_WT)

    def test_cobalt_oxide_formula_weight(self):
        self.assert_formula_weight(self.make_sample({"CoO": 100.0}), 74.932)

    def test_non_cobalt_mole_fraction_control(self):
        composition = pd.Series(MG_SI_WT, dtype="float64", name="Mg-Si mixture")
        sample = self.make_sample(composition)
        self.assert_series_result(sample, {"MgO": 0.5, "SiO2": 0.5}, units="mol_oxides")

    def test_magnesium_formula_weight_control(self):
        self.assert_formula_weight(self.make_sample({"MgO": 100.0}), 40.304)

    def test_explicit_oxide_mass_override_control(self):
        sample = self.make_sample(OVERRIDE_WT)
        masses = {"CoO": 44.01}
        before_masses = masses.copy()
        before_defaults = core.oxideMass.copy()
        try:
            self.assert_series_result(
                sample, {"CoO": 0.5, "MgO": 0.5}, units="mol_oxides", oxide_masses=masses,
            )
        finally:
            self.assertEqual(masses, before_masses)
            self.assertEqual(core.oxideMass, before_defaults)
