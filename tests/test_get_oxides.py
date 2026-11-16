import unittest

import pandas as pd

import VESIcal as v


FULL_COMPOSITION = {
    "SiO2": 50.0,
    "TiO2": 1.0,
    "Al2O3": 15.0,
    "Cr2O3": 0.1,
    "FeO": 7.0,
    "Fe2O3": 1.0,
    "MnO": 0.2,
    "MgO": 8.0,
    "CaO": 10.0,
    "Na2O": 3.0,
    "K2O": 1.0,
    "P2O5": 0.2,
    "NiO": 0.05,
    "CoO": 0.02,
    "H2O": 2.0,
    "CO2": 0.1,
    "F2O": 0.01,
}


class TestGetOxides(unittest.TestCase):
    def assert_dict_filter(self, sample, expected):
        before = sample.copy()
        try:
            actual = v.get_oxides(sample)
        finally:
            self.assertEqual(sample, before)
            self.assertEqual(list(sample), list(before))
        self.assertIs(type(actual), dict)
        self.assertEqual(actual, expected)
        self.assertEqual(list(actual), list(expected))

    def assert_series_filter(self, sample, expected):
        before = sample.copy(deep=True)
        try:
            actual = v.get_oxides(sample)
        finally:
            pd.testing.assert_series_equal(sample, before)
        self.assertIs(type(actual), pd.Series)
        pd.testing.assert_series_equal(actual, expected)

    def test_partial_dict(self):
        sample = {"sample": "MORB", "MgO": 8.0, "FeOT": 10.5,
                  "SiO2": 50.0, "T_C": 1200}
        self.assert_dict_filter(sample, {"SiO2": 50.0, "MgO": 8.0})

    def test_partial_series(self):
        sample = pd.Series({"sample": "MORB", "MgO": 8.0, "FeOT": 10.5,
                            "SiO2": 50.0, "T_C": 1200},
                           dtype=object, name="partial sample")
        expected = pd.Series({"SiO2": 50.0, "MgO": 8.0}, dtype="float64")
        self.assert_series_filter(sample, expected)

    def test_complete_dict_compatibility(self):
        sample = dict(reversed(list(FULL_COMPOSITION.items())))
        sample.update({"sample": "complete", "T_C": 1200, "FeOT": 10.5})
        self.assert_dict_filter(sample, FULL_COMPOSITION)

    def test_complete_series_compatibility(self):
        composition = dict(reversed(list(FULL_COMPOSITION.items())))
        composition.update({"sample": "complete", "T_C": 1200, "FeOT": 10.5})
        sample = pd.Series(composition, dtype=object, name="complete sample")
        expected = pd.Series(FULL_COMPOSITION, dtype="float64")
        self.assert_series_filter(sample, expected)

    def test_metadata_only_dict(self):
        self.assert_dict_filter({"sample": "metadata", "T_C": 1200, "FeOT": 10.5}, {})

    def test_metadata_only_series(self):
        sample = pd.Series({"sample": "metadata", "T_C": 1200, "FeOT": 10.5},
                           dtype=object, name="metadata only")
        self.assert_series_filter(sample, pd.Series(dtype=object))
