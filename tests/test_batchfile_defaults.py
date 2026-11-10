import unittest

import numpy as np
import pandas as pd

import VESIcal as v


# Fixed molecular masses and stoichiometry for an independent output oracle.
CHEMISTRY = {
    'SiO2': (60.083, 'Si', 1, 2),
    'MgO': (40.304, 'Mg', 1, 1),
    'Na2O': (61.979, 'Na', 2, 1),
    'Al2O3': (101.961, 'Al', 2, 3),
    'H2O': (18.015, 'H', 2, 1),
    'CO2': (44.01, 'C', 1, 2),
}
UNITS = ('wtpt_oxides', 'mol_oxides', 'mol_cations')
NORMALIZATIONS = ('none', 'standard', 'fixedvolatiles',
                  'additionalvolatiles')


def expected_composition(row, units, normalization='none',
                         exclude_volatiles=False):
    weights = {oxide: row[oxide] for oxide in CHEMISTRY
               if not (exclude_volatiles and oxide in ('H2O', 'CO2'))}
    if units == 'wtpt_oxides':
        values = weights
    else:
        moles = {oxide: weight / CHEMISTRY[oxide][0]
                 for oxide, weight in weights.items()}
        if units == 'mol_oxides':
            values = {oxide: amount / sum(moles.values())
                      for oxide, amount in moles.items()}
        else:
            cations = {CHEMISTRY[oxide][1]: amount * CHEMISTRY[oxide][2]
                       for oxide, amount in moles.items()}
            if units == 'mol_singleO':
                total = sum(amount * CHEMISTRY[oxide][3]
                            for oxide, amount in moles.items())
            else:
                total = sum(cations.values())
            values = {cation: amount / total
                      for cation, amount in cations.items()}

    if normalization == 'none':
        return values
    if normalization == 'standard':
        total = 100.0 if units == 'wtpt_oxides' else 1.0
        return {key: amount * total / sum(values.values())
                for key, amount in values.items()}

    # Volatile-preserving controls use the documented weight-percent modes.
    if units != 'wtpt_oxides':
        raise ValueError('These controls require weight-percent output')
    volatiles = {key: amount for key, amount in values.items()
                 if key in ('H2O', 'CO2')}
    dry = {key: amount for key, amount in values.items()
           if key not in volatiles}
    total = 100.0
    if normalization == 'fixedvolatiles':
        total -= sum(volatiles.values())
    result = {key: amount * total / sum(dry.values())
              for key, amount in dry.items()}
    return dict(result, **volatiles)


class TestBatchFileOutputDefaults(unittest.TestCase):
    def setUp(self):
        self.frame = pd.DataFrame({
            'SiO2': [45.0, 40.0], 'MgO': [15.0, 20.0],
            'Na2O': [10.0, 8.0], 'Al2O3': [20.0, 12.0],
            'H2O': [2.0, 3.0], 'CO2': [1.0, 0.5],
            'Notes': ['first sample', 'second sample'], 'Order': [7, 2],
        }, index=['sample-3', 'sample-1'])
        self.original = self.frame.copy(deep=True)

    def make_batch(self, units='wtpt_oxides', normalization='none',
                   setters=False):
        if setters:
            batch = v.BatchFile(None, dataframe=self.frame, label=None)
            batch.set_default_units(units)
            batch.set_default_normalization(normalization)
        else:
            batch = v.BatchFile(None, dataframe=self.frame, label=None,
                                default_units=units,
                                default_normalization=normalization)
        return batch

    def assert_output(self, batch, expected_units, expected_normalization='none',
                      **kwargs):
        stored = batch.data.copy(deep=True)
        expected = pd.DataFrame([
            expected_composition(row, expected_units, expected_normalization)
            for _, row in self.frame.iterrows()
        ], index=self.frame.index)
        for getter in ('get_composition', 'get_data'):
            with self.subTest(getter=getter):
                result = getattr(batch, getter)(asBatchFile=False, **kwargs)
                self.assertIsInstance(result, pd.DataFrame)
                self.assertTrue(result.columns.is_unique)
                self.assertEqual(result.index.tolist(), self.frame.index.tolist())
                columns = set(expected.columns)
                if getter == 'get_data':
                    columns.update(('Notes', 'Order'))
                    pd.testing.assert_frame_equal(result[['Notes', 'Order']],
                                                  self.frame[['Notes', 'Order']])
                self.assertEqual(set(result.columns), columns)
                np.testing.assert_allclose(result[expected.columns], expected,
                                           rtol=0, atol=1e-9)
                pd.testing.assert_frame_equal(batch.data, stored)
                pd.testing.assert_frame_equal(self.frame, self.original)

    def test_constructor_units_match_independent_oracle(self):
        for units in UNITS:
            with self.subTest(units=units):
                self.assert_output(self.make_batch(units=units), units)

    def test_setter_units_match_independent_oracle(self):
        for units in UNITS:
            with self.subTest(units=units):
                batch = self.make_batch(units=units, setters=True)
                self.assert_output(batch, units)

    def test_explicit_units_override_constructor_defaults(self):
        for default in UNITS:
            for units in UNITS + ('mol_singleO',):
                with self.subTest(default=default, units=units):
                    self.assert_output(self.make_batch(units=default), units,
                                       units=units)

    def test_explicit_units_override_setter_defaults(self):
        for default in UNITS:
            for units in UNITS + ('mol_singleO',):
                with self.subTest(default=default, units=units):
                    batch = self.make_batch(units=default, setters=True)
                    self.assert_output(batch, units, units=units)

    def test_explicit_none_overrides_constructor_normalization(self):
        for default in NORMALIZATIONS[1:]:
            with self.subTest(default=default):
                batch = self.make_batch(normalization=default)
                self.assert_output(batch, 'wtpt_oxides', normalization='none')

    def test_explicit_none_overrides_setter_normalization(self):
        for default in NORMALIZATIONS[1:]:
            with self.subTest(default=default):
                batch = self.make_batch(normalization=default, setters=True)
                self.assert_output(batch, 'wtpt_oxides', normalization='none')

    def test_constructor_weight_percent_normalizations(self):
        for normalization in NORMALIZATIONS:
            with self.subTest(normalization=normalization):
                batch = self.make_batch(normalization=normalization)
                self.assert_output(batch, 'wtpt_oxides', normalization)

    def test_setter_weight_percent_normalizations(self):
        for normalization in NORMALIZATIONS:
            with self.subTest(normalization=normalization):
                batch = self.make_batch(normalization=normalization,
                                        setters=True)
                self.assert_output(batch, 'wtpt_oxides', normalization)

    def test_explicit_normalization_starts_from_stored_weights(self):
        for default in NORMALIZATIONS:
            for normalization in NORMALIZATIONS:
                with self.subTest(default=default, normalization=normalization):
                    batch = self.make_batch(normalization=default)
                    self.assert_output(batch, 'wtpt_oxides', normalization,
                                       normalization=normalization)

    def test_standard_molar_defaults(self):
        for units in UNITS[1:]:
            for setters in (False, True):
                with self.subTest(units=units, setters=setters):
                    batch = self.make_batch(units, 'standard', setters=setters)
                    self.assert_output(batch, units, 'standard')

    def test_species_ignores_unrequested_defaults(self):
        expected = {
            'SiO2': self.frame['SiO2'].to_numpy(),
            'H2O': self.frame['H2O'].to_numpy(),
            'Na': np.array([
                expected_composition(row, 'mol_cations')['Na']
                for _, row in self.frame.iterrows()
            ]),
        }
        for units in UNITS:
            batch = self.make_batch(units, 'standard')
            for species, values in expected.items():
                with self.subTest(units=units, species=species):
                    result = batch.get_composition(species=species,
                                                   asBatchFile=False)
                    np.testing.assert_allclose(result[species], values,
                                               rtol=0, atol=1e-9)

    def test_exclude_volatiles_uses_stored_weights(self):
        for units in UNITS:
            for normalization in ('none', 'standard'):
                with self.subTest(units=units, normalization=normalization):
                    batch = self.make_batch(units, 'additionalvolatiles')
                    stored = batch.data.copy(deep=True)
                    expected = pd.DataFrame([
                        expected_composition(row, units, normalization,
                                             exclude_volatiles=True)
                        for _, row in self.frame.iterrows()
                    ], index=self.frame.index)
                    result = batch.get_composition(
                        normalization=normalization, exclude_volatiles=True,
                        asBatchFile=False)
                    self.assertEqual(set(result.columns), set(expected.columns))
                    np.testing.assert_allclose(result[expected.columns], expected,
                                               rtol=0, atol=1e-9)
                    pd.testing.assert_frame_equal(batch.data, stored)
                    pd.testing.assert_frame_equal(self.frame, self.original)

    def test_equal_oxide_moles_are_not_converted_twice(self):
        frame = pd.DataFrame({'SiO2': [59.851375178060906],
                              'MgO': [40.148624821939094]}, index=['sample'])
        for setters in (False, True):
            with self.subTest(setters=setters):
                batch = v.BatchFile(None, dataframe=frame, label=None,
                                    default_units='mol_oxides')
                if setters:
                    batch.set_default_units('wtpt_oxides')
                    batch.set_default_units('mol_oxides')
                result = batch.get_composition(asBatchFile=False)
                np.testing.assert_allclose(result[['SiO2', 'MgO']], [[0.5, 0.5]],
                                           rtol=0, atol=1e-9)
