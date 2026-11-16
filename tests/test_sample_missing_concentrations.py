"""Missing concentrations should behave like explicit zero concentrations."""

from copy import deepcopy
import unittest
import warnings

import numpy as np
import pandas as pd

import VESIcal as v
from VESIcal import core


INPUT_UNITS = ('wtpt_oxides', 'mol_oxides', 'mol_cations')
NORMALIZATIONS = ('none', 'standard', 'fixedvolatiles', 'additionalvolatiles')
# Published package constants, repeated here so the weight-percent oracle does
# not call the conversion routine under test.
MASSES = {'SiO2': 60.083, 'Al2O3': 101.961, 'MgO': 40.304,
          'MnO': 70.937, 'H2O': 18.02, 'CO2': 44.01}
CATIONS = {'Si': ('SiO2', 1), 'Al': ('Al2O3', 2), 'Mg': ('MgO', 1),
           'Mn': ('MnO', 1), 'H': ('H2O', 2), 'C': ('CO2', 1)}


def composition(units, manganese):
    if units == 'wtpt_oxides':
        return {'SiO2': 60.0, 'Al2O3': 15.0, 'MgO': 10.0,
                'MnO': manganese, 'H2O': 2.0, 'CO2': 0.5}
    if units == 'mol_oxides':
        return {'SiO2': 0.55, 'Al2O3': 0.2, 'MgO': 0.15,
                'MnO': manganese, 'H2O': 0.08, 'CO2': 0.02}
    return {'Si': 0.55, 'Al': 0.4, 'Mg': 0.15,
            'Mn': manganese, 'H': 0.16, 'C': 0.02}


def input_container(values, kind):
    if kind == 'series':
        return pd.Series(values, dtype='float64', name='source-sample')
    return dict(values)


def expected_weight_percent(units):
    values = composition(units, 0.0)
    if units == 'wtpt_oxides':
        return pd.Series(values)
    if units == 'mol_oxides':
        weights = {oxide: amount * MASSES[oxide] for oxide, amount in values.items()}
    else:
        weights = {CATIONS[cation][0]: amount / CATIONS[cation][1]
                   * MASSES[CATIONS[cation][0]] for cation, amount in values.items()}
    total = sum(weights.values())
    return pd.Series({oxide: 100.0 * value / total for oxide, value in weights.items()})


def normalized_weight_percent(values, normalization):
    if normalization == 'none':
        return values.copy()
    if normalization == 'standard':
        return values * 100.0 / sum(values)
    nonvolatile = [oxide for oxide in values.index if oxide not in ('H2O', 'CO2')]
    volatile_total = values.get('H2O', 0.0) + values.get('CO2', 0.0)
    scale = 100.0 - volatile_total if normalization == 'fixedvolatiles' else 100.0
    result = values.copy()
    for oxide in nonvolatile:
        result[oxide] = values[oxide] * scale / sum(values[nonvolatile])
    volatile_order = ('CO2', 'H2O') if normalization == 'fixedvolatiles' else ('H2O', 'CO2')
    return result.loc[nonvolatile + [oxide for oxide in volatile_order if oxide in values]]


class TestSampleMissingConcentrations(unittest.TestCase):
    def assert_input_unchanged(self, actual, original):
        if isinstance(actual, pd.Series):
            pd.testing.assert_series_equal(actual, original, check_exact=True)
        else:
            self.assertEqual(list(actual), list(original))
            for key in actual:
                if pd.isna(original[key]):
                    self.assertTrue(pd.isna(actual[key]))
                else:
                    self.assertEqual(actual[key], original[key])

    def assert_same_getters(self, actual, explicit_zero):
        raw = actual.get_composition(normalization='none', units='wtpt_oxides')
        self.assertFalse(raw.isna().any())
        self.assertEqual(raw['MnO'], 0.0)
        for units in INPUT_UNITS:
            for normalization in NORMALIZATIONS:
                for exclude in (False, True):
                    pd.testing.assert_series_equal(
                        actual.get_composition(units=units, normalization=normalization,
                                               exclude_volatiles=exclude),
                        explicit_zero.get_composition(units=units, normalization=normalization,
                                                      exclude_volatiles=exclude),
                        check_exact=True)
        for exclude in (False, True):
            pd.testing.assert_series_equal(
                actual.get_composition(units='mol_singleO', normalization='none',
                                       exclude_volatiles=exclude),
                explicit_zero.get_composition(units='mol_singleO', normalization='none',
                                              exclude_volatiles=exclude), check_exact=True)
            self.assertEqual(actual.get_formulaweight(exclude_volatiles=exclude),
                             explicit_zero.get_formulaweight(exclude_volatiles=exclude))
        self.assertEqual(actual.get_composition(species='MnO'), 0.0)
        self.assertEqual(actual.get_composition(species='Mn'), 0.0)
        pd.testing.assert_series_equal(
            actual.get_composition(asSampleClass=True).get_composition(),
            explicit_zero.get_composition(asSampleClass=True).get_composition(),
            check_exact=True)

    def check_constructor(self, units, kind):
        given = input_container(composition(units, np.nan), kind)
        original = deepcopy(given)
        actual = v.Sample(given, units=units)
        explicit_zero = v.Sample(input_container(composition(units, 0.0), kind), units=units)
        self.assert_input_unchanged(given, original)
        self.assert_same_getters(actual, explicit_zero)
        expected = expected_weight_percent(units)
        for normalization in NORMALIZATIONS:
            pd.testing.assert_series_equal(
                actual.get_composition(normalization=normalization, units='wtpt_oxides'),
                normalized_weight_percent(expected, normalization),
                check_names=False, check_dtype=False, rtol=1e-12, atol=1e-12)

    def test_constructor_wt_dict(self):
        self.check_constructor('wtpt_oxides', 'dict')

    def test_constructor_wt_series(self):
        self.check_constructor('wtpt_oxides', 'series')

    def test_constructor_mol_oxides_dict(self):
        self.check_constructor('mol_oxides', 'dict')

    def test_constructor_mol_oxides_series(self):
        self.check_constructor('mol_oxides', 'series')

    def test_constructor_mol_cations_dict(self):
        self.check_constructor('mol_cations', 'dict')

    def test_constructor_mol_cations_series(self):
        self.check_constructor('mol_cations', 'series')

    def check_change(self, units, kind, inplace, add_new=False):
        initial = composition('wtpt_oxides', 1.0)
        if add_new:
            del initial['MnO']
        actual = v.Sample(initial)
        explicit_zero = v.Sample(initial)
        before = actual.get_composition()
        key = 'Mn' if units == 'mol_cations' else 'MnO'
        given = input_container({key: np.nan}, kind)
        original = deepcopy(given)
        changed = actual.change_composition(given, units=units, inplace=inplace)
        expected = explicit_zero.change_composition(
            input_container({key: 0.0}, kind), units=units, inplace=inplace)
        self.assert_input_unchanged(given, original)
        if inplace:
            self.assertIs(changed, actual)
        else:
            self.assertIsNot(changed, actual)
            pd.testing.assert_series_equal(actual.get_composition(), before, check_exact=True)
        self.assert_same_getters(changed, expected)

    def test_change_wt_dict_inplace(self):
        self.check_change('wtpt_oxides', 'dict', True)

    def test_change_wt_dict_copy(self):
        self.check_change('wtpt_oxides', 'dict', False)

    def test_change_wt_series_inplace(self):
        self.check_change('wtpt_oxides', 'series', True)

    def test_change_wt_series_copy(self):
        self.check_change('wtpt_oxides', 'series', False)

    def test_change_mol_oxides_dict_inplace(self):
        self.check_change('mol_oxides', 'dict', True)

    def test_change_mol_oxides_dict_copy(self):
        self.check_change('mol_oxides', 'dict', False)

    def test_change_mol_oxides_series_inplace(self):
        self.check_change('mol_oxides', 'series', True)

    def test_change_mol_oxides_series_copy(self):
        self.check_change('mol_oxides', 'series', False)

    def test_change_mol_cations_dict_inplace(self):
        self.check_change('mol_cations', 'dict', True)

    def test_change_mol_cations_dict_copy(self):
        self.check_change('mol_cations', 'dict', False)

    def test_change_mol_cations_series_inplace(self):
        self.check_change('mol_cations', 'series', True)

    def test_change_mol_cations_series_copy(self):
        self.check_change('mol_cations', 'series', False)

    def test_reported_dataframe_example(self):
        kilauea = {'SiO2': 50.142, 'TiO2': 2.257, 'Al2O3': 12.224, 'FeO': 11.33,
                   'MnO': 0.323, 'MgO': 10.215, 'CaO': 10.476, 'Na2O': 1.754,
                   'K2O': 0.398, 'P2O5': 0.233, 'H2O': 0.3, 'CO2': 0.08}
        morb = {'SiO2': 47.4, 'TiO2': 1.01, 'Al2O3': 17.64, 'FeO': 7.98,
                'MgO': 7.63, 'CaO': 12.44, 'Na2O': 2.65, 'K2O': 0.03,
                'P2O5': 0.08, 'H2O': 0.2, 'CO2': 0.11}
        labels = ['Kilauea_unnorm', 'MORB_unnorm']
        frame = pd.DataFrame([v.Sample(kilauea).get_composition(),
                              v.Sample(morb).get_composition()], index=labels)
        original = frame.copy(deep=True)
        self.assertTrue(pd.isna(frame.loc['MORB_unnorm', 'MnO']))
        outputs = frame.copy(deep=True)
        for label, values in zip(labels, (kilauea, dict(morb, MnO=0.0))):
            sample = v.Sample(frame.loc[label])
            expected = pd.Series(values, index=frame.columns)
            for normalization in NORMALIZATIONS[1:]:
                result = sample.get_composition(normalization=normalization)
                self.assertFalse(result.isna().any())
                pd.testing.assert_series_equal(
                    result, normalized_weight_percent(expected, normalization),
                    check_names=False, check_dtype=False, rtol=1e-12, atol=1e-12)
                outputs.loc[label + '_' + normalization] = result
        self.assertEqual(len(outputs), 8)
        self.assertEqual(outputs.loc['MORB_unnorm_standard', 'MnO'], 0.0)
        pd.testing.assert_frame_equal(frame, original, check_exact=True)

    def test_missing_volatiles_and_multiple_fields(self):
        given = composition('wtpt_oxides', np.nan)
        given.update(H2O=np.nan, CO2=np.nan)
        explicit_zero = dict(given, MnO=0.0, H2O=0.0, CO2=0.0)
        self.assert_same_getters(v.Sample(given), v.Sample(explicit_zero))

    def test_nullable_series(self):
        given = pd.Series(composition('wtpt_oxides', np.nan), dtype='Float64')
        original = given.copy(deep=True)
        expected = pd.Series(composition('wtpt_oxides', 0.0), dtype='Float64')
        self.assert_same_getters(v.Sample(given), v.Sample(expected))
        pd.testing.assert_series_equal(given, original, check_exact=True)

    def test_all_missing_inputs_match_explicit_zeros(self):
        actual = v.Sample({'SiO2': np.nan, 'MgO': np.nan})
        expected = v.Sample({'SiO2': 0.0, 'MgO': 0.0})
        pd.testing.assert_series_equal(actual.get_composition(), expected.get_composition(),
                                       check_exact=True)

    def test_adding_missing_component(self):
        for units in INPUT_UNITS:
            for kind in ('dict', 'series'):
                for inplace in (False, True):
                    with self.subTest(units=units, kind=kind, inplace=inplace):
                        self.check_change(units, kind, inplace, add_new=True)

    def test_finite_constructor_controls(self):
        for units in INPUT_UNITS:
            for kind in ('dict', 'series'):
                with self.subTest(units=units, kind=kind):
                    given = input_container(composition(units, 0.0), kind)
                    original = deepcopy(given)
                    sample = v.Sample(given, units=units)
                    pd.testing.assert_series_equal(
                        sample.get_composition(), expected_weight_percent(units),
                        check_names=False, check_dtype=False, rtol=1e-12, atol=1e-12)
                    sample.set_default_units('mol_oxides')
                    sample.set_default_normalization('standard')
                    pd.testing.assert_series_equal(
                        sample.get_composition(), sample.get_composition(
                            units='mol_oxides', normalization='standard'), check_exact=True)
                    self.assert_input_unchanged(given, original)

    def test_finite_change_controls(self):
        for units in INPUT_UNITS:
            for kind in ('dict', 'series'):
                for inplace in (False, True):
                    with self.subTest(units=units, kind=kind, inplace=inplace):
                        sample = v.Sample(composition('wtpt_oxides', 1.0))
                        before = sample.get_composition()
                        key = 'Mn' if units == 'mol_cations' else 'MnO'
                        given = input_container({key: 0.0}, kind)
                        original = deepcopy(given)
                        changed = sample.change_composition(given, units=units, inplace=inplace)
                        self.assertEqual(changed.get_composition()['MnO'], 0.0)
                        self.assertFalse(changed.get_composition().isna().any())
                        self.assert_input_unchanged(given, original)
                        if inplace:
                            self.assertIs(changed, sample)
                        else:
                            self.assertIsNot(changed, sample)
                            pd.testing.assert_series_equal(sample.get_composition(), before,
                                                           check_exact=True)

    def test_zero_and_series_name_order_preserved(self):
        given = pd.Series({'MgO': 40.0, 'MnO': 0.0, 'SiO2': 60.0}, name='named-sample')
        pd.testing.assert_series_equal(v.Sample(given).get_composition(), given,
                                       check_exact=True)

    def test_invalid_inputs_keep_errors(self):
        with self.assertRaises(core.InputError):
            v.Sample([60.0, 40.0])
        with self.assertRaises(ValueError):
            v.Sample({'SiO2': 60.0, 'sample': 'metadata'})
        with self.assertRaises(core.InputError):
            v.Sample({'SiO2': 60.0}, units='invalid')
        sample = v.Sample({'SiO2': 60.0, 'MgO': 40.0})
        before = sample.get_composition()
        with self.assertRaises(core.InputError):
            sample.change_composition({'MnO': np.nan}, units='invalid')
        pd.testing.assert_series_equal(sample.get_composition(), before, check_exact=True)

    def test_infinity_is_not_treated_as_missing(self):
        sample = v.Sample({'SiO2': np.inf, 'MgO': 40.0})
        self.assertEqual(sample.get_composition()['SiO2'], np.inf)
        sample.change_composition({'MgO': -np.inf})
        self.assertEqual(sample.get_composition()['MgO'], -np.inf)

    def test_existing_series_text_is_not_dropped(self):
        given = pd.Series({'SiO2': 60.0, 'sample': 'metadata'}, name='metadata-sample')
        pd.testing.assert_series_equal(v.Sample(given).get_composition(), given,
                                       check_exact=True)
        with self.assertRaises(TypeError):
            v.Sample(given).get_composition(normalization='standard')

    def test_empty_composition_stays_empty(self):
        self.assertTrue(v.Sample({}).get_composition().empty)

    def test_finite_object_numeric_series_dtype_preserved(self):
        given = pd.Series({'SiO2': 60.0, 'MgO': 40.0, 'MnO': 0.0},
                          dtype=object, name='object-numeric')
        pd.testing.assert_series_equal(v.Sample(given).get_composition(), given,
                                       check_exact=True)

    def test_missing_object_numeric_series(self):
        given = pd.Series(composition('wtpt_oxides', np.nan),
                          dtype=object, name='object-missing')
        original = given.copy(deep=True)
        expected = pd.Series(composition('wtpt_oxides', 0.0),
                             dtype=object, name='object-missing')
        self.assert_same_getters(v.Sample(given), v.Sample(expected))
        pd.testing.assert_series_equal(given, original, check_exact=True)

    def check_derived_missing_ingestion(self, values, raw_expected, warning_count):
        sample = v.Sample(values)
        with warnings.catch_warnings(record=True) as raw_warnings:
            warnings.simplefilter('always', RuntimeWarning)
            raw = sample.get_composition(normalization='standard')
        pd.testing.assert_series_equal(raw, pd.Series(raw_expected), check_exact=True)
        with warnings.catch_warnings(record=True) as wrapped_warnings:
            warnings.simplefilter('always', RuntimeWarning)
            wrapped = sample.get_composition(normalization='standard', asSampleClass=True)
        expected_warnings = [('RuntimeWarning', 'invalid value encountered in scalar divide')
                             for _ in range(warning_count)]
        self.assertEqual([(w.category.__name__, str(w.message)) for w in raw_warnings],
                         expected_warnings)
        self.assertEqual([(w.category.__name__, str(w.message)) for w in wrapped_warnings],
                         expected_warnings)
        # The raw undefined normalization is retained. Wrapping explicitly ingests
        # that returned Series as a new Sample, whose missing values become zero.
        pd.testing.assert_series_equal(
            wrapped.get_composition(), pd.Series({'SiO2': 0.0, 'MgO': 0.0}),
            check_exact=True)

    def test_zero_total_raw_and_wrapped_missing_ingestion(self):
        self.check_derived_missing_ingestion(
            {'SiO2': 0.0, 'MgO': 0.0}, {'SiO2': np.nan, 'MgO': np.nan}, 2)

    def test_infinite_total_raw_and_wrapped_missing_ingestion(self):
        self.check_derived_missing_ingestion(
            {'SiO2': np.inf, 'MgO': 40.0}, {'SiO2': np.nan, 'MgO': 0.0}, 1)

    def check_legacy_molar_constructor(self, units, infinite):
        first = 'Si' if units == 'mol_cations' else 'SiO2'
        second = 'Mg' if units == 'mol_cations' else 'MgO'
        values = {first: np.inf if infinite else 0.0, second: 40.0 if infinite else 0.0}
        original = dict(values)
        raw = v.Sample(values, units=units).get_composition()
        expected = pd.Series({'SiO2': np.nan, 'MgO': 0.0 if infinite else np.nan})
        pd.testing.assert_series_equal(raw, expected, check_exact=True)
        self.assertEqual(values, original)

    def test_legacy_molar_oxide_zero_constructor(self):
        self.check_legacy_molar_constructor('mol_oxides', False)

    def test_legacy_molar_cation_zero_constructor(self):
        self.check_legacy_molar_constructor('mol_cations', False)

    def test_legacy_molar_oxide_infinite_constructor(self):
        self.check_legacy_molar_constructor('mol_oxides', True)

    def test_legacy_molar_cation_infinite_constructor(self):
        self.check_legacy_molar_constructor('mol_cations', True)

    def check_legacy_molar_zero_change(self, units):
        keys = ('Si', 'Mg') if units == 'mol_cations' else ('SiO2', 'MgO')
        for inplace in (False, True):
            with self.subTest(inplace=inplace):
                sample = v.Sample({'SiO2': 60.0, 'MgO': 40.0})
                before = sample.get_composition()
                given = pd.Series({key: 0.0 for key in keys})
                original = given.copy(deep=True)
                changed = sample.change_composition(given, units=units, inplace=inplace)
                pd.testing.assert_series_equal(
                    changed.get_composition(), pd.Series({'SiO2': np.nan, 'MgO': np.nan}),
                    check_exact=True)
                pd.testing.assert_series_equal(given, original, check_exact=True)
                if not inplace:
                    pd.testing.assert_series_equal(sample.get_composition(), before,
                                                   check_exact=True)

    def test_legacy_molar_oxide_zero_change(self):
        self.check_legacy_molar_zero_change('mol_oxides')

    def test_legacy_molar_cation_zero_change(self):
        self.check_legacy_molar_zero_change('mol_cations')

    def category_input(self, manganese, include_zero=False):
        categories = [0.5, 1.0, 2.0, 10.0, 15.0, 60.0]
        if include_zero:
            categories.append(0.0)
        return pd.Series(composition('wtpt_oxides', manganese),
                         dtype=pd.CategoricalDtype(categories), name='category-sample')

    def test_missing_categorical_series(self):
        given = self.category_input(np.nan)
        original = given.copy(deep=True)
        expected = self.category_input(0.0, include_zero=True)
        self.assert_same_getters(v.Sample(given), v.Sample(expected))
        pd.testing.assert_series_equal(given, original, check_exact=True)

    def test_finite_categorical_series_preserved(self):
        given = self.category_input(1.0)
        pd.testing.assert_series_equal(v.Sample(given).get_composition(), given,
                                       check_exact=True)

    def test_missing_sparse_series(self):
        given = pd.Series(composition('wtpt_oxides', np.nan),
                          dtype=pd.SparseDtype('float64', 0.0), name='sparse-sample')
        original = given.copy(deep=True)
        expected = pd.Series(composition('wtpt_oxides', 0.0),
                             dtype=pd.SparseDtype('float64', 0.0), name='sparse-sample')
        self.assert_same_getters(v.Sample(given), v.Sample(expected))
        pd.testing.assert_series_equal(given, original, check_exact=True)

    def test_finite_sparse_series_preserved(self):
        given = pd.Series(composition('wtpt_oxides', 1.0),
                          dtype=pd.SparseDtype('float64', 0.0), name='sparse-sample')
        pd.testing.assert_series_equal(v.Sample(given).get_composition(), given,
                                       check_exact=True)

    def check_missing_categorical_change(self, inplace):
        given = self.category_input(1.0)
        sample = v.Sample(given)
        before = sample.get_composition()
        changed = sample.change_composition({'MnO': np.nan}, inplace=inplace)
        expected = v.Sample(self.category_input(0.0, include_zero=True))
        self.assert_same_getters(changed, expected)
        pd.testing.assert_series_equal(given, before, check_exact=True)
        if not inplace:
            pd.testing.assert_series_equal(sample.get_composition(), before, check_exact=True)

    def test_missing_categorical_change_inplace(self):
        self.check_missing_categorical_change(True)

    def test_missing_categorical_change_copy(self):
        self.check_missing_categorical_change(False)

    def test_categorical_nonmissing_invalid_value_keeps_error(self):
        sample = v.Sample(self.category_input(1.0))
        before = sample.get_composition()
        with self.assertRaisesRegex(TypeError, 'new category'):
            sample.change_composition({'MnO': 0.0})
        pd.testing.assert_series_equal(sample.get_composition(), before, check_exact=True)


if __name__ == '__main__':
    unittest.main()
