import unittest
import VESIcal as v
import pandas as pd
import pathlib

# Allow unittest to find the file
TEST_FILE = pathlib.Path(__file__).parent.joinpath("ImportTest.xlsx")


class TestCreateBatchFileFromDataFrame(unittest.TestCase):
    def assertDataframeEqual(self, a, b, msg):
        """
        Creates a new type of unittest to assert that pd DataFrames
        are equal, inheriting from pandas testing routine
        """
        try:
            pd._testing.assert_frame_equal(a, b)
        except AssertionError as e:
            raise self.failureException(msg) from e

    def setUp(self):
        self.addTypeEqualityFunc(pd.DataFrame, self.assertDataframeEqual)

        self.myfile = v.BatchFile(TEST_FILE)

        self.data_dict = {'SiO2':  [47.95, 47.95],
                          'TiO2':  [1.67, 1.67],
                          'Al2O3': [17.32, 17.32],
                          'FeO':   [10.24, 10.24],
                          'Fe2O3': [0.1, 0.1],
                          'MgO':   [5.76, 5.76],
                          'CaO':   [10.93, 10.93],
                          'Na2O':  [3.45, 3.45],
                          'K2O':   [1.99, 1.99],
                          'P2O5':  [0.51, 0.51],
                          'MnO':   [0.1, 0.1],
                          'H2O':   [2.0, 2.0],
                          'CO2':   [0.1, 0.1],
                          'Notes': ['Normal sample', 'Duplicate sample']}

        self.df = pd.DataFrame(self.data_dict,
                               index=['test_samp',
                                      'test_samp-duplicate-1'])
        self.bf = v.BatchFile_from_DataFrame(self.df)

    def test_BatchFileFromDataFrame(self):
        self.assertEqual(self.bf.get_data(), self.myfile.get_data(),
                         'BatchFiles are different')


class TestGetData(unittest.TestCase):
    def test_function_processes_all_oxides(self):
        """Verify conversion functions can process every oxide."""
        sample_all_oxides = {ox: i+1 for i, ox in enumerate(v.core.oxides)}
        df_all_oxides = pd.DataFrame(sample_all_oxides, index=[0])
        bf_all_oxides = v.BatchFile_from_DataFrame(df_all_oxides)

        for unitname in ['wtpt_oxides', 'mol_oxides', 'mol_cations',
                         'mol_singleO']:
            result = bf_all_oxides.get_data(units=unitname)

            # Verify we got a result
            self.assertIsNotNone(result)
            self.assertGreater(len(result), 0)


class TestMolarInputDataFrame(unittest.TestCase):
    def check_input(self, units, rows):
        from tests.molar_input_fixtures import (
            assert_molar_batch, molar_frame)

        frame = molar_frame(units, rows)
        original = frame.copy(deep=True)
        try:
            batch = v.BatchFile_from_DataFrame(frame, units=units, label='Label')
        finally:
            pd.testing.assert_frame_equal(frame, original)
        assert_molar_batch(self, batch, frame, units)

    def test_mol_oxides_one_row(self):
        self.check_input('mol_oxides', 1)

    def test_mol_oxides_two_rows(self):
        self.check_input('mol_oxides', 2)

    def test_mol_cations_one_row(self):
        self.check_input('mol_cations', 1)

    def test_mol_cations_two_rows(self):
        self.check_input('mol_cations', 2)


class TestMolarInputControls(unittest.TestCase):
    def check_scaling(self, units):
        from tests.molar_input_fixtures import (
            assert_molar_batch, molar_frame)

        for scale in [1e-6, 1e6]:
            with self.subTest(scale=scale):
                frame = molar_frame(units, scale=scale)
                original = frame.copy(deep=True)
                batch = v.BatchFile_from_DataFrame(
                    frame, units=units, label='Label')
                pd.testing.assert_frame_equal(frame, original)
                assert_molar_batch(self, batch, frame, units)

    def test_mol_oxide_scaling(self):
        self.check_scaling('mol_oxides')

    def test_mol_cation_scaling(self):
        self.check_scaling('mol_cations')

    def test_numeric_labels_follow_existing_string_conversion(self):
        from tests.molar_input_fixtures import (
            assert_molar_batch, molar_frame)

        for units in ['mol_oxides', 'mol_cations']:
            with self.subTest(units=units):
                frame = molar_frame(units)
                frame['Label'] = [101, 17]
                original = frame.copy(deep=True)
                batch = v.BatchFile_from_DataFrame(
                    frame, units=units, label='Label')
                pd.testing.assert_frame_equal(frame, original)
                assert_molar_batch(self, batch, frame, units)

    def test_wtpt_input_and_metadata_are_preserved(self):
        from tests.molar_input_fixtures import molar_frame

        frame = molar_frame('mol_oxides')
        frame['SiO2'] = [80., 30.]
        frame['MgO'] = [10., 15.]
        # In wt% input, this field is metadata, not a molar cation input.
        frame['Si'] = [17., 19.]
        original = frame.copy(deep=True)
        batch = v.BatchFile_from_DataFrame(frame, label='Label')
        pd.testing.assert_frame_equal(frame, original)
        expected = frame.set_index('Label')
        pd.testing.assert_frame_equal(batch.get_data()[expected.columns],
                                      expected, check_names=False)
        self.assertEqual(batch.default_units, 'wtpt_oxides')
        self.assertEqual(batch.default_normalization, 'none')

    def test_custom_defaults_do_not_change_internal_weight_percent(self):
        from tests.molar_input_fixtures import (
            expected_wtpercent, molar_frame)
        import numpy as np

        for units in ['mol_oxides', 'mol_cations']:
            with self.subTest(units=units):
                frame = molar_frame(units)
                original = frame.copy(deep=True)
                batch = v.BatchFile(None, dataframe=frame, units=units,
                                    label='Label', default_units='mol_cations',
                                    default_normalization='additionalvolatiles')
                pd.testing.assert_frame_equal(frame, original)
                self.assertEqual(batch.default_units, 'mol_cations')
                self.assertEqual(batch.default_normalization,
                                 'additionalvolatiles')
                expected = expected_wtpercent(frame, units)
                # Check stored weight percent independently of output settings.
                np.testing.assert_allclose(batch.data[expected.columns],
                                           expected, rtol=0, atol=1e-9)
