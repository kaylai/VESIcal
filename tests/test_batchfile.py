import unittest
import VESIcal as v
import pandas as pd
import pathlib
import warnings


HEADER_PREFIX = 'Possible unrecognized oxide column headers:'
HEADER_VARIANTS = [' SiO2', 'SiO2 ', 'sio2', 'SiO2_melt', 'SiO2, wt%']


def header_frame(header='SiO2'):
    chemistry = {'SiO2': 50.142, 'TiO2': 2.257, 'Al2O3': 12.224,
                 'FeO': 11.33, 'MnO': 0.323, 'MgO': 10.215,
                 'CaO': 10.476, 'Na2O': 1.754, 'K2O': 0.398,
                 'P2O5': 0.233, 'H2O': 0.3, 'CO2': 0.08}
    frame = pd.DataFrame([chemistry], index=['sample-A'])
    frame['Notes'] = 'metadata survives'
    return frame.rename(columns={'SiO2': header})

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


class TestOxideHeaderDiagnostics(unittest.TestCase):
    def import_frame(self, frame, **kwargs):
        original = frame.copy(deep=True)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            batch = v.BatchFile_from_DataFrame(frame, **kwargs)
        pd.testing.assert_frame_equal(frame, original)
        header_warnings = [warning for warning in caught
                           if str(warning.message).startswith(HEADER_PREFIX)]
        return batch, header_warnings, caught

    def test_dataframe_variants(self):
        for header in HEADER_VARIANTS:
            with self.subTest(header=header):
                batch, header_warnings, _ = self.import_frame(header_frame(header))
                self.assertEqual(len(header_warnings), 1)
                self.assertIs(header_warnings[0].category, RuntimeWarning)
                self.assertIn(repr(header), str(header_warnings[0].message))
                self.assertAlmostEqual(batch.get_data().loc['sample-A', header], 50.142)
                composition = batch.get_composition()
                self.assertNotIn('SiO2', composition.columns)
                self.assertAlmostEqual(composition.loc['sample-A'].sum(), 49.590)

    def test_canonical_and_unrelated_metadata(self):
        frame = header_frame()
        for column in ['SiO2_sd', 'SiO2_2sigma', 'SiO2_notes', 'H2O_liq_VESIcal',
                       'SiO2_melt_sd', 'SiO2,wt%', 'SiO2 _melt', 'SiO2 (wt%)', 3, None]:
            frame[column] = 'metadata'
        batch, header_warnings, _ = self.import_frame(frame)
        self.assertEqual(header_warnings, [])
        self.assertEqual(batch.get_data().index.tolist(), ['sample-A'])
        self.assertAlmostEqual(batch.get_composition().loc['sample-A'].sum(), 99.732)
        pd.testing.assert_frame_equal(batch.get_data(), frame)

    def test_all_canonical_oxide_names(self):
        frame = pd.DataFrame([{oxide: 1.0 for oxide in v.core.oxides}], index=['sample-A'])
        batch, header_warnings, _ = self.import_frame(frame)
        self.assertEqual(header_warnings, [])
        self.assertEqual(set(batch.get_composition().columns), set(v.core.oxides))

    def test_bounded_case_and_phase_forms(self):
        for header in [' SIO2 ', 'sIO2_MELT', ' sio2_OLIVINE ', 'SIO2, WT%']:
            with self.subTest(header=header):
                _, header_warnings, _ = self.import_frame(header_frame(header))
                self.assertEqual(len(header_warnings), 1)
                self.assertIn(repr(header), str(header_warnings[0].message))

    def test_phase_and_canonical_collisions_are_grouped(self):
        for canonical_present in [False, True]:
            with self.subTest(canonical_present=canonical_present):
                frame = header_frame('SiO2_melt')
                frame['SiO2_olivine'] = 40.0
                if canonical_present:
                    frame['SiO2'] = 51.0
                batch, header_warnings, _ = self.import_frame(frame)
                self.assertEqual(len(header_warnings), 1)
                message = str(header_warnings[0].message)
                self.assertIn('SiO2_melt', message)
                self.assertIn('SiO2_olivine', message)
                self.assertIn('multiple candidates', message)
                self.assertEqual('canonical column' in message, canonical_present)
                composition = batch.get_composition()
                if canonical_present:
                    self.assertAlmostEqual(composition.loc['sample-A', 'SiO2'], 51.0)
                else:
                    self.assertNotIn('SiO2', composition.columns)
                self.assertCountEqual(batch.get_data().columns, frame.columns)
                pd.testing.assert_frame_equal(batch.get_data().loc[:, frame.columns], frame)

    def test_explicit_label_is_excluded(self):
        for label in ['sio2', ' SiO2']:
            with self.subTest(label=label):
                frame = header_frame().reset_index(names=label)
                batch, header_warnings, _ = self.import_frame(frame, label=label)
                self.assertEqual(header_warnings, [])
                self.assertEqual(batch.get_data().index.tolist(), ['sample-A'])
                self.assertAlmostEqual(batch.get_composition().loc['sample-A', 'SiO2'],
                                       50.142)

    def test_total_iron_aliases_keep_existing_warning_and_values(self):
        for alias in ['FeOt', 'FeOT', 'FeOtot', 'FeOtotal', 'FeOstar', 'FeO*']:
            with self.subTest(alias=alias):
                frame = header_frame().rename(columns={'FeO': alias})
                batch, header_warnings, caught = self.import_frame(frame)
                self.assertEqual(header_warnings, [])
                self.assertEqual(sum('Total iron column' in str(item.message)
                                     for item in caught), 1)
                self.assertAlmostEqual(batch.get_composition().loc['sample-A', 'FeO'],
                                       11.33)

    def test_warning_is_not_repeated_when_reading(self):
        batch, header_warnings, _ = self.import_frame(header_frame('sio2'))
        self.assertEqual(len(header_warnings), 1)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            for _ in range(3):
                batch.get_data()
                batch.get_composition()
        self.assertFalse(any(str(item.message).startswith(HEADER_PREFIX) for item in caught))

    def test_other_input_units_do_not_gain_header_warning(self):
        _, header_warnings, _ = self.import_frame(header_frame('sio2'), units='mol_oxides')
        self.assertEqual(header_warnings, [])

    def test_warning_can_be_filtered_without_hiding_other_runtime_warnings(self):
        with warnings.catch_warnings():
            warnings.filterwarnings('error', message=HEADER_PREFIX, category=RuntimeWarning)
            with self.assertRaisesRegex(RuntimeWarning, HEADER_PREFIX):
                v.BatchFile_from_DataFrame(header_frame('sio2'))
        with warnings.catch_warnings():
            warnings.simplefilter('error', RuntimeWarning)
            warnings.filterwarnings('ignore', message=HEADER_PREFIX, category=RuntimeWarning)
            v.BatchFile_from_DataFrame(header_frame('sio2'))
            with self.assertRaisesRegex(RuntimeWarning, 'Total iron column'):
                frame = header_frame('sio2').rename(columns={'FeO': 'FeOt'})
                v.BatchFile_from_DataFrame(frame)

    def test_dataframe_wrapper_label_defaults_are_preserved(self):
        frame = pd.DataFrame([{' SiO2': 50.142, 'MgO': 10.215}], index=['sample-A'])
        batch, header_warnings, _ = self.import_frame(frame)
        self.assertEqual(len(header_warnings), 1)
        self.assertEqual(batch.get_data().index.tolist(), ['sample-A'])
        original = frame.copy(deep=True)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            batch = v.batchfile.from_DataFrame(frame)
        pd.testing.assert_frame_equal(frame, original)
        self.assertEqual(batch.get_data().index.tolist(), ['50.142'])
        self.assertTrue(str(caught[0].message).startswith(HEADER_PREFIX))
        self.assertTrue(any('No Label column given' in str(item.message) for item in caught))
