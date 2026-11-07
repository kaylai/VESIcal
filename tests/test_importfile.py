import unittest
import VESIcal as v
import pandas as pd
import pathlib
import tempfile
import warnings

from tests.test_batchfile import HEADER_PREFIX, HEADER_VARIANTS, header_frame

# Allow unittest to find the file
TEST_FILE = pathlib.Path(__file__).parent.joinpath("ImportTest.xlsx")

class TestImportExcel(unittest.TestCase):
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
        # Add assertDataframeEqual to unittest test cases
        self.addTypeEqualityFunc(pd.DataFrame, self.assertDataframeEqual)

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

        self.myfile = v.BatchFile(TEST_FILE)

    def test_ImportExcel(self):
        self.assertEqual(self.df, self.myfile.get_data(), 
                         'DataFrames are different')


class TestHeaderDiagnosticsFileImports(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)

    def save_frame(self, frame, kind, include_label=True):
        path = pathlib.Path(self.directory.name) / ('composition.' + kind)
        kwargs = {'index_label': 'Label'} if include_label else {'index': False}
        if kind == 'csv':
            frame.to_csv(path, **kwargs)
        else:
            frame.to_excel(path, **kwargs)
        return path

    def test_header_variants_through_csv_and_xlsx(self):
        for kind in ['csv', 'xlsx']:
            for header in ['SiO2'] + HEADER_VARIANTS:
                with self.subTest(kind=kind, header=header):
                    frame = header_frame(header)
                    original = frame.copy(deep=True)
                    path = self.save_frame(frame, kind)
                    with warnings.catch_warnings(record=True) as caught:
                        warnings.simplefilter('always')
                        batch = v.BatchFile(path)
                    pd.testing.assert_frame_equal(frame, original)
                    header_warnings = [item for item in caught
                                       if str(item.message).startswith(HEADER_PREFIX)]
                    self.assertEqual(len(header_warnings), int(header != 'SiO2'))
                    self.assertAlmostEqual(batch.get_data().loc['sample-A', header], 50.142)
                    composition = batch.get_composition()
                    if header == 'SiO2':
                        self.assertAlmostEqual(composition.loc['sample-A', 'SiO2'], 50.142)
                    else:
                        self.assertNotIn('SiO2', composition.columns)

    def test_warning_precedes_fallback_index_selection(self):
        frame = pd.DataFrame([{' SiO2': 50.142, 'MgO': 10.215}])
        for kind in ['csv', 'xlsx']:
            with self.subTest(kind=kind):
                path = self.save_frame(frame, kind, include_label=False)
                with warnings.catch_warnings(record=True) as caught:
                    warnings.simplefilter('always')
                    batch = v.BatchFile(path)
                self.assertTrue(str(caught[0].message).startswith(HEADER_PREFIX))
                self.assertIn("' SiO2'", str(caught[0].message))
                self.assertTrue(any('No Label column given' in str(item.message)
                                    for item in caught))
                self.assertEqual(batch.get_data().index.tolist(), ['50.142'])
                self.assertNotIn(' SiO2', batch.get_data().columns)

    def test_custom_label_is_excluded_in_files(self):
        for kind in ['csv', 'xlsx']:
            with self.subTest(kind=kind):
                frame = header_frame().reset_index(names='sio2')
                path = self.save_frame(frame, kind, include_label=False)
                with warnings.catch_warnings(record=True) as caught:
                    warnings.simplefilter('always')
                    batch = v.BatchFile(path, label='sio2')
                self.assertFalse(any(str(item.message).startswith(HEADER_PREFIX)
                                     for item in caught))
                self.assertEqual(batch.get_data().index.tolist(), ['sample-A'])
                self.assertAlmostEqual(batch.get_composition().loc['sample-A', 'SiO2'],
                                       50.142)

    def test_canonical_columns_without_label_still_create_index(self):
        frame = pd.DataFrame([{'SiO2': 50.142, 'MgO': 10.215}])
        for kind in ['csv', 'xlsx']:
            with self.subTest(kind=kind):
                path = self.save_frame(frame, kind, include_label=False)
                with warnings.catch_warnings(record=True) as caught:
                    warnings.simplefilter('always')
                    batch = v.BatchFile(path)
                self.assertFalse(any(str(item.message).startswith(HEADER_PREFIX)
                                     for item in caught))
                self.assertEqual(batch.get_data().index.tolist(), ['0'])
                self.assertAlmostEqual(batch.get_data().loc['0', 'SiO2'], 50.142)
