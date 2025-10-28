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
        