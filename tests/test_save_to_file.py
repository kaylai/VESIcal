import unittest
import tempfile
import os
import pandas as pd

from VESIcal import sample_class, batchfile
from VESIcal import save_util

class DummyObject:
    def to_dict(self):
        return {"dummy": 1, "value": 42}


class TestSaveToFileFunction(unittest.TestCase):

    def setUp(self):
        # create temporary place to store written files
        self.temp_dir = tempfile.TemporaryDirectory()
        self.temp_path = self.temp_dir.name
        
        # Create a temporary CSV file to import to BatchFile object
        df = pd.DataFrame({
            'Sample': ['Test1'],
            'SiO2': [72.4],
            'Al2O3': [13.1],
            'FeO': [2.4]
        })
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            df.to_csv(f.name, index=False)
            temp_filename = f.name        
        
        # create objects of supported types for testing
        # VESIcal objects
        self.sample = sample_class.Sample({"SiO2": 70, "Al2O3": 15})
        self.batch = batchfile.BatchFile(temp_filename, label="Sample")
        os.remove(temp_filename)
        
        # pandas objects
        self.df = pd.DataFrame({"a": [1, 2], "b": [3, 4]})
        self.pd_series = pd.Series({"a": 1, "b": 2})
        
        # dicts
        self.dict_data = {"x": 10, "y": 20}
        
        # anything else with a `to_dict()` method
        self.object_data = DummyObject()
        
        # list of dummy objects
        self.dummy_list = [DummyObject(), DummyObject()]

    def tearDown(self):
        self.temp_dir.cleanup()

    def test_save_single_dataframe(self):
        path = os.path.join(self.temp_path, "single.csv")
        save_util.save_to_file(path, self.df)
        self.assertTrue(os.path.exists(path))

    def test_save_list_of_dicts(self):
        path = os.path.join(self.temp_path, "dicts.csv")
        save_util.save_to_file(path, [self.dict_data, self.dict_data])
        # re-import saved data to ensure proper transformation from list to
        # one file with two rows
        df = pd.read_csv(path)
        self.assertEqual(len(df), 2)

    def test_save_series(self):
        path = os.path.join(self.temp_path, "series.csv")
        save_util.save_to_file(path, self.pd_series)
        self.assertTrue(os.path.exists(path))

    def test_save_sample(self):
        path = os.path.join(self.temp_path, "sample.csv")
        save_util.save_to_file(path, self.sample)
        self.assertTrue(os.path.exists(path))

    def test_save_batchfile(self):
        path = os.path.join(self.temp_path, "batch.csv")
        save_util.save_to_file(path, self.batch)
        self.assertTrue(os.path.exists(path))

    def test_save_mixed_combined(self):
        path = os.path.join(self.temp_path, "combined.csv")
        save_util.save_to_file(
            path,
            self.df,
            self.dict_data,
            self.sample,
            self.batch,
            self.object_data,
            self.pd_series,
        )
        df = pd.read_csv(path)
        self.assertGreaterEqual(len(df), 5)

    def test_save_mixed_split(self):
        path = os.path.join(self.temp_path, "split.csv")
        save_util.save_to_file(path, self.df, self.dict_data, self.pd_series, combine=False)
        for i in range(3):
            self.assertTrue(os.path.exists(os.path.join(self.temp_path, f"split_{i}.csv")))

    # test errors
    def test_unsupported_type(self):
        path = os.path.join(self.temp_path, "bad.csv")
        with self.assertRaises(TypeError):
            save_util.save_to_file(path, object())

    def test_unknown_extension(self):
        path = os.path.join(self.temp_path, "file.unknown")
        with self.assertRaises(ValueError):
            save_util.save_to_file(path, self.df)

    def test_no_items(self):
        path = os.path.join(self.temp_path, "empty.csv")
        with self.assertRaises(ValueError):
            save_util.save_to_file(path)
    
    # additional tests for saving to excel file type
    def test_save_dataframe_to_excel(self):
        path = os.path.join(self.temp_path, "df.xlsx")
        save_util.save_to_file(path, self.df)  # filetype inferred from .xlsx
        self.assertTrue(os.path.exists(path))

        df_read = pd.read_excel(path)
        self.assertEqual(len(df_read), len(self.df))
        self.assertListEqual(list(df_read.columns), list(self.df.columns))
    
    def test_save_mixed_items_to_excel(self):
        path = os.path.join(self.temp_path, "mixed.xlsx")
        save_util.save_to_file(path, self.df, self.dict_data, self.pd_series, filetype="excel")
        self.assertTrue(os.path.exists(path))

        df_read = pd.read_excel(path)
        self.assertGreaterEqual(len(df_read), 3)
    
    def test_save_split_to_excel(self):
        path = os.path.join(self.temp_path, "split.xlsx")
        save_util.save_to_file(path, self.df, self.dict_data, combine=False, filetype="excel")
        
        for i in range(2):
            split_path = os.path.join(self.temp_path, f"split_{i}.xlsx")
            self.assertTrue(os.path.exists(split_path))
            df_read = pd.read_excel(split_path)
            self.assertGreaterEqual(len(df_read), 1)

if __name__ == "__main__":
    unittest.main()
