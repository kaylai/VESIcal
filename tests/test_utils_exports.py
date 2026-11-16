import csv
from pathlib import Path
import tempfile
import unittest

from VESIcal.utils import save_results
from VESIcal.utils.save_util import save_results as implementation


def wildcard_exports():
    namespace = {}
    exec("from VESIcal.utils import *", namespace)
    return {name: value for name, value in namespace.items() if name != "__builtins__"}


def csv_rows(path):
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.reader(stream))


class TestUtilsExports(unittest.TestCase):
    def test_named_import_preserves_callable_and_real_csv_output(self):
        self.assertIs(save_results, implementation)
        data = {"x": 1, "tag": "sample"}
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "named"
            save_results(str(path), data)
            self.assertEqual(csv_rows(path.with_suffix(".csv")),
                             [["x", "tag"], ["1", "sample"]])
        self.assertEqual(data, {"x": 1, "tag": "sample"})

    def test_wildcard_exports_only_the_complete_public_callable(self):
        namespace = wildcard_exports()
        self.assertEqual(set(namespace), {"save_results"})
        self.assertIs(namespace["save_results"], implementation)

    def test_wildcard_callable_writes_the_same_real_csv_as_named_import(self):
        namespace = wildcard_exports()
        data = {"x": 1, "tag": "sample"}
        with tempfile.TemporaryDirectory() as directory:
            named_path = Path(directory) / "named.csv"
            wildcard_path = Path(directory) / "wildcard.csv"
            save_results(str(named_path), data)
            namespace["save_results"](str(wildcard_path), data,
                                      filetype="csv", mode="single_sheet")
            self.assertEqual(csv_rows(named_path), [["x", "tag"], ["1", "sample"]])
            self.assertEqual(csv_rows(wildcard_path), [["x", "tag"], ["1", "sample"]])
            self.assertEqual(wildcard_path.read_bytes(), named_path.read_bytes())
        self.assertEqual(data, {"x": 1, "tag": "sample"})


if __name__ == "__main__":
    unittest.main()
