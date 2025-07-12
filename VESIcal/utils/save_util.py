from VESIcal import sample_class, batchfile, calculate_classes

import pandas as pd
import numpy as np
import os

def flatten(obj):
    """
    Recursively converts a wide range of nested or mixed Python objects into a 
    list of pandas DataFrames, each representing one logical data unit.

    Supported input types include:
    - Scalars (int, float, str, and NumPy scalars like np.float64)
    - dicts and lists/tuples of dicts
    - pandas.DataFrame and pandas.Series
    - VESIcal objects with data access methods:
        - Calculate: uses `.result`
        - BatchFile: uses `.get_data()`
        - Sample: uses `.get_composition()`

    Parameters
    ----------
    obj : Any
        A single object or container of objects (e.g. list of dicts, list of Samples)
        to be flattened into one or more DataFrames.

    Returns
    -------
    list of pandas.DataFrame
        A list where each item is a flattened DataFrame representation of a top-level
        object in the input. If `obj` is not a list or tuple, the result is a single-item list.

    Raises
    ------
    TypeError
        If an unsupported object type is encountered during flattening.

    Examples
    --------
    >>> flatten({"A": 1, "B": 2})
    [   A  B
     0  1  2]

    >>> flatten([Sample(...), Calculate(...), {"x": 3}])
    [ DataFrame from Sample,
      DataFrame from Calculate,
      DataFrame from dict ]
    """

    def to_dataframe(item):
        if isinstance(item, calculate_classes.Calculate):
            return to_dataframe(item.result)

        elif isinstance(item, batchfile.BatchFile):
            return item.get_data().reset_index(drop=True)

        elif isinstance(item, sample_class.Sample):
            return to_dataframe(item.get_composition())

        elif isinstance(item, pd.DataFrame):
            return item.reset_index(drop=True)

        elif isinstance(item, pd.Series):
            return pd.DataFrame([item.to_dict()])

        elif isinstance(item, dict):
            return pd.DataFrame([item])

        elif isinstance(item, (list, tuple)):
            if all(isinstance(x, dict) for x in item):
                return pd.DataFrame(item)
            else:
                return pd.DataFrame([{"value": x} for x in item])

        elif isinstance(item, (int, float, str, np.generic)):
            return pd.DataFrame([{"value": item}])

        else:
            raise TypeError(f"Unsupported type in flatten: {type(item)}")

    if isinstance(obj, (list, tuple)):
        return [to_dataframe(item) for item in obj]
    return [to_dataframe(obj)]

def save_results(filename, obj, filetype="csv", mode="single_sheet", descriptions=None):
    """
    Save one or more data objects to CSV or Excel files with flexible formatting options.

    This function supports a wide range of input types, including scalars, dictionaries,
    pandas objects, and custom scientific objects (e.g., VESIcal Calculate, Sample, and BatchFile).
    It flattens each object into a pandas DataFrame before saving, and offers options for
    saving to a single file, multiple sheets, or separate files. Metadata can optionally be
    added to track the origin or meaning of each data entry.

    Parameters
    ----------
    filename : str
        Base filename to use for saving. The appropriate file extension will be added or replaced
        based on the selected `filetype`. When `mode='multi_file'`, index-based suffixes are appended
        automatically (e.g., `output_1.csv`, `output_2.csv`).

    obj : object or list/tuple of objects
        The data to save. Each object must be one of the supported types:
        - `pandas.DataFrame`, `pandas.Series`
        - `dict`, `list` of dicts
        - `int`, `float`, `str`, or any `numpy.generic` scalar (e.g., `np.float64`)
        - `Calculate` (accesses `.result`)
        - `Sample` (accesses `.get_composition()`)
        - `BatchFile` (accesses `.get_data()`)

    filetype : str, optional
        Output file format. Options are:
        - `"csv"`: Save as CSV file(s)
        - `"excel"`: Save as Excel file(s) (`.xlsx`)
        Default is `"csv"`.

    mode : str, optional
        Controls how multiple data objects are written. Options are:
        - `"single_sheet"`: All data combined into a single file and sheet (CSV or Excel)
        - `"multi_sheet"`: One Excel file with each item on its own sheet (Excel only)
        - `"multi_file"`: Each item saved to its own file (suffixes are added automatically)
        Default is `"single_sheet"`.

    descriptions : list of str, optional
        Optional metadata labels to add to each DataFrame as a `"Description"` column.
        Must match the number of top-level items in `obj` if `obj` is a list or tuple.
        Ignored if `obj` is a single object.

    Raises
    ------
    ValueError
        If `mode` is incompatible with the chosen `filetype`, or if the number of descriptions
        does not match the number of objects being saved.

    Examples
    --------
    >>> save_results("output.csv", [df1, df2], filetype="csv", mode="multi_file")

    >>> save_results("summary.xlsx", [sample1, calc], filetype="excel",
    ...              mode="multi_sheet", descriptions=["Basalt Sample", "Solubility Result"])
    """

    dfs = flatten(obj)

    if descriptions:
        if len(descriptions) != len(dfs):
            raise ValueError("Length of descriptions must match number of items.")
        for df, desc in zip(dfs, descriptions):
            df.insert(0, "Description", desc)

    if filetype == "csv":
        if mode == "single_sheet":
            combined = pd.concat(dfs, ignore_index=True)
            outname = filename if filename.endswith(".csv") else filename + ".csv"
            combined.to_csv(outname, index=False)
        elif mode == "multi_file":
            for i, df in enumerate(dfs):
                outname = f"{os.path.splitext(filename)[0]}_{i+1}.csv"
                df.to_csv(outname, index=False)
        else:
            raise ValueError("multi_sheet mode not supported for CSV. Use 'single_sheet' or 'multi_file'.")

    elif filetype == "excel":
        outname = filename if filename.endswith(".xlsx") else filename + ".xlsx"

        if mode == "single_sheet":
            combined = pd.concat(dfs, ignore_index=True)
            with pd.ExcelWriter(outname) as writer:
                combined.to_excel(writer, index=False, sheet_name="Sheet1")

        elif mode == "multi_sheet":
            with pd.ExcelWriter(outname) as writer:
                for i, df in enumerate(dfs):
                    sheetname = f"Item_{i+1}"
                    df.to_excel(writer, index=False, sheet_name=sheetname)

        elif mode == "multi_file":
            for i, df in enumerate(dfs):
                fname = f"{os.path.splitext(filename)[0]}_{i+1}.xlsx"
                df.to_excel(fname, index=False)
        else:
            raise ValueError("Invalid mode for Excel export. Choose 'single_sheet', 'multi_sheet', or 'multi_file'.")

    else:
        raise ValueError("Unsupported filetype. Choose 'csv' or 'excel'.")
