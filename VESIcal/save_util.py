from VESIcal import sample_class, batchfile

import pandas as pd
import os

def save_to_file(filename, *items, filetype=None, combine=True, **kwargs):
    """
    Save one or more items to a CSV or Excel file.

    Parameters
    ----------
    filename : str
        Name of the output file.
        
    *items :
        One or more data objects to save. Each item may be:
        - A pandas DataFrame
        - A pandas Series
        - A dictionary
        - An object with a `to_dict()` method
        - A VESIcal Sample object
        - A VESIcal BatchFile object
        - A list containing any combination of the above
        
    filetype : str, optional
        File format to save as ("csv" or "excel"). If not provided, it will be inferred
        from the filename extension.
    
    combine (bool): If True (default), combine all items into one file. If False, indices will
        be appended to the end of the filename for each item saved out (e.g., filename1.csv,
        filename2.csv...)
        
    **kwargs :
        Additional keyword arguments passed to pandas' `to_csv()` or `to_excel()`.

    Raises
    ------
    TypeError
        If an item is of an unsupported type.
        
    ValueError
        If filetype cannot be determined and is not explicitly provided.
    """

    if not items:
        raise ValueError("No items provided to save.")

    # Infer filetype if not provided
    if filetype is None:
        ext = os.path.splitext(filename)[-1].lower().lstrip(".")
        if ext in {"csv", "tsv"}:
            filetype = "csv"
        elif ext in {"xls", "xlsx"}:
            filetype = "excel"
        else:
            raise ValueError(
                f"Cannot infer filetype from extension '.{ext}'. Please specify filetype."
            )

    if combine:
        df_list = [_coerce_to_dataframe(item) for item in items]
        combined_df = pd.concat(df_list, ignore_index=True)
        _write_df(combined_df, filename, filetype, **kwargs)
    else:
        for i, item in enumerate(items):
            this_filename = _append_index_to_filename(filename, i)
            df = _coerce_to_dataframe(item)
            _write_df(df, this_filename, filetype, **kwargs)

def _coerce_to_dataframe(item):
    """
    Converts a single item or list of items into a DataFrame.
    Returns a DataFrame directly, or raises TypeError if unsupported.
    """
    # Handle list of items recursively
    if isinstance(item, list):
        frames = [_coerce_to_dataframe(subitem) for subitem in item]
        return pd.concat(frames, ignore_index=True)
    
    # Supported single items
    if isinstance(item, batchfile.BatchFile):
        return item.get_data()
    elif isinstance(item, sample_class.Sample):
        return(pd.DataFrame([item.get_composition]))
    elif isinstance(item, pd.DataFrame):
        return item
    elif isinstance(item, pd.Series):
        return pd.DataFrame([item])
    elif isinstance(item, dict):
        return pd.DataFrame([item])
    elif hasattr(item, 'to_dict'):
        return pd.DataFrame([item.to_dict()])
    else:
        raise TypeError(f"Unsupported item type: {type(item)}")

def _write_df(df, filename, filetype, **kwargs):
    if filetype == 'csv':
        df.to_csv(filename, index=False, **kwargs)
    elif filetype == 'excel':
        df.to_excel(filename, index=False, **kwargs)
    else:
        raise ValueError(f"Unsupported filetype: {filetype}")

def _append_index_to_filename(filename, index):
    root, ext = os.path.splitext(filename)
    return f"{root}_{index}{ext}"

def _save_single_item(item, filename, filetype=None, **kwargs):
    # Handle input types
    if isinstance(item, batchfile.BatchFile):
        df = item.get_data()
    elif isinstance(item, pd.DataFrame):
        df = item
    elif isinstance(item, dict):
        df = pd.DataFrame([item])
    elif hasattr(item, 'to_dict'):
        df = pd.DataFrame([item.to_dict()])
    else:
        raise TypeError(f"Unsupported type: {type(item)}")

    # Write to disk
    if filetype == 'csv':
        df.to_csv(filename, index=False, **kwargs)
    elif filetype == 'excel':
        df.to_excel(filename, index=False, **kwargs)
    else:
        raise ValueError(f"Unsupported filetype: {filetype}")

def _append_index_to_filename(filename, index):
    root, ext = os.path.splitext(filename)
    return f"{root}_{index}{ext}"
