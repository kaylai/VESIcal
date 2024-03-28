import pandas as pd
import os
import sys
import warnings as w

from VESIcal import core
from VESIcal import sample_class

# Turn off chained assignment pandas warning
pd.options.mode.chained_assignment = None  # default='warn'


def rename_duplicates(df, suffix='-duplicate-'):
    appendents = (suffix +
                  df.groupby(level=0).cumcount().astype(str).replace('0', ''))
    appendents = appendents.replace(suffix, '')
    return df.set_index(df.index.astype(str) + appendents)


class status_bar(object):
    """Various styles of status bars that display the progress of a calculation
    within a loop
    """
    def __init__():
        pass

    def status_bar(percent, sample_name=None, btext=None, barLen=20):
        """
        Prints an updating status bar to the terminal or jupyter notebook.

        Parameters
        ----------
        percent: float
            Percent value of progress from 0 to 1

        sample_name: string
            Name of the current sample being calculated

        btext: string
            Any extra text to display next to status bar

        barLen: int
            Length of bar to print
        """
        sys.stdout.write("\r")
        sys.stdout.write("[{:<{}}] {:.0f}%".format("=" * int(barLen * percent),
                                                   barLen, percent * 100))

        sample_string = str(sample_name)
        # Set max number of characters in sample name
        max_name_length = 25
        if len(str(sample_name)) >= max_name_length:
            sample_string = str(sample_name)[0:max_name_length-1] + "..."

        # Write out sample name and trailing spaces to cover contents of
        # previous sample names left over on line
        if sample_name is not None:
            sys.stdout.write("  Working on sample " + sample_string +
                             "                            ")
        if btext is not None:
            sys.stdout.write(" " + str(btext))
        if percent == 1.0:
            sys.stdout.write("\n")
        sys.stdout.flush()


# ---------- BATCHFILE CLASS --------- #
class BatchFile(object):
    """A batch file with sample names and oxide compositions

    Attributes
    ----------
        filename: str
            Path to the batch file, e.g., "my_file.xlsx". This always needs to
            be passed, even if the user is passing a pandas DataFrame rather
            than an batch file. If passing a DataFrame, filename should be set
            to None. File can be excel file (.xlsx) or .csv.

        sheet_name: str
            OPTIONAL. For Excel files. Default value is 0 which gets the first
            sheet in the batch spreadsheet file. This implements the pandas.
            read_excel() sheet_name parameter. But functionality to read in
            more than one sheet at a time (e.g.,
            pandas.read_excel(sheet_name=None)) is not yet imlpemented in
            VESIcal. From the pandas 1.0.4 documentation:

            Available cases:
            - Defaults to 0: 1st sheet as a DataFrame
            - 1: 2nd sheet as a DataFrame
            - "Sheet1": Load sheet with name “Sheet1”

        file_type: str
            OPTIONAL. Default is 'excel', which denotes that passed file has
            extension .xlsx. Other option is 'csv', which denotes that the
            passed file has extension .csv.

        units: str
            OPTIONAL. Default is 'wtpt_oxides'. String defining whether the
            oxide composition is given in wt percent ("wtpt_oxides", which is
            the default), mole oxides (mol_oxides) or mole cations
            (mol_cations).

        default_normalization:     None or str
            The type of normalization to apply to the data by default. One of:
            - None (no normalization)
            - 'standard' (default): Normalizes an input composition to 100%.
            - 'fixedvolatiles': Normalizes major element oxides to 100 wt%,
               including volatiles. The volatile wt% will remain fixed, whilst
               the other major element oxides are reduced proportionally so
               that the total is 100 wt%.
            - 'additionalvolatiles': Normalises major element oxide wt% to
              100%, assuming it is volatile-free. If H2O or CO2 are passed to
              the function, their un-normalized values will be retained in
              addition to the normalized non-volatile oxides, summing to >100%.

        default_units     str
            The type of composition to return by default, one of:
            - wtpt_oxides (default)
            - mol_oxides
            - mol_cations

        label: str
            OPTIONAL. Default is 'Label'. Name of the column within the passed
            file referring to sample names.

        dataframe: pandas DataFrame
            OPTIONAL. Default is None in which case this argument is ignored.
            This argument is used when the user wishes to turn a pandas
            DataFrame into an BatchFile object, for example when user data is
            already in python rather than being imported from a file. In this
            case set `dataframe` equal to the dataframe object being passed in.
            If using this option, pass None to filename.
    """
    def __init__(self, filename, sheet_name=0, file_type='excel',
                 units='wtpt_oxides', label='Label',
                 default_normalization='none', default_units='wtpt_oxides',
                 dataframe=None, **kwargs):
        """Return a BatchFile object whose parameters are defined here."""
        self.units = units
        self.set_default_normalization(default_normalization)
        self.set_default_units(default_units)

        if filename is not None:
            file_name, file_extension = os.path.splitext(filename)
            if file_extension == '.xlsx' or file_extension == '.xls':
                file_type = 'excel'
            if file_extension == '.csv':
                file_type = 'csv'

        if isinstance(sheet_name, str) or isinstance(sheet_name, int):
            pass
        else:
            raise core.InputError("If sheet_name is passed, it must be of "
                                  "type str or int. Currently, VESIcal cannot "
                                  "import more than one sheet at a time.")

        # handle data if passed in as existing dataframe or as file
        if dataframe is not None:
            data = dataframe
            if label is not None:
                data = self.try_set_index(data, label)
        else:
            if file_type == 'excel':
                data = pd.read_excel(filename, sheet_name=sheet_name)
                data = self.try_set_index(data, label)
            elif file_type == 'csv':
                data = pd.read_csv(filename)
                data = self.try_set_index(data, label)
            else:
                raise core.InputError("file_type must be one of \'excel\' or "
                                      "\'csv\'.")

        # Sanitize data inputs
        data = rename_duplicates(data)  # handle any duplicated sample names
        for column in data:  # convert all oxide columns to numeric
            if column in core.oxides:
                data[column] = data[column].apply(pd.to_numeric, errors='coerce')
        data = data.dropna(how='all')  # drop any rows that are all NaNs
        data = data.fillna(0)  # fill in any missing data with 0's

        if 'model' in kwargs:
            w.warn("You don't need to pass a model here, so it will be "
                   "ignored. You can specify a model when performing "
                   "calculations on your dataset (e.g., "
                   "calculate_dissolved_volatiles())",
                   RuntimeWarning, stacklevel=2)

        if 'norm' in kwargs:
            w.warn("We noticed you passed a norm argument here. This does "
                   "nothing. You can normalize your BatchFile and save it to "
                   "a new variable name after import using "
                   "normalize(BatchFileObject). See the documentation for "
                   "more info.",
                   RuntimeWarning, stacklevel=2)

        total_iron_columns = ["FeOt", "FeOT", "FeOtot", "FeOtotal", "FeOstar",
                              "FeO*"]
        for name in total_iron_columns:
            if name in data.columns:
                if 'FeO' in data.columns:
                    for row in data.itertuples():
                        if (data.at[row.Index, "FeO"] == 0 and
                           data.at[row.Index, name] > 0):
                            w.warn("Sample " + str(row.Index) + ": " +
                                   str(name) + " value of " +
                                   str(data.at[row.Index, name]) +
                                   " used as FeO. Fe2O3 set to 0.0.",
                                   RuntimeWarning, stacklevel=2)
                            data.at[row.Index, "Fe2O3"] = 0.0
                            data.at[row.Index, "FeO"] = (
                                                      data.at[row.Index, name])
                else:
                    w.warn("Total iron column " + str(name) + " detected. " +
                           "This column will be treated as FeO. If Fe2O3 " +
                           "data are not given, Fe2O3 will be 0.0. In " +
                           "future, an option to calcualte FeO/Fe2O3 based " +
                           "on fO2 will be implemented.",
                           RuntimeWarning, stacklevel=2)
                    data['FeO'] = data[name]

        if units == "wtpt_oxides":
            pass
        if units == "mol_oxides":
            data = self._molOxides_to_wtpercentOxides(data)
        if units == "mol_cations":
            data = self._molCations_to_wtpercentOxides(data)

        for column in data:
            if column in core.oxides:
                data[column][data[column] < 0] = 0

        self.data = data

    def set_default_normalization(self, default_normalization):
        """ Set the default type of normalization to use with the
        get_composition() method.

        Parameters
        ----------
        default_normalization:    str
            The type of normalization to apply to the data. One of:
            - 'none' (no normalization)
            - 'standard' (default): Normalizes an input composition to 100%.
            - 'fixedvolatiles': Normalizes major element oxides to 100 wt%,
              including volatiles. The volatile wt% will remain fixed, whilst
              the other major element oxides are reduced proportionally so
              that the total is 100 wt%.
            - 'additionalvolatiles': Normalises major element oxide wt% to
              100%, assuming it is volatile-free. If H2O or CO2 are passed to
              the function, their un-normalized values will be retained in
              addition to the normalized non-volatile oxides, summing to >100%.

        """
        if default_normalization in ['none', 'standard', 'fixedvolatiles',
                                     'additionalvolatiles']:
            self.default_normalization = default_normalization
        else:
            raise core.InputError("The normalization method must be one of "
                                  "'none', 'standard', 'fixedvolatiles' "
                                  "or 'additionalvolatiles'.")

    def set_default_units(self, default_units):
        """ Set the default units of composition to return when using the
        get_composition() method.

        Parameters
        ----------
        default_units     str
            The type of composition to return, one of:
            - wtpt_oxides (default)
            - mol_oxides
            - mol_cations
        """
        if default_units in ['wtpt_oxides', 'mol_oxides', 'mol_cations']:
            self.default_units = default_units
        else:
            raise core.InputError("The units must be one of 'wtpt_oxides', "
                                  "'mol_oxides','mol_cations'.")

    def get_composition(self, species=None, normalization=None, units=None,
                        exclude_volatiles=False, asBatchFile=False):
        """ Returns a pandas DataFrame containing the compositional
        information for all samples in the BatchFile object

        Parameters
        ----------
        species:    NoneType or str
            The name of the oxide or cation to return the concentration of. If
            NoneType (default) the whole composition of each sample will be
            returned. If an oxide is passed, the value in wtpt will be
            returned unless units is set to 'mol_oxides', even if the default
            units for the sample object are mol_oxides. If an element is
            passed, the concentration will be returned as mol_cations, unless
            'mol_singleO' is specified as units, even if the default units for
            the sample object are mol_singleO. Unless normalization is
            specified in the method call, none will be applied.

        normalization:     NoneType or str
            The type of normalization to apply to the data. One of:
            - 'none' (no normalization)
            - 'standard' (default): Normalizes an input composition to 100%.
            - 'fixedvolatiles': Normalizes major element oxides to 100 wt%,
              including volatiles. The volatile wt% will remain fixed, whilst
              the other major element oxides are reduced proportionally so
              that the total is 100 wt%.
            - 'additionalvolatiles': Normalises major element oxide wt% to
              100%, assuming it is volatile-free. If H2O or CO2 are passed to
              the function, their un-normalized values will be retained in
              addition to the normalized non-volatile oxides, summing to >100%.
            If NoneType is passed the default normalization option will be
            used (self.default_normalization).

        units:     NoneType or str
            The units of composition to return, one of:
            - wtpt_oxides (default)
            - mol_oxides
            - mol_cations
            - mol_singleO
            If NoneType is passed the default units option will be used
            (self.default_type).

        exclude_volatiles   bool
            If True, volatiles will be excluded from the returned composition,
            prior to normalization and conversion.

        asBatchFile:    bool
            If True, returns a BatchFile object. If False, returns a
            pandas.DataFrame object.

        Returns
        -------
        pandas.DataFrame or BatchFile object
            All sample information.
        """
        data = self.data.copy()

        # Fetch the default return types if not specified in function call
        if normalization is None and species is None:
            normalization = self.default_normalization
        if units is None and species is None:
            units = self.default_units

        new_compositions = []
        sample_names = []
        for index, row in data.iterrows():
            sample_comp = self.get_sample_composition(index, units=units,
                                                      asSampleClass=True)
            new_compositions.append(sample_comp.get_composition(
                     species=species, normalization=normalization, units=units,
                     exclude_volatiles=exclude_volatiles))
            sample_names.append(index)
        if isinstance(new_compositions[0], pd.Series):
            return_frame = pd.concat(
                           [pd.DataFrame(j) for j in new_compositions], axis=1)
            return_frame = return_frame.transpose()
            return_frame["new_index"] = sample_names
            return_frame = return_frame.set_index("new_index")
            return_frame.index.name = None
        elif isinstance(new_compositions[0], float):
            species_data = {species: new_compositions}
            return_frame = pd.DataFrame(
                           species_data, index=[name for name in sample_names])
        else:
            return_frame = None

        if asBatchFile is False:
            return return_frame
        else:
            return BatchFile(filename=None, dataframe=return_frame, label=None)

    def get_data(self, normalization=None, units=None, asBatchFile=False):
        """
        Returns all data stored in a BatchFile object (both compositional and
        other data). To return only the compositional data, use
        get_composition().

        Parameters
        ----------
        normalization:     NoneType or str
            The type of normalization to apply to the data. One of:
            - 'none' (no normalization)
            - 'standard' (default): Normalizes an input composition to 100%.
            - 'fixedvolatiles': Normalizes major element oxides to 100 wt%,
              including volatiles. The volatile wt% will remain fixed, whilst
              the other major element oxides are reduced proportionally so
              that the total is 100 wt%.
            - 'additionalvolatiles': Normalises major element oxide wt% to
              100%, assuming it is volatile-free. If H2O or CO2 are passed to
              the function, their un-normalized values will be retained in
              addition to the normalized non-volatile oxides, summing to >100%.

            If NoneType is passed the default normalization option will be
            used (self.default_normalization).

        units:     NoneType or str
            The units of composition to return, one of:
            - wtpt_oxides (default)
            - mol_oxides
            - mol_cations
            - mol_singleO

            If NoneType is passed the default units option will be used
            (self.default_type).

        asBatchFile:    bool
            If True, returns a BatchFile object. If False, returns a
            pandas.DataFrame object.

        Returns
        -------
        pandas.DataFrame or BatchFile object
            All sample information.
        """
        data = self.data.copy()

        # Fetch the default return units if not specified in function call
        if units is None:
            units = self.default_units

        # Fetch the default normalization if not specified in the function call
        if normalization is None:
            normalization = self.default_normalization

        # Grab all compositional data
        compositional_data = self.get_composition(normalization=normalization,
                                                  units=units)

        # Grab all non-compositional data
        non_compositional_data = data.filter(
                       [col for col in data.columns if col not in core.oxides])

        # concatenate both compositional and non-compositional dataframes
        # into one
        return_frame = pd.concat([compositional_data, non_compositional_data],
                                 axis=1)

        if asBatchFile is False:
            return return_frame
        else:
            return BatchFile(filename=None, dataframe=return_frame, label=None)

    def get_sample_composition(self, samplename, species=None,
                               normalization=None, units=None,
                               asSampleClass=False):
        """
        Returns oxide composition of a single sample from a user-imported file
        as a dictionary

        Parameters
        ----------
        samplename: string
            Name of the desired sample

        normalization: NoneType or str
            The type of normalization to apply to the data. One of:
            - 'none' (no normalization)
            - 'standard' (default): Normalizes an input composition to 100%.
            - 'fixedvolatiles': Normalizes major element oxides to 100 wt%,
              including volatiles. The volatile wt% will remain fixed, whilst
              the other major element oxides are reduced proportionally so
              that the total is 100 wt%.
            - 'additionalvolatiles': Normalises major element oxide wt% to
              100%, assuming it is volatile-free. If H2O or CO2 are passed to
              the function, their un-normalized values will be retained in
              addition to the normalized non-volatile oxides, summing to >100%.

            If NoneType is passed the default normalization option will be
            used (self.default_normalization).

        units:     NoneType or str
            The units of composition to return, one of:
            - wtpt_oxides (default)
            - mol_oxides
            - mol_cations
            - mol_singleO
            If NoneType is passed the default units option will be used
            (self.default_type).

        asSampleClass:  bool
            If True, the sample composition will be returned as a sample
            class, with default options. In this case any normalization
            instructions will be ignored.

        Returns
        -------
        dictionary, float, or sample_class.Sample object
            Composition of the sample as oxides
        """
        # Fetch the default return types if not specified in function call
        if normalization is None and species is None:
            normalization = self.default_normalization
        if units is None and species is None:
            units = self.default_units

        # Check that normalization being chosen is one of the possible options
        if normalization in [None, 'none', 'standard', 'fixedvolatiles',
                             'additionalvolatiles']:
            pass
        else:
            raise core.InputError("The normalization method must be one of "
                                  "'none', 'standard', 'fixedvolatiles', "
                                  "or 'additionalvolatiles'.")

        data = self.data
        my_sample = pd.DataFrame(data.loc[samplename])
        sample_dict = (my_sample.to_dict()[samplename])
        sample_oxides = {}
        for item, value in sample_dict.items():
            if item in core.oxides:
                sample_oxides.update({item: value})

        _sample = sample_class.Sample(sample_oxides)

        # Get sample composition in terms of any species, units, and
        # normalization passed
        return_sample = _sample.get_composition(species=species, units=units,
                                                normalization=normalization)

        if asSampleClass:
            return sample_class.Sample(return_sample)
        else:
            if species is None:
                return dict(return_sample)
            elif isinstance(species, str):
                return return_sample

    def _molOxides_to_wtpercentOxides(self, data):
        for i, row in data.iterrows():
            sample_comp = {}
            for oxide in core.oxides:
                if oxide in data.columns:
                    sample_comp[oxide] = row[oxide]
                else:
                    sample_comp[oxide] = 0.0
            _sample = sample_class.Sample(sample_comp, units='mol_oxides')
            _sample_conv = _sample.get_composition()
            for ox in core.oxides:
                data.loc[i, oxide] = _sample_conv[oxide]
        return data

    def _molCations_to_wtpercentOxides(self, data):
        for i, row in data.iterrows():
            sample_comp = {}
            for cation in core.oxides_to_cations[core.oxides]:
                if cation in data.columns:
                    sample_comp[cation] = row[cation]
                else:
                    sample_comp[cation] = 0.0
            _sample = sample_class.Sample(sample_comp, units='mol_cations')
            _sample_conv = _sample.get_composition()
            for oxide in core.oxides:
                data.loc[i, oxide] = _sample_conv[oxide]
        return data

    def try_set_index(self, dataframe, label):
        """
        Method to handle setting the index column in an BatchFile object. If
        no column is passed that matches the default index name, then this
        method will attempt to choose the 'best' column that the user might
        want to serve as an index column.

        Parameters
        ----------
        dataframe: pandas DataFrame

        label: str
            Name of the column within the passed Excel file referring to
            sample names.
        """
        _dataframe = dataframe.copy()
        try:
            _dataframe = _dataframe.set_index(label)
        except Exception:
            label_found = False
            for col in _dataframe.columns:
                if col in core.oxides:
                    pass
                else:
                    _dataframe = _dataframe.set_index(col)
                    label_found = True
                    w.warn("No Label column given, so column '" + str(col) +
                           "' was chosen for you. To choose your own, set " +
                           "label='<column-name>'.", RuntimeWarning,
                           stacklevel=2)
                    break
            if label_found is False:
                _dataframe.index.name = 'Label'
                w.warn("No Label column given, so one was created for you. "
                       "To choose your own, set label='<column-name>'.",
                       RuntimeWarning, stacklevel=2)

        return _dataframe

    def save_excel(self, filename, calculations, sheet_names=None):
        """
        Saves data calculated by the user in batch processing mode (using the
        BatchFile class methods) to an organized Excel file, with the original
        user data plus any calculated data.

        Parameters
        ----------
        filename: string
            Name of the file. Extension (.xlsx) should be passed along with
            the name itself, all in quotes (e.g., 'myfile.xlsx').

        calculations: pandas DataFrame or list of pandas DataFrames
            A single DataFrame or list of DataFrames (e.g., calculated outputs
            from any of the core BatchFile functions:
            calculate_dissolved_volatiles, calculate_equilibrium_fluid_comp,
            and calculate_saturation_pressure). If None, only the original
            user data will be saved.

        sheet_names: None, string, or list
            OPTIONAL. Default value is None. Allows user to set the name of
            the sheet or sheets written to the Excel file.

        Returns
        -------
            Creates and saves an Excel file with data from each calculation
            saved to its own sheet.
        """
        if isinstance(calculations, list):
            if isinstance(sheet_names, list) or sheet_names is None:
                pass
            else:
                raise core.InputError("If calculations is passed as list, "
                                      "sheet_names must also be list of same "
                                      "length")
        elif calculations is None:
            pass
        else:
            calculations = [calculations]

        with pd.ExcelWriter(filename) as writer:
            self.data.to_excel(writer, 'Original_User_Data')
            if isinstance(calculations, list):
                if sheet_names is None:
                    for n, df in enumerate(calculations):
                        df.to_excel(writer, 'Calc%s' % n)
                elif isinstance(sheet_names, list):
                    pass
                else:
                    sheet_names = [sheet_names]
                if isinstance(sheet_names, list):
                    if len(sheet_names) == len(calculations):
                        pass
                    else:
                        raise core.InputError("calculations and sheet_names "
                                              "must have the same length")

                    for i in range(len(calculations)):
                        if isinstance(sheet_names[i], str):
                            calculations[i].to_excel(writer, sheet_names[i])
                        else:
                            raise core.InputError("if sheet_names is passed, "
                                                  "it must be list of strings")
            elif calculations is None:
                pass
        return print("Saved " + str(filename))

    def save_csv(self, filenames, calculations, **kwargs):
        """
        Saves data calculated by the user in batch processing mode to a
        comma-separated values (csv) file. Mirros the pandas.to_csv() method.
        Any argument that can be passed to pandas.csv() can be passed here.
        One csv file will be saved for each calculation passed.

        Parameters
        ----------
        filenames: string or list of strings
            Name of the file. Extension (.csv) should be passed along with
            the name itself, all in quotes (e.g., 'myfile.csv'). The number
            of calculations passed must match the number of filenames passed.
            If passing more than one, should be passed as a list.

        calculations: pandas DataFrame or list of pandas DataFrames
            A single variable or list of variables containing calculated
            outputs from any of the core BatchFile functions:
            calculate_dissolved_volatiles, calculate_equilibrium_fluid_comp,
            and calculate_saturation_pressure.

        Returns
        -------
            Creates and saves a CSV file or files with data from each
            calculation saved to its own file.
        """
        if isinstance(filenames, list) is False:
            filenames = [filenames]
        if isinstance(calculations, list) is False:
            calculations = [calculations]
        if len(filenames) != len(calculations):
            raise core.InputError("calculations and filenames must have the "
                                  "same length")

        for i in range(len(filenames)):
            calculations[i].to_csv(filenames[i], **kwargs)
            print("Saved " + str(filenames[i]))


def from_DataFrame(dataframe, units='wtpt_oxides', label='Label'):
    """
    Transforms any pandas DataFrame object into a VESIcal BatchFile object.

    Parameters
    ----------
    dataframe: pd.DataFrame object
        DataFrame object containing samples and oxide compositions.

    units: str
        OPTIONAL. Default is 'wtpt_oxides'. String defining whether the oxide
        composition is given in wt percent ("wtpt_oxides", which is the
        default), mole fraction oxides ("mol_oxides"), or mole fraction
        cations ("mol_cations").

    label: str
        OPTIONAL. Default is 'Label'. Name of the column within the passed
        file referring to sample names. This column will be set as the index
        column.

    Returns
    -------
    VESIcal.BatchFile object
    """
    return BatchFile(filename=None, dataframe=dataframe, units=units,
                     label=label)
