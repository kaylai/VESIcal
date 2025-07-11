# Import a Data file

```{contents}
```

For any code using the VESIcal library, the library must be imported for use. Here we import VESIcal as `v`. Any time we wish to use a function from VESIcal, that function must be preceded by v.. Specific examples of this usage follow. Here we also import some other python libraries that we will be using in the worked examples below.

```python
import VESIcal as v
```

## Import your excel file and view it

You can import an excel or csv file containing compositional data describing your samples using the `BatchFile` class. Your file should have each sample in a separate row, with data in terms of oxides. You can find the 'example_data.xlsx'\`' file in the 'manuscript' folder in the github repository.

```python
myfile = v.BatchFile('example_data.xlsx')
```

Once the BatchFile object is created and assigned to a variable, the user can then access the data loaded from their file as `variable.get_data()`. In this example, the variable corresponding to the `BatchFile` object is named `myfile` and so the data in that file can be accessed with `myfile.get_data()`. Below, `myfile.get_data()` is saved to a variable we name `data`. The variable `data` is a pandas DataFrame object, which makes displaying the data itself quite simple and aesthetically pleasing, since pandas DataFrames mimic spreadsheets.

```python
data = myfile.get_data()
data

|           Label    | SiO2  | TiO2   | Al2O3 | Fe2O3 | Cr2O3 | FeO    | MnO    | MgO    | NiO | CoO | CaO    | Na2O | K2O  | P2O5   | H2O      | CO2      | Press | Temp |
|--------------------|-------|--------|-------|-------|-------|--------|--------|--------|-----|-----|--------|------|------|--------|----------|----------|-------|------|
| BT-ex              | 77.50 | 0.0800 | 12.50 | 0.207 | 0     | 0.4730 | 0.0000 | 0.0300 | 0   | 0   | 0.4300 | 3.98 | 4.88 | 0.0000 | 5.500000 | 0.050000 | 500   | 900  |
| TVZMa-ex           | 78.37 | 0.1300 | 11.94 | 0.000 | 0     | 0.9900 | 0.0400 | 0.0500 | 0   | 0   | 0.5300 | 3.80 | 4.14 | 0.0000 | 4.060000 | 0.005000 | 600   | 800  |
| TVZOh-ex           | 77.90 | 0.0800 | 12.15 | 0.000 | 0     | 0.9500 | 0.0500 | 0.0600 | 0   | 0   | 0.5500 | 4.05 | 4.12 | 0.0000 | 4.630000 | 0.005000 | 50    | 900  |
| Oh48-FTIR1-MI1-a   | 78.27 | 0.0298 | 12.02 | 0.000 | 0     | 0.9828 | 0.0336 | 0.0515 | 0   | 0   | 0.4772 | 4.05 | 4.09 | 0.0000 | 4.214912 | 0.004566 | 250   | 950  |
| Oh48-FTIR1-MI1-b   | 78.27 | 0.0298 | 12.02 | 0.000 | 0     | 0.9828 | 0.0336 | 0.0515 | 0   | 0   | 0.4772 | 4.05 | 4.09 | 0.0000 | 4.005816 | 0.004448 | 500   | 1025 |
| Oh48-FTIR1-MI1-IRc | 78.27 | 0.0298 | 12.02 | 0.000 | 0     | 0.9828 | 0.0336 | 0.0515 | 0   | 0   | 0.4772 | 4.05 | 4.09 | 0.0000 | 3.885649 | 0.004654 | 5000  | 925  |
| Oh50-4.1           | 77.91 | 0.0984 | 12.07 | 0.000 | 0     | 1.0556 | 0.0257 | 0.0999 | 0   | 0   | 0.5216 | 4.04 | 4.18 | 0.0000 | 4.641843 | 0.004566 | 1000  | 862  |
| Oh50-4.2           | 77.91 | 0.0984 | 12.07 | 0.000 | 0     | 1.0556 | 0.0257 | 0.0999 | 0   | 0   | 0.5216 | 4.04 | 4.18 | 0.0000 | 4.402133 | 0.004448 | 100   | 770  |
| Oh49-4.1           | 77.92 | 0.0099 | 12.11 | 0.000 | 0     | 1.0020 | 0.0672 | 0.0546 | 0   | 0   | 0.5346 | 4.01 | 4.30 | 0.0000 | 4.283934 | 0.004566 | 1000  | 855  |
| Oh49-4.2           | 77.92 | 0.0099 | 12.11 | 0.000 | 0     | 1.0020 | 0.0672 | 0.0546 | 0   | 0   | 0.5346 | 4.01 | 4.30 | 0.0000 | 4.230533 | 0.004448 | 500   | 1000 |
| Ma55-5a.1          | 77.68 | 0.0096 | 12.27 | 0.000 | 0     | 1.0272 | 0.0628 | 0.0342 | 0   | 0   | 0.6064 | 3.97 | 4.35 | 0.0000 | 4.459767 | 0.004654 | 5000  | 1010 |
| Ma57-3b.2          | 77.90 | 0.0498 | 12.07 | 0.000 | 0     | 1.0844 | 0.0748 | 0.0355 | 0   | 0   | 0.4759 | 4.10 | 4.21 | 0.0000 | 3.712506 | 0.004448 | 1000  | 1012 |
| Ma57-3c.1          | 77.65 | 0.1590 | 12.28 | 0.000 | 0     | 0.9769 | 0.0597 | 0.0577 | 0   | 0   | 0.5598 | 4.08 | 4.18 | 0.0064 | 4.443973 | 0.004654 | 100   | 885  |
| Ma57-3c.2          | 77.65 | 0.1590 | 12.28 | 0.000 | 0     | 0.9769 | 0.0597 | 0.0577 | 0   | 0   | 0.5598 | 4.08 | 4.18 | 0.0064 | 4.283171 | 0.004645 | 1000  | 885  |
```

## Extract and view a single sample

Defined within the BatchFile() class, the method get_sample_composition() allows for the extraction of a melt composition from a loaded file.

```python
SampleName = 'BT-ex'
extracted_bulk_comp = myfile.get_sample_composition(SampleName)
extracted_bulk_comp

{'SiO2': 77.5,
 'TiO2': 0.08,
 'Al2O3': 12.5,
 'Fe2O3': 0.207,
 'Cr2O3': 0.0,
 'FeO': 0.473,
 'MnO': 0.0,
 'MgO': 0.03,
 'NiO': 0.0,
 'CoO': 0.0,
 'CaO': 0.43,
 'Na2O': 3.98,
 'K2O': 4.88,
 'P2O5': 0.0,
 'H2O': 5.5,
 'CO2': 0.05}
```

In order to use this sample for VESIcal calculations, you'll need to create a Sample object to hold this compositional data. You can extract a single sample as a Sample object by changing the above example to:

```python
SampleName = 'BT-ex'
extracted_bulk_comp = myfile.get_sample_composition(SampleName, asSampleClass=True)
```

You can create any Sample object from a dict or pandas Series as:

```python
mysample = v.Sample({'SiO2': 77.5,
 'TiO2': 0.08,
 'Al2O3': 12.5,
 'Fe2O3': 0.207,
 'Cr2O3': 0.0,
 'FeO': 0.473,
 'MnO': 0.0,
 'MgO': 0.03,
 'NiO': 0.0,
 'CoO': 0.0,
 'CaO': 0.43,
 'Na2O': 3.98,
 'K2O': 4.88,
 'P2O5': 0.0,
 'H2O': 5.5,
 'CO2': 0.05})
```

or simply

```python
mysample = v.Sample(extracted_bulk_comp)
```

if extracted_bulk_comp is a dict or pandas Series.
