# Normalizing Compositions

```{contents}
```

By default, VESIcal assumes your data are input in terms of wt% oxides and applies no normalization to your data. You may wish to normalize your dataset (using one of VESIcal's three normalization routines) after import, translate your wt% oxide data into some other units (mol fraction oxides or cations), or you may with to import data already in terms of mol fraction oxides or cations (in which case, you need to inform VESIcal of this, otherwise it will assume the values are in wt% oxides).

To normalize a dataset upon import, use the `default_normalizaion` argument when creating your Sample or BatchFile object:

```python
my_sample = v.Sample({'SiO2':  77.3,
 'TiO2':   0.08,
 'Al2O3': 12.6,
 'Fe2O3':  0.207,
 'Cr2O3':  0.0,
 'FeO':    0.473,
 'MnO':    0.0,
 'MgO':    0.03,
 'NiO':    0.0,
 'CoO':    0.0,
 'CaO':    0.43,
 'Na2O':   3.98,
 'K2O':    4.88,
 'P2O5':   0.0,
 'H2O':    6.5,
 'CO2':    0.05},
 default_normalization='standard')
```

```python
myfile = v.BatchFile('path/to/your/file.xlsx.xlsx', default_normalization='standard')
```

To convert units from wt% oxides to something else (in this example, mol fraction oxides) upon import, use the `default_units` argument when creating your Sample or BatchFile object:

```python
my_sample = v.Sample({'SiO2':  77.3,
 'TiO2':   0.08,
 'Al2O3': 12.6,
 'Fe2O3':  0.207,
 'Cr2O3':  0.0,
 'FeO':    0.473,
 'MnO':    0.0,
 'MgO':    0.03,
 'NiO':    0.0,
 'CoO':    0.0,
 'CaO':    0.43,
 'Na2O':   3.98,
 'K2O':    4.88,
 'P2O5':   0.0,
 'H2O':    6.5,
 'CO2':    0.05},
 default_units='mol_oxides')
```

```python
myfile = v.BatchFile('path/to/your/file.xlsx.xlsx', default_units='mol_oxides')
```

To instruct VESIcal that you are inputting your data in terms of units other than wt% oxides (here mol fraction oxidxes), use the `units` argument when creating your Sample object:

```python
my_sample = v.Sample({'SiO2':  0.67,
 'TiO2':   0.00053,
 'Al2O3':  0.065,
 'Fe2O3':  0.00068,
 'Cr2O3':  0.0,
 'FeO':    0.0035,
 'MnO':    0.0,
 'MgO':    0.00039,
 'NiO':    0.0,
 'CoO':    0.0,
 'CaO':    0.0040,
 'Na2O':   0.0337,
 'K2O':    0.0272,
 'P2O5':   0.0,
 'H2O':    0.189,
 'CO2':    0.0006},
 units='mol_oxides')
```

```python
myfile = v.BatchFile('path/to/your/file.xlsx.xlsx', units='mol_oxides')
```

Note that, by default, your sample composition(s) will be returned to you in wt% oxides unless you also specify `default_units='moloxides'`.

## Specific Normalization Types

Before performing model calculations on a dataset, it may be desired to normalize the input composition(s) to a total of 100%. VESIcal has multiple built-in methods for doing so. It should be noted that this procedure is by no means required and not necessarily advised depending on what the user intends to model.

In some cases, data transformations internal to model calculations (e.g., converting between wt% and mol fraction) in effect cause normalization of the input bulk composition anyways, and so normalizing ahead of time will make no difference in the final modeled result. For example, `calculate_dissolved_volatiles` is agnostic to any a priori normalization of the data since the volatiles are handled separately from the dry bulk. On the other hand, `calculate_saturation_pressure` depends very much on any normalization performed, since the calculated pressure depends directly and strongly on the proportion of volatiles in the bulk composition.

To normalize your dataset upon import, please see the section above. This section will cover working with already imported data in VESIcal.

### Standard normalization

Returns the composition normalized to 100%, including any volatiles.

```python
standard = mysample.get_composition(normalization="standard")
```

If you wish to update the composition in mysample to the normalized one, you can then do:

```python
mysample.change_composition(standard)
```

### FixedVolatiles Normalization

Normalizes the oxides to 100%, but volatiles remain fixed while other major element oxides are reduced proporitonally so that the total is 100 wt%.

```python
fixed = mysample.get_composition(normalization="fixedvolatiles")
mysample.change_composition(fixed)
```

### AdditionalVolatiles Normalization

Normalizes oxides to 100% assuming the sample is volatile-free. If H2O or CO2 concentrations are passed to the function, their un-normalized values will be retained in addition to the normalized non-volatile oxides, summing to >100%.

```python
additional = mysample.get_composition(normalization="additionalvolatiles")
mysample.change_composition(additional)
```

### Normalize a BatchFile object

One might wish to normalize all samples within a BatchFile object. To do so, you can extract and normalize all of the data from your BatchFile object and then create a new BatchFile object with the now normalized data:

```python
my_normed_data = myfile.get_data(normalization="standard")
myNewData = v.BatchFile(filname=None, dataframe=my_normed_data)
```

The value for normalization can be any of "standard", "fixedvolatiles", or "additionalvolatiles".
