# Calculating liquid densities

```{contents}
```

The {py:meth}`VESIcal.calculate_liquid_density()` function calculates the density of the silicate liquid given composition, temperature, and pressure. The function uses the DensityX model of Iacovino and Till (2019). No other models are currently available for this calculation.

**Method Structure**

Single sample:

```python
def calculate_liquid_density(self, sample, pressure, temperature).result
```

BatchFile process:

```python
def calculate_liquid_density(self, pressure, temperature)
```

**Required inputs:**

`sample`: *Only for single-sample calculations*. The composition of a sample as Sample class.

`pressure`: The pressure in bars. For BatchFile calculations, if pressure information is present in the file (e.g., as a column with unique pressure values for each sample), this can be accessed by passing the column name in quotes to the pressure variable.

`temperature`: The temperature in degres C. For BatchFile calculations, if temperature information is present in the file (e.g., as a column with unique temperature values for each sample), this can be accessed by passing the column name in quotes to the temperature variable.

**Calculated outputs:**
The density of the liquid in grams per liter, rounded to 3 dp.

## For an entire dataset

### Import a data file

```python
myfile = v.BatchFile('example_data.xlsx')
myfile.get_data()
```

```{eval-rst}
.. csv-table:: Output
   :file: tables/example_data.csv
   :header-rows: 1
```

### Do the calculation

```python
densities = myfile.calculate_liquid_density(pressure=1000, temperature=900)
densities
```

```{eval-rst}
.. csv-table:: Output
   :file: tables/densities.csv
   :header-rows: 1
```

## For a single sample

### Extract a single sample from your dataset

```python
SampleName = 'BT-ex'
extracted_bulk_comp = myfile.get_sample_composition(SampleName, asSampleClass=True)
```

### Do the calculation

```python
v.calculate_liquid_density(sample=extracted_bulk_comp, pressure=1000, temperature=900).result
```

```python
2142.827
```
