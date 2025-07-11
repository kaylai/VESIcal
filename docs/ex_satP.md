# Calculating saturation presures

```{contents}
```

The {py:meth}`VESIcal.calculate_saturation_pressure()` function calculates the minimum pressure at which a given silicate melt with known temperature and H2O and CO2 concentrations would be saturated with fluid. This is calcualted by finding the pressure at which the smallest amount of vapor is present. This function also calculates the composition of the vapor in equilibrium with the melt at those conditions.

The function works by calculating the equilibrium state of the given melt at very high pressure (2,000 MPa) and then decreasing the pressure in steps of 100 MPa until the mass of vapor is >0 grams. At this point, the pressure space is narrowed and searched in steps of 10 MPa and then in steps of 1 MPa until the saturation pressure is found.

**Method structure:**

Single sample:

```python
def calculate_saturation_pressure(self, sample, temperature, verbose=False).result
```

BatchFile process:

```python
def calculate_saturation_pressure(self, temperature, print_status=False)
```

**Required inputs:**

`sample`: *Only for single-sample calculations*. The composition of a sample as Sample class.

`temperature`: The temperature in degres C. For BatchFile calculations, if temperature information is present in the file (e.g., as a column with unique temperature values for each sample), this can be accessed by passing the column name in quotes to the temperature variable.

**Optional inputs:**

`verbose`: *Only for single-sample calculations.* Default value is False. If set to True, additional parameters are returned in a dictionary: saturation pressure in bars, H2O and CO2 concentrations in the fluid, mass of the fluid in grams, and proportion of the fluid in the system in wt%.

`print_status`: *Only for ExcelFile batch calcualtions*. The default value is False. If True is passed, the progress of the calculation will be printed to the terminal.

**Calculated outputs:**
If a single sample is passed to sample, the saturation pressure in bars is returned as a numerical value (float) (plus additional variables ‘XH2O_fl’, ‘XCO2_fl’, ‘FluidMass_grams’, and ‘FluidProportion_wtper’ if verbose is set to True).

If mutliple samples are passed as a BatchFile object, a pandas DataFrame is returned with sample information plus calculated saturation pressures, equilibrium fluid compositions, mass of the fluid in grams, and proportion of the fluid in the system in wt%. Temperature (in degrees C) is always returned.

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
satPs = myfile.calculate_saturation_pressure(temperature=925.0)
satPs
```

```{eval-rst}
.. csv-table:: Output
   :file: tables/satP.csv
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
v.calculate_saturation_pressure(sample=extracted_bulk_comp, temperature=925.0).result
```

```python
2490.0
```
