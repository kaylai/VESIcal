# Calculating dissolved volatile concentrations

```{contents}
```

The {py:meth}`VESIcal.calculate_dissolved_volatiles()` function calcutions the concentration of dissolved H2O and CO2 in the liquid at a given pressure-temperature condition and with a given H2O-CO2 fluid composition, defined as the mole fraction of H2O in an H2O-CO2 fluid (XH2Ofluid). The default MagmaSat model relies on the underlying functionatlity of MELTS, whose basic function is to calculate the equilibrium phase assemblage given the bulk composition of the system and pressure-temperature conditions. To calculate dissolved volatile concentrations thus requires computing the equilibrium state of a system at fixed pressure and temperature over a range of bulk volatile concentrations until a solution is found that satisfies the user defined fluid composition.

First, the function makes an initial guess at the appropriate bulk volatile concentrations by finding the minimum dissolved volatile concentrations in the liquid at saturation, while asserting that the weight fraction of H2O/(total volatiles) in the system is equal to the user input mole fraction of H2O/(total volatiles) in the fluid. This is done by increasing the H2O and CO2 concentrations appropriately until a fluid phase is stable. Once fluid saturation is determined, the code then performs directional, iterative, and progressively more refined searches, increasing the proportion of H2O or CO2 in the system if the mole fraction of H2O calculated in the fluid is greater than or less than that defined by the user, respectively. Four iterative searches are performed; the precision of the match between the calculated and defined XH2Ofluid increases from 0.1 in the first iteration to 0.01, 0.001, and finally to 0.0001. Thus, the calculated dissolved volatile concentrations correspond to a system with XH2Ofluid within 0.0001 of the user defined value.

**Method structure:**

Single sample:

```python
def calculate_dissolved_volatiles(self, sample, temperature, pressure, X_fluid=1, verbose=False).result
```

ExcelFile batch process:

```python
def calculate_dissolved_volatiles(self, temperature, pressure, X_fluid=1, print_status=False)
```

**Required inputs:**

`sample`: *Only for single-sample calculations*. The composition of a sample. A single sample may be passed as a dictionary of values, with compositions of oxides in wt%.

`temperature`, `pressure`, and `X_fluid`: the temperature in degrees C, the pressure in bars, and the mole fraction of H2O in the H2O-CO2 fluid, XH2Ofluid. Temperature and pressure of the sample or samples must be passed unless an ExcelFile object with a column for temperature and/or pressure is passed to sample. XH2Ofluid is optional, with a default value of 1 (pure H2O fluid). If a numerical (float) value is passed for either temperature, pressure, or X_fluid, that will be the value used for one or all samples. If, alternatively, the user wishes to use temperature, pressure, and/or X_fluid information in their ExcelFile object, the title of the column containing temperature, pressure, or X_fluid data should be passed in quotes (as a string) to temperature, pressure, and/or X_fluid, respectively. Note for batch calculations that if temperature, pressure, or XH2Ofluid information exists in the ExcelFile but a single numerical value is defined for one or both of these variables, both the original information plus the values used for the calculations will be returned.

**Optional inputs:**

`verbose`: *Only for single-sample calculations.* Default value is False. If set to True, additional parameters are returned in a dictionary: H2O and CO2 concentrations in the fluid in mole fraction, temperature, pressure, and proportion of the fluid in the system in wt%.

`print_status`: *Only for ExcelFile batch calcualtions*. The default value is False. If True is passed, the progress of the calculation will be printed to the terminal. The user may desire to see the status of the calculation, as this particular function can be quite slow, averaging between 3-5 seconds per sample.

**Calculated outputs:**
If a single sample is passed to sample, a dictionary with keys ‘H2O’ and ‘CO2’ corresponding to the calculated dissolved H2O and CO2 concentrations in the liquid is returned (plus additional variables ‘temperature’ in degrees C, ‘pressure’ in bars, ‘XH2O_fl’, ‘XCO2_fl’, and ‘FluidProportion_wtper’ (the proportion of the fluid in the system in wt%) if verbose is set to True).

If mutliple samples are passed as an ExcelFile object, a pandas DataFrame is returned with sample information plus calculated dissolved H2O and CO2 concentrations in the liquid, the fluid composition in mole fraction, and the proportion of the fluid in the system in wt%. Pressure (in bars) and Temperature (in degrees C) columns are always returned.

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
dissolved = myfile.calculate_dissolved_volatiles(temperature=900.0, pressure=1000.0, X_fluid=0.5, print_status=True)
dissolved
```

```{eval-rst}
.. csv-table:: Output
   :file: tables/dissolved.csv
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
v.calculate_dissolved_volatiles(sample=extracted_bulk_comp, temperature=900.0, pressure=2000.0, X_fluid=0.5).result
```

```python
{'CO2': 0.0704089917125897, 'H2O': 3.40549411877139}
```
