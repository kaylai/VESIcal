# API Documentation

```{contents}
```

## Modules

### VESIcal

#### Model()

```{eval-rst}
.. autoclass:: VESIcal.model_classes.Model
        :members:
```

#### FugacityModel()

```{eval-rst}
.. autoclass:: VESIcal.fugacity_models.FugacityModel
        :members:
```

#### activity_model()

```{eval-rst}
.. autoclass:: VESIcal.activity_models.activity_model
        :members:
```

#### Calculate()

```{eval-rst}
.. autoclass:: VESIcal.calculate_classes.Calculate
        :members:
```

##### calculate_dissolved_volatiles(Calculate)

```{eval-rst}
.. autoclass:: VESIcal.calculate_classes.calculate_dissolved_volatiles
        :members:
```

##### calculate_equilibrium_fluid_comp(Calculate)

```{eval-rst}
.. autoclass:: VESIcal.calculate_classes.calculate_equilibrium_fluid_comp
        :members:
```

##### calculate_saturation_pressure(Calculate)

```{eval-rst}
.. autoclass:: VESIcal.calculate_classes.calculate_saturation_pressure
        :members:
```

##### calculate_isobars_and_isopleths(Calculate)

```{eval-rst}
.. autoclass:: VESIcal.calculate_classes.calculate_isobars_and_isopleths
        :members:
```

##### calculate_degassing_path(Calculate)

```{eval-rst}
.. autoclass:: VESIcal.calculate_classes.calculate_degassing_path
        :members:
```

##### calculate_liquid_density(Calculate)

```{eval-rst}
.. autoclass:: VESIcal.thermo.thermo_calculate_classes.calculate_liquid_density
        :members:
```

##### calculate_liquid_viscosity(Calculate)

```{eval-rst}
.. autoclass:: VESIcal.thermo.thermo_calculate_classes.calculate_liquid_viscosity
        :members:
```

### batchfile module

Functions defined in VESIcal.batchfile

#### BatchFile()

```{eval-rst}
.. autoclass:: VESIcal.batchfile.BatchFile
        :members:
```

#### status_bar()

```{eval-rst}
.. autofunction:: VESIcal.batchfile.status_bar
```

#### BatchFile_from_DataFrame()

```{eval-rst}
.. autoclass:: VESIcal.batchmodel.BatchFile_from_DataFrame
        :members:
```

### Sample

Functions defined in VESIcal.sample_class

#### Sample()

```{eval-rst}
.. autoclass:: VESIcal.sample_class.Sample
        :members:
```

### Fugacity Models

#### fugacity_idealgas(FugacityModel)

```{eval-rst}
.. autoclass:: VESIcal.fugacity_models.fugacity_idealgas
        :members:
```

#### fugacity_KJ81_co2(FugacityModel)

```{eval-rst}
.. autoclass:: VESIcal.fugacity_models.fugacity_KJ81_co2
        :members:
```

#### fugacity_KJ81_h2o(FugacityModel)

```{eval-rst}
.. autoclass:: VESIcal.fugacity_models.fugacity_KJ81_h2o
        :members:
```

#### fugacity_ZD09_co2(FugacityModel)

```{eval-rst}
.. autoclass:: VESIcal.fugacity_models.fugacity_ZD09_co2
        :members:
```

#### fugacity_RedlichKwong(FugacityModel)

```{eval-rst}
.. autoclass:: VESIcal.fugacity_models.fugacity_RedlichKwong
        :members:
```

#### fugacity_HollowayBlank(FugacityModel)

```{eval-rst}
.. autoclass:: VESIcal.fugacity_models.fugacity_HollowayBlank
        :members:
```

#### fugacity_HB_co2(FugacityModel)

```{eval-rst}
.. autoclass:: VESIcal.fugacity_models.fugacity_HB_co2
        :members:
```

#### fugacity_HB_h2o(FugacityModel)

```{eval-rst}
.. autoclass:: VESIcal.fugacity_models.fugacity_HB_h2o
        :members:
```

### Activity Models

#### activity_idealsolution(activity_model)

```{eval-rst}
.. autoclass:: VESIcal.activity_models.activity_idealsolution
        :members:

```

### Pure Fluid Models

#### ShishkinaCarbon(Model)

```{eval-rst}
.. autoclass:: VESIcal.models.shishkina.carbon
        :members:
```

#### ShishkinaWater(Model)

```{eval-rst}
.. autoclass:: VESIcal.models.shishkina.water
        :members:
```

#### DixonCarbon(Model)

```{eval-rst}
.. autoclass:: VESIcal.models.dixon.carbon
        :members:
```

#### DixonWater(Model)

```{eval-rst}
.. autoclass:: VESIcal.models.dixon.water
        :members:
```

#### IaconoMarzianoCarbon(Model)

```{eval-rst}
.. autoclass:: VESIcal.models.iaconomarziano.carbon
        :members:
```

#### IaconoMarzianoWater(Model)

```{eval-rst}
.. autoclass:: VESIcal.models.iaconomarziano.water
        :members:
```

#### LiuWater(Model)

```{eval-rst}
.. autoclass:: VESIcal.models.liu.water
        :members:
```

#### LiuCarbon(Model)

```{eval-rst}
.. autoclass:: VESIcal.models.liu.carbon
        :members:
```

#### MooreWater(Model)

```{eval-rst}
.. autoclass:: VESIcal.models.moore.water
        :members:
```

#### AllisonCarbon(Model)

```{eval-rst}
.. autoclass:: VESIcal.models.allison.carbon
        :members:

```

### Mixed Fluid Models

#### MixedFluid(Model)

```{eval-rst}
.. autoclass:: VESIcal.model_classes.MixedFluid
        :members:
```

#### MagmaSat(Model)

```{eval-rst}
.. autoclass:: VESIcal.models.magmasat.MagmaSat
        :members:
```

### VESIcal Plotting Functions

Functions defined in VESIcal.vplot

```{eval-rst}
.. autofunction:: VESIcal.vplot.plot
```

```{eval-rst}
.. autofunction:: VESIcal.vplot.scatterplot
```

```{eval-rst}
.. autofunction:: VESIcal.vplot.smooth_isobars_and_isopleths
```

```{eval-rst}
.. autofunction:: VESIcal.vplot.calib_plot
```

### Data Transformation Functions

```{eval-rst}
.. autofunction:: VESIcal.fluid_molfrac_to_wt
```

```{eval-rst}
.. autofunction:: VESIcal.fluid_wt_to_molfrac
```

```{eval-rst}
.. autofunction:: VESIcal.get_oxides
```

### Universal Informative Functions

```{eval-rst}
.. autofunction:: VESIcal.get_model_names


















```
