# Updates

```{contents}
```

## VESIcal is now hosted on Curvenote

Through AGU's Notebooks Now! initiative in partnership with Curvenote, VESIcal Part I was among the first manuscripts to be transformed into a Curvenote style interactive manuscript along with a binderized version of the manuscript as a jupyter notebook also hosted on the Curvenote servers. This allows for a more robust binderized environment that users can go to to get coding right away without any installations on their local machine. Since VESIcal and the ENKI ThermoEngine (required to run MagmaSat) are installed on Curvenote's servers, you can not only execute the code in the manuscript, but you can get to writing your own.

- Check out the [Curvenote implementation of VESIcal Part I](https://agu.curve.space/articles/NN0001)
- Go directly to the binderized [jupyter notebook of VESIcal Part I](https://agu-binder.curvenote.dev/user/2be900e9-fb5d-4-9778d16a48c.zip-4wrcztow/lab/tree/Manuscript.ipynb?token=EzBUfh6US4qFq4UW0MSkYA)

## VESIcal can now be used without needing to install thermoengine!

Without thermoengine installed, MagmaSat (the default model) cannot be used, however all other models, all plotting capability, and the thermo package can be used. Simply make sure you are using version 1.2.0 or higher following the simple installation instructions below. Then be sure to explicitly tell VESIcal which model you want to use (any model other than MagmaSat. See {doc}`models` for a complete list) by passing `model="some-model-name"` when performing a calculation.

For example:

```python
v.calculate_saturation_pressure(sample=<your-sample-here>, model="<some-model-name-here>")
```

If thermoengine is not installed and you import VESIcal, you will be warned that you won't be able to use MagmaSat, but everything else will work as expected. Remember that if you do not pass a model name, VESIcal will default to MagmaSat, and an error will be generated telling you that you need to pass a model name.
