# FAQ's

```{contents}
```

These are still a work in progress. Have a Q? Get in touch by email: <mailto:kayla.iacovino@nasa.gov>. Or submit an issue on the [github repository](https://github.com/kaylai/VESIcal).

## I'm running a BatchFile with MagmaSat, and the calculation hangs or crashes

Sometimes MagmaSat can fail to converge on an answer for a specific sample. In this case, your batch calculation will either hang on this sample, never spitting out a result, or a segmentation fault will occur within the ENKI thermoengine module (which MagmaSat runs on top of). We cannot catch these errors in VESIcal, since they are occuring within thermoengine. First, run `.get_data()` on your BatchFile object and visually inspect your data, ensuring values are present for all species as you expect them to be. Next, determine which sample is causing an issue (use the status bar as a guide). Try removing that sample from your BatchFile and running the calculate command again.

### Video on how to troubleshoot a MagmaSat crash

**Useful files:**

> - Dataset: [Example_Crash.xlsx](https://github.com/kaylai/VESIcal/raw/master/docs/jupyter_notebooks/Example_Crash.xlsx)
> - Jupyter notebook (right click and choose 'Save Link As...'): [MagmaSat_Crashing.ipynb](https://github.com/kaylai/VESIcal/blob/master/docs/jupyter_notebooks/MagmaSat_Crashing.ipynb)

```{raw} html
<iframe width="560" height="315" src="https://www.youtube.com/embed/NTasRQsAmHA" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
```
