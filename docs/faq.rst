#####
FAQ's
#####
.. contents::

These are still a work in progress. Have a Q? Get in touch by email: kiacovino@seti.org Or submit an issue on the `github repository <https://github.com/kaylai/VESIcal>`_.

A calculation hangs or crashes
==============================
- Possible Causes:
	- **Your column headers weren't recognized as oxide names**
		- Try:
			- Inspect both the raw data (``.data()`` on a Sample or ``.get_data()`` on a BatchFile) and ``.get_composition()``. A value visible in raw data may not be recognized as a chemical component. Check that the expected oxide values are present in the composition.
			- For ``wtpt_oxides`` BatchFile input, the warning beginning ``Possible unrecognized oxide column headers:`` identifies selected whitespace, case, phase and unit variants before sample-index selection. It does not rename columns or select a phase. Check the sample label, phase, units and possible canonical-name collisions before explicitly renaming intended concentration columns; see :doc:`ex_import`.
	- **MagmaSat fails to converge on an answer or crashes entirely**
		- In this case, your batch calculation will either hang on this sample, never spitting out a result, or a segmentation fault will occur within the ENKI thermoengine module (which MagmaSat runs on top of). We cannot catch these errors in VESIcal, since they are occuring within thermoengine. 
		- Try:
			- Delete any MnO concentrations in your dataset. Concentrations are low, and they have no effect on the result of a solubility calculation.
			- Determine which sample is causing an issue (use the status bar as a guide). Remove that sample from your BatchFile and running the calculate command again.

Video on how to troubleshoot a MagmaSat crash
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Useful files:**

	* Dataset: `Example_Crash.xlsx <https://github.com/kaylai/VESIcal/raw/master/docs/jupyter_notebooks/Example_Crash.xlsx>`_
	* Jupyter notebook (right click and choose 'Save Link As...'): `MagmaSat_Crashing.ipynb <https://github.com/kaylai/VESIcal/blob/master/docs/jupyter_notebooks/MagmaSat_Crashing.ipynb>`_

.. raw:: html

	<iframe width="560" height="315" src="https://www.youtube.com/embed/NTasRQsAmHA" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
