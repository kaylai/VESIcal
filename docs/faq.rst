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
			- Run `.data()` on a Sample object or `.get_data()` on a BatchFile object and visually inspect your data, ensuring values are present for all species as you expect them to be.
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
