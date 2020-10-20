#############
About VESIcal
#############

VESIcal is an open-source thermodynamic model engine for mixed volatile (H2O-CO2) solubility in silicate melts. It was designed by Kayla Iacovino, Simon Matthews, Penny Wieser, Gordon Moore, and Florence Begue in 2020.

Current Version 0.1.2 (Pre-Review)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
At present, a manuscript, which will serve as a user's guide and introduction to VESIcal, is in preparation. VESIcal is not yet peer-reviewed but will be submitted for peer-review soon. The current stable release of VESIcal is version 0.1.2. This version is archived on zenodo and has a citable doi. 

How to Cite VESIcal
^^^^^^^^^^^^^^^^^^^
To cite computations done using VESIcal, please cite the in prep manuscript, the VESIcal version number, as well as the model(s) used. For example: “Calculations were performed using VESIcal (v. 0.1.2; Iacovino et al., in prep) with the models of Shishkina et al. (2014) and Dixon (1997; “VolatileCalc”).” The web-app always runs on the most up-to-date version of the VESIcal code, but it is best practice to note if the web-app was used (“Calculations were perforumed using the VESIcal web-app (v. 0.1.2; Iacovino et al., 2020)...”). We also encourage users to be as explicit as possible as to the conditions used for modelling. This includes stating the pressure, temperature, volatile concentration, and bulk magma composition used in modelling. In the best case, VESIcal users will provide their code (e.g., as a jupyter notebook or .py file) along with their publication such that it can be easily replicated.

The full citation of the (in prep) manuscript is:

	Iacovino K., Matthews S., Wieser P.E., Moore G.M., and Begue F. (in prep) VESIcal: An open-source thermodynamic model engine for mixed volatile (H2O-CO2) solubility in silicate melts, Earth and Space Sciences.

To cite the zenodo archive of the current stable release:

	Iacovino, Kayla, Matthews, Simon, Wieser, Penny E., Moore, Gordon M., & Begue, Florence. (2020, October 16). VESIcal v. 0.1.2 (Review) (Version 0.1.2). Zenodo. `http://doi.org/10.5281/zenodo.4096463 <http://doi.org/10.5281/zenodo.4096463>`_

Please also cite whichever model(s) within VESIcal were used. If using the default model (if no model was specified, then you were using the default model), cite MagmaSat as:

	Ghiorso, M. and Gualda, G. (2015) An H2O-CO2 mixed fluid saturation model compatible with rhyolite-MELTS, Contributions to Mineralogy and Petrology, 169, 1–30.

Many Ways to Use VESIcal
^^^^^^^^^^^^^^^^^^^^^^^^
VESIcal can be accessed and used in a variety of ways. From most flexible (advanced) to least flexible (novice), these are:

	- local installation of the VESIcal library
	- through a jupyter notebook hosted online (the VESIcal manuscript will be in the form of a jupyter notebook once published)
	- via the web-app (`https://vesical.anvil.app/ <https://vesical.anvil.app/>`_)

This documentation mainly serves to help users interact with VESIcal as a python library. That is, the user will instal VESIcal onto their machine (or run VESIcal in a jupyter notebook) and type python code to execute VESIcal.