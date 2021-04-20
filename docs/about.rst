#############
About VESIcal
#############

VESIcal is an open-source thermodynamic model engine for mixed volatile (H2O-CO2) solubility in silicate melts. It was designed by Kayla Iacovino, Simon Matthews, Penny Wieser, Gordon Moore, and Florence Begue in 2021.

Current Version 0.9.11 (version in press)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The VESIcal manuscript (currently in press) serves, along with this documentation, as a user's guide and introduction to VESIcal. The manuscript will be available as an executable jupyter notebook or as a static PDF published in the journal Earth and Space Science.

	Iacovino K., Matthews S., Wieser P.E., Moore G.M., and Begue F. (in press) VESIcal Part I: An open-source thermodynamic model engine for mixed volatile (H2O-CO2) solubility in silicate melts, Earth and Space Sciences. (a link to the manuscript will be here once it is public)

This version of VESIcal (v 0.9.11) is archived on zenodo with a citable doi at `https://zenodo.org/record/4652839 <https://zenodo.org/record/4652839>`_


How to Cite VESIcal
^^^^^^^^^^^^^^^^^^^
To cite computations done using VESIcal, please cite the in prep manuscript, the VESIcal version number, as well as the model(s) used. For example: “Calculations were performed using VESIcal (v. 0.9.11; Iacovino et al., 2021) with the models of Shishkina et al. (2014) and Dixon (1997; “VolatileCalc”).” The web-app always runs on the most up-to-date version of the VESIcal code, but it is best practice to note if the web-app was used (“Calculations were perforumed using the VESIcal web-app (v. 0.9.11; Iacovino et al., 2021)...”). We also encourage users to be as explicit as possible as to the conditions used for modelling. This includes stating the pressure, temperature, volatile concentration, and bulk magma composition used in modelling. In the best case, VESIcal users will provide their code (e.g., as a jupyter notebook or .py file) along with their publication such that it can be easily replicated.

The full citation of the (in press) manuscript is:

	Iacovino K., Matthews S., Wieser P.E., Moore G.M., and Begue F. (in press) VESIcal Part I: An open-source thermodynamic model engine for mixed volatile (H2O-CO2) solubility in silicate melts, Earth and Space Sciences.

To cite the zenodo archive of the current stable release:

	Iacovino, Kayla, Matthews, Simon, Wieser, Penny E., Moore, Gordon M., & Begue, Florence. (2021, March 31). VESIcal v. 0.9.11. Zenodo. `https://doi.org/10.5281/zenodo.4291043 <https://doi.org/10.5281/zenodo.4652839>`_

Please also cite whichever model(s) within VESIcal were used. If using the default model (if no model was specified, then you were using the default model), cite MagmaSat as:

	Ghiorso, M. and Gualda, G. (2015) An H2O-CO2 mixed fluid saturation model compatible with rhyolite-MELTS, Contributions to Mineralogy and Petrology, 169, 1–30.

Many Ways to Use VESIcal
^^^^^^^^^^^^^^^^^^^^^^^^
VESIcal requires installation of not only the VESIcal library but also some other python libraries, one of which is a bit tricky to install (ENKI/thermoengine aka the engine behind MELTS). But, we have a solution! All dependencies and the latest version of VESIcal are all installed on the ENKI server, within a Jupyter Notebook Hub. Steps to use VESIcal on the ENKI server are:

	1. Create a (free) GitLab account, which you'll use to sign into ENKI here: (`https://gitlab.com/users/sign_up <https://gitlab.com/users/sign_up>`_)
	2. Access the ENKI Production Server by going to `http://enki-portal.org/ <http://enki-portal.org/>`_ and clicking "SERVERS" > "PRODUCTION SERVER"
	3. Sign in with your GitLab credentials: You are now in your own jupyter notebook workspace! You can upload and create files here. They won't be accessible to anyone else. 
	4. Click the green "CLOSE THIS SCREEN" button
	5. Create a new notebook by clicking the blue plus button and then selecting Python3 under Notebook. Or select from the menu File > New > Notebook
	6. Be sure to import VESIcal as v at the top of your file, and now you are ready to get to work!

In general, VESIcal can be accessed and used in a variety of ways. From most flexible (advanced) to least flexible (novice), these are:

	- local installation of the VESIcal library
	- through a jupyter notebook hosted online (the VESIcal manuscript will be in the form of a jupyter notebook once published)
	- via the web-app (`https://vesical.anvil.app/ <https://vesical.anvil.app/>`_)

This documentation mainly serves to help users interact with VESIcal as a python library. That is, the user will instal VESIcal onto their machine (or run VESIcal in a jupyter notebook) and type python code to execute VESIcal.