#############
VESIcal Manuscripts
#############

Two VESIcal manuscripts serve, along with this documentation, as a user's guide and introduction to VESIcal. The Part I manuscript is available as an executable jupyter notebook or as a static PDF published in the journal Earth and Space Science. Part II is available as a PDF with jupyter notebooks as supplementary files.

The manuscript version of VESIcal Part I (v 1.0.1) is archived on zenodo with a citable doi at `https://zenodo.org/record/5095382 <https://zenodo.org/record/5095382>`_



Interactive Manuscript on Curvenote
-----------------------------------
Through AGU's Notebooks Now! initiative in partnership with Curvenote, VESIcal Part I was among the first manuscripts to be transformed into a Curvenote style interactive manuscript along with a binderized version of the manuscript as a jupyter notebook also hosted on the Curvenote servers. This allows for a more robust binderized environment that users can go to to get coding right away without any installations on their local machine. Since VESIcal and the ENKI ThermoEngine (required to run MagmaSat) are installed on Curvenote's servers, you can not only execute the code in the manuscript, but you can get to writing your own.

- Check out the `Curvenote implementation of VESIcal Part I <https://agu.curve.space/articles/NN0001>`_
- Go directly to the binderized `jupyter notebook of VESIcal Part I <https://agu-binder.curvenote.dev/user/2be900e9-fb5d-4-9778d16a48c.zip-4wrcztow/lab/tree/Manuscript.ipynb?token=EzBUfh6US4qFq4UW0MSkYA>`_


Older Interactive versions of manuscripts
-----------------------------------------

   - Direct link to interactive VESIcal Part I manuscript: `Manuscript on Binder <https://mybinder.org/v2/gh/kaylai/vesical-binder/HEAD?filepath=Manuscript.ipynb>`_
   - Jupyter Notebook hub with VESIcal: `Jupyter hub on Binder <https://mybinder.org/v2/gh/kaylai/vesical-binder/HEAD>`_

PDF versions of manuscripts
---------------------------

   - `VESIcal Part I: An Open-Source Thermodynamic Model Engine for Mixed Volatile (H2O-CO2) Solubility in Silicate Melts <https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2020EA001584>`_

	Iacovino K., Matthews S., Wieser P.E., Moore G.M., and Begue F. (2021) VESIcal Part I: An open-source thermodynamic model engine for mixed volatile (H2O-CO2) solubility in silicate melts, Earth and Space Science, 8, e2020EA001584. https://doi.org/10.1029/2020EA001584.
	
   - `VESIcal Part II: A critical approach to volatile solubility modelling using an open-source Python3 engine <https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2021EA001932>`_
	
	Wieser P.E., Iacovino K., Matthews S., Moore G.M., Allison C.M. (2022) VESIcal Part II: A critical approach to volatile solubility modelling using an open-source Python3 engine, Earth and Space Science. https://doi.org/10.1029/2021EA001932


How to Cite VESIcal
^^^^^^^^^^^^^^^^^^^
To cite computations done using VESIcal, please cite the 2021 manuscript, the VESIcal version number, and as the model(s) used. For example: “Calculations were performed using VESIcal (v. 1.2.6; Iacovino et al., 2021) with the models of Shishkina et al. (2014) and Dixon (1997; “VolatileCalc”).” The web-app always runs on the most up-to-date version of the VESIcal code, but it is best practice to note if the web-app was used (“Calculations were performed using the VESIcal web-app (v. 1.2.5; Iacovino et al., 2021)...”). We also encourage users to be as explicit as possible as to the conditions used for modelling. This includes stating the pressure, temperature, volatile concentration, and bulk magma composition used in modelling. In the best case, VESIcal users will provide their code (e.g., as a jupyter notebook or .py file) along with their publication such that it can be easily replicated.

The full citation of the manuscript is:

	Iacovino K., Matthews S., Wieser P.E., Moore G.M., and Begue F. (2021) VESIcal Part I: An open-source thermodynamic model engine for mixed volatile (H2O-CO2) solubility in silicate melts, Earth and Space Science, 8, e2020EA001584. https://doi.org/10.1029/2020EA001584

Please also cite whichever model(s) within VESIcal were used. If using the default model (if no model was specified, then you were using the default model), cite MagmaSat as:

	Ghiorso, M. and Gualda, G. (2015) An H2O-CO2 mixed fluid saturation model compatible with rhyolite-MELTS, Contributions to Mineralogy and Petrology, 169, 1–30.

To cite the zenodo archive of the current stable release (rare that you will need to do this):

	Iacovino, Kayla, Matthews, Simon, Wieser, Penny E., Moore, Gordon M., & Begue, Florence. (2021, June 8). VESIcal v. 1.0.1. Zenodo. `https://zenodo.org/record/5095382 <https://zenodo.org/record/5095382>`_
