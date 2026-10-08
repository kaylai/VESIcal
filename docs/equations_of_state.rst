##################
Equations of State
##################
.. contents::

VESIcal implements multiple equations of state that are used in the VESIcal submodels, which are accessible to the user. This modular implementation of the components of the submodels allows them to be adjusted or even hybridized.

Ideal Gas Law
-------------
.. list-table::
   :header-rows: 1
   :widths: 30 10 5

   * - Description
     - Phase
     - Function
   * - Returns the fugacity of an ideal gas, i.e., the partial pressure.
     - H2O, CO2
     - ``fugacity_idealgas()``

Kerrick and Jacobs (1981)
-------------------------
.. list-table::
   :header-rows: 1
   :widths: 30 10 5

   * - Description
     - Phase
     - Function
   * - Implementation of the Kerrick and Jacobs (1981) EOS for mixed fluids.
     - H2O, CO2
     - ``fugacity_KJ81_h2o()``, ``fugacity_KJ81_co2()``

Zhang and Duan (2009)
---------------------
.. list-table::
   :header-rows: 1
   :widths: 30 10 5

   * - Description
     - Phase
     - Function
   * - Implementation of the Zhang and Duan (2009) fugacity model for pure CO2 fluids.
     - CO2
     - ``fugacity_ZD09_co2()``

Modified Redlich Kwong
-------------------------
.. list-table::
   :header-rows: 1
   :widths: 30 10 5

   * - Description
     - Phase
     - Function
   * - Fugacity model as used by VolatileCalc. Python implementation by D. J Rasmussen (github.com/DJRgeoscience/VolatileCalcForPython), based on VB code by Newman & Lowenstern (2002).
     - H2O, CO2
     - ``fugacity_MRK_h2o()``, ``fugacity_MRK_co2()``
   * - Modified Redlich Kwong by Holloway and Blank (1994)
     - H2O, CO2 
     - ``fugacity_HB_h2o()``, ``fugacity_HB_co2()``
   * - Implementation of the Modified Redlich Kwong presented in Holloway and Blank (1994) Reviews in Mineralogy and Geochemistry vol. 30. Originally written in Quickbasic. CO2 calculations translated to Matlab by Chelsea Allison and translated to python by K. Iacovino for VESIcal. H2O calculations translated to VisualBasic by Gordon M. Moore and translated to python by K. Iacovino for VESIcal. 
     - CO2
     - ``fugacity_HollowayBlank()``

Redlich Kwong
-------------
.. list-table::
   :header-rows: 1
   :widths: 30 10 5

   * - Description
     - Phase
     - Function
   * - Implementation of the Redlich Kwong EoS by Patrick J. Barrie 30 October 2003. Code derived from http://people.ds.cam.ac.uk/pjb10/thermo/pure.html 
     - H2O, CO2, mixed H2O-CO2
     - ``fugacity_RK_h2o()``, ``fugacity_RK_co2()``, ``fugacity_RedlichKwong``

Example Use Cases
=================

- :doc:`Model hybridization </adv_newcalcs>`
- Custom Equations of State
- Extracting thermodynamic values from a VESIcal calculation