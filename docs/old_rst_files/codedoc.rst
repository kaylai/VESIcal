#################
API Documentation
#################
.. contents::

*******
Modules
*******

VESIcal
=======

Model()
-------
.. autoclass:: VESIcal.model_classes.Model
	:members:

FugacityModel()
---------------
.. autoclass:: VESIcal.fugacity_models.FugacityModel
	:members:

activity_model()
----------------
.. autoclass:: VESIcal.activity_models.activity_model
	:members:

Calculate()
-----------
.. autoclass:: VESIcal.calculate_classes.Calculate
	:members:

calculate_dissolved_volatiles(Calculate)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: VESIcal.calculate_classes.calculate_dissolved_volatiles
	:members:

calculate_equilibrium_fluid_comp(Calculate)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: VESIcal.calculate_classes.calculate_equilibrium_fluid_comp
	:members:

calculate_saturation_pressure(Calculate)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: VESIcal.calculate_classes.calculate_saturation_pressure
	:members:

calculate_isobars_and_isopleths(Calculate)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: VESIcal.calculate_classes.calculate_isobars_and_isopleths
	:members:

calculate_degassing_path(Calculate)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: VESIcal.calculate_classes.calculate_degassing_path
	:members:

calculate_liquid_density(Calculate)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: VESIcal.thermo.thermo_calculate_classes.calculate_liquid_density
	:members:

calculate_liquid_viscosity(Calculate)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: VESIcal.thermo.thermo_calculate_classes.calculate_liquid_viscosity
	:members:

batchfile module
================
Functions defined in VESIcal.batchfile

BatchFile()
-----------
.. autoclass:: VESIcal.batchfile.BatchFile
	:members:

status_bar()
------------
.. autofunction:: VESIcal.batchfile.status_bar

BatchFile_from_DataFrame()
--------------------------
.. autoclass:: VESIcal.batchmodel.BatchFile_from_DataFrame
	:members:

Sample
======
Functions defined in VESIcal.sample_class

Sample()
--------
.. autoclass:: VESIcal.sample_class.Sample
	:members:

Fugacity Models
===============

fugacity_idealgas(FugacityModel)
--------------------------------
.. autoclass:: VESIcal.fugacity_models.fugacity_idealgas
	:members:

fugacity_KJ81_co2(FugacityModel)
--------------------------------
.. autoclass:: VESIcal.fugacity_models.fugacity_KJ81_co2
	:members:

fugacity_KJ81_h2o(FugacityModel)
--------------------------------
.. autoclass:: VESIcal.fugacity_models.fugacity_KJ81_h2o
	:members:

fugacity_ZD09_co2(FugacityModel)
--------------------------------
.. autoclass:: VESIcal.fugacity_models.fugacity_ZD09_co2
	:members:

fugacity_RedlichKwong(FugacityModel)
------------------------------------
.. autoclass:: VESIcal.fugacity_models.fugacity_RedlichKwong
	:members:

fugacity_HollowayBlank(FugacityModel)
-------------------------------------
.. autoclass:: VESIcal.fugacity_models.fugacity_HollowayBlank
	:members:

fugacity_HB_co2(FugacityModel)
------------------------------------
.. autoclass:: VESIcal.fugacity_models.fugacity_HB_co2
	:members:

fugacity_HB_h2o(FugacityModel)
------------------------------------
.. autoclass:: VESIcal.fugacity_models.fugacity_HB_h2o
	:members:

Activity Models
===============

activity_idealsolution(activity_model)
--------------------------------------
.. autoclass:: VESIcal.activity_models.activity_idealsolution
	:members:


Pure Fluid Models
=================

ShishkinaCarbon(Model)
----------------------
.. autoclass:: VESIcal.models.shishkina.carbon
	:members:

ShishkinaWater(Model)
---------------------
.. autoclass:: VESIcal.models.shishkina.water
	:members:

DixonCarbon(Model)
------------------
.. autoclass:: VESIcal.models.dixon.carbon
	:members:

DixonWater(Model)
-----------------
.. autoclass:: VESIcal.models.dixon.water
	:members:

IaconoMarzianoCarbon(Model)
---------------------------
.. autoclass:: VESIcal.models.iaconomarziano.carbon
	:members:

IaconoMarzianoWater(Model)
--------------------------
.. autoclass:: VESIcal.models.iaconomarziano.water
	:members:

LiuWater(Model)
-------------------
.. autoclass:: VESIcal.models.liu.water
	:members:

LiuCarbon(Model)
-------------------
.. autoclass:: VESIcal.models.liu.carbon
	:members:

MooreWater(Model)
-----------------
.. autoclass:: VESIcal.models.moore.water
	:members:

AllisonCarbon(Model)
--------------------
.. autoclass:: VESIcal.models.allison.carbon
	:members:


Mixed Fluid Models
==================

MixedFluid(Model)
-----------------
.. autoclass:: VESIcal.model_classes.MixedFluid
	:members:

MagmaSat(Model)
-----------------
.. autoclass:: VESIcal.models.magmasat.MagmaSat
	:members:

VESIcal Plotting Functions
==========================
Functions defined in VESIcal.vplot

.. autofunction:: VESIcal.vplot.plot

.. autofunction:: VESIcal.vplot.scatterplot

.. autofunction:: VESIcal.vplot.smooth_isobars_and_isopleths

.. autofunction:: VESIcal.vplot.calib_plot

Data Transformation Functions
=============================

.. autofunction:: VESIcal.fluid_molfrac_to_wt

.. autofunction:: VESIcal.fluid_wt_to_molfrac

.. autofunction:: VESIcal.get_oxides

Universal Informative Functions
===============================

.. autofunction:: VESIcal.get_model_names



















