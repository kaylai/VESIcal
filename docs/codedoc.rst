==========================
VESIcal Code Documentation
==========================
.. contents::

Major Classes
=============

Model()
-------
.. autoclass:: VESIcal.Model
	:members:

ExcelFile()
-----------
.. autoclass:: VESIcal.ExcelFile
	:members:

FugacityModel()
---------------
.. autoclass:: VESIcal.FugacityModel
	:members:

activity_model()
----------------
.. autoclass:: VESIcal.activity_model
	:members:

Calculate()
-----------
.. autoclass:: VESIcal.Calculate
	:members:

calculate_dissolved_volatiles(Calculate)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: VESIcal.calculate_dissolved_volatiles
	:members:

calculate_equilibrium_fluid_comp(Calculate)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: VESIcal.calculate_equilibrium_fluid_comp
	:members:

calculate_saturation_pressure(Calculate)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: VESIcal.calculate_saturation_pressure
	:members:

calculate_isobars_and_isopleths(Calculate)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: VESIcal.calculate_isobars_and_isopleths
	:members:

calculate_degassing_path(Calculate)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: VESIcal.calculate_degassing_path
	:members:


Fugacity Models
===============

fugacity_idealgas(FugacityModel)
--------------------------------
.. autoclass:: VESIcal.fugacity_idealgas
	:members:

fugacity_KJ81_co2(FugacityModel)
--------------------------------
.. autoclass:: VESIcal.fugacity_KJ81_co2
	:members:

fugacity_KJ81_h2o(FugacityModel)
--------------------------------
.. autoclass:: VESIcal.fugacity_KJ81_h2o
	:members:

fugacity_ZD09_co2(FugacityModel)
--------------------------------
.. autoclass:: VESIcal.fugacity_ZD09_co2
	:members:

fugacity_RedlichKwong(FugacityModel)
------------------------------------
.. autoclass:: VESIcal.fugacity_RedlichKwong
	:members:


Activity Models
===============

activity_idealsolution(activity_model)
--------------------------------------
.. autoclass:: VESIcal.activity_idealsolution
	:members:


Pure Fluid Models
=================

ShishkinaCarbon(Model)
----------------------
.. autoclass:: VESIcal.ShishkinaCarbon
	:members:

ShishkinaWater(Model)
---------------------
.. autoclass:: VESIcal.ShishkinaWater
	:members:

DixonCarbon(Model)
------------------
.. autoclass:: VESIcal.DixonCarbon
	:members:

DixonWater(Model)
-----------------
.. autoclass:: VESIcal.DixonWater
	:members:

IaconoMarzianoWater(Model)
--------------------------
.. autoclass:: VESIcal.IaconoMarzianoWater
	:members:

IaconoMarzianoCarbon(Model)
---------------------------
.. autoclass:: VESIcal.IaconoMarzianoCarbon
	:members:

EguchiCarbon(Model)
-------------------
.. autoclass:: VESIcal.EguchiCarbon
	:members:

MooreWater(Model)
-----------------
.. autoclass:: VESIcal.MooreWater
	:members:

AllisonCarbon(Model)
--------------------
.. autoclass:: VESIcal.AllisonCarbon
	:members:


Mixed Fluid Models
==================

MixedFluid(Model)
-----------------
.. autoclass:: VESIcal.MixedFluid
	:members:

MagmaSat(Model)
-----------------
.. autoclass:: VESIcal.MagmaSat
	:members:

MagmaSat Plotting Functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Plot Isobars and Isopleths
++++++++++++++++++++++++++
.. autofunction:: VESIcal.plot_isobars_and_isopleths

Plot Degassing Paths
++++++++++++++++++++
.. autofunction:: VESIcal.plot_degassing_paths
















