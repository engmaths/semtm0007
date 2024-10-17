Wind turbine SCADA data
=======================

Data is taken from <https://www.kaggle.com/datasets/pythonafroz/wind-turbine-structural-health-monitoring-data>

It is licensed under the CC BY-NC-SA 4.0 and should be cited as follows:
  
> Sara Fogelström, Håkan Johansson, Ola Carlson, Martin Hofsäß, Oliver
> Bischoff, Yuriy Marykovskiy, & Imad Abdallah. (2023). Björkö Wind Turbine
> Version 1 (45kW) high frequency Structural Health Monitoring (SHM) data
> [Data set]. Zenodo. <https://doi.org/10.5281/zenodo.8230330>

The data is a subset of the B1_CL4_20.csv file from the dataset. The subsets are
contiguous in time, though occasional single sample dropouts are present. The
subsets have been chosen such that SysMode == 12 (running) and Dig_IO_States ==
8063.

The sampling frequency is 20 Hz.

Columns of particular interest are:

- (0) Time: time in seconds since 00:00:00 01/01/1904 in UTC (s)
- (9) DCV: DC voltage output from the generator (A)
- (10) DCC: DC current output from the generator (V)
- (59) WDN: wind direction at the nacelle; 180deg means fully facing the wind direction (degrees)
- (60) WSN: wind speed at the nacelle (m/s)
