"""
This script helps you make a profile from ERA5 data

The profile can then be input to LBLRTM

Written by sophie but reckon I got inspiration from sanjee/laura

"""
import os
import numpy as np
import subprocess as sub
import glob as glob
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt

"Here is your INPUT raw data net cdf format"
ERA5_data_whole = xr.open_dataset('/disk1/sm4219/Andoya_key/ERA5_raw_data/ERA5_20_02_23.nc')

output_file = '/disk1/sm4219/Andoya_key/ERA5_13pm_02_20th.txt'
# time = '2023-03-08 09:00:00'
time = np.datetime64('2023-02-20T13:00:00')

ERA5_data_4 = ERA5_data_whole.sel(latitude=69.16, longitude=16.0, method='nearest')

# Check the available dimensions and coordinates first
print(ERA5_data_whole.dims)
print(ERA5_data_whole.coords)

# If 'time' is indeed the dimension/coordinate, use this
ERA5_data = ERA5_data_4.sel(valid_time=np.datetime64('2023-02-20T13:00:00'))

"GETTING PRESSURE LEVELS"
# Constants
Re = 6371.  # Radius of the Earth in meters
g = 9.81
# Function to convert geopotential to geo height to altitude using Equation 7
def geopotential_to_altitude(Zg):
    "Converts from geopotential to geo height to altitude, using earth radius and g"
    # print('GEO ITSELF', Zg.values)
    geo_height = (Zg/g)/1000.
    # print('GEOPOT HEIGHT', geo_height)
    real_z = geo_height/(1 - (geo_height/Re))
    return real_z


"""
All the code below here is for fitting the surface to surface measurements by FINESSE!

You will need a CO2 value, I assume a constant CO2 profile as it is a well mixed gas, but I don't get this 
value from ERA5...
"""
PTH_file = pd.read_csv('/disk1/Andoya/sp1016/Raw_Data/20230220/PTH_all.txt', sep = '\t')
PTH_file['Unnamed: 0'] = pd.to_datetime(PTH_file['Unnamed: 0'])
PTH_file2 = PTH_file.loc[PTH_file['Unnamed: 0']<pd.Timestamp(time)+pd.Timedelta(seconds=30)]
PTH_file3 = PTH_file2.loc[pd.Timestamp(time)-pd.Timedelta(seconds=30)<PTH_file2['Unnamed: 0']]
# print('time:', time)
# print('HEY PRINTING HERE pth file1, 2, 3:', PTH_file, PTH_file2, PTH_file3)
surface_pressure = PTH_file3['P, hPa'].mean()
surface_temp = PTH_file3['T, C'].mean() + 273.14
surface_RH = PTH_file3['RH, %'].mean()

print('SURFACE VALUES ARE PRINTED HERE % HPA K ',surface_RH, surface_pressure, surface_temp)
"GET CO2 VALUES"
CO2_file = pd.read_csv('/disk1/Andoya/sp1016/Raw_Data/20230220/CO2_all.txt', sep = '\t')
print(CO2_file)
CO2_file['Unnamed: 0'] = pd.to_datetime(CO2_file['Unnamed: 0'], errors='coerce')
CO2_file2 = CO2_file.loc[CO2_file['Unnamed: 0']<pd.Timestamp(time)+pd.Timedelta(seconds=30)]
CO2_file3 = CO2_file2.loc[pd.Timestamp(time)-pd.Timedelta(seconds=30)<CO2_file2['Unnamed: 0']]
CO2_file4 = CO2_file3.astype({'CO2, ppm':'float'})
# print('values:', (CO2_file4['CO2, ppm'].values), type((CO2_file4['CO2, ppm'].values)))
surface_co2  = np.mean((CO2_file4['CO2, ppm'].values))
print('mean here', surface_co2)
# surface_co2_new = CO2_file3['CO2, ppm'].mean()
# print('co2 here:', surface_co2)

pressure_levels = ERA5_data['pressure_level']
ERA5_data_cut = ERA5_data.sel(pressure_level=pressure_levels[pressure_levels < (surface_pressure)])

# ERA5_data_cut = ERA5_data_cut2.sel(pressure_level=pressure_levels[pressure_levels > 1])
# print('old', pressure_levels)
# print('new', ERA5_data_cut['pressure_level'])
# ERA5_data_cut = ERA5_data.interp(pressure_level=surface_pressure)
# - 5 to get rid of 953 value
# print('WORKING')

Zg_geos = ERA5_data_cut['z']
# print('ZG GEOPOTENTIAL HERE:',Zg_geos.values)
altitude = geopotential_to_altitude(Zg_geos)

altitude5 = np.zeros_like(altitude)

surface_co2_2 = 3.956166666666667311e+02
# surface_co2_2 = 394.783
co2= np.full_like(altitude5, surface_co2)
print(co2)

"""
Here is where you are putting everything together
Building a stack of arrays to make your final txt file :)

"""

arrays = [
    np.flip(altitude),
    (ERA5_data_cut['pressure_level'].values),
    (ERA5_data_cut['t'].values),
    (ERA5_data_cut['r'].values),
    # (ERA5_data_cut['q'].values)*1000,
    co2,  # Ensure correct length
    (ERA5_data_cut['o3'].values * 1E3)  # Convert O3 to g/kg foR lblrtm
    ]

"Quick check everything is the same size here"
# Print the shapes of each array for verification
# for z in arrays:
#     print(np.shape(z))



# Combine the arrays into a final profile
final_profile = np.vstack(arrays).T

# Prepare the final data to save, including the surface values at the beginning
# THE OZONE VALUE IS THE 950 VALUE 2.630400000000000000e+02
final_data = np.vstack(([399, surface_pressure, surface_temp, surface_RH, surface_co2, 9.193959704134613276e-05], final_profile))
# final_data = np.vstack(([399, surface_pressure, surface_temp, 1.283994, surface_co2_2, 9.193959704134613276e-05], final_profile))

# Save to the specified output file fmt='%.6f'
np.savetxt(output_file, final_data, header="# Contains: Altitude (km), Pressure (hPa), Temperature (K), RH (%), CO2 (ppmv), O3 (g/kg)", delimiter='\t') 
# %%
