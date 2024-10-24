# Python_LBLRTM_apodise
The files you need to process LBLRTM using python including apodisation etc for FINESSE

# TO RUN LBLRTM:
LBLRTM is a line by line radiative transfer model. It needs an input of an atmospheric profile, it can assume a standard atmosphere for various gases but it needs to know temperature, h2o .. (etc - need to check).

We use .txt files. It outputs radiance over a wavenumber range (that you can specify).

Use the script run_lbl_example_pres_ERA5.py
-You will need to change the path input
-Depending on your type of profile you will need to make sure its reading in each variable correctly and specify the number of layers etc
-The flags of the various gases (flags let the code know which units you are using)
-The path line of saving of the TAPE files (these are the outputs)

To understand the code you wil want to familiarise yourself with:
write_tape_5_CFC_pres_ERA5.py

This code refers to various records which are explained on this website:
https://www2.atmos.umd.edu/~ytma/lblrtm/lblrtm_instructions.html

# Fitting to FINESSE
magic all in one file.py is for applying fft and apodisation to your LBLRTM output such that you can accurately compare simulations to observations!! 
