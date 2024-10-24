"""
Run LBLRTM 
Authour: Laura Warwick
Written: 25/05/2019
WORKING I THINK?THIS IS THE ONE WE USE MAY 2024
/disk1/sm4219/lblrtm_12.13_sgl/lblrtm_v12.13_linux_intel_sgl
"""
import shutil
import os
import numpy as np
import subprocess as sub
import write_tape_5_CFC_pres_ERA5 as wr
import glob as glob

# This class allows you to change the current working
# directory as you would in linux. IT is needed to run 
# LBLRTM later on. You coud probably
# do this yourself using os if you prefered
class cd:
    """Context manager for changing the current working directory
    from https://stackoverflow.com/questions/431684/how-do-i-change-directory-cd-in-python/13197763#13197763"""
    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)


# Specify location to save the output files
save_location = "/net/sirocco/disk1/sm4219/LBLRTM_FOR_SOPHIE/"

# profile_name = "/disk1/sm4219/GIT/ANDOYA_NEW/2023-03-08T09:00:00.000000000_ERA5_profile.txt"
# profile_name ="/disk1/sm4219/GIT/ANDOYA_NEW/2023-03-08T08:00:00.000000000_ERA5_profile_NEW.txt"
# profile_name= "/disk1/sm4219/GIT/ANDOYA_NEW/2023-03-08T08:00:00.000000000_ERA5_profile_edited.txt"
# profile_name= "/disk1/sm4219/GIT/ANDOYA_NEW/updated_data.txt"
# profile_name = "/disk1/sm4219/GIT/ANDOYA_NEW/updated_data_8am.txt"
# profile_name = "/disk1/sm4219/GIT/ANDOYA_NEW/RADIOSONDE_profile_8am.txt"
# profile_name = "/disk1/sm4219/GIT/ANDOYA_NEW/2023-11-08T08:00:00.000000000_ERA5_profile_fix_sh.txt"

# 8 am path
# profile_name = "/disk1/sm4219/GIT/ANDOYA_NEW/2023-03-08T08:00:00.000000000_ERA5_profile_fix_sh_08.txt"
# profile_name = "/disk1/sm4219/GIT/ANDOYA_NEW/2023-03-08T08:00:00.000000000_ERA5_profile_fix_sh_08_fixfirstlineahh.txt"
# profile_name = "/disk1/sm4219/GIT/ANDOYA_NEW/2023-03-08T08:00:00.000000000_ERA5_SPH_FIXED.txt"
# profile_name = "/disk1/sm4219/GIT/ANDOYA_NEW/2023-03-08T08:00:00.000000000_ERA5_profile_fix_sh_new_k.txt"
# profile_name = "/disk1/sm4219/GIT/ANDOYA_NEW/2023-03-08T07:00:00.000000000_ERA5_profile_fix_sh_new_k.txt"
# /disk1/sm4219/GIT/ANDOYA_NEW/2023-03-08T08:00:00.000000000_ERA5_profile_NEW.txt
# /disk1/sm4219/GIT/ANDOYA_NEW/2023-03-08T08:00:00.000000000_ERA5_profile.txt

# profile_name = "/disk1/sm4219/GIT/ANDOYA_NEW/RADIOSONDE_profile_8am_25april_cutup.txt"
# profile_name = "/disk1/sm4219/GIT/ANDOYA_NEW/2023-03-08T08:00:00.000000000_ERA5_profile_edited.txt"
# profile_name = "/disk1/sm4219/GIT/ANDOYA_NEW/2023-03-08T08:00:00.000000000_ERA5_profile_SANjee.txt"
# profile_name = "/disk1/sm4219/GIT/ANDOYA_NEW/2023-11-08T08:00:00.000000000_ERA5_profile_SH_MAY_new.txt"
# profile_name = "/disk1/sm4219/GIT/ANDOYA_NEW/2023-03-08T07:00:00.000000000_ERA5_profile_RH_MAY_new.txt"
# Contains: Altitude (km), Pressure (hPa), Temperature (K), h2o (RH), co2 (ppmv), o3 (gkg-1)
"HERE IS MOST CURRENT PROFILE USED"
# profile_name  = "/disk1/sm4219/GIT/ANDOYA_NEW/RADIOSONDE_profile_8am_PURE_NEW_2_edit.txt"
# profile_name = "/disk1/sm4219/Andoya_key/ERA5_8am_rh_pure.txt"
profile_name = "/disk1/sm4219/Andoya_key/ERA5_8am_RH_24oct.txt"
# profile_name = "/disk1/sm4219/Andoya_key/ERA5_8am_sh_TRIAL.txt"
# profile_name = "/disk1/sm4219/Andoya_key/ERA5_8am_PLUS_hatpro_2ndoct.txt"
# profile_name = "/disk1/sm4219/Andoya_key/ERA5_8am_PLUS_hatpro_3rdoct.txt"

# profile_name = "/disk1/sm4219/Andoya_key/ERA5_8am_sh_cut.txt"
# profile_name = "/disk1/sm4219/Andoya_key/ERA5_8am_SH_15oct.txt"
# profile_name = "/disk1/sm4219/Andoya_key/ERA5_8am_v1_cut.txt"
# profile_name = "/disk1/sm4219/Andoya_key/ERA5_8am_rh_cut.txt"

# profile_name = "/disk1/sm4219/GIT/ANDOYA_NEW/2023-03-08T07:00:00.000000000_ERA5_profile_RH_test.txt"
# profile_name = "/disk1/sm4219/GIT/ERA5_ENSEMBLE/2023-03-08T09:00:00.000000000_ERA5_ensemble_means_profile.txt"
# profile_name = "/disk1/sm4219/GIT/ANDOYA_NEW/RADIOSONDE_profile_8am_25april_PURE.txt"
# profile_name = "/disk1/sm4219/GIT/ERA5_ENSEMBLE/ENSEMBLE__profile_9_VAS.txt"
# profile_name = "/disk1/sm4219/GIT/ANDOYA_NEW/RADIOSONDE_profile_8am_PURE_NEW_2_edit.txt"

# profile_name = "/disk1/sm4219/GIT/ANDOYA_NEW/2023-03-08T09:00:00.000000000_ERA5_profile_reanalyis_MEANcoords.txt"
# profile_name = "/disk1/sm4219/GIT/ANDOYA_NEW/2023-03-08T09:00:00.000000000_ERA5_profile_RH_newcoords.txt"
# profile_name = "/disk1/sm4219/GIT/ANDOYA_NEW/RADIOSONDE_profile_8am_PURE_aug13.txt"

# profile_name = "/disk1/sm4219/GIT/ANDOYA_NEW/2023-03-08T09:00:00.000000000_ERA5_profile_RH_newcoords.txt"
# profile_name = "/disk1/sm4219/GIT/ANDOYA_NEW/_ERA5_point_DIFFERENT_TIMES_1.txt"
# profile_name = "/disk1/sm4219/GIT/ANDOYA_NEW/2023-03-08_ERA5_DIFFERENTTIMES_composite_trackingrs.txt"


# /disk1/sm4219/GIT/ANDOYA_NEW/RADIOSONDE_profile_8am_PURE_NEW.txt
# profile_name = "/disk1/sm4219/GIT/ANDOYA_NEW/2023-03-08_ERA5_DIFFERENTTIMES_composite_edit.txt"
z, pressure, temp, h2o,co2,o3 = np.loadtxt(profile_name, unpack=True)
# Specify location of LBLRTM executable
# lbl_location = "/net/sirocco/disk1/sm4219/lblrtm_12.13/"
lbl_location = "/net/sirocco/disk1/sm4219/lblrtm_12.13/"
# /disk1/sm4219/lblrtm_12.13/lblrtm_v12.13_linux_intel_dbl
# Specify start position of run, angle, wn range,
# supplementary atmosphere, mode, blackbody




# h_start = pressure[0]  # Radiation calculation starts from altitude in km
# h_obs = pressure[-1] # 9.479191666666669107e+02 #961.386  # Km height of observation (in this case aircraft)

h_start =  pressure[-1] # Radiation calculation starts from altitude in km
h_obs = pressure[0] # 9.479191666666669107e+02 #961.386  # Km height of observation (in this case aircraft)

# ANGLE WAS 9 BEFORE / SM
angle = 0 
wn_range = [300,1600]
atm = 5  # 5 Sub Arctic winter.
        # See write_tape_5.py instructions for types of atmosphere
mode = 1  # 1 = radiance, 0 = transmission
OD = 0  # 0 = no optical depth files
blackbody_surface = True  # Assuming black body surface
surface_temp = 2  # Surface temp K it was 2 before
print(h_start, h_obs, angle)
# n_levels = 33
 
# print('HELLO', len(z))
# n.b. the number of levels must be a multiple of 8!!

 # This is the step where the TAPE5 file is written
 # The t5 keyword MUST save the TAPE5 file to the folder where the
# h2o_g = h2o*1000
# print('altitude:', np.flip(z))
# print('pressure', pressure)
# Select every 100th point starting from the first point
# num_samples = 100
# indices = np.unique(np.logspace(0, np.log10(len(z) - 1), num_samples, endpoint=True, base=10, dtype=int))

# # Sort the indices to keep them in order
# indices.sort()
# num_selected_indices = min(num_samples, len(indices))

# # Select elements from each variable using the generated indices
# z = z[::20]
# z = z[indices]

# num_skip = 30 
# z = z[::num_skip]
# pressure = pressure[::num_skip]
# temp = temp[::num_skip]
# h2o = h2o[::num_skip]
# co2 = co2[::num_skip]
# o3 = o3[::num_skip]

n_levels = 35
# continuum2 = [0, 0, 0, 0, 0, 0, 0] 
# specify 'H' for relative humidity, 'C' for mass mixing ratio (g/kg), else ppmv selected
 # LBLRTM executable is located with the name TAPE5 and nothing else
wr.write_tape5(
    (z/1000),
    (pressure),
    temp,
    h2o,
    surface_temp,
    h_obs,
    h_start,
    wn_range,
    angle,
    n_levels,
    atm,
    mode,
    pro_o3 = o3,
    # CHECK NO. of ZEROS 1000 too big
    pro_oco = co2,
    od=OD,
    t5="/net/sirocco/disk1/sm4219/lblrtm_12.13/TAPE5",
    # /disk1/sm4219/lblrtm_v12.13_linux_intel_dbl
    h2o_flag = 'C',
    o3_flag = 'C',
    blackbody=blackbody_surface,
    continuum=False)
print("Starting LBLRTM")
# Change folder to where the LBLRTM executable is located
# save_location = "/disk1/sm4219/LBLRTM_FOR_SOPHIE/from_jon_READODINT"
with cd(lbl_location):
    # Run LBLRTM (the string should be the name of the LBLRTM executable)
    print('got here')
    sub.call("lblrtm_v12.13_linux_intel_dbl")
    # sub.call("lblrtm_v12.13_linux_intel_sgl")
    # print('call ex')
    # Move all files somewhere to be saved
    # IF this step is skipped, LBLRTM will overwrite these
    # files when it is run again
    os.rename("TAPE5", save_location + "/TAPE5")
    print('written tap5')
    os.rename("TAPE6", save_location + "/TAPE6")
    os.rename("TAPE7", save_location + "/TAPE7")
    os.rename("TAPE11", save_location + "/TAPE11")
    os.rename("TAPE12", save_location + "/TAPE12")
    print('written tap12 now')
    # Also move OD files if created
    # if OD == 1:
    #     save_location = "/disk1/sm4219/LBLRTM_FOR_SOPHIE/D10/"
    #     optical_files = glob.glob("OD*")
    #     for o_file in optical_files:
    #         os.rename(o_file, save_location + o_file)


    if OD == 1:
        save_location = "/disk1/sm4219/LBLRTM_FOR_SOPHIE/from_jon_READODINT/9am_newcoords_OD/"
        optical_files = glob.glob("OD*")
        for o_file in optical_files:
            shutil.move(o_file, save_location + o_file)



"""
idl -e read_lbl_file_example

"""