"Made by sophie with help from lots of people"
import numpy as np
from scipy.fftpack import fft, ifft
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import panel_file as panpy
from scipy.io import readsav
import numpy as np
import xarray
import pandas as pd



"STEP 1 LETS LOAD THE TAPE 12 FILE"
panel_data = panpy.panel_file('/net/sirocco/disk1/sm4219/LBLRTM_FOR_SOPHIE/TAPE12', do_load_data = True)

"your final output name"
path_new_out_final = '/net/sirocco/disk1/sm4219/LBLRTM_FOR_SOPHIE/LBLRTM_SIMGA_BOUNDS_40/aposided_RS_542layer_aug13.txt'

# Handy print statements 
# print(panel_data)
# print(panel_data.hdr)
# # print(panel_data.hdr.secant)

"""
If the data is loaded, then object attributes v, data1 are defined for a single panel file, 
and in addition, data2 for a double panel file. 
That is data1 will be radiance and data2 will be the tranmission. 
"""
# what does single panel and double panel mean??? - this bit above is from Laura

# WRITE CODE HERE TO ACCES THE RADIANCES
# print("Radiance values, data 1:", panel_data.data1)
# print("Transmission values, data 2:", panel_data.data2)

rad_in_raw = panel_data.data1
wn_in_raw = panel_data.v


"""
STEP 2 applying FFT to the thing
the original file for this code is: /disk1/sm4219/GIT/fft_for_sophie/fft_for_sophie.py

"""


def apodise_spectrum(
    frequency,
    radiance,
    fre_grid,
    st,
    ed,
    new_pd,
    apodisation_func=False,
    test_delta=False,
):
    """
    Apodise a high resolution spectrum using a boxcar or
    triangle function

    Adapted from apodise_spectra_boxcar_v1.pro
    ;
    ; Original Author: J Murray (14-Oct-2020)
    ;
    ; Additional comments by R Bantges
    ; Version 1: Original
    ;
    Params
    ------
    frequency array
        Original wavenumber scale (cm^-1)
    radiance array
        Original spectrum
    fre_grid float
        The frequency of the  output grid for the apodised spectra (cm^-1)
    st float
        Wavenumber to start apodised spectrum (cm^-1)
    ed float
        Wavenumber to end apodised spectrum (cm^-1)
    new_pd float
        Optical path difference i.e. width of boxcar to apodise (cm)
    apodisation_func string
        deafult=False
        Function to use in addition to boxcar to apodise the spectrum
        Options
        -------
        "triangle" - Triangle function, running from 1 at centre of interferogram
        to zero at edge of interferogram
    test_delta bool
        deafult=False
        If True, the spectrum is taken to be a delta function, can be
        used to test the apodisation. This should return the ILS which is a sinc
        function in the case of a boxcar
        If False input spectrum is used

    Returns
    -------
    wn array
        Wavenumber of apodised spectrum (cm^-1)
    radiance array
        Radiance or transmission of apodised spectrum
        (same units as input)
    """
    # Determine the number of samples making up the output spectra
    samples = int(np.round((ed - st) / fre_grid))

    # Define the wavenumber grid resolution (Fixed high resolution grid.
    # The Monochromatic spectra will be interpolated onto this grid for
    # convenience and potentially reduce time taken for the FFT, the arbitrary
    # number of points in the spectra can be such that it significantly slows
    # the FFT.
    # NB: 0.0001 cm-1 was chosen to resolve the spectral features in the
    # high resolution simulation
    dum_new_res = 0.0001
    dum_samples = int(np.round((ed - st) / dum_new_res))
    # The number of samples in the high res frequency scale

    # ********** Define the arrays for the re-interpolated radiance files **********
    # generate a wavenumber scale running from st - ed wavenumbers
    # at new_res cm-1
    new_fre = np.arange(st, ed, fre_grid)
    # generate a wavenumber scale running from st - ed wavenumbers at 0.001 cm-1
    dum_new_fre = np.arange(st, ed, dum_new_res)
    # ******************************************************************************

    # ********** Interpolate the high res radiance to new array scales **********
    f_dum_spec = interp1d(frequency, radiance)
    dum_spec = f_dum_spec(dum_new_fre)
    if test_delta:
        dum_spec = np.zeros_like(
            dum_spec
        )  # These can be set to produce a delta function to check the sinc
        dum_spec[int(15000000 / 2) : int(15000000 / 2) + 101] = 100.0
    # *****************************************************************************

    # FFT the interpolated LBLRTM spectrum
    int_2 = fft(dum_spec)
    # sampling=1./(2*0.01)/samples/100.   # Sampling interval of the interferogram in cm these are the same for the 0.001 and 0.01 spectra
    sampling = 1.0 / (2 * fre_grid) / samples / 100.0
    # Sampling interval of the interferogram in cm these are the same for the 0.001 and 0.01 spectra

    # ********** Apodise the LBLRTM sim and transform **********
    Q = int(
        round(new_pd / 100.0 / sampling / 2.0)
    )  # number of samples required to extend the path difference to 1.26cm
    # *****************************************************************************

    # Define an array to hold the folded out inteferogram
    int_1 = np.zeros(samples, dtype=np.cdouble)

    # 'int_2' - this interferogram is equivalent to a sampling grid of 0.001 cm-1
    # in the spectral domain, this statement applies a boxcar apodisation over +/-1.26 cm
    int_2[Q:-Q] = 0.0

    # The following two lines reduce the output spectra to a sampling grid of 0.01 cm-1
    # while copying in the truncated interferogram from the high resolution interferogram
    int_1[0 : int(round((samples / 2)))] = int_2[0 : int(round((samples / 2)))]
    int_1[int(round((samples / 2))) : samples] = int_2[
        (dum_samples) - int(round((samples / 2))) : dum_samples
    ]

    if apodisation_func == "triangle":
        print("Apodising with triangle")
        int_1_unapodised = np.copy(int_1)
        triangle_left = [1, 0]
        triangle_left_x = [0, Q]
        triangle_left_x_all = np.arange(len(int_1[0:Q]) + 1)
        f_triangle_left = interp1d(triangle_left_x, triangle_left)
        triangle_right = [0, 1]
        triangle_right_x = [len(int_1) - Q - 1, len(int_1)]
        triangle_right_x_all = np.arange(len(int_1) - Q - 1, len(int_1), 1)
        f_triangle_right = interp1d(triangle_right_x, triangle_right)

        int_1[0 : Q + 1] = int_1[0 : Q + 1] * f_triangle_left(triangle_left_x_all)
        int_1[-Q - 2 : -1] = int_1[-Q - 2 : -1] * f_triangle_right(triangle_right_x_all)

    elif not apodisation_func:
        print("Apodising with boxcar")

    else:
        print("No recognised function selected, defaulting to boxcar")

    new_lbl_spec = ifft(int_1)

    # ***********************************************************************
    apodised_spectra = np.real(new_lbl_spec / (fre_grid / dum_new_res))
    return new_fre, apodised_spectra


# Example with test spectrum
# wn_in, rad_in = np.loadtxt("/disk1/sm4219/GIT/fft_for_sophie/test_spectrum.txt", unpack=True)
# wn_in, rad_in  = np.loadtxt("/disk1/sm4219/LBLRTM_FOR_SOPHIE/LBLRTM_SIMGA_BOUNDS_40/TAPE12_andoya_9am.txt", unpack = True, dtype=np.float64) #137

# wn_in, rad_in  = np.loadtxt("/disk1/sm4219/LBLRTM_FOR_SOPHIE/LBLRTM_SIMGA_BOUNDS_40/TAPE12_andoya_8am_NEW.txt", unpack = True, dtype=np.float64) 

# wn_in, rad_in  = np.loadtxt("/disk1/sm4219/LBLRTM_FOR_SOPHIE/wavenumbers_radiance.txt",  unpack = True, dtype=np.float64) 
wn_in = wn_in_raw
rad_in = rad_in_raw

wn_out, rad_out = apodise_spectrum(wn_in, rad_in, 0.2, 300, 1600, 1.21)
# print(wn_out,rad_out)
# np.savetxt('/disk1/sm4219/LBLRTM_FOR_SOPHIE/apodized_wn_IDL_aug.txt',wn_out)
# np.savetxt('/disk1/sm4219/LBLRTM_FOR_SOPHIE/apodized_rad_IDL_aug.txt',rad_out)
print('DONE step 1 FFT')


"""
STEP 3: APPLYING APODIATION ils?
the original code is found at /disk1/sm4219/LBLRTM_FOR_SOPHIE/apodise_finesse.py

"""


# Here I have redefined the start and end fre such that it is an array, to get the bins 
# the function requires the inputs start_fre and end_fre to be arrays of wavenumber rather than float values
start_fre=np.array([400,450,560,630,730,850,950,1050,1150,1250,1360,1450,1550,1650,1750,1800,1900])
end_fre=np.array([450,560,630,730,850,950,1050,1150,1250,1360,1450,1550,1650,1750,1800,1900,1950])



# Didn't change anything within the def apply_ILS_sav
def apply_ILS_sav(ILS, start_fre, end_fre, wn, spectrum, pad_length=10):
    """Apply ILS to a spectrum

    Args:
        ILS (array (ILS, frequency bin)): ILS axis 0 is the
            ILS axis 1 is the frequency bin as defined by
            start_fre, end_fre
        start_fre (array): Start frequency for wn bin (cm-1)
        end_fre (array): end frequency for wn bin (cm-1)
        wn (array): wn scale of spectrum (cm-1)
        spectrum (array): spectrum
        padlength (int): amount to add to end of each wavenumber
            section to remove edge effects. Expressed in units of
            wavenumber
    """
    # Specify frequency scale of ILS
    ILS_frequency_scale = np.linspace(-5, 5, np.shape(ILS)[0])

    # Loop through each chunk of spectrum and apply the ILS
    # to that chunk
    for i in range(len(start_fre)):
        # , END_WN, ils_now in zip(
        #     ILS["start_fre_range"],
        #     ILS["end_fre_range"],
        #     ILS["ils_function"],
        # comment out lines 50-54 in Sophie's file as the wavenumber interpolation would no longer be necessary
        # Trim to correct chunk of spectrum
        # Add extra for convolution overlap
        index = np.where(
            np.logical_and(
                wn >= start_fre[i] - pad_length,
                wn <= end_fre[i] + pad_length,
            )
        )
        wn_now = wn[index]
        spectrum_now = spectrum[index]
        # Interpolate ILS onto frequency of signal COMMENT OUT LINES BELOW 
        ILS_frequency_scale_interp = np.arange(-5, 5, np.average(np.diff(wn_now)))
        ils_now_interp = np.interp(
            ILS_frequency_scale_interp, ILS_frequency_scale, ILS[:, i]
        )[::-1]
        # convolution below, REPLACED ILS_INTERP WITH WN_NOW

        spectrum_interp = np.convolve(
            spectrum_now,
            ils_now_interp,
            mode="same",
        ) / sum(ils_now_interp)


        # Trim so only spectrum in area of interest is
        # retained
        index_out = index = np.where(
            np.logical_and(
                wn_now >= start_fre[i],
                wn_now < end_fre[i],
            )
        )

        wn_out = wn_now[index_out]
        spectrum_out = spectrum_interp[index_out]
        # plt.cla()
        # plt.plot(wn_now, spectrum_now)
        # plt.plot(wn_now, spectrum_interp)
        # plt.plot(wn_out, spectrum_out)
        # plt.savefig(
        #     "/net/shamal/disk1/lw2214/Water emissivity results/"
        #     + "Water 3/test_convolve_1.png")
        if i == 0:
            wn_all = wn_out
            spectrum_all = spectrum_out
        else:
            wn_all = np.append(wn_all, wn_out)
            spectrum_all = np.append(spectrum_all, spectrum_out)
    return wn_all, spectrum_all

# ILS_LOCATION = '/net/fennec-store/disk2/lw2214/Bruker_Data/EM27_ILS.sav'
ILS_LOCATION = '/disk1/sm4219/EM27_ILS.sav'
ils = readsav(ILS_LOCATION)
ILS = ils['em27_ils'][:]

# HERE WE LOAD IN THE VARAIBLES FROM EARLIER FROM STEP 2 
wnT = wn_out
spectrumT = rad_out

print(wnT)

# WHY IS THIS START FRE AND END FRE DEFINED TWICE... 
start_fre=np.array([300,350,400,450,560,630,730,850,950,1050,1150,1250,1360,1450,1550])
end_fre=np.array([350,400,450,560,630,730,850,950,1050,1150,1250,1360,1450,1550,1600])

# EDIT IT HERE!!!

# print(start_fre[0],1950, np.shape(spectrumT)[0])

wnT = np.linspace(300,1600,np.shape(spectrumT)[0]) 
for hour in ['0800_only_radiosonde']: #['0200', '0300', '0400','0500','0600','0700']:
    path = '/net/sirocco/disk1/sm4219/LBLRTM_FOR_SOPHIE/LBLRTM_SIMGA_BOUNDS_40/'
    # wn_B,spectrum_B = np.loadtxt(path+'example_spectrum_joncopy.txt', unpack = True) 
    # /net/sirocco/disk1/sm4219/LBLRTM_FOR_SOPHIE/LBLRTM_SIMGA_BOUNDS_40/apodised_AFTERHAT_9am_comp_new.txt
    # spectrum2 = xarray.DataArray(data = spectrum_B*1E4, dims = ['wn'], 
    #  coords=dict(wn=wn_B))
    # spectrum = rad_jon
    spectrum = spectrumT
    wn = wnT
    # np.linspace(start_fre[0],1950,np.shape(spectrum)[0])     # changed end freq to the final number 
    apodised_wn, apodised_spectrum = apply_ILS_sav(ILS, start_fre, end_fre, wn,
                    spectrum, pad_length=10)
    
    # changed the number here to 8500, in the .interp
    # apodised_interp = xarray.DataArray(data = apodised_spectrum*1E4, dims = ['wn'], coords=dict(wn=apodised_wn)).interp(wn=np.linspace(300,2000,8500))


    np.savetxt(path_new_out_final,np.vstack([apodised_wn, apodised_spectrum]).T)
    # np.savetxt(path+'/apodised_spectra_NEW_old.txt',np.vstack([wn,apodised_interp]).T)

# podised_AFTERHAT_ANDoya.txt", unpack=True)
# print(len(apodised_interp), len(wn))
print('DONE step 2 APODISE')

print('All done and ready to plot', 'saved:', path_new_out_final)