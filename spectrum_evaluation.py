# -*- coding: utf-8 -*-
"""
Created on Fri Mar 29 10:26:57 2019
Handbook of X-Ray Spectrometry, 2nd Ed., ISBN:  0-8247-0600-5
Chapter 2, Spectrum Evalutaion, P. van Espen
Fortran code converted into Python
@author: windover
"""

import numpy as np
import scipy as sp
import scipy.signal as signal
from os import walk
import bruker_io as bruker_io
import hyperspy as hs
import copy
import pandas as pd

# Savitsky and Golay Poly Smoothing (pg 315 in Fortran)
#
# Input:  Y          Original Spectrum
#         NCHAN      Number of channels in the spectrum
#         ICH1,ICH2  First and last channel number to be smoothed
#         IWID       Width of the filter (2m+1), IWDI<52
# Output: S          Smoothed spectrum, only defined between ICH1 & ICH2

def SGSMITH(Y, NCHAN, ICH1, ICH2, IWID):
    #calculate filter coefficients
    #
    IW = np.min([IWID, 51])
    C = np.zeros(IW)   #changed this to the width of the filter
    M = np.int((IW-1)/2)   # needed to add a -1 term from Fortran
    SUM = (2*M-1)*(2*M+1)*(2*M+3)  #no change
    for j in np.arange(IW): #index from 0 to 2M instead of -M to M
        C[j] = 3*(3*M**2 + 3*M-1-5*(j-M)**2)  #need to subtract M
    #convolute spectrum with filter
    #
    JCH1 = np.max([ICH1, M])
    JCH2 = np.min([ICH2, NCHAN-1-M])
    S = np.zeros(NCHAN)   #initialize and empty filtered spectrum
    for i in np.arange(JCH1, JCH2): #this is a subset of the full spectra
        for j in np.arange(IW):  # this is 0 to 2m
            S[i] = S[i] + C[j]*Y[i+(j-M)]  # need to subtract M
        S[i] = S[i]/SUM
    return S

# Peak stripping - SNIP algorithm (pg 319 in Fortran)
#
# Input:  Y          Spectrum
#         NCHAN      Number of channels
#         ICH1,ICH2  First and Last channels of region for continuum calc
#         FWHM       Width parame for smoothing &stripping algorithm
#                    set to ave. FWHM of peaks (typical 8)
#         NITER      Number of iterations for SNIP algorithm
# Output: YBACK      Caltculated continuum over ICH1 to ICH2
# comment:  uses subroutine SGSMITH



def SNIPBG(Y, NCHAN, ICH1, ICH2, FWHM, NREDUC, NITER):
    #Smooth spectrum
    IW = np.int(FWHM)
    I1 = np.max([ICH1-IW, 0])
    I2 = np.min([ICH2+IW, NCHAN-1])
    YBACK = SGSMITH(Y, NCHAN, I1, I2, IW)
    zeros = np.zeros(NCHAN)
    test = np.zeros(NCHAN)
    #Square root transformation over region
    YBACK = np.sqrt(np.maximum(YBACK, zeros))
    #for i in np.arange(I1,I2):
    #    YBACK[i] = np.sqrt(np.max([YBACK[i],0]))
        #YBACK[i] = np.log(np.log(np.sqrt(np.max([YBACK[i],0]) +1)+1)+1)
    #Peak stripping
    REDFAC = 1
    for n in np.arange(0, NITER):
        if n+1 > NITER-NREDUC:
            REDFAC = REDFAC/np.sqrt(2)
        IW = np.max([np.int(REDFAC*FWHM), 1])
        for i in np.arange(ICH1, ICH2):
            I1 = np.int(np.max([i-IW, 0]))
            I2 = np.int(np.min([i+IW, NCHAN-1]))
            test[i] = 0.5*(YBACK[I1]+YBACK[I2])
            YBACK[i] = np.minimum.reduce([YBACK[i], test[i]])
    YBACK = np.square(YBACK)
    #for i in np.arange(I1,I2):
    #    YBACK[i] = np.square(YBACK[i])
        #YBACK[i] = np.square(np.exp(np.exp(YBACK[i])-1)-1)-1
        #YBACK[i] = np.square(np.max([YBACK[i],0]))
    return YBACK


# TOPHAT TOPHAT filtering protram (pg 323)
#
# Input:  IN        Spectrum
#         NCHAN     Number of channels
#         IFIRST    First channel of region for top hat filter
#         ILAST     Last channel of region for top hat filter
#         IWIDTH    Width parameter for tophat region
#         MODE      0 or !0 to set filtered spectra of weights
# OUTPUT  OUT       Output spectrum or weights
# Has two modes :
#    Mode = 0 (calculates the filtered spectrum)
#    Mode ne 0 (calculates weights)

def TOPHAT(IN, NCHAN, IFIRST, ILAST, IWIDTH, MODE):
    # Mode = 0 (calculate filtered spctrum)
    # Mode != 0 (calculate weights)
    # Calculate filter constants
    #print("mode: ", MODE)
    OUT = np.zeros(NCHAN)
    IW = IWIDTH
    # makes sure IW is odd by adding 1 to any even IW
    if np.mod(IW, 2) == 0: IW = IW + 1
    FPOS = 1./np.float(IW)
    #print("FPOS: ", FPOS)
    KPOS = np.int(IW/2 -0.5)
    IV = 2*np.int(IW/2 - 0.5)
    FNEG = -1./np.float(2*IV)
    #print("FNEG: ", FNEG)
    KNEG1 = np.int(IW/2 +1)
    KNEG2 = np.int(IW/2 +IV)
    N = 0
    #loop over all requested channels
    for channel in np.arange(IFIRST, ILAST, 1):
        #central positive part
        YPOS = 0
        YNEG = 0
        for xPOS  in np.arange(-KPOS, KPOS +1, 1):
            IK = np.min([np.max([channel+xPOS, 1]), NCHAN])
            YPOS = YPOS + IN[IK]
        # left and right negagive part
        for xNEG in np.arange(KNEG1, KNEG2+1, 1):
            IK = np.min([np.max([channel-xNEG, 1]), NCHAN])
            YNEG = YNEG + IN[IK]
            IK = np.min([np.max([channel+xNEG, 1]), NCHAN])
            YNEG = YNEG + IN[IK]
        #calc filtered spectra
        if MODE == 0:
            OUT[channel] = FPOS*YPOS + FNEG*YNEG
        #calc variance of the spectra
        else:
            VAR = FPOS*FPOS*YPOS + FNEG*FNEG*YNEG
            OUT[channel] = 1/np.max([VAR, 1])
        N = N +1
    return OUT




# TOPHATFAST TOPHAT filtering protram (Mofiied by DW 20190403)
#
# Input:  IN        Spectrum
#         NCHAN     Number of channels
#         IFIRST    First channel of region for top hat filter
#         ILAST     Last channel of region for top hat filter
#         IWIDTH    Width parameter for tophat region
#         MODE      0 or !0 to set filtered spectra of weights
# OUTPUT  OUT       Output spectrum or weights
# Has two modes :
#    Mode = 0 (calculates the filtered spectrum)
#    Mode ne 0 (calculates weights)


def TOPHATFAST(IN, NCHAN, IWIDTH, MODE):
    if np.mod(IWIDTH, 2) == 0: IWIDTH = IWIDTH + 1
    ones = np.ones(NCHAN)
    KNEG = np.int(IWIDTH/2 -0.5)
    #print('KNEG: ', KNEG)
    KPOS = np.int(IWIDTH -1)
    #print('KPOS: ', KPOS)
    hatpos = np.zeros(KPOS) +1/KPOS
    hatneg = np.zeros(KNEG) -1/(2*KNEG)
    tophat = np.hstack((hatneg, hatpos, hatneg))
    #print('tophat: ',tophat)
    hatpos_var = np.zeros(KPOS) +(1/KPOS)**2
    hatneg_var = np.zeros(KNEG) + (-1/(2*KNEG))**2
    tophat_var = np.hstack((hatneg_var, hatpos_var, hatneg_var))
    #print('tophat_var: ',tophat_var)
    #print('tophat: ', np.sum(tophat))
    #print(tophat)
    # call to scipy.signal for S.G. filter
    IN = signal.savgol_filter(IN, IWIDTH, 2)
    if MODE == 0:
        OUT = signal.convolve(IN, tophat, mode='same', method = 'direct')
    #calc variance of the spectra
    else:
        OUT = signal.convolve(IN, tophat_var, mode='same', method = 'direct')
        OUT = 1/np.maximum(OUT, ones)
    return OUT





# Peak stripping - SNIP algorithm (Modified by DW 20190403)
# NOTE:  instead of altering YBACK as it is used in each loop, we convolve
#        a +/-IW 1 function with the YBACK and find the min of this with YBACK
#        this removes the asymmetric bias in the Pg 319 implementation of acting
#        on YBACK during its used in a NCHAN for loop.  This method requires
#        more loops for stripping, but is much faster 10ms verus 3s
#
# Input:  Y          Spectrum
#         NCHAN      Number of channels
#         ICH1,ICH2  First and Last channels of region for continuum calc
#         FWHM       Width parame for smoothing &stripping algorithm
#                    set to ave. FWHM of peaks (typical 8)
#         NITER      Number of iterations for SNIP algorithm
# Output: YBACK      Caltculated continuum over ICH1 to ICH2
# comment:  uses subroutine SGSMITH


def SNIPFAST (Y, NCHAN, FWHM, NREDUC, NITER):
    REDFAC = 1
    if np.mod(FWHM,2) == 0: FWHM = FWHM + 1
    #Smooth spectrum using scipy.signal function of S.G.
    YBACK = Y
    #YBACK = signal.savgol_filter(Y, FWHM, 2 )
    #initialize two NCHAN arrays for calulation of sums
    zeros = np.zeros(NCHAN)
    YBACK_sum = np.zeros(NCHAN)

    # we are using a non linear square/square root scaling circa Van espen    
    YBACK = np.sqrt(np.maximum(YBACK, zeros))
    #alternative scaling not used
    #YBACK = np.log(np.log(np.sqrt(np.maximum(YBACK,zeros) +1) +1) +1) 
    # for loop for number of times we wish to strip background
    for n in np.arange(0, NITER):
        #after a suffiicient number of loops, reduce 'FWHM' by sqrt(2) until = 1
        if n+1 > NITER-NREDUC:
            REDFAC = REDFAC/np.sqrt(2)
        #Allow for a reduction in 'FWHM' over loops with a minimum IW of 1
        IW = np.max([np.int(REDFAC*FWHM),1])
        #make a function [1,0,0,0,...1] width 2W+1 for convolve
        straddle = np.zeros(2*IW+1)
        straddle[0]=1
        straddle[-1]=1
        # use scipy.signal.convolve to determine +/- average for test
        YBACK_sum = 0.5*signal.convolve(YBACK, straddle, mode='same')
        YBACK = np.minimum(YBACK, YBACK_sum) 
    # we are using a non linear square/square root scaling circa Van espen
    YBACK = np.square(YBACK)
    
    #alternative scaling not used
    #YBACK = np.square(np.exp(np.exp(YBACK)-1)-1)-1
    return YBACK




def pulse_pileup_removal(fittingdata):
    """Removal tool for first-order pulse-pileups in XRF data.
    
    Parameters
    ----------
    
    energy_scale : array of MCA energy at each channel (in eVs)
    channels : array from MCA giving counts in each channel
    shaping_time : events per second that can be processed
    

    Examples
    --------
    
    """
    pos_channels = \
    fittingdata.channels[np.nonzero(fittingdata.energy_scale>0)[0][0]:-1]
    pos_channels_per_s = pos_channels/(fittingdata.life_time_in_ms/1000)
    pileup_sum = np.zeros(len(pos_channels))   
    for i in np.arange(100,len(pos_channels)):
        forward = pos_channels_per_s[0:i]
        reverse = np.flip(forward, axis=0)
        shape_factor = (0.006/fittingdata.shaping_time)
        pileup = shape_factor*(forward)*(reverse)
        pileup_sum[i] = sum(pileup)
    fittingdata.channels[np.nonzero(fittingdata.energy_scale>0)[0][0]:-1] = \
    (pos_channels_per_s - pileup_sum) * (fittingdata.life_time_in_ms/1000)
    return 
    
       
def SCALEDSNIP(fittingdata):
    pos_energy_scale = \
    fittingdata.energy_scale[np.nonzero(fittingdata.energy_scale>0)[0][0]:-1]
    pos_channels = \
    fittingdata.channels[np.nonzero(fittingdata.energy_scale>0)[0][0]:-1]
    spectrum_function = \
    sp.interpolate.interp1d(np.sqrt(pos_energy_scale),
                            pos_channels, kind = 'linear',
                            fill_value = (0,0), bounds_error = False )
    energy_scale_sqrt = np.arange(0,np.sqrt(max(pos_energy_scale)),
                                  np.sqrt(max(pos_energy_scale))
                                  /len(pos_channels))
    channels_sqrt = spectrum_function(energy_scale_sqrt)
    data_bg = SNIPFAST(channels_sqrt, len(channels_sqrt), 13, 10, 1000)
    #data_bg = SNIPBG(channels_sqrt, len(channels_sqrt), 0, len(channels_sqrt)-1, 13, 10, 1000)
    new_data_corr = channels_sqrt - data_bg
    spectrum_function_squared = \
    sp.interpolate.interp1d(np.square(energy_scale_sqrt), 
                            new_data_corr, kind = 'linear', 
                            fill_value = (0,0), bounds_error = False )
    channels_corr = spectrum_function_squared(pos_energy_scale)
    channels_corr = channels_corr.clip(min=0)
    fittingdata.channels[np.nonzero(fittingdata.energy_scale>0)[0][0]:-1] = \
    channels_corr
    return


def polycap_remove(fittingdata):
    number_of_points = 40 #changed from 40 on 20190725 dw
    scaling = int(4000/number_of_points)
    bg_energy_scale = np.zeros(number_of_points)
    bg_channels = np.zeros(number_of_points)
    for i in np.arange(bg_energy_scale.shape[0]):
        bg_energy_scale[i] = fittingdata.energy_scale[i*scaling + 100]
        bg_channels[i] = np.min(fittingdata.channels[i*scaling + 50:(i+1)*scaling+50])
    bg_function = \
    sp.interpolate.interp1d(bg_energy_scale, bg_channels, kind = 'cubic',
                            fill_value = (0,0), bounds_error = False)
    bg_intensity = bg_function(fittingdata.energy_scale)
    bg_intensity = bg_intensity.clip(min=0)
    bg_corrected = fittingdata.channels - bg_intensity
    fittingdata.channels = bg_corrected.clip(min=0)
    return
        



def spectra_fit(directory_path, fitter, method, elements):
    files = []
    spx_files =[]
    roi_data = elements
    model_data = elements
    print(elements)
    life_time_in_ms = []
    for (dirpath, dirnames, filenames) in walk(directory_path):
        files.extend(filenames)
        break
    for file in files:
        if '.spx' in file:
            spx_files.append(file) 
    for file in spx_files:
        print(file)
        spx_file = file
        spx = bruker_io.FittingData(directory_path + '/' + spx_file)
        bruker_io.bruker_spx_import(spx)
        pulse_pileup_removal(spx)
        SCALEDSNIP(spx)
        polycap_remove(spx)
        hsEDS = hs.signals.EDSSEMSpectrum(spx.channels)
        hsEDS.set_microscope_parameters(50000)
        hsEDS.axes_manager[0].name = 'XRF spectra'
        hsEDS.axes_manager[0].offset = spx.calibration_abs
        hsEDS.axes_manager[0].scale = spx.calibration_lin
        hsEDS.axes_manager[0].units = 'eV'
        hsEDS.add_elements(elements)
        hsEDS.add_lines()
        line_names = hsEDS.metadata.Sample.xray_lines
        #print(line_names)
        mod = hsEDS.create_model()
        mod.remove('background_order_6')
        new_roi = np.zeros(len(elements))
        new_model = np.zeros(len(elements))
        mod.fit(fitter= fitter, method= method)
        for i in np.arange(len(elements)):
            new_roi[i] = np.float(hsEDS.get_lines_intensity([line_names[i]])[0].data[0])
            model_call = ''.join(['mod.components.',line_names[i],'.A.value'])
            new_model[i] = eval(model_call)
            #test_param[i] = ''.join(['mod.components.', line_names[i], '.A.value'])
            #new_model[i]= np.float(test_param)
        #Hf_La_model.append(np.float())
        #Si_Ka_model.append(np.float(mod.components.Si_Ka.A.value))
        #Hf_Si_ratio.append(np.float(mod.components.Hf_La.A.value)/np.float(mod.components.Si_Ka.A.value))
        life_time_in_ms.append(spx.life_time_in_ms)
        roi_data = np.vstack((roi_data,new_roi))
        model_data = np.vstack((model_data,new_model))
    roi_data = roi_data[1:,:]
    model_data = model_data[1:,:]    
    #spx_files = np.array(spx_files)
    #print(spx_files)
    #spx_files = spx_files[np.newaxis]    
    #spx_files = spx_files.transpose
    #spx_files.shape
    #element_data.shape
    #columns_df = ['filename'].append(line_names)
    #data_df = np.hstack((spx_files, element_data))   
    roi_df = pd.DataFrame(data=roi_data, columns = line_names)
    roi_df.insert(0,'filename', spx_files)
    roi_df['life time in ms'] = life_time_in_ms
    model_df = pd.DataFrame(data=model_data, columns = line_names)
    model_df.insert(0,'filename', spx_files)
    model_df['life time in ms'] = life_time_in_ms
    return roi_df, model_df

#def model_lookup(mod, line_name):
#    model_call = {'N_Ka': mod.components.N_Ka.A.value,
#                  'O_Ka': mod.components.O_Ka.A.value,
#                  'F_Ka': mod.components.F_Ka.A.value,
#                  'Ne_Ka': mod.components.Ne_Ka.A.value,
#                  'Na_Ka': mod.components.Na_Ka.A.value, 
#                  'Mg_Ka': mod.components.Mg_Ka.A.value,
#                  'Al_Ka': mod.components.Al_Ka.A.value,
#                  'Si_Ka': mod.components.Si_Ka.A.value,
#                  'P_Ka': mod.components.P_Ka.A.value,
#                  'S_Ka': mod.components.S_Ka.A.value,
#                  'Cl_Ka': mod.components.Cl_Ka.A.value,
#                  'Ar_Ka': mod.components.Ar_Ka.A.value,
#                  'K_Ka': mod.components.K_Ka.A.value,
#                  'Ca_Ka': mod.components.Ca_Ka.A.value,
#                  'Sc_Ka': mod.components.Sc_Ka.A.value,
#                  'Ti_Ka': mod.components.Ti_Ka.A.value,
#                  'V_Ka': mod.components.V_Ka.A.value,
#                  'Cr_Ka': mod.components.Cr_Ka.A.value,
#                  'Mn_Ka': mod.components.Mn_Ka.A.value,
#                  'Fe_Ka': mod.components.Fe_Ka.A.value,
#                  'Co_Ka': mod.components.Co_Ka.A.value,
#                  'Ni_Ka': mod.components.Ni_Ka.A.value,
#                  'Cu_Ka': mod.components.Cu_Ka.A.value,
#                  'Zn_Ka': mod.components.Zn_Ka.A.value,
#                  'Ga_Ka': mod.components.Ga_Ka.A.value,
#                  'Ge_Ka': mod.components.Ge_Ka.A.value,
#                  'As_Ka': mod.components.As_Ka.A.value,
#                  'Se_Ka': mod.components.Se_Ka.A.value,
#                  'Br_Ka': mod.components.Br_Ka.A.value,
#                  'Kr_Ka': mod.components.Kr_Ka.A.value,
#                  'Rb_Ka': mod.components.Rb_Ka.A.value,
#                  'Sr_Ka': mod.components.Sr_Ka.A.value,
#                  'Y_Ka': mod.components.Y_Ka.A.value,
#                  'Zr_Ka': mod.components.Zr_Ka.A.value,
#                  'Nb_Ka': mod.components.Nb_Ka.A.value,
#                  'Mo_Ka': mod.components.Mo_Ka.A.value,
#                  'Tc_Ka': mod.components.Tc_Ka.A.value,
#                  'Ru_Ka': mod.components.Ru_Ka.A.value,
#                  'Pd_Ka': mod.components.Pd_Ka.A.value,
#                  'Ag_Ka': mod.components.Ag_Ka.A.value,
#                  'Cd_Ka': mod.components.Cd_Ka.A.value,
#                  'In_Ka': mod.components.In_Ka.A.value,
#                  'Sn_Ka': mod.components.Sn_Ka.A.value,
#                  'Sb_Ka': mod.components.Sb_Ka.A.value,
#                  'Te_Ka': mod.components.Te_Ka.A.value,
#                  'I_Ka': mod.components.I_Ka.A.value,
#                  'Xe_Ka': mod.components.Xe_Ka.A.value,
#                  'Cs_Ka': mod.components.Cs_Ka.A.value,
#                  'Ba_Ka': mod.components.Ba_Ka.A.value,
#                  'La_Ka': mod.components.La_Ka.A.value,
#                  'Hf_La': mod.components.Hf_La.A.value,
#                  'Ta_La': mod.components.Ta_La.A.value,
#                  'W_La': mod.components.W_La.A.value,
#                  'Re_La': mod.components.Re_La.A.value,
#                  'Os_La': mod.components.Os_La.A.value,
#                  'Ir_La': mod.components.Ir_La.A.value,
#                  'Pt_La': mod.components.Pt_La.A.value,
#                  'Au_La': mod.components.Au_La.A.value,
#                  'Hg_La': mod.components.Hg_La.A.value,
#                  'Tl_La': mod.components.Tl_La.A.value,
#                  'Pb_La': mod.components.Pb_La.A.value,
#                  'Bi_La': mod.components.Bi_La.A.value,
#                  'U_La': mod.components.U_La.A.value,
#                  'xx_La': mod.components.xx_La.A.value}
#    call = model_call[line_name]
#    return call





#Work in progress - it performs a tophat on the fly, and is very
#confusing to read through 20190401 daw    
## Peak Search - LOCPEAKS uses tophyhat filter
##
## Input:  Y           Spectrum
##         NCHAN       Number of channels
##         R           Peak search sensitivity factor (typcially 2 to 4)
##         IWID        Width of the filter, approx, FWHM of the peaks
##         MAXP        Maxiumum number of peaks allowed (sets max arrays)
## OUtput: NPEAK       Number of peaks found
##         IPOS        Array of peak positions
#
#def LOCPEAKS (Y, NCHAN, IWID, R, MAXP):
#    # Width of filter (number of channels in the tophat)
#    # must be odd and at least 3
#    NP = np.int(np.max([(IWID/2)*2 + 1, 3]))
#    print ("NP: ", NP)
#    NPEAKS = 0 #initiallize the number of peaks 
#    # Calculate the half width and start and stop channel
#    N = np.int(NP/2)
#    I1 = NP
#    I2 = NCHAN - NP
#    # INisitalize the running sums
#    I = I1
#    TOTAL = 0
#    TOP = 0
#    for i in np.arange(2*NP+1):
#        TOTAL = TOTAL + Y[i]
#    print("TOTAL: ", TOTAL)
#    for i in np.arange(2*N+1):
#        TOP = TOP + (i)
#    print("TOP: ", TOP)
#    #Loop ofver all channels
#    LASTPOS = 0
#    SENS = R**2
#    FI = 0
#    FNEXT = 0
#    SNEXT = 0
#    for i in np.arange(I2-I1):
#        TOP = TOP - Y[i+I1 - N] + Y[i+I1 + N]
#        TOTAL = TOTAL - Y[i+I1-NP] + Y[i+I1+NP]
#    print("TOP: ", TOP)
#    print("TOTAL: ", TOTAL)
#    return