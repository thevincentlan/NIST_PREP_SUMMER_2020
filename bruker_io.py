# -*- coding: utf-8 -*-
"""
.. module:: bruker_io
   :platform: Windows
   :synopsis: input, output, and conversion of Bruker XRF spectra files

.. moduleauthor:: Donald Windover (windover@nist.gov)

Created on Fri Apr 26 09:49:43 2019

Edited by Vincent Lan (vlan@umd.edu) June/July 2020

Changelog:
    General:
        -Fixed spelling errors in comments
        -Changed 'fittingdata' to 'fitting_data' for greater legibility
        -Converted duplicated code fragments into helper functions for modularity
        -Replaced np.zeros with np.empty
        -Minor aesthetic chances in exported .txt files
    FittingData:
    Functions:
        bruker_txt_test:
            -Changed .readlines() to .readline() since only checking first line
        bruker_txt_import:
            -Changed for loop to use enumerate
            -Created a constant data_file_lines to use in array creation
        bruker_spx_import:
            -Changed 'level' to 'sublevel' to better show parent-child relationship
            -Rearranged loop structure and added provisions
                to account for data with channels != 4096
        write_converted_file:
            -Changes file endings
            -Writes to file from a list, adding a newline after each element

"""
#

############################
#  20190426 Donald Windover
#  20200706 Vincent Lan
#
#  These are are series of functions used to read in data from the Bruker M4
#  and to allow editing of the spectrum contained in the file.  They also
#  us to read in fitting information, and write out to a txt format similar
#  to the Bruker .txt conversion as well as the .msa file format
#
#  we use these functions to read in data for analysis in later packages
#
###########################
#
import xml.etree.ElementTree as ET
import re as re
from datetime import datetime
import numpy as np


class FittingData:
    """ all parameters imported or exported from Bruker Spectra Files

    def __init__(self, file_name):

        **self.file_name:** str [?]
            file to be read (always input for an instance)

    **calibration_abs:** float [-955.1]
        set value for zero of Bruker M4 at 40keV setting

    **calibration_lin:** float [10]
        set value for linearity of Bruker M4 spectr at 40keV setting

    **channels:** np.array [4096,]
        assuming 4096 MCA channels (Bruker m4)  EDAX uses 4000 instead

    **data_lines:** str [""]
        dummy string for the parsed lines of an ascii readable file

    **date_measure:** str [""]
        date read from Bruker M4 .spx files by ETREE

    **detector_thickness:** float [0]
        detector thickness in μm

    **detector_type:** str [""]
        detector type from Bruker M4 .spx

    **energy_scale:** np.array [4096,]
        numpy array of 4096 energy spectra scalings  (initialized as empty)

    **file_content:** str [""]
        entire ascii file read for import of data

    **file_line:** str ['']
        the current single line of a data file

    **file_lines:** list of str ['']
        python list of strings; each list being a line of an imported file

    **file_status:** bool [False]
        this boolean is used in test of .txt file = Bruker spectra

    **header_lines:** list of str ['']
        lines preceding the start of spectra data in import file

    **life_time_in_ms:** int [0]
        live time for MCA spectra collection

    **line_count:** int [0]
        dummy variable used in counting for start of spectra data

    **mn_fwhm:** float [143.796]
        manganese fwhm in eV for the detector used in spectra collection

    **modification:** str ['_modified.txt']
        name to be appended to .txt files after modified by functions

    **no_channels:** int [4096]
        number of channels in MCA default to Bruker value

    **pulse_density:** str ['']
        left as a string, as we are not using this for any calc

    **real_time_in_ms:** int [0]
        real time for MCA spectra collection

    **replace_lines:** list of str []
        list of string lines containing modified spectra data

    **shaping_time:** float [0]
        detector count rate shaping time

    **si_dead_layer:** float [0]
        dead layer of detector used in modeling quantitative data

    **start_count:** int [21]
        line directly after header info in spectra file (start of spectra)

    **time_measure:** str ['']
        time read from Bruker M4 .spx files by ETREE

    **window_type:** str ['']
        material used in detector (used in quant modeling)

    """

    def __init__(self, file_name):
        """ name of spectrum file imported or exported"""
        self.file_name = file_name

    calibration_abs = -955.1
    calibration_lin = 10
    channels = np.empty(4096)
    data_lines = ''
    date_measure = ''
    detector_thickness = 0
    detector_type = ''  #
    energy_scale = np.empty(4096)
    file_content = ''
    file_line = ''
    file_lines = ''
    file_status = False
    header_lines = ''
    modification = '_modified.txt'
    life_time_in_ms = 0
    line_count = 0
    mn_fwhm = 143.796
    no_channels = 0
    pulse_density = ''
    real_time_in_ms = 0
    replace_lines = []
    shaping_time = 0
    si_dead_layer = ''
    start_count = 21
    time_measure = ''
    window_type = ''


###########################
#  20190426 Donald Windover
#  This function tests if a .txt file present in the directory is a Bruker
#  spectrum txt file by reading the first line and comparing to expected value
#
def bruker_txt_test(fitting_data):
    """ Function to test if *.txt* is a Bruker spectra file

    We test if it is a Bruker *.txt* file by reading the first line and
    comparing the result against the known header line from the Bruker txt
    files (see the following line of code:)

    Parameters
    ----------

    fitting_data : see Class FittingData


    Example
    -------

    >>>> fitting_data.file_status = bool(r'Bruker Nano GmbH Berlin, Germany\\n'
    >>>>                                in fitting_data.file_lines[0]) 

    See Also
    --------

    FittingData

    """
    #    .. code-block:: python
    #
    #        fitting_data.file_status = bool(r'Bruker Nano GmbH Berlin, Germany\\n'
    #                                       in fitting_data.file_lines[0]) """Build and fit a model of an EDS Signal1D.
    #
    fitting_data.file_content = open(fitting_data.file_name)
    fitting_data.file_line = fitting_data.file_content.readline()
    fitting_data.file_content.close()
    fitting_data.file_status = bool('Bruker Nano GmbH Berlin, Germany\n'
                                    == fitting_data.file_line)
    reset_to_default_values(fitting_data)
    return


###########################
#  20190426 Donald Windover
#  20200706 Vincent Lan
#  This function reads in the .txt file to provide data for fitting routines
#
def bruker_txt_import(fitting_data):
    """ Function to import the Bruker *.txt* spectra data

    The working portion of this import breaks the 2 column (energy,counts)
    spectral data into two numpy arrays via an inefficient *for* loop:

    .. code-block:: python

        fitting_data.energy_scale[i] = np.float(split_line[0])
        fitting_data.channels[i] = np.float(split_line[1])

    """
    # open Bruker .txt file and read lines
    txt_read_lines(fitting_data)
    # determine the start position of the energy, counts data
    txt_start_count(fitting_data)
    # extract the data beginning at the start count
    fitting_data.data_lines = fitting_data.file_lines[fitting_data.start_count:]

    # a constant for use in future arrays and loops
    data_lines_range = len(fitting_data.data_lines)

    # keeps only the lines of spectral data
    fitting_data.energy_scale = np.empty(data_lines_range)
    fitting_data.channels = np.empty(data_lines_range)

    # split the data into energy and channels

    for index, file_line in enumerate(fitting_data.data_lines):
        split_line = file_line.split()
        fitting_data.energy_scale[index] = np.float(split_line[0])
        fitting_data.channels[index] = np.float(split_line[1])

    # for index, file_line in enumerate(fitting_data.data_lines):
    #     split_line = [np.float(item) for item in file_line.split()]
    #     fitting_data.energy_scale[index] = split_line[0]
    #     fitting_data.channels[index] = split_line[1]

    # provides two 1D arrays with the energy and counts data
    print('import size: ', fitting_data.channels.shape)
    reset_to_default_values(fitting_data)
    return  # these counts have been pulse pile up modified


###########################
#  20190426 Donald Windover
#  This function reads in the .txt file, passes the channels and energy
#  for modification, and then resaves channels with a "modified" name change
#
def bruker_txt_mod(fitting_data):
    """ function to export a Bruker *.txt* with modified spectra data

    This function opens the Bruker *.txt* file, separates the header data
    from the spectra data section, and rewrites the content of the spectra
    section from *energy_scale* and *channels* currently in the fitting_data
    instance.  The data is then recombined and saved under a new name
    combined from the original plus the *modification* string.

    """

    print('size into string on export: ', fitting_data.channels.shape)
    txt_read_lines(fitting_data)
    txt_start_count(fitting_data)

    # extracts lines with information and lines with data
    fitting_data.header_lines = fitting_data.file_lines[:fitting_data.start_count]
    fitting_data.data_lines = fitting_data.file_lines[fitting_data.start_count:]

    fitting_data.replace_lines = []

    # keeps only the lines of spectral data
    for data_line in fitting_data.data_lines:
        split_line = data_line.split()
        replace_line = '    '.join(split_line) + '\n'
        fitting_data.replace_lines.append(replace_line)

    text_list = fitting_data.header_lines + fitting_data.replace_lines

    write_converted_file(fitting_data, text_list, 'modification', False)
    reset_to_default_values(fitting_data)
    return


###########################
#  20190426 Donald Windover
#  This function reads in the *.txt file, passes the channels and energy
#  for modification
#

def bruker_msa_import(fitting_data):
    """function to open Bruker MSA format spectra files"""
    print(fitting_data.file_name)
    fitting_data.file_content = open(fitting_data.file_name)
    fitting_data.file_lines = fitting_data.file_content.readlines()
    fitting_data.file_content.close()
    fitting_data.line_count = 0
    for fitting_data.file_line in fitting_data.file_lines:
        # finds the start of the error data
        fitting_data.line_count = fitting_data.line_count + 1
        # looks for the word 'Spectrum' in each line
        if fitting_data.file_line.find('XPERCHAN') != -1:
            splitline = fitting_data.file_line.split(':')
            fitting_data.calibration_lin = 1000 * float(splitline[1])
        if fitting_data.file_line.find('OFFSET') != -1:
            splitline = fitting_data.file_line.split(':')
            fitting_data.calibration_abs = -10 * float(splitline[1])
        if fitting_data.file_line.find('SPECTRUM') != -1:
            # startMSA local variable only for indexing end of MSA header
            start_msa = fitting_data.line_count
    fitting_data.data_lines = fitting_data.file_lines[start_msa:-1]
    # keeps only the lines of error data
    string = ''.join(fitting_data.data_lines)
    new_string = re.sub("\n", '', string)
    fitting_data.channels = np.fromstring(new_string, sep=',')
    fitting_data.energy_scale = (fitting_data.calibration_abs +
                                 np.arange(4096) * fitting_data.calibration_lin)
    reset_to_default_values(fitting_data)
    return


###########################
#  20190426 Donald Windover
#  This function reads in the *.SPX file, passes the channels and energy
#  for modification
#
def bruker_spx_import(fitting_data):
    """function to import channels and energy info from Bruker *.spx* file"""
    #
    # establish 4096 array to take the channel data from an spx file
    # for the Bruker spx data (assumes all 4096 channels present)
    #
    # counts_spx = np.empty((1,4096))
    # prints which file is being converted
    #
    try:
        # opens the XML file
        tree = ET.parse(open(fitting_data.file_name, "r"))
        root = tree.getroot()
        # print(r'SPXFile: ', fitting_data.file_name)
    except TypeError:
        # fails gracefully, if filename or format is not XML.
        print("Unable to open and parse input definition file: "
              + fitting_data.FileName)
    # pulls in the channel data
    for sublevel_two in root:
        # pulls in the parameters needed for the txt file
        for sublevel_three in sublevel_two:
            if sublevel_three.tag == 'TRTHeaderedClass':
                for sublevel_four in sublevel_three:
                    # pulls in the collection time information
                    if {'Type': 'TRTSpectrumHardwareHeader'} == sublevel_four.attrib:
                        for sublevel_five in sublevel_four:
                            if sublevel_five.tag == 'RealTime':
                                fitting_data.real_time_in_ms = np.float(
                                    sublevel_five.text)
                            if sublevel_five.tag == 'LifeTime':
                                fitting_data.life_time_in_ms = np.float(
                                    sublevel_five.text)
                            if sublevel_five.tag == 'PulseDensity':
                                fitting_data.pulse_density = sublevel_five.text
                            if sublevel_five.tag == 'ShapingTime':
                                fitting_data.shaping_time = np.float(
                                    sublevel_five.text)
                    if {'Type': 'TRTDetectorHeader'} == sublevel_four.attrib:
                        for sublevel_five in sublevel_four:
                            # pulls is in the detector info
                            if sublevel_five.tag == 'Type':
                                fitting_data.detector_type = sublevel_five.text
                            if sublevel_five.tag == 'DetectorThickness':
                                fitting_data.detector_thickness = float(sublevel_five.text)
                            if sublevel_five.tag == 'SiDeadLayerThickness':
                                fitting_data.si_dead_layer = sublevel_five.text
                            if sublevel_five.tag == 'WindowType':
                                fitting_data.window_type = sublevel_five.text
            if {'Type': 'TRTSpectrumHeader'} == sublevel_three.attrib:
                for sublevel_four in sublevel_three:
                    # pulls in the energy calibration info
                    if sublevel_four.tag == 'Date':
                        date = sublevel_four.text
                    if sublevel_four.tag == 'Time':
                        time = sublevel_four.text
                    if sublevel_four.tag == 'ChannelCount':
                        fitting_data.no_channels = sublevel_four.text
                        # formats channels and energy arrays to match
                        # the number of channels
                        fitting_data.channels = np.empty(int(sublevel_four.text))
                        fitting_data.energy_scale = np.empty(int(sublevel_four.text))
                        if sublevel_four.text != "4096":
                            print("NOTE: Number of channels is " + sublevel_four.text
                                  + ", instead of the default 4096.")
                    if sublevel_four.tag == 'CalibAbs':
                        calibration_abs = np.float(sublevel_four.text)
                    if sublevel_four.tag == 'CalibLin':
                        calibration_lin = np.float(sublevel_four.text)
                    if sublevel_four.tag == 'SigmaAbs':
                        sigma_abs = np.float(sublevel_four.text)
                    if sublevel_four.tag == 'SigmaLin':
                        sigma_lin = np.float(sublevel_four.text)
        if sublevel_two.find('Channels') is not None:
            channels = sublevel_two.find('Channels')
            fitting_data.channels = np.asarray(channels.text.split(','), dtype=int)
            # print('import size: ', fitting_data.channels.shape)
    # converts the time to the correct format
    time = datetime.strptime(time, "%H:%M:%S")
    fitting_data.time_measure = time.strftime("%I:%M:%S %p")
    # converts the date to the correct format
    date = datetime.strptime(date, "%d.%m.%Y")
    fitting_data.date_measure = date.strftime("%m/%d/%Y")
    # rescales energy calibration factors for the txt format
    fitting_data.calibration_abs = 1000 * calibration_abs
    fitting_data.calibration_lin = 1000 * calibration_lin
    # Energy used in the calucation of Mn FWHM (approximated on 2017/10/19)
    mn_energy = 5.900
    # Formula given by Bruker (Falk Reinhardt) on 2017/10/19
    sigma = np.sqrt(sigma_abs + mn_energy * sigma_lin)
    fwhm_factor = 1000 * np.sqrt(8 * np.log(2)) * sigma
    fitting_data.mn_fwhm = float(fwhm_factor)  # we now know the calc rather than needing a const.
    # print(fitting_data.mn_fwhm)
    # Energy scale calculation
    fitting_data.energy_scale = np.empty(int(fitting_data.no_channels))
    for i in np.arange(int(fitting_data.no_channels)):
        fitting_data.energy_scale[i] = (fitting_data.calibration_abs +
                                        fitting_data.calibration_lin * i) / 1000
    return  # provides the comma delimited list of channel intensity


###########################
#  20190426 Donald Windover
#  20200706 Vincent Lan
#  This function reads in the *.spx file, passes the channels and energy
#  for modification, and generates a .txt file in the Bruker output format
#
def bruker_spx_to_txt_convert(fitting_data):
    """function converting Bruker *.spx* to *.txt* """
    # prints which file is being converted
    print(fitting_data.file_name)
    #
    bruker_spx_import(fitting_data)
    # text file data formatting
    text_header = []
    text_header.append(r'Bruker Nano GmbH Berlin, Germany')
    text_header.append(r'esprit 1.9')
    text_header.append(r'')
    text_header.append(r'Date: ' + fitting_data.date_measure + ' '
                       + fitting_data.time_measure)
    text_header.append(r'Real time: ' + '%.0f' % fitting_data.real_time_in_ms)
    text_header.append(r'Life time: ' + '%.0f' % fitting_data.life_time_in_ms)
    text_header.append(r'Pulse density: ' + fitting_data.pulse_density)
    text_header.append(r'')
    text_header.append(r'')
    text_header.append(r'Detector type: ' + fitting_data.detector_type)
    text_header.append(r'Window type: ' + fitting_data.window_type)
    text_header.append(r'Detector thickness: ' + str(fitting_data.detector_thickness))
    text_header.append(r'Si dead layer: ' + fitting_data.si_dead_layer)
    text_header.append(r'')
    text_header.append(r'Calibration, lin.: ' + str(fitting_data.calibration_lin))
    text_header.append(r'Calibration, abs.: ' + '%.3f' % fitting_data.calibration_abs)
    text_header.append(r'Mn FWHM: ' + '%.3f' % fitting_data.mn_fwhm)
    text_header.append(r'Fano factor: 0.116')
    text_header.append(r'Channels: ' + fitting_data.no_channels)
    text_header.append(r'')
    text_header.append(r'Energy     Counts')
    # including energy and counts
    for index in np.arange(int(fitting_data.no_channels)):
        text_header.append('%.4f' % fitting_data.energy_scale[index] +
                           '    ' + '%.0f' % fitting_data.channels[index])
    write_converted_file(fitting_data, text_header, 'spx_to_txt')
    return


def test_io():
    """Function tests the spx, msa, and txt readers using sample files

    Parameters
    ----------

    txt_file : name of test Bruker XRF txt file
    msa_file : name of test Bruker XRF msa file
    spx_file : name of txt Bruker XRF spx file
    txt, msa, spx : instances of Class FittingData


    Example
    -------

    >>>> txt_file = r'test.txt'
    >>>> txt = FittingData(txt_file)
    >>>> bruker_txt_test(txt)
    >>>> bruker_txt_import(txt)
    >>>> bruker_txt_mod(txt)
    >>>> msa_file = r'test.msa'
    >>>> msa = FittingData(msa_file)
    >>>> bruker_msa_import(msa)
    >>>> spx_file = r'test2.spx'
    >>>> spx = FittingData(spx_file)
    >>>> bruker_spx_import(spx)
    >>>> bruker_spx_to_txt_convert(spx)
    """
    txt_file = r'test.txt'
    txt = FittingData(txt_file)
    bruker_txt_test(txt)
    bruker_txt_import(txt)
    bruker_txt_mod(txt)
    msa_file = r'test.msa'
    msa = FittingData(msa_file)
    bruker_msa_import(msa)
    spx_file = r'test2.spx'
    spx = FittingData(spx_file)
    bruker_spx_import(spx)
    bruker_spx_to_txt_convert(spx)
    return

#
# HELPER FUNCTIONS
#


# helper function that opens a .txt and reads the lines
def txt_read_lines(fitting_data):
    fitting_data.file_content = open(fitting_data.file_name)
    fitting_data.file_lines = fitting_data.file_content.readlines()
    fitting_data.file_content.close()


# helper function to determine the start line of the energy, counts data
def txt_start_count(fitting_data):
    # checks whether default start_count line is where
    # the energy, counts data begins
    # bypasses start count verification if default values hold true
    if fitting_data.file_lines[fitting_data.start_count].find('Counts') != -1:
        pass
    # else count the line where energy, counts data begins
    else:
        fitting_data.line_count = 0
        # finds the start of the error data
        for fitting_data.file_line in fitting_data.file_lines:
            fitting_data.line_count += 1
            # looks for the word 'Counts' in each line
            if fitting_data.file_line.find('Counts') != -1:
                start_count = fitting_data.line_count
                # print a warning if start count is different than default
                if fitting_data.start_count != start_count:
                    print('warning: start of channels != normal value')
                    fitting_data.start_count = start_count


# helper function to clear values
def reset_to_default_values(fitting_data):
    fitting_data.file_content = ''
    fitting_data.file_line = ''
    fitting_data.file_lines = ''
    fitting_data.header_lines = ''
    fitting_data.data_lines = ''
    fitting_data.replace_lines = []


# helper function to write lines from the text source into
# a new file with a new file ending
def write_converted_file(fitting_data, text_source, operation, separator='\n'):
    # determines whether to use a newline separator
    if separator != '\n':
        separator = ''
    # determines the type of export operation
    if operation == 'spx_to_txt':
        file_name_mod = fitting_data.file_name.replace('.spx', '.txt')
    elif operation == 'modification':
        file_name_mod = fitting_data.file_name.replace('.txt', fitting_data.modification)
    # writing the file
    file = open(file_name_mod, "w")
    # writes text_source to a txt file, adding a newline after each line
    file.writelines(["%s" % line + separator for line in text_source])
    file.close()
