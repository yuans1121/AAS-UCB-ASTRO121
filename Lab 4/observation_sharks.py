from astropy.table import QTable, Table
from astropy.io import ascii
from astropy import units as u
import numpy as np
import glob
from astropy.time import Time
from astropy.coordinates import SkyCoord, Galactic, AltAz, EarthLocation
import os
import errno
import ugradrio as ugr
import Leuschner as leu
import time

# 'Constants'

leus = EarthLocation(lat = 37.9183 * u.deg, lon = -122.1067 * u.deg, height = 304.0 * u.m)

MIN_TEL_ALT = 15 * u.deg
MAX_TEL_ALT = 85 * u.deg

MIN_TEL_AZ = -5 * u.deg
MAX_TEL_AZ = 365 * u.deg

INTIAL_POINTING_FILENAME = 'target_tables/target_table_000.csv' 


def load_target_points_table(path = None):
    """Loads a table of target points. Columns ra, dec, and collected (true/false).
    
    The file with the greatest number will be loaded. Files must be in a directory 'target_tables/'. 
    
    Path can be set if you have some particular table you want to load.
    """

    if path == None:

        #glob returns a list of filenames that match the format given
        # max returns the filename that is the 'largest'.
        # assumes all filenames are EXACTLY the same except for the last few characters
        path = max(glob.glob('target_tables/target_table_*'))
        

    # the column names for the table we will make from the file
    names = ('gal_long_ell',
             'gal_lat_bee',
             'ra',
             'dec',
             'collected',
             'pointing_number')


    # the table has units where approprtaite (see save_table)
    return QTable(ascii.read(path), names = names)


def save_table(table):
    """
    Saves a table of target points.
    The filename is automatically created to be
    one number higher than the previous highest number in the directory 'target_tables/'. 
    
    The file target_table_test_000.csv MUST be manually created and present in the directory.
    """

    # for verbosity

    # file the highest number filename and increment that number by 1
    current_max_filename = max(glob.glob('target_tables/target_table_*'))
    current_max_filename = current_max_filename.replace("_", " ")
    current_max_filename = current_max_filename.replace(".", " ")
    current_max_number = [int(s) for s in current_max_filename.split() if s.isdigit()][0]
    new_number_string = "{0:0=3d}".format(current_max_number + 1)
    
    new_filename = 'target_tables/target_table_' + new_number_string + '.csv'

    table.write(new_filename, format='ascii.ecsv')


def can_target(target_coords):
    
    alt_OK = False
    
    if target_coords.alt.value >= MIN_TEL_ALT.value:
        if target_coords.alt.value <= MAX_TEL_ALT.value:
            alt_OK == True
            
    az_OK = False
    
    if target_coords.az.value >= MIN_TEL_AZ.value:
        if target_coords.az.value <= MAX_TEL_AZ.value:
            az_OK == True
            
    return alt_OK and az_OK


def target_alt_az(target):

	galactic = SkyCoord(target['gal_long_ell'],target['gal_lat_bee'], frame='galactic', unit = "deg")
	return galactic.transform_to(AltAz(obstime=Time.now(),location=leus))



def julian_date_str():    ### ADDED: allows calling this to obtain the current JD as a string with '_' instead of '.'
    
    current_time = Time.now()
    julian_date = str(current_time.jd)

    return julian_date.replace('.', '_')


def observe(lo_main_value = 635, 
            lo_alternate_value = 639, 
            forwards = True,
            num_spectra = 20, 
            num_noise_spectra = 5):

    ### Setup

    telescope = ugr.leusch.LeuschTelescope()
    spec = leu.Spectrometer('10.0.1.2')
    LO = ugr.agilent.SynthDirect() 
    noise_maker = ugr.leuschLeuschNoise()

    lo_main_value_interger_part_str = str(int(lo_main_value))

    lo_main_value_fractional_part_str = str(int((lo_main_value % 1) * 1000))

    lo_alternate_value_interger_part_str = str(int(lo_alternate_value))

    lo_alternate_value_fractional_part_str = str(int((lo_alternate_value % 1) * 1000))

    number_collected_this_run = 0

	# what are all the targets we need
    target_points = load_target_points_table()


    ### Pointing and collection

    # what direction are we scanning through points in?
    # forwards = -10 to 130
    direction = 1 if forwards == True else -1

	# go through al the needed targets
	for i, target in enumerate(target_points[::direction]):

		print('Target %d' %target['pointing_number'])

		if target['collected'] == True:
			print('Target %d collected.' %target['pointing_number'])

		else: # still need to collect this target

			print('Checking availability.')

			target_coords = target_alt_az(target)

			# if we can collect this target
			if can_target(target_coords) == True:

                # setLO _Main 
                LO.set_frequency(lo_main_value)

                #LO-main, Noise-off
                noise_maker.off()
                directory = 'LO_MAIN_NOISE_OFF'
                file = 'LO_' + lo_main_value_interger_part_str + '_' + lo_main_value_fractional_part_str + '_NOISE_OFF_JD_' + julian_date_str() + '.fits'
                filename = directory + '/' + file
                telescope.point(target_coords.alt, target_coords.az) # move the scope
                spec.read_spec(filename, num_spectra, (target['gal_long_ell'].value, target['gal_lat_bee'].value), 'ga')# get data

                #LO-main, Noise-on
                noise_maker.on()
                directory = 'LO_MAIN_NOISE_ONN'
                file = 'LO_' + lo_main_value_interger_part_str + '_' + lo_main_value_fractional_part_str + '_NOISE_ONN_JD_' + julian_date_str() + '.fits'
                filename = directory + '/' + file
                spec.read_spec(filename, num_spectra, (target['gal_long_ell'].value, target['gal_lat_bee'].value), 'ga')

                #LO-alternate, Noise-on
                LO.set_frequency(lo_alternate_value)
                directory = 'LO_ALTERNATE_NOISE_ONN'
                file = 'LO_' + lo_alternate_value_interger_part_str + '_' + lo_alternate_value_fractional_part_str + '_NOISE_ONN_JD_' + julian_date_str() + '.fits'
                filename = directory + '/' + file
                spec.read_spec(filename, num_spectra, (target['gal_long_ell'].value, target['gal_lat_bee'].value), 'ga')
      
                #LO-alternate, Noise-off
                noise_maker.off()
                directory = 'LO_ALTERNATE_NOISE_OFF'
                file = 'LO_' + lo_alternate_value_interger_part_str + '_' + lo_alternate_value_fractional_part_str + '_NOISE_OFF_JD_' + julian_date_str() + '.fits'
                filename = directory + '/' + file
                telescope.point(target_coords.alt, target_coords.az)
                spec.read_spec(filename, num_spectra, (target['gal_long_ell'].value, target['gal_lat_bee'].value), 'ga')

                
                row_number = i if forwards == True else (len(target_points) - 1) - i
				target_points[row_number]['collected'] = True

				save_table(target_points)

                number_collected_this_run = number_collected_this_run + 1
                print('Target %d collected.' %target['pointing_number'])
			else:
				print('Target %d not availbale.' %target['pointing_number'])
				pass # couldnt get this target

    print('Number collected this run %d' %number_collected_this_run)