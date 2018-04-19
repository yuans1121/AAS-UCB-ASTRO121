
from astropy.io import ascii
from astropy.table import QTable, Table
from astropy import units as u
from astropy.table import vstack
from astropy.time import Time as astro_time
import glob
import os
from numpy.random import random
from time import sleep
import numpy as np

def create_FAKE_volt_data(num_files = 1, time_between_files = 5):
    
    """ # example usuage:
        create_FAKE_volt_data(num_files = 5)
        # creates fils in the folder HPM_Data/"""
    
    fake_volts = random(600)
    fake_times = np.linspace(0,5,600)
    
    for i in range(num_files):
        
        JD = astro_time.now().jd
        table_to_write = Table(data = [fake_volts, fake_times],
                               names = ('volts', 'times'),
                               meta = {'JD'   : astro_time.now().jd})
        
        folder_volt_table   = 'HPM_Data/FAKEDATA_observation_'          
        filename_volt_table = folder_volt_table + str(JD).replace(".", "_") + '_volt_table.csv'
        table_to_write.write(filename_volt_table, format='ascii.csv')
        
        fake_volts = np.append(fake_volts,fake_volts)
        fake_times = np.append(fake_times,fake_times)
        sleep(time_between_files)


def create_FAKE_scope_settings(num_files = 1, time_between_files = 5):
    """
    # example usuage:
    create_FAKE_scope_settings(num_files = 5)
    # creates files in the folder Scope_Settings/"""

    
    for i in range(num_files):
            
        east_alt = random(6)*90
        east_az = random(6)*360
        
        west_alt = random(6)*90
        west_az = random(6)*360
        
        JD = np.linspace(1,6,6)*astro_time.now().jd
    
        table_to_write = Table(data = [JD,
                                    east_alt, east_az,
                                    west_alt, west_az],
                               names = ('JD', 'east_alt', 'east_az', 'west_alt', 'west_az'))
        
        folder_scope_table   = 'Scope_Settings/FAKEDATA_observation_'          
        filename_scope_table = folder_scope_table + str(astro_time.now().jd).replace(".", "_") + '_scope_settings.csv'
        table_to_write.write(filename_scope_table, format='ascii.csv')
        
        sleep(time_between_files)


def load_collection_of_scope_settings(foldername = 'Scope_Settings/'):
    """All the scope settings in one folder will be loaded into a table
        example usage:
        scope_settings = load_collection_of_scope_settings()"""
    
    list_of_qtables = []
    list_of_files = glob.glob('Scope_Settings/*')
    
    for filepath in sorted(list_of_files):
        
        names = ['obs_JD', 'east_alt', 'east_az', 'west_alt', 'west_az']
        
        table = QTable(ascii.read(filepath), names = names)
        
        table['east_alt'].unit = u.deg
        table['east_az' ].unit = u.deg
        table['west_alt'].unit = u.deg
        table['west_az' ].unit = u.deg
        
        list_of_qtables.append(table)

    return  QTable(vstack([Table(qtable) for qtable in list_of_qtables]))

def load_specific_scope_settings(path):
    
    """One scope settings file will be loaded into a table
    example usage:
    PATH = 'Scope_Settings/FAKEDATA_observation_2458180_7314443593_scope_settings.csv'
    scope_settings = load_specific_scope_settings(PATH)"""
    
    
    names = ['obs_JD', 'east_alt', 'east_az', 'west_alt', 'west_az']
    
    table = QTable(ascii.read(path), names = names)
    
    table['east_alt'].unit = u.deg
    table['east_az' ].unit = u.deg
    table['west_alt'].unit = u.deg
    table['west_az' ].unit = u.deg
    
    return table


def load_specific_volt_data(path, LST = False):
    """One volt data file will be loaded into a table
    example usage:
    PATH = 'HPM_Data/FAKEDATA_observation_2458180_7314443593_volt_table.csv'
    volt_data = load_specific_volt_data(PATH)"""
    
    names = ['volts', 'times']
    
    table = QTable(ascii.read(path), names = names)
    
    table['volts'].unit = u.V
    table['times' ].unit = u.s if LST == False else u.rad
    
    return table




def load_most_recent_volt_data(foldername = 'HPM_Data/'):
    """Most recent volt data file will be loaded into a table
    example usage:
    volt_data = load_most_recent_volt_data()"""
    
    list_of_files = glob.glob(foldername + '*')
    path_to_latest_file = max(list_of_files, key=os.path.getctime)

    names = ['volts', 'times']
    table = QTable(ascii.read(path_to_latest_file), names = names)
    
    table['volts'].unit = u.V
    table['times' ].unit = u.s
    
    return table