# Observation Main Function
from astropy import units as u
from astropy.time import Time as astro_time # use for JD
from time import sleep, time # use sleep to pause code
import csv   # for writring scope settings to csv file
import ugradio as ugr
from astropy.table import Table
from astropy.io import ascii

def observe(target_RA = None,
            target_DEC = None,
            total_observation_duration = 3600.,
            sample_collect_interval = 1,
            write_data_interval = 300.,
            scope_repoint_interval = 10,
            verbose = False,
            ugr_verbose = False,
            sun = False,
            moon = False):

    """
    Main data collection function.
    
    Parameters
    ----------
    target_RA                  (float)    : Right ascension in degrees for target (J2000). 
    target_DEC                 (float)    : Declination in degrees for target (J2000).
    total_observation_duration (float)    : Total time in seconds.
    sample_collect_interval    (float)    : Time in seconds between voltage readings.
    write_data_interval        (float)    : Time in seconds between writing data to disk.
    scope_repoint_interval     (float)    : Time in seconds between repointing the scope.
    verbose                    (bool)     : Be verbose, for our functions only.
    ugr_verbose                (bool)     : Be verbose, for ugr functions only.
    sun                        (bool)     : Target is the sun.  Overrides target_RA and target_DEC.
    moon                       (bool)     : Target is the moon. Overrides target_RA and target_DEC, and sun.

    Return
    ------
    Returns 0 if the az or alt is out of range.
    """

    # verbosity function
    verboseprint = print if verbose else lambda *a, **k: None

    MIN_IFM_ALT = 5
    MIN_IFM_AZ = 95
    MAX_IFM_AZ = 295
        
    if total_observation_duration % write_data_interval == 0:

        verboseprint('Observing for           %d seconds.' %total_observation_duration)
        verboseprint('Repointing IFM every    %d seconds.' %scope_repoint_interval)
        verboseprint('HPM samples every       %d seconds.' %sample_collect_interval)
        verboseprint('Writing volt data every %d seconds.' %write_data_interval)


        # ######
        #
        # Before loop set the IFM to the target coordinates and start recording data
        #
        # #####

        # get the current location
        nch_longitude = ugr.nch.lon
        nch_latittude = ugr.nch.lat
        nch_altitude  = ugr.nch.alt

        # get the current julian date
        observation_date = astro_time.now().jd

        #set the equinox
        equinox = 'J2000'
        
        # override RA DEC if target is sun or moon
        if sun == True:
            target_RA, target_DEC = ugr.coord.sunpos(observation_date)
            equinox = 'J2018'
            
        if moon == True:
            target_RA, target_DEC = ugr.coord.moonpos(observation_date)
            equinox = 'J2018'


        # convert target coordinates to AltAz using observetory location and observation_date
        target_altitude, target_azimuth = ugr.coord.get_altaz(target_RA,
                                                                  target_DEC,
                                                                  jd = observation_date,
                                                                  equinox = equinox)


        if target_altitude < MIN_IFM_ALT:
            verboseprint('target_altitude = %f' %target_altitude)
            print('Target altitude out of range. Ending observation before it began.' %MIN_IFM_ALT)
            return 0

        if target_azimuth <= MIN_IFM_AZ:

            target_azimuth  = target_azimuth + 180
            target_altitude = 180 - target_altitude

            if target_azimuth < MIN_IFM_AZ or target_azimuth > MAX_IFM_AZ:
                verboseprint('target_azimuth = %f' %target_azimuth)
                print('Target azimuth out of range. Ending observation before it began.')
                return 0

        if target_azimuth > MAX_IFM_AZ:

            target_azimuth  = target_azimuth - 180
            target_altitude = 180 - target_altitude

            if target_azimuth < MIN_IFM_AZ or target_azimuth > MAX_IFM_AZ:
                verboseprint('target_azimuth = %f' %target_azimuth)
                print('Target azimuth out of range. Ending observation before it began.')
                return 0

        # intialize interferometer
        ifm = ugr.interf.Interferometer()
        
        # point the IFM to the target coordinates
        ifm.point(target_altitude, target_azimuth, verbose = verbose)

        # initialize multimeter
        hpm = ugr.hp_multi.HP_Multimeter()

        # tell multimeter to being recording at desired rate
        hpm.start_recording(sample_collect_interval)

        # ######
        #
        # Collect data
        #
        # #####

        # compute number of reocording intervals
        number_of_collection_intervals =  int(total_observation_duration / write_data_interval)

        print('Total number of collection intervals = %d.' %number_of_collection_intervals)

        # repeat data collection number_of_collection_intervals times.
        for collection_interval in range(number_of_collection_intervals):

            print('Collection interval: %d.\n' %collection_interval)

            # julian day at start of collection_interval
            # used for a few things in the collection_interval loop
            observation_date = astro_time.now().jd

            # filename for scope settings file.
            # Include julian day of start of collection_interval
            folder_scope_settings = 'Scope_Settings/observation_'
            folder_scope_settings = 'Scope_Settings/sun_observation_' if sun == True else folder_scope_settings
            folder_scope_settings = 'Scope_Settings/moon_observation_' if moon == True else folder_scope_settings
            filename_scope_settings = folder_scope_settings + str(observation_date).replace(".", "_") + '_scope_settings.csv'

            # #####
            #
            # While loop will run for write_data_interval seconds.
            #
            # #####

            # capture current time in unix seconds.
            # this time will be updated each run through the while loop.
            # The end condition for the while loop will be based on this.
            now_unix = astro_time.now().unix
            # print('now_unix = %f' %now_unix)

            # capture current time in unix seconds.
            # this remains fixed.
            # The end condition for the while loop will be based on this.
            collection_interval_start_time = astro_time.now().unix

            # For the end condition for the while loop
            collection_interval_end_time = collection_interval_start_time + write_data_interval
            # print('collection_interval_end_time = %f' %collection_interval_end_time)


            # Run continuously until write_data_interval have passed
            while(now_unix < collection_interval_end_time):
                

                # #####
                #
                # Point the IFM
                #
                # #####

                # Capture the julian day at the start of the loop
                current_JD = astro_time.now().jd

                # set the equinox
                equinox = 'J2000'

                # override RA DEC if target is sun or moon
                if sun == True:
                    target_RA, target_DEC = ugr.coord.sunpos(current_JD)
                    equinox = 'J2018'

                if moon == True:
                    target_RA, target_DEC = ugr.coord.moonpos(current_JD)
                    equinox = 'J2018'

        
                # Recompute the target coordinates based on the new julian day
                # convert target coordinates to AltAz using observetory location and observation_date
                target_altitude, target_azimuth = ugr.coord.get_altaz(target_RA,
                                                                      target_DEC,
                                                                       jd = current_JD,
                                                                       equinox = equinox)

                # target_altitude, target_azimuth = transfrom_RA_DEC_to_ALT_AZ(ra = target_RA, dec = target_DEC, jd = observation_date)
                # target_altitude = target_altitude.value
                # target_azimuth = target_azimuth.value

                if target_altitude < MIN_IFM_ALT:

                    print('target_azimuth = %f' %target_azimuth)
                    print('Target altitude out of range. Ending observation and stowing IFM. The most recent voltage collection has been DISCARDED.')
                    hpm.end_recording()
                    ifm.stow(verbose = verbose)
                    return 0
                if target_azimuth <= MIN_IFM_AZ:

                    target_azimuth = target_azimuth + 180
                    target_altitude = 180 - target_altitude

                    if target_azimuth < MIN_IFM_AZ or target_azimuth > MAX_IFM_AZ:
                        print('Target azimuth out of range. Ending observation and stowing IFM. The most recent voltage collection has been DISCARDED.')
                        hpm.end_recording()
                        ifm.stow(verbose = verbose)

                        return 0
                if target_azimuth > MAX_IFM_AZ:

                    target_azimuth  = target_azimuth - 180
                    target_altitude = 180 - target_altitude

                    if target_azimuth < MIN_IFM_AZ or target_azimuth > MAX_IFM_AZ:
                        print('Target azimuth out of range. Ending observation and stowing IFM. The most recent voltage collection has been DISCARDED.')
                        hpm.end_recording()
                        ifm.stow(verbose = verbose)
                        return 0
        
                # point IFM to target coordinates
                # This takes a nonzero amount of time
                # Do not repoint more than once every 5 seconds
                ifm.point(target_altitude, target_azimuth, verbose = verbose)

                # #####
                #
                # Write the scope settings to the scope settings file
                #
                # #####
                
                # Get the positions of the scopes from the IFM
                ifm_actual_settings = ifm.get_pointing(verbose = verbose)
                
                # Package the positions with the julian day
                line_for_scope_settings_file = [astro_time.now().jd,
                                                ifm_actual_settings['ant_e'][0], ifm_actual_settings['ant_e'][1],
                                                ifm_actual_settings['ant_w'][0], ifm_actual_settings['ant_w'][1]]
                
                # write the positions to disk
                with open(filename_scope_settings, 'a') as f:
                    
                    writer = csv.writer(f)
                    writer.writerow(line_for_scope_settings_file)
                    
                # sleep for scope_repoint_interval seconds
                # This controls the rate at which we repoint the scope
                # Do not set this below 5 seconds!
                sleep(scope_repoint_interval)
                    
                # Capture the current time in unix seconds.
                # If this is write_data_interval more than the start time, the loop will end
                now_unix = astro_time.now().unix

            # #####
            #
            # At the end of the write_data_interval second loop, collect the volt data
            #
            # #####

            # Ask the HPM for the volt data
            volts, times = hpm.get_recording_data()
            print('%d volts captured.' %len(volts))

            # package the volt data in a table
            table_to_write = Table(data = [volts, times],
                                   names = ('volts', 'times'),
                                   meta = {'JD'   : observation_date})

            # filename for scope settings file.
            # Include julian day of start of collection_interval
            folder_volt_table   = 'HPM_Data/observation_'
            folder_volt_table   = 'HPM_Data/sun_observation_' if sun == True else folder_volt_table
            folder_volt_table   = 'HPM_Data/moon_observation_' if moon == True else folder_volt_table          
            filename_volt_table = folder_volt_table + str(observation_date).replace(".", "_") + '_volt_table.csv'

            # write the tabke to disk
            table_to_write.write(filename_volt_table, format='ascii.csv')  

           
    else:
        print('total_observation_duration % write_data_interval != 0')

    print('Observation & Recording Ending.')
    hpm.end_recording()
    ifm.stow(verbose = verbose)


