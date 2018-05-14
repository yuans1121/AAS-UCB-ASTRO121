import glob

def _paths():
    """
    internal.
    
    returns a list of lists of files in Data folders and their (hardcoded) subdirectories.
    """
    return [sorted(glob.glob('Data/LO_MAIN_NOISE_OFF/*')),
            sorted(glob.glob('Data/LO_MAIN_NOISE_ONN/*')), 
            sorted(glob.glob('Data/LO_ALTERNATE_NOISE_OFF/*')), 
            sorted(glob.glob('Data/LO_ALTERNATE_NOISE_ONN/*'))]

path_collection = _paths()

from numpy import empty, array
from astropy.io import fits

def _ell_order():
    """
    internal.
    
    returns a list of lists of the ell-order of the files.
    """
    ell_order = empty((len(path_collection), len(path_collection[0])))
    for i, paths in enumerate(path_collection):
        for j, path in enumerate(paths):
            with fits.open(path) as file:
                ell_order[i][j] =  int(file[0].header['L'])
                
    return ell_order

from numpy import argsort

ell_order = _ell_order()

def _ell_sorted_paths():
    """
    internal.
    
    returns a list of lists of the ell-sorted paths.
    
    """

    
    

    return [array(path_collection[0])[argsort(ell_order[0])],
            array(path_collection[1])[argsort(ell_order[1])],
            array(path_collection[2])[argsort(ell_order[2])],
            array(path_collection[3])[argsort(ell_order[3])]]

ell_sorted_paths = _ell_sorted_paths()

def fits_ell(ell, LO_main = True, noise_off = True):
    """Get the fits file for a specfic ell value, with certain parameters.
    
    Parameters
    ----------
    ell (int) : The galacitc ell value of the required fits file
    LO_main (bool) default = True : Main LO value or alternate?
    noise_off (bool) default = True : Noise off or on?
    
    Returns
    -------
    fits file
    
    """


    ell_to_index = lambda x:  int((x + 8)/2)

    if LO_main == True:
        if noise_off == True:
            path_i = 0
        if noise_off == False:
            path_i = 1
    else:
        if noise_off == True:
            path_i = 2
        if noise_off == False:
            path_i = 3
    
    return fits.open(ell_sorted_paths[path_i][ell_to_index(ell)])

