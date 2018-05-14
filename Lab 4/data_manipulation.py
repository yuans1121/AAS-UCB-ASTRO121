from lab_4_data_loading import fits_ell

from freqs import freqs
freqs = freqs().value

from numpy import empty, mean

def mean_spectra(f_file, spectra_type = 'auto0_real', combine = False):
    """The mean of the spectra for a fits file. If the file is a noise file, the first of the 5 spectra is discarded.
    
    Parameters
    ----------
    f_file (fits) : The fits file
    spectra_type (string) default = 'auto0_real' : Which of the 4 types of spectra
    combine (bool) default = False : if true, combine the two spectra.
    Returns
    -------
    spectra (ndarray)
    
    """
        
    num_spectra = f_file[0].header['NSPEC']
    num_channels = f_file[0].header['NCHAN']
    
    if num_spectra == 5:
        num_spectra = 4
        i = 1
    else:
        i = 0
    
    spectra_2D = empty((num_spectra, num_channels))
    if combine == False:   
        
        for j in range(i, num_spectra + i):
            spectra_2D[j - 1, :] =  f_file[j + 1].data[spectra_type]
        
        return mean(spectra_2D, axis = 0)
    
    if combine == True:

        spectra_2D_2 = empty((num_spectra, num_channels))
        
        for j in range(i, num_spectra + i):
            spectra_2D[j - 1, :] =  f_file[j + 1].data['auto0_real']
            spectra_2D_2[j - 1, :] =  f_file[j + 1].data['auto1_real']

        return mean(spectra_2D + spectra_2D_2, axis = 0)
        
        

def spectra_in_window(y, L = 1419.405751786, R = 1421.405751786 ):
    """Truncate the spectra data to be in a window determined by frequency values.
    
    Parameters
    ----------
    y (ndarray) : The spectra to truncate
    L (float) default = 1419.5 : Left bound on the frequency
    R (float) default = 1421.5 : Right bound on the frequency
    
    Returns
    -------
    spectra (ndarray)

    """
    x = freqs
    x0 = x[x < R]
    y0 = y[x < R]
    y0 = y0[x0 > L]
    x0 = x0[x0 > L]
    
    return y0

def freqs_in_window(L = 1419.405751786, R = 1421.405751786 ):
    """Truncate the frequency data to be in a window determined by frequency values. 
    Make sure the window matches the spectra window.
    
    Parameters
    ----------
    L (float) default = 1419.5 : Left bound on the frequency
    R (float) default = 1421.5 : Right bound on the frequency
    
    Returns
    -------
    freqs (ndarray)

    """
    x = freqs
    x0 = x[x < R]
    x0 = x0[x0 > L]
    
    return x0

freqs_window = freqs_in_window()

from numpy import median

def data_under(x, y, median_y = True):
    """Find all y data below the median (mean if median is false) of y and return the corresponding x and y data."""
    threshold = median(y) if median_y == True else mean(y)
    return x[y < threshold], y[y < threshold]

from numpy import load

def loadnpz(filename):
    """Load a .npz file."""
    a = load(filename)
    d = dict(zip(("data1{}".format(k) for k in a), (a[k] for k in a)))
    
    return d['data1arr_0']

def gain(ell, T_N = 70, combine = False, mean_g = True):
    """Compute the mean gain for a given ell value"""
    P_prime = spectra_in_window(mean_spectra(fits_ell(ell, noise_off = False), combine = combine))
    P       = spectra_in_window(mean_spectra(fits_ell(ell, noise_off = True ), combine = combine))
    
    if mean_g == True : return mean((P_prime - P) / T_N)
    if mean_g == False: return (P_prime - P) / T_N

from numpy import polyfit, polyval

def shift(spectra, median_y = True):
    """Subtract a linear fit from the spectra. The linear fit is found using only the data under the mean of the data."""
    x, y = data_under(freqs_window, spectra, median_y = median_y)
    p = polyfit(x, y, 2)
    fit = polyval(p, freqs_window)
    return spectra - fit

def intensity(ell, combine = False,  mean_g = True):
    """Convert from measured power to intensity. Uses the mean gain of the ell."""
    ms = spectra_in_window(mean_spectra(fits_ell(ell), combine = combine))
    g = gain(ell, combine = combine, mean_g = mean_g)
    
    return ms / g

def shifted_intensity(ell, combine= False, mean_g  = True, median_y = True):
    """For a given ell value, determine the intensity of the spectra and shift its base to near 0K."""
    return shift(intensity(ell, combine = combine, mean_g = mean_g), median_y = median_y)

from astropy.units import MHz, km, s
from astropy.constants import c

doppler_corrections = loadnpz('doppler_corrections.npz')

def corrected_velocity(ell):
    """Convert frequency to velocity and subtract off the doppler correction for a given ell."""
    ell_to_index = lambda x:  int((x + 8)/2)

    f0 = 1420.405751786 * MHz
    f = freqs_window * MHz

    deltaf = f - f0

    v = ((-c * deltaf / f0).to(km / s)).value


    return v + doppler_corrections[ell_to_index(ell)] / 1000

from numpy import array 

def perpendicular_vector(v):
    """Find the perindicular vector (2D)"""
    return array([v[1], -v[0]])

def mag(v):
    """Find the magnitude of a vector"""
    return sum(v**2)**0.5

def unit_vector(v):
    """Find the unit vector. Careful about zero vector.."""
    return v / mag(v)

from numpy import cos, sin, arcsin
from astropy.units import deg
def point(gamma, R_0 = -8.5, R = 15.):
    """
    Finds the endpoint of a line starting at (0,R_0) and ending on a circle of radius R.
    Gamma is the angle of the line of vertical.
    """
    alpha = lambda gamma: arcsin(sin(gamma * deg) * abs(R_0) / R)
    beta = lambda gamma: 180 * deg - alpha(gamma) - gamma * deg
    beta_prime = lambda gamma: 180 * deg - beta(gamma)
    
    return R * array([ sin(beta_prime(gamma)), cos(beta_prime(gamma))])

def radial_line(gamma, R = 15, sol = [0, -8.5]):
    """Returns [x[0], x[1]], [y[0], y[1]]. 0 is start 1 is end"""
    return array(list(zip(*(sol, point(gamma, R = R)))))

from numpy import c_, round, int32, linspace, abs, diff
def connect(ends):
    """
    In a 2D space of INTEGERS, this returns all the points between two endpoints.
    
    ends = np.array([ [x[0], y[0]]
                    , [x[1], y[1]] ] )
                    
    credit: Paul Panzer, https://stackoverflow.com/a/47705495
    """
    d0, d1 = abs(diff(ends, axis=0))[0]
    if d0 > d1: 
        return c_[linspace(ends[0, 0], ends[1, 0], d0+1, dtype = int32),
                     round(linspace(ends[0, 1], ends[1, 1], d0+1)).astype(int32)]
    else:
        return c_[round(linspace(ends[0, 0], ends[1, 0], d1+1)).astype(int32),
                     linspace(ends[0, 1], ends[1, 1], d1+1, dtype = int32)]
    
def close(num, arr):
    """Returns the value from the array that is closest to the given number"""
    return min(arr, key = lambda x: abs(x - num))

from numpy import where
def close_index(num, arr):
    """Returns the index in the array that corresponds to the value closest to the given number"""
    return where(arr == close(num, arr))[0][0]

from numpy import dot
import sys

def make_vel_grid( x_arr, y_arr,R = 30):
    
    V_0 = 220 * array([-1, 0])
    pos_0 = array([0, -8.5])
    
    N_rows = len(y_arr)
    N_cols = len(x_arr)
    
    vel_grid = empty((N_rows,N_cols))
    
    for i, row in enumerate(range(N_rows)[::-1]): 
        for j, col in enumerate(range(N_cols)): 
            
            # start at the highest row and work down
            # start at the leftmost column and work right
            # i goes forwards, row goes backwards.
            
            x = x_arr[col]
            y = y_arr[row]
            
            # position and velocity vectors for the point in question
            pos_1 = array([x,y])
            V_1 = 220 * unit_vector(perpendicular_vector(pos_1))
            
            # difference in position vectors
            r_01 = (pos_1 - pos_0)
            
            
            vel_grid[i][j] = dot((V_1 - V_0), unit_vector(r_01))
            
            sys.stdout.write('\r'+str('%d, %d' %(i, j) ))    
    return vel_grid



def find_zeros(spectra):
    left_zeros = []
    right_zeros = []
    for i, power in enumerate(spectra):
    
        if power == 0:
            if i < 4096:
                left_zeros.extend([i])
            else:
                right_zeros.extend([i])
                
    return left_zeros, right_zeros

def fix_zeros(spectra):
   
    left_zeros, right_zeros = find_zeros(spectra)
    right_of_left = left_zeros[-1] + 1
    spectra[:right_of_left] = spectra[right_of_left]
    
    left_of_right = right_zeros[0] - 1
    spectra[left_of_right + 1:] = spectra[left_of_right]
    
    return spectra

def x_y_indicies_along_line(ell, x, y, x_arr, y_arr, R = 30):
    """
    Get the x and y indices of a line in R2 integer space that starts at (0,-8.5)
    and ends a distance R away, where the line is at an angle ell from the vertical.
    
    Parameters
    ----------
    x, y = [x0, x1], [y0,y1] (0 is start, 1 is end)
    """
    
    sol_index = close_index(x[0], x_arr), close_index(y[0], y_arr)
    end_index = close_index(x[1], x_arr), close_index(y[1], y_arr)
    
    line_coordinates = connect(array([[sol_index[0], sol_index[1]], [ end_index[0], end_index[1]]])).T # in index space

    
    return line_coordinates[0], line_coordinates[1]

def expected_velocities(ell, x, y, x_arr, y_arr, R = 30, vel_grid = []):
    """
    Get the expected velocities along a line at an angle ell and length R
    """

    x_line_indices, y_line_indices = x_y_indicies_along_line(ell, x, y, x_arr, y_arr, R = R)
    
    red_line_vels = empty(len(y_line_indices))
    
    for i in range(len(y_line_indices)):
        red_line_vels[i]  = vel_grid[y_line_indices[i]][x_line_indices[i]]
        
    return red_line_vels

def expected_distances(num_ds, x, y):    
    return linspace(0,1,num_ds) * mag(array([abs(x[1] - x[0]), abs(y[1] - y[0])]))    


def map_intensities_to_distance(vels, distances, ell, combine = True, mean_g = True, median_y = True):
    
    corr_vells = corrected_velocity(ell)
    intensities = shifted_intensity(ell,combine = combine, mean_g = mean_g, median_y = median_y)
    
    intensity_along_ell = empty(len(distances))
    
    for i in range(len(distances)):
        intensity_along_ell[i] = intensities[close_index(vels[i], corr_vells)]    
    return intensity_along_ell