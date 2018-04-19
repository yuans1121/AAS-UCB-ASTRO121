from astropy import units as u, constants as c
import numpy as np
from astropy.units.quantity import Quantity
import sys

def is_Quantity(x):
    assert (type(x) == Quantity), 'type != astropy.units.quantity.Quantity'

def is_angle(x):
    is_Quantity(x)
    assert (x.unit == u.deg or x.unit == u.rad), 'unit != degree or radian'

def can_multiply(A,x):
    assert (A.shape[1] == x.shape[0]), "Matrix shape mismatch."

def can_solve(A,x):
    can_multiply(A.T,x)
    
def constrain_angle_to_one_circle(angle):
    """Constrain an angle to be between 0 and 360 degrees (or equivalent radians).
    Ex: -10degrees -> 350degrees.
    
    
    Parameters
    ----------
    angle : (Quantity <float, degrees or radians>)
    
    Returns
    -------
    (Quantity <float, degrees or radians>)
    return has same units as input
    """
    
    is_angle(angle)
    
    circle = 360 * u.deg if angle.unit == u.deg else 2 * np.pi * u.rad
    
    constrain = lambda x: x % circle
    
    if angle.isscalar:
        return constrain(angle)
                
    elif not angle.isscalar:
        for i in range(len(angle)):
            angle[i] = constrain(angle[i])
        return angle  

def hour_angle(LST, alpha):
    """Get the local hour angle of a target.
    
    Parameters
    ----------
    LST : (Quantity <float, degrees or radians>)
    
    Local Sideral Time
    
    alpha : (Quantity <float, degrees or radians>)
    Right Ascnesion of the target
    
    Note: Parameters can be arrays of the same size.
    
    Returns
    -------
    
    array of Quantity <float, degrees>
    Hour angle
    
    """
    is_angle(LST)
    is_angle(alpha)

    hs = (LST - alpha).to(u.deg)
    
    return constrain_angle_to_one_circle(hs)

def wavelength(frequency):
    """Convert frequency to wavelength.

    Parameters
    ----------
    frequency (Quantity < float, Hertz>)

    Returns
    -------
    Quantity <float, meters>
     """
    is_Quantity(frequency)
    return (c.c / frequency).to(u.m)

def Q_ew(dec, B_ew = 20 * u.m, wavelength = wavelength(10.7 * u.GHz)):
    """Compute the constants Q_ew for some parameters
    
    Parameters
    ----------
    dec Quantity <float, degrees or radians>
    Declination of target

    B_ew Quantity <meters>, default = 20 * u.m
    Basline distance in east west direction

    Wavelength Quantity <meters>, default = wavelength(10.7 * u.GHz)
    Wavelenth of the signal from the target

    Returns
    -------
    float

    """
    is_Quantity(dec)
    is_Quantity(B_ew)
    is_Quantity(wavelength)
    
    return (B_ew * np.cos(dec) / wavelength).value

def Q_ns(dec, B_ns = 20 * u.m, wavelength = wavelength(10.7 * u.GHz), L = 37.873199 * u.deg):
    """Compute the constants Q_ns for some parameters.
    
    Parameters
    ----------
    dec Quantity <float, degrees or radians>
    Declination of target

    B_ns Quantity <meters>, default = 20 * u.m
    Basline distance in north south direction

    Wavelength Quantity <meters>, default = wavelength(10.7 * u.GHz)
    Wavelenth of the signal from the target

    L Quantity <float, degrees or radians> default = 37.873199 * u.deg
    Terrestrial lattitude

    Returns
    -------
    float

    """
    is_Quantity(dec)
    is_Quantity(B_ns)
    is_Quantity(wavelength)
    is_angle(L)

    return (B_ns * np.cos(dec) * np.sin(L) / wavelength).value

def B_ew(dec, Q_ew_value, wavelength = wavelength(10.7 * u.GHz)):
    """Compute the baseline distance for some parameters.
    
    Parameters
    ----------
    dec Quantity <float, degrees or radians>
    Declination of target

    Q_ew float
    A given constant value.

    Wavelength Quantity <meters>, default = wavelength(10.7 * u.GHz)
    Wavelenth of the signal from the target

    Returns
    -------
    Quantity <meters>

    """
    is_Quantity(dec)
    is_Quantity(wavelength)
    
    return ( wavelength * Q_ew_value / np.cos(dec) ).to(u.m)

def B_ns(dec, Q_ns_value, wavelength = wavelength(10.7 * u.GHz), L = 37.873199 * u.deg):
    """Compute the baseline distance for some parameters.
    
    Parameters
    ----------
    dec Quantity <float, degrees or radians>
    Declination of target

    Q_ns float
    A given constant value.

    Wavelength Quantity <meters>, default = wavelength(10.7 * u.GHz)
    Wavelenth of the signal from the target

    L Quantity <float, degrees or radians> default = 37.873199 * u.deg
    Terrestrial lattitude

    Returns
    -------
    Quantity <meters>

    """
    is_Quantity(dec)
    is_angle(L)
    is_Quantity(wavelength)
    
    return ( wavelength * Q_ns_value / ( np.sin(L) * np.cos(dec) ) ).to(u.m)

def nu_tau_of_Q(h_s, Q_ew, Q_ns):
    """Express the delay in number of turns.
    
    equation 11, page 9 (with the substitution constants -> Q (see page 10).
    https://github.com/AaronParsons/ugradio/blob/master/lab_interf/interf.pdf
    
    Parameters
    -------
    h_s (Quantity <float, degrees or radians>)
    Hour angle of observation
    
    Q_ew (float)
    A constant value.
    
    Q_ns (float)
    A constant value.   
    
    Returns
    -------
    np.array() of Quantity <float, radians>
    
    """
    assert (np.isscalar(Q_ew) and np.isscalar(Q_ns)), "Both Q values must be scalars."
            
    is_angle(h_s)
    
    return (Q_ew * np.sin(h_s) + Q_ns * np.cos(h_s) ) * u.rad

def solve_x(A,b):
    """Solve Ax = b for x, where A does not have to be a square matrix.
    
    Parameters
    -------
    M (np.matrix(), shape = (n,m))
    The matrix of the coefficients of the system of equations.
    
    b (np.matrix(), shape = (n,1))
    The matrix (column vector) of the constants of the system of equations.
    
    
    Returns
    -------
    (A^T A)^(-1) * A^T * b
    
    """
    
    can_solve(A,b)
    
    return (A.T * A).I * A.T * b
            
def _fringe_term_A(hs, Q_ew, Q_ns):
    """First term of the equation, not including the constant multiple, A"""
    return np.cos(2 * np.pi * nu_tau_of_Q(hs, Q_ew, Q_ns)).value

def _fringe_term_B(hs,Q_ew, Q_ns):
    """Second term of the equation, not including the constant multiple, B"""
    return np.sin(2 * np.pi * nu_tau_of_Q(hs, Q_ew, Q_ns)).value  
                   
def fringe_amplitude_fit_matrix(h_s, Q_ew, Q_ns):
    """Compute the fit matrix for the given parameters.
    The fit matrix A must be multiplied by the fit constants (vector x)
    to get the fit.
        
    Parameters
    -------
    h_s : array:(Quantity <float, degrees or radians>) shape = (m,1)
    Hour angle of observation
    
    A (float)
    A constant value

    B (float)
    A constant value
    
    Q_ew (float)
    A constant value.
    
    Q_ns (float)
    A constant value.   
    
    
    Returns
    -------
    np.matrix() : array: shape = (m,2)
    Matrix of modeled fringe amplitude coefficients.
    """
    return np.matrix([_fringe_term_A(h_s, Q_ew, Q_ns),_fringe_term_B(h_s, Q_ew, Q_ns)]).T
            
def fringe_amplitude_fit(A, x):
    """returns a matrix

    Parameters
    ----------
    A :np.matrix(), shape = (m,n)

    x: np.matrix(), shape = (n,1)


    Returns
    -------
    np.matrix(), shape = (m,1)
    """
    can_multiply(A,x)
    return A * x

def sum_squared_residuals(actual, model, sigma = 1):
    """Compute the sum of the squared residuals
    Parameters
    ----------
    actual: np.matrix(), shape = (m,1)

    model: np.matrix(), shape = (m,1)

    sigma: scalar or np.array, shape = (m, ) (optional)

    Returns
    -------
    np.array(), shape = (m, ) 
    """
    return (((model.A1 - actual.A1)/(sigma))**2).sum()




def argmin_2D(A):
    """Get the indices of the minimum value of a 2D array

    Parameters
    ----------
    A np.matrix() (or np.array()) ndim = 2

    Returns
    -------
    (i, j)
    The minimum indicies.
    """
    return np.unravel_index(A.argmin(), A.shape)

def min_2D(A):
    """Get the minimum value of a 2D array

    Parameters
    ----------
    A np.matrix() (or np.array()) ndim = 2

    Returns
    -------
    The minimum value
    """
    x = min_2D(A)  
    return A[x[0],x[1]]

def get_polyfit(x,y,degree = 1):
    """TODO: DOCSTRING"""
    p = np.polyfit(x, y, degree)
    return  np.polyval(p, x)


def fringe_freqs(ha, dec, B_ew = 14.74438572777332 * u.m, B_ns = 1.4358938157675671 * u.m):
    return Q_ew(dec, B_ew = B_ew) * np.cos(ha) - Q_ns(dec, B_ns = B_ns) * np.sin(ha)

# def S_mins_image(Q_ns_range, Q_ew_range, ha, b, filename = None):
    
#     num_Q = len(Q_ns_range)
#     S_2D = np.zeros((num_Q, num_Q))
    
#     # compute minimums
#     for i in range(num_Q):#ns
#         for j in range(num_Q):#ew
#             M = fringe_amplitude_fit_matrix(ha, Q_ew_range[j], Q_ns_range[i])
#             x = solve_x(M,b)
#             fit = fringe_amplitude_fit(M,x)
#             S_2D[i, j] = sum_squared_residuals(b,fit)
#             sys.stdout.write('\r'+str('%d, %d' %(i,j) ))
#     # show plot       
#     plt.imshow(np.log10(S_2D), cmap = 'gray')
#     plt.xlabel('Q_ew indices')
#     plt.ylabel('Q_ns indices')
#     plt.title('S(Q_ns,Q_ew,ha)')
#     plt.colorbar()
#     plt.show()
    
#     if filename != None:
#         np.savez(filename, S_2D)
    
#     # return Q values of minimium
#     return argmin_2D(S_2D) 

# def Q_to_baseline(Q_mins, Q_ns_range, Q_ew_range):
#     Q_min_ns = Q_ns_range[Q_mins[0]]
#     Q_min_ew = Q_ew_range[Q_mins[1]]
    
#     print(B_ns(orion_dec, Q_min_ns))
#     print(B_ew(orion_dec, Q_min_ew))