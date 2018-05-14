

def multiple_peak_finder(signal, threshold = None, width = 10, separation = 10):
    """
    Finds many peaks in a set of data
    
    Parameters
    ----------
    
    signal (numpy array) : The data.
    threshold (float) default = np.mean(sigal):  Data below this value is not considered to be a peak. 
    
    width (int) default = 1: How many data points on each side to compare a possible peak to.
                             Increasing this reduces the number of found peaks.
                             
    separation (int) default = 1: How far apart the peaks must be ( in indices)
                                  Increasing this reduces the number of found peaks.
    
    Returns
    -------
    returns a list of the INDEX locations of the peaks from the signal
    """

    peak_indicies = []

    for i in range(width,len(signal) - width):
        
        left_wing = signal[i - width: i]
        right_wing = signal[i + 1: i + 1 + width]
        
        if threshold == None:
            threshold = np.mean(signal)
        
        if signal[i] > threshold:
            
            if np.all(left_wing < signal[i]) == True:
                if np.all(signal[i] > right_wing) == True:
                    
                    if len(peak_indicies) == 0:
                        peak_indicies.append(i)
                        
                    if len(peak_indicies) > 0 and abs(peak_indicies[-1] - i) >= separation:
                        peak_indicies.append(i)
    return peak_indicies