def peak_finder(signal, threshold):

    peaks = [] #x positions of the peaks, or rather, their index

    for i in range(2,len(signal)-2):

        if signal[i - 2] < signal[i] and signal[i - 1] < signal[i]and signal[i] > signal[i + 1] and signal[i] > signal[i + 2]: #four
                                    #conditions to be a peak (see description)
            if signal[i] > threshold: #is the value of the spectrum at i higher than our
                            #threshold?

                peaks.append(i)
    return peaks