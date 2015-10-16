from scipy import signal

def hanning(x,N):
    """ Filter a time series x with a Hanning window of length N '
    Inputs:
    x - a numpy array to be filtered
    N - width of window
    Output: numpy array of filtered time series
    """
    win = signal.hann(N)
    filtered = signal.convolve(x,win,mode='same')/sum(win)
    
