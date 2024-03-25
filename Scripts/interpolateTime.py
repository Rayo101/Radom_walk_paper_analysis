import numpy as np

def interpolateUnique(msTime: np.ndarray) -> np.ndarray:
    """
    Function that interpolates the LOR times linearly.  

    Input: 
        msTime: The raw output time
    
    Output:
        newTime: The interpolated time
    
    """

    uniqueTimes, uniqueTimesIndex, counts = np.unique(msTime, return_index = True, return_counts = True) # Find the unique times and the indices where they occur
    newTime = [] # list to store new times

    for i in range(len(uniqueTimes)-1): # loop over unique times
        timeInterval = msTime[uniqueTimesIndex[i]: uniqueTimesIndex[i+1]+1] # find interval
        newTimeInterval = np.linspace(timeInterval[0], timeInterval[-1], len(timeInterval))[0:-1] # interpolate linearly betwen the times
        newTime.extend(newTimeInterval.tolist()) # add new time to the newTime list

    tNew = np.linspace(uniqueTimes[-1], uniqueTimes[-1] + 1, counts[-1]) # This is to handle the edge case of the final timestamp
    newTime.extend(tNew) # extend with edge case

    return newTime
