import pandas as pd
from scipy.stats import f, skew, kurtosis
import numpy as np
from sklearn.cluster import DBSCAN
from sklearn.neighbors import NearestNeighbors
from scipy.stats import mode
from kneed import DataGenerator, KneeLocator
from scipy.signal import find_peaks

# ------------------------------------------------------------------------------
# Methods for Step 1: Identify all jump events
# ------------------------------------------------------------------------------

def read_data(*args, **kwargs):
    return pd.read_csv(*args, **kwargs)

def calculate_moving_average(data, window_size):
    """
    Method designed to emulate final F-stat calculation

    Parameters:
    -----------
    parameter1 : type
        Description of Paramter

    # TODO: Remove usage of this function
    """
    moving_avg = data.rolling(window=window_size).mean()
    moving_avg = moving_avg.fillna(0) 
    moving_avg = moving_avg.mean(axis=1)  # average across all dimensions
    pd.DataFrame(moving_avg)
    return moving_avg


def find_peaks_in_data(data, height = 300):
    """
    Given a 1D array, find all local maxima. 

    Parameters:
    -----------
    parameter1 : type
        Description of parameter
    height : int
        The vertical threshold for which points are considered to be peaks

    #TODO Update peak finding to match matlab implementation
    """
    peaks, _ = find_peaks(data, height=height)
    return peaks

def segment_data(window_size, gap_size, peaks, data):
    """
    For each peak passed through, capture a window around the peak. The following 
    is a rough sketch of what is stored as a segment

            Window 1            Gap              Window 2     
        [ ----- n1 ----- ) [ ----- g -----) [ ------ n2 -------)
        |                |                  |                  |
        peak-window-gap/2          peak                 peak + window + gap/2
    
    Parameters:
    -----------
    parameter1 : type
        Description of parameter
    """
    segments = []
    for peak in peaks:
        start = max(0, int(peak - window_size - gap_size/2))
        end = min(len(data), int(peak + window_size + gap_size/2))
        segments.append(data[start:end])
    return segments

# ------------------------------------------------------------------------------
# Methods for Step 2: Identifying Single-Event Jumps
# ------------------------------------------------------------------------------


def compute_fwhm(values):
    """
    A helper method designed to calculate the full width at half maximum (FWHM).

    Parameters:
    -----------
    parameter1 : type
        The 1D array whose FWHM will be calculated.
    """
    # Compute the maximum value and its index
    max_value = np.max(values)
    max_index = np.argmax(values)

    # Compute the half maximum value
    half_max = max_value / 2

    # Find the indices where the values is above the half maximum
    above_half_max = np.where(values > half_max)[0]

    # Find the indices of the first and last points above half maximum
    first_index = above_half_max[0]
    last_index = above_half_max[-1]

    # Calculate the FWHM
    fwhm = last_index - first_index

    return fwhm

def get_peak_features(data):
    """
    Given an ndarray array, this method will compute various summary features that are 
    used to reduce the dimensionality of the array.

    Parameters:
    -----------
    data : numpy.ndarray
        An array with shape ()
    """
    # TODO: Check FWHM Code
    FWHM = compute_fwhm(data)
    M1 = np.mean(data)
    M2 = np.var(data)
    M3 = skew(data)
    M4 = kurtosis(data)
    return ((FWHM, M1, M2, M3, M4))


def normalize_features(features, mode = 'normalize'):
    """
    Given a  ndarray array of features, this method will normalize each feature such that they 
    are given equal weight for clustering. 

    
    Parameters:
    -----------
    features : numpy.ndarray
        An array with shape (n_events, n_features). Each feature is to be normalized such that it lies between 0 and 1. 
    mode : string
        A string indicating which type of normalization scheme should be used when normalizing
            - 'normalize'
                Force each feature to be in the range [0,1]
            - 'quartile' 
                A generalization of normalize in which 0 and 1 are reserved for the 25th and 75th
            - 'stanrdize'
                Use z-score normalization to ensure each feature has mean 0 and standard deviation 1.

    """
    if mode == 'normalize':
        data_min = np.min(features, axis=0) 
        data_max = np.max(features, axis=0)  
        # Subtract min and divide by range to normalize to [0, 1]
        normalized = (features - data_min) / (data_max - data_min)

    elif mode == "quartile":
        data_low_quartile = np.percentile(features, 25, axis=0)
        data_high_quartile = np.percentile(features, 75, axis=0)
        # Subtract low quartile and divide by interquartile range to normalize
        normalized = (features - data_low_quartile) / (data_high_quartile - data_low_quartile)
        
    elif mode == "standardize":
        normalized = (features - features.mean(axis = 0)) / features.std(axis = 0)

    return normalized

def get_eps(features, mode = 'knee', proportion = None):
    """
    From a given an ndarray of features, compute an estimated value for eps based on noise

    Parameters:
    -----------
    features : numpy.ndarray
        An array with shape (n_events, n_features). Each feature is to be normalized such that it lies between 0 and 1. 
    mode : string
        A string indicated which method should be used for determing value for eps
            - 'knee'
                Find the knee point for the sorted distances and report as estimated eps
            - 'proportion'
                Find the value of eps that produces a cluster based on the proportion parameter
    proportion : float
        A decimal between 0 and 1 representing the proportion of the data to be normalized.
        For instance, when in proportion mode, proportion = 0.5 indicates the largest cluster 
        will contain around 50% of the data.
    
    TODO: Implement proportion mode
    """

    # Compute the distance matrix based on NN
    neigh = NearestNeighbors(n_neighbors=2*features.shape[1])
    nbrs = neigh.fit(features)
    distances, _ = nbrs.kneighbors(features)

    # Sort the distances
    sorted_distances = np.sort(distances, axis=None)

    # Use knee finder function to estimate knee_eps
    kneedle = KneeLocator(np.arange(len(sorted_distances)), sorted_distances, S=1.0, curve="convex", direction="increasing")
    knee_eps = sorted_distances[kneedle.knee]

    return knee_eps

def get_class_labels(features, eps):
    """
    Use DBSCAN to asign the labels based on the previously determined eps. Note this method is different than traditional
    DBSCAN as it will assign 1 to the largest cluster and 0 to all other clusters and/or noise.

    Parameters:
    -----------
    features : numpy.ndarray
        An array with shape (n_events, n_features). Each feature is to be normalized such that it lies between 0 and 1. 
    eps : float
        A value indicated the eps used for the DBSCAN algorithm
    """
    # Perfrom traditional DBSCAN
    eps = eps  
    dbscan = DBSCAN(eps = eps, min_samples = 2 * features.shape[1])
    dbscan.fit(features)
    labels = dbscan.labels_

    # Exclude the value -1 from consideration
    valid_labels = labels[labels != -1]

    # Find the most populous value (largest cluster)
    most_populous_value = np.bincount(valid_labels).argmax()

    # Create the binary representation
    binary_representation = np.where(labels == most_populous_value, 1, 0)

    return(binary_representation)

# ------------------------------------------------------------------------------
# Methods for Step 3: Measure Events
# ------------------------------------------------------------------------------

def normalize_jump(timeseries, a = 0, b = 1):
    """
    Helper function designed to normalize a time series segment to be between the range of a and b.
    Each mode present in the time series should be normalized independently. 

    Parameters:
    -----------
    timeseries : ndarray
        An ndarray in the shape (n_sampes, n_modes)
    a : float
        Maximum value of new time series
    b : float
        Minimum value of new time series
    """
    #TODO: Allow for different normalization schemes similar to normalize jump

    # To get values between a and b:
    # Normalize Y = a + (X - X_min) * (b - a) / (X_max - X_min)
    array_min = timeseries.min(axis=0)
    array_max = timeseries.max(axis=0)
    normalized_array = a + ((timeseries - array_min) * (b - a)) / (array_max - array_min)
    return normalized_array

def calculate_median_jump(segment_list):
    """
    Given a list of segments, this function will compute the "median jump" othewise known as the "jump signature" for all the segments

    Parameters:
    -----------
    segment_list : list of ndarrays
        Description of parameter
    """
    jump_shape = segment_list[int(len(segment_list)/2)].original.shape
    segment_arr = np.array([normalize_jump(segment.original) for segment in segment_list if segment.original.shape == jump_shape])
    median_jump = np.median(segment_arr, axis=0)

    
    return median_jump

