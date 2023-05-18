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


def find_peaks_in_data(data):
    """
    Description of method.

    Parameters:
    -----------
    parameter1 : type
        Description of parameter

    #TODO Update peakinding to match matlab implementation
    """
    peaks, _ = find_peaks(data, prominence=0.01)
    return peaks

def segment_data(window_size, peaks, data):
    """
    Description of method.

    Parameters:
    -----------
    parameter1 : type
        Description of parameter
    """
    segments = []
    for peak in peaks:
        start = max(0, peak - window_size)
        end = min(len(data), peak + window_size)
        segments.append(data[start:end])
    return segments

# ------------------------------------------------------------------------------
# Methods for Step 2: Identifying Single-Event Jumps
# ------------------------------------------------------------------------------


def compute_fwhm(values):
    """
    Description of method.

    Parameters:
    -----------
    parameter1 : type
        Description of parameter
    """
    # Compute the maximum value and its index
    max_value = np.max(values)
    max_index = np.argmax(values)

    # Compute the half maximum value
    half_max_value = max_value / 2

    # Find the indices where the values cross the half maximum value
    left_index = np.argmin(np.abs(values[:max_index] - half_max_value))
    right_index = max_index + np.argmin(np.abs(values[max_index:] - half_max_value))

    # Compute the FWHM
    fwhm = right_index - left_index

    return fwhm

def get_peak_stats(data):
    """
    Description of method.

    Parameters:
    -----------
    parameter1 : type
        Description of parameter
    """
    # TODO: Check FWHM Code
    # FWHM = compute_fwhm(df['FStat'])
    FWHM = 0
    M1 = np.mean(data)
    M2 = np.var(data)
    M3 = skew(data)
    M4 = kurtosis(data)
    return ((FWHM, M1, M2, M3, M4))


def normalize_features(df, features):
    """
    Description of method.

    Parameters:
    -----------
    parameter1 : type
        Description of parameter
    """

    features = ['M1', 'M2', 'M3']
    subset = df[features]
    normalized = (subset - subset.mean()) / subset.std()
    return normalized

def get_eps(normalized):
    """
    From a normalized df, compute use knee point to determine eps

    Parameters:
    -----------
    parameter1 : type
        Description of parameter
    
    
    TODO: implement proportion based eps similar to matlab
    """

    # Compute the distance matrix based on NN
    neigh = NearestNeighbors(n_neighbors=2*normalized.shape[1])
    nbrs = neigh.fit(normalized)
    distances, _ = nbrs.kneighbors(normalized)

    # Sort the distances
    sorted_distances = np.sort(distances, axis=None)

    # Use knee finder function to estimate knee_eps
    kneedle = KneeLocator(np.arange(len(sorted_distances)), sorted_distances, S=1.0, curve="convex", direction="increasing")
    knee_eps = sorted_distances[kneedle.knee]

    return knee_eps

def get_class_labels(normalized, eps):
    """
    Use DB scan to find labeling based on previously determined eps

    Parameters:
    -----------
    parameter1 : type
        Description of parameter
    
    """
    eps = eps  
    dbscan = DBSCAN(eps=eps, min_samples=normalized.shape[1])
    dbscan.fit(normalized)

    return(dbscan.labels_)

# ------------------------------------------------------------------------------
# Methods for Step 3: Measure Events
# ------------------------------------------------------------------------------

