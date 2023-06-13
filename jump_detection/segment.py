"""
Module containing the definition of Segment.

This module defines the Segment class and its associated methods.
It provides functionality for storing each individual jump segments
"""

import numpy as np
from .utils import get_peak_features

class Segment:
    """
    Class representing a time series segment associated with a peak.

    Attributes:
    -----------
    original : pd.DataFrame
        The original time series segment.
    Fstats : pd.DataFrame
        The moving Fstats time series segment.
    """

    def __init__(self, original, f_stats):
        """
        Initialize a Segment.

        Parameters:
        -----------
        original : pd.DataFrame
            The original time series segment.
        covariance : pd.DataFrame
            The moving covariance time series segment.
        """
        self.original = original
        self.f_stats = f_stats
        self.features = get_peak_features(self.f_stats, 'normalized')
        self.diff = 0

    def calculate_freq_shift(self, window_size):
        """
        For a single segment, calculate corresonding relative frequency shift

        Parameters:
        -----------
        parameter1 : type
        """
        x_1 = np.mean(self.original[0:window_size], axis=0)
        x_2 = np.mean(self.original[-window_size:], axis=0)
        self.diff = (x_2 - x_1) / x_1
        
    def get_features(self):
        '''
        A getter method for retrieving the features of a particular jump.
        '''
        return self.features