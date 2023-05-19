import numpy as np
import pandas as pd
from .utils import *
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

    def __init__(self, original, Fstats):
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
        self.Fstats = Fstats
        self.features = []
        self.diff = 0

    def calculate_features(self):
        """
        Calculate set of reduced features and store

        Parameters:
        -----------
        parameter1 : type
        """
        self.features = get_peak_features(self.Fstats)

    def calculate_freq_shift(self, window_size):
        """
        For a single segment, calculate corresonding relative frequency shift

        Parameters:
        -----------
        parameter1 : type
        """
        x1 = np.mean(self.original[0:window_size], axis = 0)
        x2 = np.mean(self.original[-window_size :], axis = 0)
        self.diff = (x2 - x1) / x1
        pass