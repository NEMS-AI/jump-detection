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

    def calculate_features(self):
        """
        Calculate set of reduced features and store

        Parameters:
        -----------
        parameter1 : type
            Desired p-value for bootstrap calculation.
        """
        self.features = get_peak_features(self.Fstats)