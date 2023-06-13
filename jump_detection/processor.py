"""
Module containing the definition of TimeSeriesProcessor.

This module defines the TimeSeriesProcessor class and its associated methods.
It provides functionality for processing entire time-series datasets
and storing each time-series segment around detected jumps.
"""

import numpy as np
import pandas as pd
from .segment import Segment
from .utils import find_peaks_in_data, normalize_features, segment_data
from .rolling_Ftest import rolling_F_statistic


class TimeSeriesProcessor:
    """
    Class for processing time series data.

    Attributes:
    -----------
    window_size : int
        The window size for moving calculations.
    """

    def __init__(self, window_size=10, gap_size=0) -> None:
        """
        Initialize a TimeSeriesProcessor.

        Parameters:
        -----------
        window_size : int
            The window size for moving calculations; represents number of 
            time steps within each side of window
        """
        self.window_size = window_size
        self.gap_size = gap_size
        self.data = None
        self.moving_fstats = None
        self.segments = []
        self.valid_jumps = []


    def load_data(self, filename):
        """
        Load multidimensional time series data from a file.

        Parameters:
        -----------
        filename : str
            The name of the file to load.
        """
        self.data = pd.read_csv(filename).values

    def process_data(self):
        """
        Process the time series data, creating segments around the peaks.
        """
        moving_fstat = rolling_F_statistic(
            self.data, self.window_size, self.window_size, self.gap_size)
        peaks = find_peaks_in_data(moving_fstat)
        self.moving_fstats = moving_fstat

        # Segment out
        # Offset time segments duw to different in original timeseries and
        # fstat
        t_offset = int(self.window_size + self.gap_size / 2)
        original_segments = segment_data(
            self.window_size, self.gap_size, peaks + t_offset, self.data)
        fstat_segments = segment_data(
            self.window_size, self.gap_size, peaks, moving_fstat)

        self.segments = [Segment(original, Fstats) for original, Fstats in zip(
            original_segments, fstat_segments)]

    def get_all_features(self):
        """
        Iterate through all collected event segments and return a ndarray with all the

        Returns:
        --------
        pd.DataFrame
            The loaded data.
        """
        jump_features = []
        for _, segment in enumerate(self.segments):
            jump_features.append(segment.features)
        jump_features = np.array(jump_features)
        normalized_features, self.valid_jumps = normalize_features(
            jump_features)
        return normalized_features

    def get_all_diffs(self):
        '''
        Iterate through all the collected events, calculate the relative frequency shift
        and then return all of the frequency shifts.
        '''
        diffs = []
        for segment in self.segments:
            segment.calculate_freq_shift(self.window_size)
            diffs.append(segment.diff)
        diffs = np.array(diffs)
        return diffs[self.valid_jumps]
