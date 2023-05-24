import numpy as np
import pandas as pd
from .segment import Segment
from .utils import *
from .rolling_Ftest import *

class TimeSeriesProcessor:
    """
    Class for processing time series data.

    Attributes:
    -----------
    window_size : int
        The window size for moving calculations.
    """

    def __init__(self, window_size = 10, gap_size = 0) -> None:
        """
        Initialize a TimeSeriesProcessor.

        Parameters:
        -----------
        window_size : int
            The window size for moving calculations; represents number of time steps within each side of window
        """
        self.window_size = window_size
        self.gap_size = gap_size

    def load_data(self, filename):
        """
        Load multidimensional time series data from a file.

        Parameters:
        -----------
        filename : str
            The name of the file to load.

        Returns:
        --------
        pd.DataFrame
            The loaded data.
        """
        self.data = pd.read_csv(filename).values
        return self.data


    def process_data(self):
        """
        Process the time series data, creating segments around the peaks.

        Parameters:
        -----------
        filename : str
            The name of the file to load.

        Returns:
        --------
        list of Segment
            The time series and Fstat segments associated with each peak.
        """
        moving_fstat = rolling_F_statistic(self.data, self.window_size, self.window_size, self.gap_size)
        peaks = find_peaks_in_data(moving_fstat)
        self.moving_fstats = moving_fstat

        # Segment out
        # Offset time segments duw to different in original timeseries and fstat
        t_offset = int(self.window_size+self.gap_size/2)
        original_segments = segment_data(self.window_size, self.gap_size, peaks + t_offset, self.data)
        fstat_segments = segment_data(self.window_size,self.gap_size, peaks,  moving_fstat)
        
        segments = [Segment(original, Fstats) for original, Fstats in zip(original_segments, fstat_segments)]
        return segments
    
