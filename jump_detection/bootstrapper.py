import numpy as np
import pandas as pd
from .segment import Segment
from .utils import *
from .rolling_Ftest import *
class NoiseBootstrapper:
    """
    Class for processing noise data .

    Attributes:
    -----------
    p-value : double
        Desired p-value for bootstrap calculation.
    """
    def __init__(self, p_value) -> None:
        """
        Initialize a TimeSeriesProcessor.

        Parameters:
        -----------
        p-value : double
        Desired p-value for bootstrap calculation.
        """
        self.p_value = 0.003


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
        self.data = pd.read_csv(filename)
        return self.data

    def process_data(self, filename):
        """
        Process the time series data, returning values of .

        Parameters:
        -----------
        filename : str
            The name of the file to load.

        Returns:
        --------
        (window_sizes, thresholds): A tuple containing two lists containing the window sizes and 
        associated F-stat thresholds.
        """
        
        # Initialize Paramters for moving block bootstarpping
        #TODO: Move window_size range to paramters
        min_window_size = 10
        max_window_size = 100
        F_stat_tresh = 0
        window_sizes = []
        F_stat_tresh_vals = []

        data = self.load_data(filename).values
        # Iterate through possible window sizes, and get corresponding F_stat threshold
        for window_size in range(min_window_size, max_window_size, 10):
            F_stats = rolling_F_statistic(data, window_size, window_size, 1)
            F_stat_tresh = np.percentile(F_stats, (1-self.p_value)*100)
            window_sizes.append(window_size)
            F_stat_tresh_vals.append(F_stat_tresh)
        

        return (window_sizes, F_stat_tresh_vals)