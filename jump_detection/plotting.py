"""
Module containing code for plotting a few key plots.
"""

import matplotlib.pyplot as plt
from .utils import calculate_median_jump


def plot_segment(segment):
    """
    Plots a histogram.

    Parameters:

    """
    _, (ax1, ax3) = plt.subplots(2, 1)

    column1 = segment.original[:, 0]
    column2 = segment.original[:, 1]

    # Plot the first column on the left y-axis
    ax1.plot(column1, 'b-', label='Column 1')
    ax1.set_xlabel('Index')
    ax1.set_ylabel('Mode 1', color='b')
    ax1.tick_params('y', colors='b')

    # Create a twin axes object and plot the second column on the right y-axis
    ax2 = ax1.twinx()
    ax2.plot(column2, 'r-', label='Column 2')
    ax2.set_ylabel('Mode 2', color='r')
    ax2.tick_params('y', colors='r')

    ax3.plot(segment.f_stats)
    plt.show()


def plot_jump_signature(segment_list, labels):
    """
    Plots the jump signature based on median.

    Parameters:

    """
    median_jump = calculate_median_jump(segment_list, labels)

    _, ax1 = plt.subplots(1, 1)

    column1 = median_jump[:, 0]
    column2 = median_jump[:, 1]

    ax1.set_title("Jump Signature")
    # Plot the first column on the left y-axis
    ax1.plot(column1, 'b-', label='Column 1')
    ax1.set_xlabel('Index')
    ax1.set_ylabel('Mode 1', color='b')
    ax1.tick_params('y', colors='b')

    # Create a twin axes object and plot the second column on the right y-axis
    ax2 = ax1.twinx()
    ax2.plot(column2, 'r-', label='Column 2')
    ax2.set_ylabel('Mode 2', color='r')
    ax2.tick_params('y', colors='r')

    plt.show()
