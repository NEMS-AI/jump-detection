"""
Rolling F-test

Methods compute a rolling F-test for a pair of moving windows.
"""

import numpy as np
import scipy as sc
from numba import njit

# ------------------------------------------------------------------------------
# Convenience methods for testing purposes
# ------------------------------------------------------------------------------
def t2_statistic(X_samples, Y_samples):
    """Return the multivariate t2 statistic for collections X & Y
    
    Parameters
    ----------
    X_samples : (nx, d) ndarray
        samples of d-dimensional random vector X
    Y_samples : (ny, d) ndarray
        samples of d-dimensional random vector Y
        
    Returns
    -------
    float
        t2 statistic
    """
    X_samples = cast_to_matrix(X_samples)
    Y_samples = cast_to_matrix(Y_samples)
    
    nx, d = X_samples.shape
    ny, d = Y_samples.shape
    
    X_bar = np.mean(X_samples, axis=0)
    Y_bar = np.mean(Y_samples, axis=0)
    
    sigma_X = sample_cov_multivariate(X_samples)
    sigma_Y = sample_cov_multivariate(Y_samples)
    sigma = 0.5 * (sigma_X + sigma_Y)
    
    
    xmy = X_bar - Y_bar
    rhs_vec = np.linalg.solve(sigma, xmy)
    
    t2_stat = (nx * ny) / (nx + ny) * np.dot(xmy, rhs_vec)
    return t2_stat
    
    
def F_statistic(X_samples, Y_samples):
    """Return the multivariate F statistic for collections X & Y
    
    Parameters
    ----------
    X_samples : (nx, d) ndarray
        samples of d-dimensional random vector X
    Y_samples : (ny, d) ndarray
        samples of d-dimensional random vector Y
        
    Returns
    -------
    float
        F statistic
    """
    X_samples = cast_to_matrix(X_samples)
    Y_samples = cast_to_matrix(Y_samples)
    
    nx, d = X_samples.shape
    ny, d = Y_samples.shape
    
    t2_stat = t2_statistic(X_samples, Y_samples)
    
    F_stat = (nx + ny - d - 1) / (d * (nx + ny - 2)) * t2_stat
    return F_stat
    

def sample_cov_multivariate(X_samples):
    """Return unbiased estimator of multivariate sample covariance
    
    Parameters
    ----------
    X_samples : (nx, d) ndarray
        samples of d-dimensional random vector X
    
    Returns
    -------
    (d, d) ndarray
        unbiased estimator of sample covariance
    """
    X_samples = cast_to_matrix(X_samples)
    
    nx, d = X_samples.shape
    X_bar = np.mean(X_samples, axis=0)
    
    acc = np.zeros((d, d))
    
    for i in range(nx):
        diff = X_samples[i, :] - X_bar
        acc += diff[:, np.newaxis] * diff
    
    cov = acc / (nx - 1)
    return cov
    

def cast_to_matrix(X):
    """Return matrix form of X by casting to ndarray"""
    X = np.asanyarray(X)
    if X.ndim == 1:
        X = X.reshape(-1, 1)
    return X


def rolling_F_statistic_slow(X_samples, n1, n2, g):
    X_samples = cast_to_matrix(X_samples)
    N, d = X_samples.shape 
    N_f_stats = N - n1 - n2 - g + 1
  
    f_stats = np.zeros(N_f_stats)
    for i in range(N_f_stats):
        W1 = X_samples[i:i + n1]
        W2 = X_samples[i + n1 + g:i + n1 + n2 + g]
        f_stats[i] = F_statistic(W1, W2)
        
    return f_stats


# ------------------------------------------------------------------------------
# Performant methods for production
# ------------------------------------------------------------------------------
@njit
def self_outer(x):
    """Return outer product of vector x"""
    return np.outer(x, x)


@njit
def _sigma(n, W_bar, W_sq):
    return W_sq - n / (n - 1) * self_outer(W_bar)


@njit
def _F_stat(d, n1, n2, W1_bar, W2_bar, W1_sq, W2_sq):
    W1_sigma = _sigma(n1, W1_bar, W1_sq)
    W2_sigma = _sigma(n2, W2_bar, W2_sq)
    
    sigma = 0.5 * (W1_sigma + W2_sigma)

    W_bar = W1_bar - W2_bar
    rhs_vec = np.linalg.solve(sigma, W_bar)
    
    F_stat_unnormed = np.dot(W_bar, rhs_vec)
    return F_stat_unnormed


@njit('float64[:, :](float64[:, :])')
def _expectation_sq(X):
    n, d = X.shape
    exp_sq = np.zeros((d, d))
    for i in range(n):
        for j in range(d):
            for k in range(d):
                exp_sq[j, k] += X[i, j] * X[i, k]
    exp_sq /= n - 1
    return exp_sq


@njit #('(float64[:, :], int64, int64, int64)')
def rolling_F_statistic(X_samples, n1, n2, g):
    """Return F-statistic of rolling window over sequence of random samples
    
    To detect changes in noisy multivariate time-series data it is useful
    to consider a null hypothesis for a pair of moving windows, namely,
    their is no difference between their average values. The statistical
    approach of analysis of variance (ANOVA) provides a formal framework
    for asking this question by conducting an F-test comparing the pair
    of windows.
    
    For this filter we consider a pair of windows with lengths `n1` and `n2`
    with a gap of length `g` in between.
    
    
              Window 1            Gap              Window 2     
        [ ----- n1 ----- ) [ ----- g -----) [ ------ n2 -------)
        |                |                  |                  |
        i           i + n1 - 1         i + n1 + g         i + n1 + n2 + g - 1
    
    
    For a sequence of d-dimensional random samples in array (N, d) in `X_samples`
    at step `i = 0, ..., N - n1 - n2 - g` the windows range over indices,
    
        window 1:  [i,          ..., i + n1 - 1]            (inclusive)
        window 2:  [i + n1 + g, ..., i + n1 + n2 + g - 1]   (inclusive)
        
    Parameters
    ----------
    X_samples : (N, d) ndarray
        random d-dimensional random samples
    n1 : int
        number of samples in window 1
    n2 : int
        number of samples in window 2
    g : int
        number of samples between windows
    
    Returns
    -------
    (N - n1 - n2 - g + 1) ndarray
        F-statistic of moving window
    """
    N, d = X_samples.shape
    C = (n1 + n2 - d - 1) / (d * (n1 + n2 - 2)) * n1 * n2 / (n1 + n2)
    
    N_f_stats = N - n1 - n2 - g + 1
    
    # Initialise Window means and square expectation
    W1_init = X_samples[:n1, :]
    W2_init = X_samples[n1 + g:n1 + n2 + g, :]
    
    W1_bar = np.sum(W1_init, axis=0) / n1
    W2_bar = np.sum(W2_init, axis=0) / n2
    
    W1_sq = _expectation_sq(W1_init)
    W2_sq = _expectation_sq(W2_init)
    
    
    # Compute rolling filters
    f_stats = np.zeros(N_f_stats)
    f_stats[0] = C * _F_stat(d, n1, n2, W1_bar, W2_bar, W1_sq, W2_sq)
    
    for i in range(0, N_f_stats - 1):
        W1_bar += (X_samples[i + n1] - X_samples[i]) / n1
        W2_bar += (X_samples[i + n1 + n2 + g] - X_samples[i + n1 + g]) / n2

        W1_sq += (self_outer(X_samples[i + n1]) - self_outer(X_samples[i])) / (n1 - 1)
        W2_sq += (
            self_outer(X_samples[i + n1 + n2 + g]) - self_outer(X_samples[i + n1 + g])
        ) / (n2 - 1)
        
        f_stats[i + 1] = C * _F_stat(d, n1, n2, W1_bar, W2_bar, W1_sq, W2_sq)
        
    return f_stats

# ------------------------------------------------------------------------------
# Implementation Notes
# ------------------------------------------------------------------------------
# Code above almost realises theoretical speed up possible with rolling F-test
# implementation but could be improved by moving to C implementation. This would
# allow control over memory allocation in the main loop.
#
# The code above does not take advantage of the symmetry of the square expectation
# and covariance matrices. This could have the numbe of multiplications in the
# main loop.
# ------------------------------------------------------------------------------


def rolling_times(time_stamps, n1, n2, g):
    """Return the window center times for rolling F-test
    
    WARNING: Non-uniformly spaced time-samples will ruin your day!!!
    
    The center time of the window is defined as the time-stamp associated with the 
    sample in the mid-point of the gap. For the i-th step the gap has indices
    
        gap indices : [i + n1, ..., i + n1 + g - 1]   (inclusive)
        
        
        ts_i = [ time_stamps[t0]                         : if g odd
               [ (time_stamps[t0] + time_stamps[t0 + 1]) : else
                            
    where t0 = i + n1 + (g - 1) // 2 and  i = 0, ..., N - n1 - n2 - g.
    
    Parameters
    ----------
    time_stamps : (N,) ndarray
        uniformly spaced time-stamps of random samples
    n1 : int
        number of samples in window 1
    n2 : int
        number of samples in window 2
    g : int
        number of samples between windows
    
    Returns
    -------
    (N - n1 - n2 - g + 1) ndarray
        time-stamps of center of moving window
    """
    N = len(time_stamps)
    t0 = n1 + (g - 1) // 2
    t_end = N - n1 - n2 - g + t0
    
    if g % 2 == 0:
        ts = 0.5 * (time_stamps[t0:t_end + 1] + time_stamps[t0 + 1:t_end + 2])
    else:
        ts = time_stamps[t0:t_end + 1]
    
    return ts


def F_statistic_conf(F_stats, d, n1, n2):
    """Return confidence in rejecting null hypothesis of F-test
    
    Parameters
    ----------
    F_stats : (N_fstats, d) ndarray
        random d-dimensional random samples
    n1 : int
        number of samples in window 1
    n2 : int
        number of samples in window 2
    g : int
        number of samples between windows
    
    Returns
    -------
    (N_fstats) ndarray
        confidence in rejecting null hypothesis of F-test
    """
    dist = sc.stats.f(d, n1 + n2 - d - 1)
    return dist.sf(F_stats)
