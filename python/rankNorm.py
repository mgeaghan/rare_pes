# Rank normalisation function
# Author: Michael Geaghan
# Date: 2020-05-29
import numpy as np
import pandas as pd
from scipy.stats import norm

def rankNorm(x, method = "general", alpha = np.pi/8,
         complete = False, na_option = "keep", skipna = True,
         adjustN = True, min_val = 1, max_val = 10, **kwargs):
    x = pd.Series(x)
    if complete is True:
        x.dropna(inplace=True)
    elif complete is not False:
        raise ValueError("Argument to 'complete' must be of type bool.")
    Ranks = x.rank(na_option = na_option, **kwargs)
    if adjustN is True:
        N = len(x.dropna())
    elif adjustN is False:
        N = len(x)
    else:
        raise ValueError("Argument to 'adjustN' must be of type bool.")
    if method == "blom":
        Score = norm.ppf((Ranks - 0.375)/(N + 0.25))
    elif method == "vdw":
        Score = norm.ppf((Ranks)/(N + 1))
    elif method == "tukey":
        Score = norm.ppf((Ranks - 1/3)/(N + 1/3))
    elif method == "rankit":
        Score = norm.ppf((Ranks - 1/2)/(N))
    elif method == "elfving":
        Score = norm.ppf((Ranks - np.pi/8)/(N - np.pi/4 + 1))
    elif method == "general":
        Score = norm.ppf((Ranks - alpha)/(N - 2 * alpha + 1))
    elif method == "zscore":
        x_mean = x.mean(skipna=skipna)
        x_std = x.std(skipna=skipna)
        Score = (x - x_mean)/x_std
    elif method == "scale":
        x_min = x.min(skipna=skipna)
        x_max = x.max(skipna=skipna)
        Score = (((x - x_min)*(max_val - min_val))/(x_max - x_min) + min_val)
    else:
        raise ValueError("Invalid argument to 'method'.")
    return(Score)
