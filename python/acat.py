# Author: Michael Geaghan
# Date: 2020-06-16
# Adapted from: ACAT.R from https://github.com/yaowuliu/ACAT/blob/master/R/ACAT.R

import numpy as np
import pandas as pd
from scipy.stats import cauchy

def acat(pvals, weights=None):
    """Performs the aggregated Cauchy association test."""
    p = pd.Series(pvals)
    if weights is not None:
        w = pd.Series(weights)
    # Check for null values
    if p.isnull().sum() > 0:
        raise ValueError("Cannot have NaNs in the p-values.")
    # Check that p-values are between 0 and 1.
    if ((p < 0).sum() + (p > 1).sum()) > 0:
        raise ValueError("P-values must be between 0 and 1.")
    # Check if any p-values are either 0 or 1 exactly.
    if (p == 0).sum() > 0:
        if (p == 1).sum() > 0:
            raise ValueError("Cannot have both 0 and 1 p-values.")
        else:
            return(0)
    elif (p == 1).sum() > 0:
        print("WARNING: There are p-values exactly equal to 1.")
        return(1)
    # Default: equal weights.
    # If weights have been supplied, check if valid and standardise.
    if weights is None:
        w = pd.Series([(1 / len(p))] * len(p))
    elif len(w) != len(p):
        raise ValueError("The length of weights should be the same as that of the p-values.")
    elif (w < 0).sum() > 0:
        raise ValueError("All weights must be positive.")
    else:
        w = w / w.sum()
    # Calculate test statistic. Also check for very small non-zero p-values.
    isSmall = p < 1e-16
    isNotSmall = (isSmall == False)
    if isSmall.sum() == 0:
        acatStat = (w * np.tan((0.5 - p) * np.pi)).sum()
    else:
        acatStat = ((w[isSmall] / p[isSmall]) / np.pi).sum()
        acatStat = acatStat + (w[isNotSmall] * np.tan((0.5 - p[isNotSmall]) * np.pi)).sum()
    # Check if statistic is very large
    if acatStat > 1e15:
        acatP = (1 / acatStat) / np.pi
    else:
        acatP = 1 - cauchy.cdf(acatStat)
    return(acatP)
