#!/usr/bin/env python
# coding: utf-8
#
# Linear regression class for gene set association
#
# Author: Michael Geaghan
# Date: 2020-04-23

import statsmodels.api as sm
import pandas as pd
import numpy as np
from scipy.stats import norm
import argparse
import sys
import csv
#import seaborn as sns
#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D  # for 3D plots

class geneSetLinearRegression:
    gene_sets = {}
    gene_sets_significant = None
    def reg(self, df, sets, covars):
        """Perform gene set association using a linear regression model."""
        for s in sets:
            Y = df['z']
            X = df[[s] + covars]
            X = sm.add_constant(X)  # explicitly specify a constant (intercept, beta_0)
            model = sm.OLS(Y, X).fit()
            self.gene_sets[s] = model

    def sig_sets(self, correction, alpha=0.05):
        """Return the significant gene sets."""
        supported_methods = ['bonferroni', 'fdr_bh']
        sets = list(self.gene_sets.keys())
        pvalues = [self.gene_sets[s].pvalues[s] for s in sets]
        if correction in supported_methods:
            mult_test = sm.stats.multipletests(pvalues, alpha=alpha, method=correction)[0:2]
            reject_null = mult_test[0]
            pvalues_corrected = mult_test[1]
            sig_sets = [sets[i] for i,j in enumerate(reject_null) if j]
            sig_results = [self.gene_sets[s] for s in sig_sets]
            self.gene_sets_significant = (sig_sets, sig_results, pvalues_corrected[reject_null])
        elif correction == 'none':
            reject_null = [p < alpha for p in pvalues]
            pvalues_alpha = [pvalues[i] for i,j in enumerate(reject_null) if j]
            sig_sets = [sets[i] for i,j in enumerate(reject_null) if j]
            sig_results = [self.gene_sets[s] for s in sig_sets]
            self.gene_sets_significant = (sig_sets, sig_results, pvalues_alpha)
        else:
            raise ValueError("Invalid correction method.")
