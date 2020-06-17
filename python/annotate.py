#!/usr/bin/env python
# coding: utf-8
#
# Author: Michael Geaghan
#
# Date: 2020-05-19

import statsmodels.api as sm
import pandas as pd
import numpy as np
from scipy.stats import beta
import argparse
import sys
import csv
#from pandas_plink import read_plink1_bin
#import seaborn as sns
#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D  # for 3D plots


def import_annot(variant_info, annot_file, chr_col, pos_col, ref_col, alt_col, af_col, cadd_col, na_symbol, p_alpha, p_beta):
    variant_annotations_cols = ["var", "chr", "pos", "ref", "alt", "af", "score", "vf", "vs", "weight"]
    variant_annotations_data = {c: [] for c in variant_annotations_cols}
    print("Importing annoation file...")
    var_annot = pd.read_table(annot_file, sep='\t', na_values=na_symbol)
    var_annot_cols = {chr_col: "chr", pos_col: "pos", ref_col: "ref", alt_col: "alt", af_col: "af", cadd_col: "score"}
    cols = [c for c in var_annot_cols]
    var_annot = var_annot[cols]
    print("Renaming columns...")
    var_annot.rename(columns={c: var_annot_cols[c] for c in cols}, inplace=True)
    print("Renaming indices...")
    vc = var_annot["chr"].astype(str)
    vp = var_annot["pos"].astype(str)
    vr = var_annot["ref"].astype(str)
    va = var_annot["alt"].astype(str)
    new_index = vc + ":" + vp + ":" + vr + ":" + va
    var_annot.rename(index=new_index, inplace=True)
    # Drop variants not in variant_info
    print("Dropping variants not in plink data...")
    var_annot = var_annot.loc[[v for v in variant_info.values() if v in var_annot.index]]
    # Add variant name to dataframe
    print("Adding variant names...")
    variant_info_2 = {variant_info[v]: v for v in variant_info}
    var_annot = var_annot.assign(var=var_annot.apply(axis=1, func=lambda x: variant_info_2[x.name]))
    # Replace NaNs in af column with 0
    print("Replacing AF NaNs with 0...")
    var_annot.fillna(value={"af": 0}, inplace=True)
    # Drop variants where af is > 0.5; they will have no CADD score for the minor allele because REF is the minor allele
    print("Dropping variants with AF > 0.5 (no CADD available)...")
    var_annot = var_annot[var_annot.af <= 0.5]
    # Drop variants with NaN for score column
    print("Dropping variants with missing CADD...")
    var_annot.dropna(axis=0, subset=["score"], inplace=True)
    # Sort by variant name, then by decreasing AF, then by increasing CADD; drop duplicate variants, keeping the first
    print("Dropping duplicate annotations...")
    var_annot.sort_values(axis=0, by=["chr", "pos", "ref", "alt", "af", "score"], ascending=[True, True, True, True, False, True], inplace=True)
    var_annot.drop_duplicates(subset=["chr", "pos", "ref", "alt"], keep="first", inplace=True)
    # Create vf, vs and weight columns
    print("Calculating vf, vs and weights...")
    var_annot = var_annot.assign(vf=lambda x: beta.pdf(x.af, p_alpha, p_beta)/p_beta)
    var_annot = var_annot.assign(vs=lambda x: 1 - (10**(x.score/(-10))))
    var_annot = var_annot.assign(weight=lambda x: x.vf * x.vs)
    return(var_annot)
