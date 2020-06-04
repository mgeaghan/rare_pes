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
from pandas_plink import read_plink1_bin
import seaborn as sns
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # for 3D plots
import argparse
import sys
import csv


def import_annot(variant_info, ref_allele, annot_file, chr_col, pos_col, ref_col, alt_col, af_col, cadd_col, na_symbol, p_alpha, p_beta):
    variant_annotations_cols = ["var", "chr", "pos", "ref", "alt", "af", "score", "vf", "vs", "weight"]
    variant_annotations_data = {c: [] for c in variant_annotations_cols}
    var_annot = pd.read_table(annot_file, sep='\t', na_values=na_symbol)
    var_annot_cols = {chr_col: "chr", pos_col: "pos", ref_col: "ref", alt_col: "alt", af_col: "af", cadd_col: "score"}
    cols = [c for c in var_annot_cols]
    var_annot = var_annot[cols]
    # var_annot.dropna(inplace=True)
    var_annot.drop_duplicates(inplace=True)
    var_annot.rename(columns={c: var_annot_cols[c] for c in cols}, inplace=True)
    vc = var_annot["chr"].astype(str)
    vp = var_annot["pos"].astype(int)
    vr = var_annot["ref"]
    va = var_annot["alt"]
    for v in variant_info:
        v_i = variant_info[v]
        v_c = str(v_i[0])
        v_p = int(v_i[1])
        v_a0 = v_i[2]
        v_a1 = v_i[3]
        if ref_allele == 0:
            ref = v_a0
            alt = v_a1
            annot = var_annot[(vc == v_c) & (vp == v_p) & (vr == v_a0) & (va == v_a1)]
        elif ref_allele == 1:
            ref = v_a1
            alt = v_a0
            annot = var_annot[(vc == v_c) & (vp == v_p) & (vr == v_a1) & (va == v_a0)]
        else:
            raise ValueError("ref_allele argument must be either '0' or '1'")
        # MAF
        af = annot["af"].dropna()
        if len(af) > 1:
            unique_af = af.unique()
            if len(unique_af) > 1:
                af = max(unique_af)  # This condition seems unlikely, but if it happens, be conservative and take the larger AF value
                vf = beta.pdf(af, p_alpha, p_beta)/p_beta
            else:
                af = unique_af[0]
                vf = beta.pdf(af, p_alpha, p_beta)/p_beta
        elif len(af) == 0:
            af = 0  # a missing allele frequency should mean it has a very low frequency - so set to 0
            vf = beta.pdf(af, p_alpha, p_beta)/p_beta
        else:
            af = af.values[0]
            vf = beta.pdf(af, p_alpha, p_beta)/p_beta
        # Calculate MAF from AF
        if af > 0.5:
            # if ALT is the major allele, then REF will be the minor allele and there will be no CADD score for the minor allele. Therefore, discard this variant.
            continue
        # CADD
        cadd = annot["score"].dropna()
        if len(cadd) > 1:
            unique_cadd = cadd.unique()
            if len(unique_cadd) > 1:
                cadd = min(unique_cadd)  # This condition seems unlikely, but if it happens, be conservative and take the smaller CADD value
                vs = 1 - (10**(cadd/(-10)))
            else:
                cadd = unique_cadd[0]
                vs = 1 - (10**(cadd/(-10)))
        elif len(cadd) == 0:
            continue  # missing CADD can't be used for scoring - discard variant
        else:
            cadd = cadd.values[0]
            vs = 1 - (10**(cadd/(-10)))
        weight = (vs * vf)
        variant_annotations_data["var"].append(v)
        variant_annotations_data["chr"].append(v_c)
        variant_annotations_data["pos"].append(v_p)
        variant_annotations_data["ref"].append(ref)
        variant_annotations_data["alt"].append(alt)
        variant_annotations_data["af"].append(af)
        variant_annotations_data["score"].append(cadd)
        variant_annotations_data["vf"].append(vf)
        variant_annotations_data["vs"].append(vs)
        variant_annotations_data["weight"].append(weight)
    return(pd.DataFrame(variant_annotations_data, index=variant_annotations_data["var"]))
