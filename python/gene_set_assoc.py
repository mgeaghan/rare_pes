#!/usr/bin/env python
# coding: utf-8
#
# Author: Michael Geaghan
#
# Date: 2020-04-23

# TODO:
# Output summary data in a more concise table

import statsmodels.api as sm
import pandas as pd
import numpy as np
from scipy.stats import norm
import argparse
import sys
import csv
from linregression import *
from logregression import *
from rankNorm import rankNorm
#import seaborn as sns
#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D  # for 3D plots

model_list = ["linear", "logistic"]

def check_args(args=None):
    """Parse arguments"""
    parser = argparse.ArgumentParser(description="Perform gene set association analysis using genic p-values.", formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-i', '--input', help="Comma-delimited (CSV) input file containing genic summary statistics. Must contain at least three columns: gene name, minor allele count and genic p-value. Must contain a header.", required=True)
    parser.add_argument('-g', '--gene_sets', help="GMT-formatted gene set file.", required=True)
    parser.add_argument('-l', '--gene_lengths', help="Tab-delimited file containing gene length information. Must contain two columns - gene name and gene length. Must be headerless.", required=True)
    parser.add_argument('-n', '--gene_name', help="Header name for gene name column. Default: gene_name.", default="gene_name")
    parser.add_argument('-s', '--stat', help="Header name for test statistic column.", default="")
    parser.add_argument('-p', '--pvalue', help="Header name for p-value column. Overrides -s/--stat", default="")
    parser.add_argument('-z', '--zvalue', help="Header name for pre-computed z-value column. Overrides -s/--stat and -p/--pvalue", default="")
    parser.add_argument('-c', '--mac', help="Header name for minor allele count column.", default="")
    parser.add_argument('-a', '--mac_con', help="Header name for minor allele count in controls column. Used if -c/--mac not found. Must be used in conjunction with -A/--mac_cas.", default="")
    parser.add_argument('-A', '--mac_cas', help="Header name for minor allele count in cases column. Used if -c/--mac not found. Must be used in conjunction with -a/--mac_con.", default="")
    parser.add_argument('-m', '--model', help=("Model to use for gene set association. Possible options:\n\t"
                                              + ", ".join(model_list)
                                              + "\n\tDefault: " + model_list[0]), default=model_list[0])
    parser.add_argument('-N', '--min_num_genes', help="Minimum number of genes for a gene set to be considered valid. Default: 11.", default=11)
    parser.add_argument('-t', '--alpha', help="Alpha for multiple test correction. Default = 0.05.", default=0.05)
    parser.add_argument('-C', '--correction_method', help="Method for multiple test correction. For no correction, use 'none'. Default = 'bonferroni'.", default="bonferroni")
    parser.add_argument('-o', '--output', help="Prefix for output files. Default = 'output'.", default="output")
    return(parser.parse_args(args))

def main(args):
    """Main function"""
    # Import p-values
    df = pd.read_csv(args.input, index_col=args.gene_name)
    # drop unnecessary columns and duplicate rows
    colNames = [args.mac, args.mac_con, args.mac_cas, args.zvalue, args.pvalue, args.stat]
    colNames = [i for i in colNames if not i == ""]
    colSelect = [i for i in colNames if i in df.columns]
    df_filt = df[colSelect]
    df_filt = df_filt.loc[~df_filt.index.duplicated(keep="first")]
    # re-name MAC column(s) and remove null/non-finite rows
    if (not args.mac == "") and (args.mac in df_filt.columns):
        df_filt = df_filt[df_filt[args.mac].notnull()]
        df_filt.rename({args.mac: "mac"}, axis=1, inplace=True)
    elif ((not args.mac_con == "") and (args.mac_con in df_filt.columns)) and ((not args.mac_cas == "") and (args.mac_cas in df_filt.columns)):
        df_filt = df_filt[np.isfinite(df_filt[[args.mac_con, args.mac_cas]]).any(1)]
        df_filt.rename({args.mac_con: "mac_con", args.mac_cas: "mac_cas"}, axis=1, inplace=True)
    else:
        raise ValueError("Minor allele column(s) missing.")
    # re-name stat/p/z column
    if (not args.zvalue == "") and (args.zvalue in df_filt.columns):
        df_filt.rename({args.zvalue: "z"}, axis=1, inplace=True)
        # filter out genes with null z values
        df_filt = df_filt[df_filt["z"].notnull()]
    elif (not args.pvalue == "") and (args.pvalue in df_filt.columns):
        df_filt.rename({args.pvalue: "pval"}, axis=1, inplace=True)
        # filter out genes with null p values
        df_filt = df_filt[df_filt["pval"].notnull()]
        # calculate z-scores
        # convert 0 and 1 p-values to the machine epsilon (eps) and 1 - eps
        # this appears to be the approach used in MAGMA (v1.07).
        p_values = np.copy(df_filt["pval"])
        min_p = np.finfo(type(p_values[0])).eps
        max_p = 1 - min_p
        p_values[p_values > max_p] = max_p
        p_values[p_values < min_p] = min_p
        df_filt = df_filt.assign(pval=p_values)
        z_scores = norm.ppf(1 - p_values)
        df_filt = df_filt.assign(z=z_scores)
    elif (not args.stat == "") and (args.stat in df_filt.columns):
        # filter out genes with null stat values
        df_filt = df_filt[df_filt[args.stat].notnull()]
        # calculate z-scores using rank-based inverse normalisation (Blom method)
        stat_scores = np.copy(df_filt[args.stat])
        z_scores = rankNorm(stat_scores, method="blom")
        df_filt = df_filt.assign(z=z_scores)
        df_filt.rename({args.stat: ("test_stat_" + args.stat)}, axis=1, inplace=True)
    else:
        raise ValueError("No test statistic/p-value/z-value column provided.")
    df_filt.index.name = "gene_name"
    # import gene lengths
    gene_lengths = pd.read_csv(args.gene_lengths, index_col=0, sep='\t')
    gene_lengths.index.name = "gene_name"
    if len(gene_lengths.columns) == 1:
        gene_lengths.columns = ["len"]
    else:
        raise ValueError("Incorrectly-formatted gene lengths file.")
    # remove genes for which len = 0
    gene_lengths = gene_lengths[gene_lengths["len"] != 0]
    # get log_len
    gene_lengths = gene_lengths.assign(log_len = np.log(gene_lengths["len"]))
    # merge gene lengths with summary statistics
    df_filt_len = df_filt.merge(gene_lengths, left_index=True, right_index=True)
    # add minor allele count (case_lof + control_lof) and log minor allele count
    if 'mac' not in df_filt_len.columns:
        df_filt_len = df_filt_len.assign(mac=(df_filt_len['mac_con'] + df_filt_len['mac_cas']))
        # remove genes for which MAC = 0
        df_filt_len = df_filt_len[df_filt_len["mac"] != 0]
    # get log_mac
    df_filt_len = df_filt_len.assign(log_mac=(np.log(df_filt_len['mac'])))
    # import gene sets
    gene_sets = {}
    with open(args.gene_sets, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        for line in reader:
            set = line[0]
            gene_list = line[2:]
            gene_sets[set] = gene_list
    # convert gene sets to a matrix (rows = genes, columns = sets)
    n_genes = int(args.min_num_genes)
    gene_sets_names = list(gene_sets.keys())
    gene_sets_names_final = []
    df_filt_len_set = df_filt_len.copy()
    for set in gene_sets_names:
        in_set = df_filt_len_set.index.isin(gene_sets[set]).astype(int)
        # only include sets with at least n_genes in the current data
        if np.sum(in_set) >= n_genes:
            kwargs = {set: in_set}
            df_filt_len_set = df_filt_len_set.assign(**kwargs)
            gene_sets_names_final.append(set)
    # construct the model for gene set association
    if args.model == "linear":
        regression = geneSetLinearRegression()
    elif args.model == "logistic":
        regression = geneSetLogisticRegression()
    else:
        raise ValueError("Invalid model selection.")
    regression.reg(df_filt_len_set, gene_sets_names_final, ["log_len", "log_mac"])
    regression.sig_sets(args.correction_method, args.alpha)
    # Export summary data
    summary_file = args.output + ".summary.txt"
    with open(summary_file, 'w') as f:
        for set in regression.gene_sets:
            s = regression.gene_sets[set]
            f.write(set + '\n\n' + str(s.summary()) + '\n\n' + "==========" + '\n\n')
    # Export list of significant sets
    sig_sets_list_file = args.output + ".sig.sets.txt"
    with open(sig_sets_list_file, 'w') as f:
        f.write('\n'.join(regression.gene_sets_significant[0]))
    # Export main data frame
    df_file = args.output + ".genes.txt"
    df_filt_len_set.to_csv(df_file)

if __name__ == "__main__":
    arguments = check_args(sys.argv[1:])
    arguments.alpha = float(arguments.alpha)
    arguments.min_num_genes = int(arguments.min_num_genes)
    main(arguments)
