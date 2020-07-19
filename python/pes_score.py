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
from annotate import *
import subprocess
import re
#from pandas_plink import read_plink1_bin
#import seaborn as sns
#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D  # for 3D plots

def check_args(args=None):
    """Parse arguments"""
    parser = argparse.ArgumentParser(description="Calculate rare variant PES scores for individuals in multiple pathways.", formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-g', '--genes', help="Gene/pathway data file. Must be a CSV file with a column for gene IDs and a column for every pathway to be scored (0/1 encoded).", required=True)
    parser.add_argument('-c', '--genes_column', help="Column name for gene IDs in --genes file. Default = 'gene_name'.", default="gene_name")
    parser.add_argument('-l', '--loc', help="Gene location data file. Must be a headerless, tab-delimited file, with four columns: gene ID, chromosome, start position and end position (1-based, inclusive). Chromosome codes must match PLINK .bim file for X/Y/MT.", required=True)
    parser.add_argument('-p', '--pathways', help="File containing pathways to be scored. One pathway ID per line. Each ID must match a column header present in the --genes file.", required=True)
    parser.add_argument('-B', '--bfile', help="PLINK file prefix.", required=True)
    parser.add_argument('-r', '--ref_allele', help="Reference allele in PLINK data. Must be either '1' or '2', corresponding to the A1 or A2 allele in the PLINK .bim file, i.e. PLINK data must be calculated relative to the reference genome.", required=True)
    parser.add_argument('-v', '--annot', help="Variant annotation file. Must contain a header line and at least four columns: chromosome, position, REF allele and ALT allele. Chromosome codes must match PLINK .bim file for X/Y/MT.", required=True)
    parser.add_argument('-f', '--af', help="Allele frequency column header for variant annotation file.", required=True)
    parser.add_argument('-s', '--cadd', help="CADD score column header for variant annotation file.", required=True)
    parser.add_argument('-C', '--chr', help="Chromosome column header for variant annotation file.", required=True)
    parser.add_argument('-P', '--pos', help="Position column header for variant annotation file.", required=True)
    parser.add_argument('-R', '--ref', help="REF allele column header for variant annotation file.", required=True)
    parser.add_argument('-A', '--alt', help="ALT allele column header for variant annotation file.", required=True)
    parser.add_argument('-n', '--na', help="NA values for annotation file. Default = '.'.", default=".")
    parser.add_argument('-a', '--alpha', help="Alpha parameter for beta density function. Default = 1.", default=1)
    parser.add_argument('-b', '--beta', help="Beta parameter for beta density function. Default = 25.", default=25)
    parser.add_argument('-o', '--output', help="Prefix for output files. Default = 'output'.", default="output")
    return(parser.parse_args(args))

def parse_sets(set_file):
    """Parse a .set file from PLINK."""
    sets = {}
    with open(set_file, 'r') as f:
        lineIsSetName = True
        setName = ""
        setVars = []
        for line in f:
            line = line.strip('\n')
            if line == "":
                continue
            elif lineIsSetName:
                setName = line
                lineIsSetName = False
            elif line == "END":
                if setVars != []:
                    try:
                        sets[setName].append(setVars)
                    except KeyError:
                        sets[setName] = setVars
                    sets[setName] = list(set(sets[setName]))
                setName = ""
                setVars = []
                lineIsSetName = True
            else:
                setVars.append(line)
    return(sets)

def plink_sets(bfile, locfile, newbfile):
    """Takes a plink binary file prefix and a gene location file (headerless, four columns: CHR, START, END, GENE_NAME) and creates a new plink binary file set subset by the gene list and a plink .set file. Function returns the .set file name."""
    subprocess.run(["plink", "--bfile", bfile, "--keep-allele-order", "--make-set", locfile, "--gene-all", "--write-set", "--make-bed", "--out", newbfile])
    setFile = newbfile + ".set"
    return(setFile)

def get_var_info(bfile, ref_allele):
    """Takes a plink binary file prefix and a ref_allele argument (either 1 or 2 corresponding to either A1 or A2 being the REF allele) and returns a dictionary of information for each variant ({variant: [CHR, POS, REF, ALT]})."""
    bimFile = bfile + ".bim"
    var_info = {}
    with open(bimFile, 'r') as f:
        reader = csv.reader(f, delimiter="\t")
        for line in reader:
            var_name = str(line[1])
            chr = str(line[0])
            pos = str(line[3])
            if ref_allele == 1:
                ref = str(line[4])
                alt = str(line[5])
            elif ref_allele == 2:
                ref = str(line[5])
                alt = str(line[4])
            else:
                raise ValueError("get_var_info: ref_allele argument must be an integer 1 or 2.")
            var_info[var_name] = [chr, pos, ref, alt]
    return(var_info)

def get_var_info_concat(bfile, ref_allele):
    """Takes a plink binary file prefix and a ref_allele argument (either 1 or 2 corresponding to either A1 or A2 being the REF allele) and returns a dictionary of information for each variant ({variant: CHR:POS:REF:ALT})."""
    var_info = get_var_info(bfile, ref_allele)
    var_info_concat = {v: ":".join(var_info[v]) for v in var_info}
    return(var_info_concat)

def get_dosages(bfile, ref_allele):
    """Get the dosages for each variant and each individual."""
    if ref_allele == 1:
        bimFile = bfile + ".bim"
        newbfile = bfile + ".swap_alleles"
        print("Swapping A1 and A2 alleles in PLINK...")
        subprocess.run(["plink", "--bfile", bfile, "--a2-allele", bimFile, 5, 2, "--make-bed", "--out", newbfile])
        plinkFile = newbfile
    else:
        plinkFile = bfile
    print("Getting dosages from PLINK...")
    newbfile = plinkFile + ".dosages"
    subprocess.run(["plink", "--bfile", plinkFile, "--keep-allele-order", "--recode", "A", "tab", "--out", newbfile])
    dosageFile = newbfile + ".raw"
    dosages_df = pd.read_table(dosageFile, sep="\t", na_values="NA", index_col="IID")
    cols = [c for c in list(dosages_df.columns) if c not in ["FID", "PAT", "MAT", "SEX", "PHENOTYPE"]]
    dosages_df = dosages_df[cols]/2
    dosages_df.fillna(value=0, inplace=True)  # missing values should be set to 0
    regex = re.compile("^(.*)_[^_]+$")
    cols = {c: regex.match(c).group(1) for c in cols}
    dosages_df.rename(columns=cols, inplace=True)
    return(dosages_df)

def weight_dosages(dosages_df, annotations, weights):
    """Weights dosages dataframe by the weights column in the annotations dataframe."""
    d_cols = list(dosages_df.columns)
    a_idx = list(annotations.index)
    keep_cols = [c for c in d_cols if c in a_idx]
    return(dosages_df[keep_cols].apply(axis=0, func=lambda x: x * annotations.loc[x.name, weights]))

def get_set_variant_info(variant_list, gene_sets, gene_info, gene_variants, gene_weights):
    """Returns a tuple of two dataframes of variants and gene sets; the first with the respective gene weight * Z-score for each pair; the second with just the gene weights for each pair."""
    set_var_z = pd.DataFrame(index=variant_list, columns=list(gene_sets.keys())).fillna(0)
    set_var_w = pd.DataFrame(index=variant_list, columns=list(gene_sets.keys())).fillna(0)
    for s in gene_sets:
        seen_variants = []
        for g in gene_sets[s]:
            if g in gene_variants:
                for v in gene_variants[g]:
                    if v in variant_list:
                        if v not in seen_variants:
                            # get gene z value
                            z = gene_info.loc[g, 'z']
                            # get gene weight
                            w = gene_weights[g]
                            set_var_z.loc[v, s] = z
                            set_var_w.loc[v, s] = w
                            seen_variants.append(v)
                        else:
                            # check if new z value is higher or lower than the previous z value - update info if higher
                            z_prev = set_var_z.loc[v, s]
                            z_new = gene_info.loc[g, 'z']
                            if z_new > z_prev:
                                w_new = gene_weights[g]
                                set_var_z.loc[v, s] = z_new
                                set_var_w.loc[v, s] = w_new
                            elif z_new == z_prev:
                                w_prev = set_var_w.loc[v, s]
                                w_new = gene_weights[g]
                                if w_new > w_prev:
                                    set_var_z.loc[v, s] = z_new
                                    set_var_w.loc[v, s] = w_new
                                else:
                                    continue
                            else:
                                continue
                    else:
                        continue
            else:
                continue
    return((set_var_z * set_var_w, set_var_w))

def main(args):
    """Main function."""
    # Intersect variants and genes using PLINK
    print("Intersecting variants and genes in PLINK...")
    newBFile = args.bfile + ".pes_score.subset"
    setsFile = plink_sets(args.bfile, args.loc, newBFile)
    gene_variants = parse_sets(setsFile)
    # Import gene and pathway data file
    print("Importing gene and pathway data...")
    genes = pd.read_csv(args.genes, index_col=args.genes_column)
    genes.drop_duplicates(inplace=True)
    # Import list of gene sets to score
    print("Importing gene sets to score...")
    gene_sets_to_score = []
    with open(args.pathways, 'r') as f:
        reader = csv.reader(f)
        for line in reader:
            set_name = line[0]
            gene_sets_to_score.append(set_name)
    # Gather gene set information
    gene_sets = {s: [] for s in gene_sets_to_score}
    gene_names = list(genes.index)
    genes_columns = genes.columns
    for g in gene_names:
        for s in gene_sets_to_score:
            if (s in genes_columns) and (genes.loc[g, s] == 1):
                gene_sets[s].append(g)
    # Import variant information
    print("Importing variant information...")
    variant_info = get_var_info_concat(newBFile, args.ref_allele)
    # Import variant annotations and check for issues
    print("Importing variant annotations...")
    variant_annotations = import_annot(variant_info, args.annot, args.chr, args.pos, args.ref, args.alt, args.af, args.cadd, args.na, args.alpha, args.beta)
    # Remove un-annotated variants
    var_list = list(variant_annotations["var"])
    gene_variants = {g: [v for v in gene_variants[g] if v in var_list] for g in gene_variants}
    gene_variants = {g: gene_variants[g] for g in gene_variants if not gene_variants[g] == []}
    # Calculate weights for genes
    ## Get min and max variants per gene and calculate the relative fraction of variants per gene for the cohort
    print("Calculating gene weights...")
    variants_per_gene = {g: len(gene_variants[g]) for g in gene_variants}
    max_variants_per_gene = max(list(variants_per_gene.values()))
    min_variants_per_gene = min(list(variants_per_gene.values()))
    rel_variants_per_gene_weights = {g: 1 - ((variants_per_gene[g] - min_variants_per_gene)/(1 + (max_variants_per_gene - min_variants_per_gene))) for g in variants_per_gene}
    # For each pathway, get the variants and unique gene for each variant; if a variant overlaps multiple genes, use the higher Z-value, then the higher weight value to select the unique gene
    print("Gathering gene information for each variant in each pathway...")
    gene_sets_variants_info = get_set_variant_info(var_list, gene_sets, genes, gene_variants, rel_variants_per_gene_weights)
    # Calculate PES for each pathway for each individual
    print("Calculating PES...")
    dosages = weight_dosages(get_dosages(newBFile, args.ref_allele), variant_annotations.set_index("var"), "weight")
    pes = dosages.dot(gene_sets_variants_info[0])/np.sqrt((dosages**2).dot(gene_sets_variants_info[1]**2))
    # write to file
    print("Writing to output...")
    output_file = args.output + ".pes.txt"
    pes.to_csv(output_file)
    print("DONE!")

if __name__ == "__main__":
    arguments = check_args(sys.argv[1:])
    # Check ref allele code (0, 1 or None)
    if arguments.ref_allele not in ['1', '2', 1, 2]:
        raise ValueError("Argument --ref_allele must be either '1', or '2'")
    else:
        arguments.ref_allele = int(arguments.ref_allele)
    # Ensure alpha and beta are floats
    arguments.alpha = float(arguments.alpha)
    arguments.beta = float(arguments.beta)
    main(arguments)
