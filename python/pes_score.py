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
from annotate import *

def check_args(args=None):
    """Parse arguments"""
    parser = argparse.ArgumentParser(description="Calculate rare variant PES scores for individuals in multiple pathways.", formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-g', '--genes', help="Gene/pathway data file. Must be a CSV file with a column for gene IDs and a column for every pathway to be scored (0/1 encoded).", required=True)
    parser.add_argument('-c', '--genes_column', help="Column name for gene IDs in --genes file. Default = 'gene_name'.", default="gene_name")
    parser.add_argument('-l', '--loc', help="Gene location data file. Must be a headerless, tab-delimited file, with four columns: gene ID, chromosome, start position and end position (1-based, inclusive). Chromosome codes must match PLINK .bim file for X/Y/MT.", required=True)
    parser.add_argument('-p', '--pathways', help="File containing pathways to be scored. One pathway ID per line. Each ID must match a column header present in the --genes file.", required=True)
    parser.add_argument('-B', '--bfile', help="PLINK file prefix.", required=True)
    parser.add_argument('-r', '--ref_allele', help="Reference allele in PLINK data. Must be either '0' or '1', i.e. PLINK data must be calculated relative to the reference genome.", required=True)
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

def get_gene_vars(gene_loc_df, plink_data):
    """Takes a gene location dataframe and PLINK data. Returns a tuple of two dictionaries: {gene: [variant_list]} and {variant: [gene_list]}"""
    # Organise data, re-name chr codes, sort by chr and pos/start
    plink_df = plink_data.variant.to_dataframe()
    new_chr = plink_df['chrom'].str.replace('X', '23').str.replace('Y', '24').str.replace('XY', '25').str.replace('MT', '26').str.replace('M', '26')
    plink_df = plink_df.assign(chrom=new_chr)
    plink_df = plink_df[["snp", "chrom", "pos"]].astype({"chrom": "int32", "pos": "int32"})
    plink_df.sort_values(by=['chrom', 'pos'], inplace=True)
    genes_df = gene_loc_df.copy()
    new_chr = genes_df['chr'].str.replace('X', '23').str.replace('Y', '24').str.replace('XY', '25').str.replace('MT', '26').str.replace('M', '26')
    genes_df = genes_df.assign(chr=new_chr)
    genes_df = genes_df.astype({"chr": "int32", "start": "int32", "end": "int32"})
    genes_df.sort_values(by=['chr', 'start', 'end'], inplace=True)
    genes_df.reset_index(drop=True, inplace=True)
    # Iterate over rows of both dataframes
    variant_genes = {}
    gene_variants = {}
    plink_iter = plink_df.iterrows()
    genes_iter = genes_df.iterrows()
    plink_idx, plink_val = next(plink_iter)
    genes_idx, genes_val = next(genes_iter)
    i = 0
    while True:
        try:
            var_chr = plink_val[1]
            gen_chr = genes_val[1]
            if var_chr < gen_chr:
                plink_idx, plink_val = next(plink_iter)  # next variant
            elif var_chr > gen_chr:
                genes_idx, genes_val = next(genes_iter)  # next gene
            else:
                var_pos = plink_val[2]
                gen_start = genes_val[2]
                gen_end = genes_val[3]
                if var_pos < gen_start:
                    plink_idx, plink_val = next(plink_iter)  # next variant
                elif var_pos > gen_end:
                    genes_idx, genes_val = next(genes_iter)  # next gene
                else:
                    var_name = plink_val[0]
                    genes_iter_2 = genes_df.iloc[genes_idx:].iterrows()
                    for genes_2_idx, genes_2_val in genes_iter_2:
                        gen_2_start = genes_2_val[2]
                        gen_2_end = genes_2_val[3]
                        if var_pos >= gen_2_start and var_pos <= gen_2_end:
                            gen_2_name = genes_2_val[0]
                            try:
                                if gen_2_name not in variant_genes[var_name]:
                                    variant_genes[var_name].append(gen_2_name)
                            except KeyError:
                                variant_genes[var_name] = [gen_2_name]
                            try:
                                if var_name not in gene_variants[gen_2_name]:
                                    gene_variants[gen_2_name].append(var_name)
                            except KeyError:
                                gene_variants[gen_2_name] = [var_name]
                        else:
                            plink_idx, plink_val = next(plink_iter)  # next variant
                            break
        except StopIteration:
            break
    return((gene_variants, variant_genes))

def get_var_chr_pos_a0_a1(snp, plink_data):
    """Returns variant chromosome, position, a0 and a1 alleles."""
    var = plink_data.sel(variant=list(plink_data.variant.values[plink_data.snp.values == snp]))
    chr = var.chrom.values[0]
    pos = var.pos.values[0]
    a0 = var.a0.values[0]
    a1 = var.a1.values[0]
    return((chr, pos, a0, a1))

def get_cohort_af(snp, plink_data):
    """Calculate SNP allele frequency based on PLINK data."""
    var = plink_data.sel(variant=list(plink_data.variant.values[plink_data.snp.values == snp]))
    return(np.sum(var.values)/(len(plink_data.sample.values)*2))

def get_a1_dosage(snp, sample, plink_data):
    """Returns the dosage of the a1 allele for the given SNP in the given individual for the given PLINK data set"""
    dosage = int(plink_data.sel(variant=list(plink_data.variant.values[plink_data.snp.values==snp]), sample=sample).values)/2
    return(dosage)

def get_a0_dosage(snp, sample, plink_data):
    """Returns the dosage of the a0 allele for the given SNP in the given individual for the given PLINK data set"""
    return(1 - get_a1_dosage(snp, sample, plink_data))

def get_set_variant_info(gene_info, gene_set_info, variant_list, gene_variant_info, gene_weights):
    gene_sets_variants_info = {s: {} for s in gene_set_info}
    for s in gene_set_info:
        seen_variants = []
        for g in gene_set_info[s]:
            if g in gene_variant_info:
                for v in gene_variant_info[g]:
                    if v in variant_list:
                        if v not in seen_variants:
                            # get gene z value
                            z = gene_info.loc[g, 'z']
                            # get gene weight
                            w = gene_weights[g]
                            gene_sets_variants_info[s][v] = (g, z, w)
                            seen_variants.append(v)
                        else:
                            # check if new z value is higher or lower than the previous z value - update info if higher
                            z_prev = gene_sets_variants_info[s][v][1]
                            z_new = gene_info.loc[g, 'z']
                            if z_new > z_prev:
                                w_new = gene_weights[g]
                                gene_sets_variants_info[s][v] = (g, z_new, w_new)
                            elif z_new == z_prev:
                                w_prev = gene_sets_variants_info[s][v][2]
                                w_new = gene_weights[g]
                                if w_new > w_prev:
                                    gene_sets_variants_info[s][v] = (g, z_new, w_new)
                                else:
                                    continue
                            else:
                                continue
                    else:
                        continue
            else:
                continue
    gene_sets_variants_info_short = {s: gene_sets_variants_info[s] for s in gene_sets_variants_info if not gene_sets_variants_info[s] == {}}
    return(gene_sets_variants_info_short)

def pes_df_init(sample_list, gene_set_list):
    len_samples = len(sample_list)
    df = pd.DataFrame({s: [None]*len_samples for s in gene_set_list}, index=sample_list)
    return(df)

def score_pes(sample_list, gene_set_list, gene_sets_variants_info, var_annotations, ref_allele, plink_data):
    if ref_allele == 1:
        dosage_df = 1 - pd.DataFrame(plink_data.sel().values, index=plink_data.sample.values, columns=plink_data.snp.values)/2  # dosage of a0 ALT allele
    else:
        dosage_df = pd.DataFrame(plink_data.sel().values, index=plink_data.sample.values, columns=plink_data.snp.values)/2  # dosage of a1 ALT allele
    pes_df = pes_df_init(sample_list, gene_set_list)
    for i in sample_list:
        for s in gene_set_list:
            pes_n = 0
            pes_d = 0
            for v in gene_sets_variants_info[s]:
                vd = dosage_df.loc[i, v]
                vsvf = var_annotations.loc[v, 'weight']  # score weight * frequency weight
                # vs = var_annotations.loc[v, 'vs']  # score weight
                # vf = var_annotations.loc[v, 'vf']  # frequency weight
                wg = gene_sets_variants_info[s][v][2]
                zg = gene_sets_variants_info[s][v][1]
                # weight = (vd * vs * vf * wg)
                weight = (vd * vsvf * wg)
                pes_n += (weight * zg)
                pes_d += (weight ** 2)
            if not pes_d == 0:
                pes = pes_n / np.sqrt(pes_d)
            else:
                pes = None
            pes_df.loc[i, s] = pes
    return(pes_df)

def main(args):
    """Main function"""
    # Import PLINK data
    bed = args.bfile + ".bed"
    bim = args.bfile + ".bim"
    fam = args.bfile + ".fam"
    G = read_plink1_bin(bed, bim, fam, verbose = False)
    # Import gene location information
    gene_loc = pd.read_table(args.loc, sep='\t', header=None)
    gene_loc.columns = ["gene", "chr", "start", "end"]
    # Create dictionaries of {gene: [variant_list]} and {variant: [gene_list]}
    gene_variants, variant_genes = get_gene_vars(gene_loc, G)
    # Create a dictionary of variant information
    variant_info = {v: get_var_chr_pos_a0_a1(v, G) for v in variant_genes}
    # Import gene/pathway data
    genes = pd.read_csv(args.genes, index_col=args.genes_column)
    genes.drop_duplicates(inplace=True)
    # Import list of gene sets to score
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
    # Gather variant annotation information
        # Import variant annotations
    variant_annotations = import_annot(variant_info, args.ref_allele, args.annot, args.chr, args.pos, args.ref, args.alt, args.af, args.cadd, args.na, args.alpha, args.beta)
    # Update gene_variants, variant_genes and variant_info to remove un-annotated variants
    variant_genes = {v: variant_genes[v] for v in variant_genes if v in variant_annotations.index}
    variant_info = {v: variant_info[v] for v in variant_info if v in variant_annotations.index}
    gene_variants = {g: [v for v in gene_variants[g] if v in variant_annotations.index] for g in gene_variants}
    gene_variants = {g: gene_variants[g] for g in gene_variants if not gene_variants[g] == []}
    # get min and max variants per gene and calculate the relative fraction of variants per gene for the cohort
    variants_per_gene = {g: len(gene_variants[g]) for g in gene_variants}
    max_variants_per_gene = max(list(variants_per_gene.values()))
    min_variants_per_gene = min(list(variants_per_gene.values()))
    rel_variants_per_gene_weights = {g: 1 - ((variants_per_gene[g] - min_variants_per_gene)/(1 + (max_variants_per_gene - min_variants_per_gene))) for g in variants_per_gene}
    # get the unique gene for each variant in each pathway. If a variant overlaps two genes in a single pathway, take the gene with the higher Z score.
    variant_list = list(variant_annotations.index)
    gene_sets_variants_info = get_set_variant_info(genes, gene_sets, variant_list, gene_variants, rel_variants_per_gene_weights)
    # calculate PES for each pathway for each individual
    sample_list = list(G.sample.values)
    gene_set_list = list(gene_sets_variants_info.keys())
    pes_df = score_pes(sample_list, gene_set_list, gene_sets_variants_info, variant_annotations, args.ref_allele, G)
    # write to file
    output_file = args.output + ".pes.txt"
    pes_df.to_csv(output_file)

if __name__ == "__main__":
    arguments = check_args(sys.argv[1:])
    # Check ref allele code (0, 1 or None)
    if arguments.ref_allele not in ['0', '1', 0, 1]:
        raise ValueError("Argument --ref_allele must be either '0', or '1'")
    else:
        arguments.ref_allele = int(arguments.ref_allele)
    # Ensure alpha and beta are floats
    arguments.alpha = float(arguments.alpha)
    arguments.beta = float(arguments.beta)
    main(arguments)
