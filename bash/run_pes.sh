#!/bin/bash
#PBS -l select=1:ncpus=1:mem=32GB
#PBS -l walltime=24:00:00
cd $PBS_O_WORKDIR

eval "$(/home/mpg993/anaconda3/anaconda/bin/conda shell.bash hook)"
conda activate rpes

cd /home/mpg993/rpes/analysis/out/pes

python3 ../../../rare_pes/python/pes_score.py -g ../DISORDER/DISORDER.SUBSET.gsa.GENESETS.genes.txt -c GENEHEADER -l ../../db/GENERANGES -p ../DISORDER/DISORDER.SUBSET.gsa.GENESETS.SIGSETS.txt -B ../../db/plink/wgs.eur.rare_MAF -r 2 -v ../../db/annot/wgs.eur.snv.indel.g_indel.annot.rare_MAF.txt -f AF_nfe_exome -s PHRED -C Chr -P Start -R Ref -A Alt -n . -o ./DISORDER.SUBSET.GENESETS.SIGSETS.MAF.pes

