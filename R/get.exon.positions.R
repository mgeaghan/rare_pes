library(biomaRt)
options(scipen=999)
ensembl <- useMart('ENSEMBL_MART_ENSEMBL', 'hsapiens_gene_ensembl', host = 'http://grch37.ensembl.org')
# get gene and transcript IDs, trancript type, gene name, chromosome and genomic coding start and end
bm <- getBM(attributes = c('ensembl_gene_id', 'ensembl_transcript_id', 'ensembl_gene_id_version', 'ensembl_transcript_id_version', 'transcript_biotype', 'external_gene_name', 'chromosome_name', 'genomic_coding_start', 'genomic_coding_end'), mart = ensembl)
chrs <- c(seq(1, 22), 'X', 'Y', 'MT')
chrs
bm.filt <- bm[bm$chromosome_name %in% chrs & !is.na(bm$genomic_coding_start) & !is.na(bm$genomic_coding_end) & bm$transcript_biotype == 'protein_coding', ]
# Also remove two erroneous gene entries (gene names are duplicated and are on different chromosomes):
# CKS1B on chr5 (real gene is on chr1)
# LSP1 on chr13 (real gene is on chr11)
bm.filt <- bm.filt[!(bm.filt$external_gene_name %in% c("CKS1B", "LSP1") & bm.filt$chromosome_name %in% c("5", "13")),]
bm.filt.sort <- bm.filt[order(bm.filt$external_gene_name, bm.filt$genomic_coding_start), ]
head(bm.filt.sort)
# Merge overlapping exons
prevRow <- bm.filt.sort[1,]
bm.filt.sort.merge <- list()
listIdx = 1
for (i in 2:dim(bm.filt.sort)[1]) {
  row <- bm.filt.sort[i,]
  gene <- as.character(row[6])
  prevGene <- as.character(prevRow[6])
  chr <- as.character(row[7])
  prevChr <- as.character(prevRow[7])
  start <- as.numeric(row[8])
  prevStart <- as.numeric(prevRow[8])
  end <- as.numeric(row[9])
  prevEnd <- as.numeric(prevRow[9])
  if ((chr == prevChr) && (gene == prevGene) && (start >= prevStart) && (start <= prevEnd) && (end > prevEnd)) {
    prevRow$genomic_coding_end <- end
    next
  } else if ((chr == prevChr) && (gene == prevGene) && (start >= prevStart) && (start <= prevEnd) && (end <= prevEnd)) {
    next
  } else {
    bm.filt.sort.merge[[listIdx]] <- prevRow
    listIdx = listIdx + 1
    prevRow <- row
  }
}
bm.filt.sort.merge[[listIdx]] <- prevRow
bm.filt.sort.merge <- do.call('rbind', bm.filt.sort.merge)
head(bm.filt.sort.merge)

# Make a data frame of gene name, chromosome, start and end positions
gene.genomic_coding <- data.frame(gene_name = bm.filt.sort.merge$external_gene_name, chr = bm.filt.sort.merge$chromosome_name, start = bm.filt.sort.merge$genomic_coding_start, end = bm.filt.sort.merge$genomic_coding_end)
head(gene.genomic_coding)

plink.gene.genomic_coding <- gene.genomic_coding[c(2, 3, 4, 1)]
head(plink.gene.genomic_coding)

# Write to file
#write.table(gene.genomic_coding, "bm.exon.pos.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = '\t')
write.table(plink.gene.genomic_coding, "bm.exon.pos.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = '\t')
