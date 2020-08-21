library(biomaRt)
ensembl <- useMart('ENSEMBL_MART_ENSEMBL', 'hsapiens_gene_ensembl', host = 'http://grch37.ensembl.org')
# get gene and transcript IDs, trancript type, gene name, chromosome and genomic coding start and end
bm <- getBM(attributes = c('ensembl_gene_id', 'ensembl_transcript_id', 'ensembl_gene_id_version', 'ensembl_transcript_id_version', 'transcript_biotype', 'external_gene_name', 'chromosome_name', 'genomic_coding_start', 'genomic_coding_end'), mart = ensembl)
bm2 <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'), mart = ensembl)
bm2 <- unique(bm2)
colnames(bm2) <- c("gene.ensembl", "gene.name")
bm2 <- bm2[bm2$gene.name != "",]
bm2 <- bm2[(!is.na(bm2$gene.ensembl)) & (!is.na(bm2$gene.name)),]
bm2 <- bm2[!(bm2$gene.ensembl %in% bm2$gene.ensembl[duplicated(bm2$gene.ensembl)]),]
rownames(bm2) <- bm2$gene.ensembl
chrs <- c(seq(1, 22), 'X', 'Y', 'MT')
chrs
bm.filt <- bm[bm$chromosome_name %in% chrs & !is.na(bm$genomic_coding_start) & !is.na(bm$genomic_coding_end) & bm$transcript_biotype == 'protein_coding', ]
# Also remove two erroneous gene entries (gene names are duplicated and are on different chromosomes):
# CKS1B on chr5 (real gene is on chr1)
# LSP1 on chr13 (real gene is on chr11)
bm.filt <- bm.filt[!(bm.filt$external_gene_name %in% c("CKS1B", "LSP1") & bm.filt$chromosome_name %in% c("5", "13")),]
bm.filt$gene.name <- bm2[bm.filt$ensembl_gene_id,]$gene.name
bm.filt <- bm.filt[!is.na(bm.filt$gene.name),]
bm.filt.sort <- bm.filt[order(bm.filt$gene.name, bm.filt$genomic_coding_start), ]
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

# Get exon union lengths
bm.filt.sort.merge$genomic_coding_len <- bm.filt.sort.merge$genomic_coding_end - (bm.filt.sort.merge$genomic_coding_start - 1)
gene.list <- unique(bm.filt.sort.merge$gene.name)
gene.genomic_coding.len <- data.frame(gene_name = gene.list)
gene.genomic_coding.len$len <- do.call('c', lapply(gene.list, function(x) sum(bm.filt.sort.merge$genomic_coding_len[bm.filt.sort.merge$gene.name == x])))
head(gene.genomic_coding.len)

# Write to file
write.table(gene.genomic_coding.len, "bm.exon.union.length.txt", quote = FALSE, row.names = FALSE, col.names = TRUE, sep = '\t')
