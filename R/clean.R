# Download BioMart table to convert gene names
library(reshape2)
library(biomaRt)
#ensembl <- useMart('ENSEMBL_MART_ENSEMBL', 'hsapiens_gene_ensembl', host = 'http://grch37.ensembl.org')
ensembl <- useMart('ENSEMBL_MART_ENSEMBL', 'hsapiens_gene_ensembl')
bm <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'), mart = ensembl)
bm2 <- getBM(attributes = c('ensembl_gene_id', 'external_synonym'), mart = ensembl)
bm <- unique(bm)
bm2 <- unique(bm2)
colnames(bm) <- c("gene.ensembl", "gene.name")
colnames(bm2) <- c("gene.ensembl", "gene.name")
bm <- bm[bm$gene.name != "",]
bm2 <- bm2[bm2$gene.name != "",]
bm <- bm[(!is.na(bm$gene.ensembl)) & (!is.na(bm$gene.name)),]
bm2 <- bm2[(!is.na(bm2$gene.ensembl)) & (!is.na(bm2$gene.name)),]
bm2 <- unique(rbind(bm2, bm))
convert <- merge(bm2, bm, by = "gene.ensembl", all = TRUE)
colnames(convert) <- c("gene.ensembl", "gene.name", "gene.name.ensembl")
convert2 <- unique(convert[c("gene.name", "gene.name.ensembl")])
convert2 <- unique(do.call(rbind, apply(convert2, 1, function(x) {
  gene <- as.character(x[1])
  gene.ens <- as.character(x[2])
  tmp <- convert2[convert2$gene.name==gene,2]
  if (gene %in% tmp) {
    return(data.frame(gene.name = gene, gene.name.ensembl = gene))
  } else {
    return(c(gene.name = gene, gene.name.ensembl = gene.ens))
  }
})))
# Function to return a vector of the number of synonyms for each of a vector of gene symbols
numSyn <- function(x) {
  n <- sapply(x, function(y) {
    return(dim(convert2[convert2$gene.name==y,])[1])
  })
  return(n)
}
# Function to return a vector of Ensembl gene names for a given vector of gene symbols
ensGene <- function(x) {
  id <- sapply(x, function(y) {
    ids <- convert2$gene.name.ensembl[convert2$gene.name == y]
    if (length(ids) == 1) {
      return(ids)
    } else {
      return(NA)
    }
  })
  return(id)
}
# Function to return a vector of Ensembl gene names for a given vector of Ensembl gene IDs
ensId2Gene <- function(x) {
  id <- sapply(x, function(y) {
    ids <- bm$gene.name[bm$gene.ensembl == y]
    if (length(ids) == 1) {
      return(ids)
    } else {
      return(NA)
    }
  })
  return(id)
}
########## CLEAN DATA ##########
# SCZ Leonenko
df <- read.csv("original_csv/scz.leonenko.csv", header = TRUE, stringsAsFactors = FALSE)
colnames(df) <- c("gene", "chr", "gene.start", "gene.end", "p.skat.o", "p.burden", "or.burden", "n")
df$gene.ensembl <- ensGene(df$gene)
df <- df[c("gene.ensembl", "gene", "chr", "gene.start", "gene.end", "p.skat.o", "p.burden", "or.burden", "n")]
df_na <- df
df <- df[!is.na(df$gene.ensembl),]
dfDupNames <- df$gene.ensembl[duplicated(df$gene.ensembl)]
dfDups <- do.call(rbind, lapply(dfDupNames, function(x) {
  dup <- df[df$gene.ensembl == x,]
  if (x %in% dup$gene) {
    dup <- dup[dup$gene == x,]
    if (dim(dup)[1] == 1) {
      return(dup)
    } else {
      return(dup[dup$p.skat.o == max(dup$p.skat.o),][1,])
    }
  } else {
    return(dup[dup$p.skat.o == max(dup$p.skat.o),][1,])
  }
}))
df <- rbind(df[!(df$gene.ensembl %in% dfDupNames),],
            dfDups)
write.csv(df, file = "scz.leonenko.csv", row.names = FALSE, col.names = TRUE, na = "")

# SCZ exTADA
df <- read.csv("original_csv/scz.extada.csv", header = TRUE, stringsAsFactors = FALSE)
colnames(df) <- c("gene", "denovos", "mutation.rates",
                  "cases.dbs1", "cases.dbs3", "cases.clin1", "cases.clin2", "cases.clin3",
                  "cases.swe1", "cases.swe2", "cases.swe3",
                  "controls.dbs1", "controls.dbs3", "controls.clin1", "controls.clin2", "controls.clin3",
                  "controls.swe1", "controls.swe2", "controls.swe3",
                  "bf", "pp", "q", "p")
df$n <- rowSums(df[grepl("^(cases|controls)\\.", colnames(df), perl = TRUE)], na.rm = TRUE)
df$n[is.na(df$n)] <- 0
df$gene.ensembl <- ensGene(df$gene)
df <- df[c("gene.ensembl", "gene", "denovos", "mutation.rates",
           "cases.dbs1", "cases.dbs3", "cases.clin1", "cases.clin2", "cases.clin3",
           "cases.swe1", "cases.swe2", "cases.swe3",
           "controls.dbs1", "controls.dbs3", "controls.clin1", "controls.clin2", "controls.clin3",
           "controls.swe1", "controls.swe2", "controls.swe3",
           "bf", "pp", "q", "p", "n")]
df_na <- df
df <- df[!is.na(df$gene.ensembl),]
dfDupNames <- df$gene.ensembl[duplicated(df$gene.ensembl)]
dfDups <- do.call(rbind, lapply(dfDupNames, function(x) {
  dup <- df[df$gene.ensembl == x,]
  if (x %in% dup$gene) {
    dup <- dup[dup$gene == x,]
    if (dim(dup)[1] == 1) {
      return(dup)
    } else {
      return(dup[dup$p == max(dup$p),][1,])
    }
  } else {
    return(dup[dup$p == max(dup$p),][1,])
  }
}))
df <- rbind(df[!(df$gene.ensembl %in% dfDupNames),],
            dfDups)
write.csv(df, file = "scz.extada.csv", row.names = FALSE, col.names = TRUE, na = "")

# EPI: ASC
df <- read.csv("original_csv/asc.csv", header = TRUE, stringsAsFactors = FALSE)
colnames(df) <- c("gene", "description",
                  "cases.ptv.denovo", "controls.ptv.denovo",
                  "cases.mis_a.denovo", "controls.mis_a.denovo",
                  "cases.mis_b.denovo", "controls.mis_b.denovo",
                  "cases.ptv.dbs", "controls.ptv.dbs",
                  "cases.ptv.swe", "controls.ptv.swe",
                  "transmitted", "untransmitted", "q")
asc.p <- read.delim("original_csv/asc.pvalues.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
df$p <- sapply(df$gene, function(x) {
  return(asc.p$p[asc.p$gene == x])
})
df$gene.ensemblId <- sapply(df$gene, function(x) {
  return(asc.p$ensembl_gene_id[asc.p$gene == x])
})
df$gene.ensembl <- ensId2Gene(df$gene.ensemblId)
df$n <- rowSums(df[grepl("^(cases|controls)\\.", colnames(df), perl = TRUE)], na.rm = TRUE)
df$n[is.na(df$n)] <- 0
df <- df[c("gene.ensembl", "gene", "gene.ensemblId", "description",
           "cases.ptv.denovo", "controls.ptv.denovo",
           "cases.mis_a.denovo", "controls.mis_a.denovo",
           "cases.mis_b.denovo", "controls.mis_b.denovo",
           "cases.ptv.dbs", "controls.ptv.dbs",
           "cases.ptv.swe", "controls.ptv.swe",
           "transmitted", "untransmitted", "q", "p", "n")]
df_na <- df
df <- df[!is.na(df$gene.ensembl),]
dfDupNames <- df$gene.ensembl[duplicated(df$gene.ensembl)]
dfDups <- do.call(rbind, lapply(dfDupNames, function(x) {
  dup <- df[df$gene.ensembl == x,]
  if (x %in% dup$gene) {
    dup <- dup[dup$gene == x,]
    if (dim(dup)[1] == 1) {
      return(dup)
    } else {
      return(dup[dup$p == max(dup$p),][1,])
    }
  } else {
    return(dup[dup$p == max(dup$p),][1,])
  }
}))
df <- rbind(df[!(df$gene.ensembl %in% dfDupNames),],
            dfDups)
write.csv(df, file = "asc.csv", row.names = FALSE, col.names = TRUE, na = "")

# EPI: DEE
df <- read.csv("original_csv/dee.csv", header = TRUE, stringsAsFactors = FALSE)
colnames(df) <- c("gene", "description",
                  "cases.lof", "controls.lof", "p.lof",
                  "cases.mpc", "controls.mpc", "p.mpc",
                  "cases.indel.inframe", "controls.indel.inframe", "p.indel.inframe",
                  "p")
df$n <- rowSums(df[grepl("^(cases|controls)\\.", colnames(df), perl = TRUE)], na.rm = TRUE)
df$n[is.na(df$n)] <- 0
df$gene.ensembl <- ensGene(df$gene)
df <- df[c("gene.ensembl", "gene", "description",
           "cases.lof", "controls.lof", "p.lof",
           "cases.mpc", "controls.mpc", "p.mpc",
           "cases.indel.inframe", "controls.indel.inframe", "p.indel.inframe",
           "p", "n")]
df_na <- df
df <- df[!is.na(df$gene.ensembl),]
dfDupNames <- df$gene.ensembl[duplicated(df$gene.ensembl)]
dfDups <- do.call(rbind, lapply(dfDupNames, function(x) {
  dup <- df[df$gene.ensembl == x,]
  if (x %in% dup$gene) {
    dup <- dup[dup$gene == x,]
    if (dim(dup)[1] == 1) {
      return(dup)
    } else {
      return(dup[dup$p == max(dup$p),][1,])
    }
  } else {
    return(dup[dup$p == max(dup$p),][1,])
  }
}))
df <- rbind(df[!(df$gene.ensembl %in% dfDupNames),],
            dfDups)
write.csv(df, file = "dee.csv", row.names = FALSE, col.names = TRUE, na = "")

# EPI: EPI
df <- read.csv("original_csv/epi.csv", header = TRUE, stringsAsFactors = FALSE)
colnames(df) <- c("gene", "description",
                  "cases.lof", "controls.lof", "p.lof",
                  "cases.mpc", "controls.mpc", "p.mpc",
                  "cases.indel.inframe", "controls.indel.inframe", "p.indel.inframe",
                  "p")
df$n <- rowSums(df[grepl("^(cases|controls)\\.", colnames(df), perl = TRUE)], na.rm = TRUE)
df$n[is.na(df$n)] <- 0
df$gene.ensembl <- ensGene(df$gene)
df <- df[c("gene.ensembl", "gene", "description",
           "cases.lof", "controls.lof", "p.lof",
           "cases.mpc", "controls.mpc", "p.mpc",
           "cases.indel.inframe", "controls.indel.inframe", "p.indel.inframe",
           "p", "n")]
df_na <- df
df <- df[!is.na(df$gene.ensembl),]
dfDupNames <- df$gene.ensembl[duplicated(df$gene.ensembl)]
dfDups <- do.call(rbind, lapply(dfDupNames, function(x) {
  dup <- df[df$gene.ensembl == x,]
  if (x %in% dup$gene) {
    dup <- dup[dup$gene == x,]
    if (dim(dup)[1] == 1) {
      return(dup)
    } else {
      return(dup[dup$p == max(dup$p),][1,])
    }
  } else {
    return(dup[dup$p == max(dup$p),][1,])
  }
}))
df <- rbind(df[!(df$gene.ensembl %in% dfDupNames),],
            dfDups)
write.csv(df, file = "epi.csv", row.names = FALSE, col.names = TRUE, na = "")

# EPI: GGE
df <- read.csv("original_csv/gge.csv", header = TRUE, stringsAsFactors = FALSE)
colnames(df) <- c("gene", "description",
                  "cases.lof", "controls.lof", "p.lof",
                  "cases.mpc", "controls.mpc", "p.mpc",
                  "cases.indel.inframe", "controls.indel.inframe", "p.indel.inframe",
                  "p")
df$n <- rowSums(df[grepl("^(cases|controls)\\.", colnames(df), perl = TRUE)], na.rm = TRUE)
df$n[is.na(df$n)] <- 0
df$gene.ensembl <- ensGene(df$gene)
df <- df[c("gene.ensembl", "gene", "description",
           "cases.lof", "controls.lof", "p.lof",
           "cases.mpc", "controls.mpc", "p.mpc",
           "cases.indel.inframe", "controls.indel.inframe", "p.indel.inframe",
           "p", "n")]
df_na <- df
df <- df[!is.na(df$gene.ensembl),]
dfDupNames <- df$gene.ensembl[duplicated(df$gene.ensembl)]
dfDups <- do.call(rbind, lapply(dfDupNames, function(x) {
  dup <- df[df$gene.ensembl == x,]
  if (x %in% dup$gene) {
    dup <- dup[dup$gene == x,]
    if (dim(dup)[1] == 1) {
      return(dup)
    } else {
      return(dup[dup$p == max(dup$p),][1,])
    }
  } else {
    return(dup[dup$p == max(dup$p),][1,])
  }
}))
df <- rbind(df[!(df$gene.ensembl %in% dfDupNames),],
            dfDups)
write.csv(df, file = "gge.csv", row.names = FALSE, col.names = TRUE, na = "")

# EPI: NAFE
df <- read.csv("original_csv/nafe.csv", header = TRUE, stringsAsFactors = FALSE)
colnames(df) <- c("gene", "description",
                  "cases.lof", "controls.lof", "p.lof",
                  "cases.mpc", "controls.mpc", "p.mpc",
                  "cases.indel.inframe", "controls.indel.inframe", "p.indel.inframe",
                  "p")
df$n <- rowSums(df[grepl("^(cases|controls)\\.", colnames(df), perl = TRUE)], na.rm = TRUE)
df$n[is.na(df$n)] <- 0
df$gene.ensembl <- ensGene(df$gene)
df <- df[c("gene.ensembl", "gene", "description",
           "cases.lof", "controls.lof", "p.lof",
           "cases.mpc", "controls.mpc", "p.mpc",
           "cases.indel.inframe", "controls.indel.inframe", "p.indel.inframe",
           "p", "n")]
df_na <- df
df <- df[!is.na(df$gene.ensembl),]
dfDupNames <- df$gene.ensembl[duplicated(df$gene.ensembl)]
dfDups <- do.call(rbind, lapply(dfDupNames, function(x) {
  dup <- df[df$gene.ensembl == x,]
  if (x %in% dup$gene) {
    dup <- dup[dup$gene == x,]
    if (dim(dup)[1] == 1) {
      return(dup)
    } else {
      return(dup[dup$p == max(dup$p),][1,])
    }
  } else {
    return(dup[dup$p == max(dup$p),][1,])
  }
}))
df <- rbind(df[!(df$gene.ensembl %in% dfDupNames),],
            dfDups)
write.csv(df, file = "nafe.csv", row.names = FALSE, col.names = TRUE, na = "")
