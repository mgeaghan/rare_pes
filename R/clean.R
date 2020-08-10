# clean exome sequencing gene-level summary statistics for SCZ, ASD, EPI, and BIP
# Author: Michael Geaghan
# Date: August 2020
source("../../rare_pes/R/gsa.R")
convert <- get_bm_convert()
meta.fisher <- function(x) {
  x2 <- (-2) * sum(log(x), na.rm = TRUE)
  k <- sum(is.finite(x))
  p <- pchisq(q = x2, df = 2*k, lower.tail = FALSE)
  return(p)
}
source("../../rare_pes/R/acat.R")
########## CLEAN DATA ##########
# SCZ Leonenko
df <- read.csv("original_csv/scz.leonenko.csv", header = TRUE, stringsAsFactors = FALSE)
colnames(df) <- c("gene", "chr", "gene.start", "gene.end", "p.skat.o", "p.burden", "or.burden", "n")
df$gene.ensembl <- gene2EnsGene(df$gene, convert)
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
df_na_2 <- df
df <- df[!is.na(df$gene.ensembl) & !is.na(df$p.skat.o),]
df$p.skat.o[df$p.skat.o==1] <- (1 - .Machine$double.eps)
df$p.burden[df$p.burden==1] <- (1 - .Machine$double.eps)
write.csv(df, file = "scz.leonenko.csv", row.names = FALSE, na = "")

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
df$gene.ensembl <- gene2EnsGene(df$gene, convert)
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
df_na_2 <- df
df <- df[!is.na(df$gene.ensembl) & !is.na(df$p),]
df$p[df$p==1] <- (1 - .Machine$double.eps)
write.csv(df, file = "scz.extada.csv", row.names = FALSE, na = "")

# ASC
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
df$gene.ensembl <- ensId2EnsGene(df$gene.ensemblId, convert)
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
df_na_2 <- df
df <- df[!is.na(df$gene.ensembl) & !is.na(df$p),]
df_1 <- df
df_no1 <- df
df_1$p[df_1$p==1] <- (1 - .Machine$double.eps)
df_no1 <- df_no1[df_no1$p != 1,]
write.csv(df_1, file = "asc.csv", row.names = FALSE, na = "")
write.csv(df_no1, file = "asc.no1.csv", row.names = FALSE, na = "")

# EPI: DEE
df <- read.csv("original_csv/dee.csv", header = TRUE, stringsAsFactors = FALSE)
colnames(df) <- c("gene", "description",
                  "cases.lof", "controls.lof", "p.lof",
                  "cases.mpc", "controls.mpc", "p.mpc",
                  "cases.indel.inframe", "controls.indel.inframe", "p.indel.inframe",
                  "p")
df$n <- rowSums(df[grepl("^(cases|controls)\\.", colnames(df), perl = TRUE)], na.rm = TRUE)
df$n[is.na(df$n)] <- 0
df$gene.ensembl <- gene2EnsGene(df$gene, convert)
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
df_na_2 <- df
df <- df[!is.na(df$gene.ensembl) & !is.na(df$p),]
df_1 <- df
df_no1 <- df
df_1$p[df_1$p==1] <- (1 - .Machine$double.eps)
df_no1 <- df_no1[df_no1$p != 1,]
write.csv(df_1, file = "dee.csv", row.names = FALSE, na = "")
write.csv(df_no1, file = "dee.no1.csv", row.names = FALSE, na = "")

# EPI: EPI
df <- read.csv("original_csv/epi.csv", header = TRUE, stringsAsFactors = FALSE)
colnames(df) <- c("gene", "description",
                  "cases.lof", "controls.lof", "p.lof",
                  "cases.mpc", "controls.mpc", "p.mpc",
                  "cases.indel.inframe", "controls.indel.inframe", "p.indel.inframe",
                  "p")
df$n <- rowSums(df[grepl("^(cases|controls)\\.", colnames(df), perl = TRUE)], na.rm = TRUE)
df$n[is.na(df$n)] <- 0
df$gene.ensembl <- gene2EnsGene(df$gene, convert)
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
df_na_2 <- df
df <- df[!is.na(df$gene.ensembl) & !is.na(df$p),]
df_1 <- df
df_no1 <- df
df_1$p[df_1$p==1] <- (1 - .Machine$double.eps)
df_no1 <- df_no1[df_no1$p != 1,]
write.csv(df_1, file = "epi.csv", row.names = FALSE, na = "")
write.csv(df_no1, file = "epi.no1.csv", row.names = FALSE, na = "")

# EPI: GGE
df <- read.csv("original_csv/gge.csv", header = TRUE, stringsAsFactors = FALSE)
colnames(df) <- c("gene", "description",
                  "cases.lof", "controls.lof", "p.lof",
                  "cases.mpc", "controls.mpc", "p.mpc",
                  "cases.indel.inframe", "controls.indel.inframe", "p.indel.inframe",
                  "p")
df$n <- rowSums(df[grepl("^(cases|controls)\\.", colnames(df), perl = TRUE)], na.rm = TRUE)
df$n[is.na(df$n)] <- 0
df$gene.ensembl <- gene2EnsGene(df$gene, convert)
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
df_na_2 <- df
df <- df[!is.na(df$gene.ensembl) & !is.na(df$p),]
df_1 <- df
df_no1 <- df
df_1$p[df_1$p==1] <- (1 - .Machine$double.eps)
df_no1 <- df_no1[df_no1$p != 1,]
write.csv(df_1, file = "gge.csv", row.names = FALSE, na = "")
write.csv(df_no1, file = "gge.no1.csv", row.names = FALSE, na = "")

# EPI: NAFE
df <- read.csv("original_csv/nafe.csv", header = TRUE, stringsAsFactors = FALSE)
colnames(df) <- c("gene", "description",
                  "cases.lof", "controls.lof", "p.lof",
                  "cases.mpc", "controls.mpc", "p.mpc",
                  "cases.indel.inframe", "controls.indel.inframe", "p.indel.inframe",
                  "p")
df$n <- rowSums(df[grepl("^(cases|controls)\\.", colnames(df), perl = TRUE)], na.rm = TRUE)
df$n[is.na(df$n)] <- 0
df$gene.ensembl <- gene2EnsGene(df$gene, convert)
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
df_na_2 <- df
df <- df[!is.na(df$gene.ensembl) & !is.na(df$p),]
df_1 <- df
df_no1 <- df
df_1$p[df_1$p==1] <- (1 - .Machine$double.eps)
df_no1 <- df_no1[df_no1$p != 1,]
write.csv(df_1, file = "nafe.csv", row.names = FALSE, na = "")
write.csv(df_no1, file = "nafe.no1.csv", row.names = FALSE, na = "")

# SCZ SCHEMA
df <- read.csv("original_csv/scz.schema.csv", header = TRUE, stringsAsFactors = FALSE)
colnames(df) <- c("gene", "description", "case.lof", "control.lof", "case.miss.gte.3", "control.miss.gte.3", "case.miss.gte.2.lt.3", "control.miss.gte.2.lt.3", "denovo.lof", "denovo.miss", "p.meta")
df$n <- rowSums(df[c("case.lof", "control.lof", "case.miss.gte.3", "control.miss.gte.3", "case.miss.gte.2.lt.3", "control.miss.gte.2.lt.3", "denovo.lof", "denovo.miss")])
df$gene.ensembl <- ensId2EnsGene(df$gene, convert)
df_na <- df
df <- df[!is.na(df$gene) & !is.na(df$p.meta),]
df_1 <- df
df_no1 <- df
df_1$p.meta[df_1$p.meta==1] <- (1 - .Machine$double.eps)
df_no1 <- df_no1[df_no1$p.meta != 1,]
write.csv(df_1, file = "scz.schema.csv", row.names = FALSE, na = "")
write.csv(df_no1, file = "scz.schema.no1.csv", row.names = FALSE, na = "")

# BIP
df <- read.csv("original_csv/bip.csv", header = TRUE, stringsAsFactors = FALSE)
colnames(df) <- c("gene", "description", "cases", "controls", "cases.ptv", "controls.ptv", "p.ptv", "or.ptv", "cases.d.miss", "controls.d.miss", "p.d.miss", "or.d.miss")
df$n <- rowSums(df[grepl("^(cases|controls)\\.", colnames(df), perl = TRUE)], na.rm = TRUE)
df$n[is.na(df$n)] <- 0
df_na <- df
df <- df[!(is.na(df$p.d.miss) & is.na(df$p.ptv)),]
df$p.meta <- apply(df, 1, function(x) { meta.fisher(c(as.numeric(x[7]), as.numeric(x[11]))) })
minP <- .Machine$double.xmin
maxP <- 1 - .Machine$double.eps
df$p.acat <- apply(df, 1, function(x) {
  p.ptv <- as.numeric(x[7])
  p.d.miss <- as.numeric(x[11])
  if(is.na(p.ptv) & is.na(p.d.miss)) {
    return(NA)
  } else if(is.na(p.ptv)) {
    return(p.d.miss)
  } else if(is.na(p.d.miss)) {
    return(p.ptv)
  } else {
    if(p.ptv > maxP) { p.ptv <- maxP }
    if(p.d.miss > maxP) { p.d.miss <- maxP }
    return(ACAT(c(p.ptv, p.d.miss)))
  }
})
df$gene.ensembl <- ensId2EnsGene(df$gene, convert)
df_na_2 <- df
df <- df[!is.na(df$gene) & (!is.na(df$p.meta) | !is.na(df$p.acat)),]
df_1 <- df
df_no1meta <- df
df_no1acat <- df
df_1$p.meta[df_1$p.meta==1] <- (1 - .Machine$double.eps)
df_1$p.acat[df_1$p.acat==1] <- (1 - .Machine$double.eps)
df_no1meta <- df_no1meta[df_no1meta$p.meta != 1,]
df_no1acat <- df_no1acat[df_no1acat$p.acat != 1,]
write.csv(df_1, file = "bip.csv", row.names = FALSE, na = "")
write.csv(df_no1meta, file = "bip.no1meta.csv", row.names = FALSE, na = "")
write.csv(df_no1acat, file = "bip.no1acat.csv", row.names = FALSE, na = "")

# BIP1
df <- read.csv("original_csv/bip.1.csv", header = TRUE, stringsAsFactors = FALSE)
colnames(df) <- c("gene", "description", "cases", "controls", "cases.ptv", "controls.ptv", "p.ptv", "or.ptv", "cases.d.miss", "controls.d.miss", "p.d.miss", "or.d.miss")
df$n <- rowSums(df[grepl("^(cases|controls)\\.", colnames(df), perl = TRUE)], na.rm = TRUE)
df$n[is.na(df$n)] <- 0
df_na <- df
df <- df[!(is.na(df$p.d.miss) & is.na(df$p.ptv)),]
df$p.meta <- apply(df, 1, function(x) { meta.fisher(c(as.numeric(x[7]), as.numeric(x[11]))) })
minP <- .Machine$double.xmin
maxP <- 1 - .Machine$double.eps
df$p.acat <- apply(df, 1, function(x) {
  p.ptv <- as.numeric(x[7])
  p.d.miss <- as.numeric(x[11])
  if(is.na(p.ptv) & is.na(p.d.miss)) {
    return(NA)
  } else if(is.na(p.ptv)) {
    return(p.d.miss)
  } else if(is.na(p.d.miss)) {
    return(p.ptv)
  } else {
    if(p.ptv > maxP) { p.ptv <- maxP }
    if(p.d.miss > maxP) { p.d.miss <- maxP }
    return(ACAT(c(p.ptv, p.d.miss)))
  }
})
df$gene.ensembl <- ensId2EnsGene(df$gene, convert)
df_na_2 <- df
df <- df[!is.na(df$gene) & (!is.na(df$p.meta) | !is.na(df$p.acat)),]
df_1 <- df
df_no1meta <- df
df_no1acat <- df
df_1$p.meta[df_1$p.meta==1] <- (1 - .Machine$double.eps)
df_1$p.acat[df_1$p.acat==1] <- (1 - .Machine$double.eps)
df_no1meta <- df_no1meta[df_no1meta$p.meta != 1,]
df_no1acat <- df_no1acat[df_no1acat$p.acat != 1,]
write.csv(df_1, file = "bip.1.csv", row.names = FALSE, na = "")
write.csv(df_no1meta, file = "bip.1.no1meta.csv", row.names = FALSE, na = "")
write.csv(df_no1acat, file = "bip.1.no1acat.csv", row.names = FALSE, na = "")

# BIP2
df <- read.csv("original_csv/bip.2.csv", header = TRUE, stringsAsFactors = FALSE)
colnames(df) <- c("gene", "description", "cases", "controls", "cases.ptv", "controls.ptv", "p.ptv", "or.ptv", "cases.d.miss", "controls.d.miss", "p.d.miss", "or.d.miss")
df$n <- rowSums(df[grepl("^(cases|controls)\\.", colnames(df), perl = TRUE)], na.rm = TRUE)
df$n[is.na(df$n)] <- 0
df_na <- df
df <- df[!(is.na(df$p.d.miss) & is.na(df$p.ptv)),]
df$p.meta <- apply(df, 1, function(x) { meta.fisher(c(as.numeric(x[7]), as.numeric(x[11]))) })
minP <- .Machine$double.xmin
maxP <- 1 - .Machine$double.eps
df$p.acat <- apply(df, 1, function(x) {
  p.ptv <- as.numeric(x[7])
  p.d.miss <- as.numeric(x[11])
  if(is.na(p.ptv) & is.na(p.d.miss)) {
    return(NA)
  } else if(is.na(p.ptv)) {
    return(p.d.miss)
  } else if(is.na(p.d.miss)) {
    return(p.ptv)
  } else {
    if(p.ptv > maxP) { p.ptv <- maxP }
    if(p.d.miss > maxP) { p.d.miss <- maxP }
    return(ACAT(c(p.ptv, p.d.miss)))
  }
})
df$gene.ensembl <- ensId2EnsGene(df$gene, convert)
df_na_2 <- df
df <- df[!is.na(df$gene) & (!is.na(df$p.meta) | !is.na(df$p.acat)),]
df_1 <- df
df_no1meta <- df
df_no1acat <- df
df_1$p.meta[df_1$p.meta==1] <- (1 - .Machine$double.eps)
df_1$p.acat[df_1$p.acat==1] <- (1 - .Machine$double.eps)
df_no1meta <- df_no1meta[df_no1meta$p.meta != 1,]
df_no1acat <- df_no1acat[df_no1acat$p.acat != 1,]
write.csv(df_1, file = "bip.2.csv", row.names = FALSE, na = "")
write.csv(df_no1meta, file = "bip.2.no1meta.csv", row.names = FALSE, na = "")
write.csv(df_no1acat, file = "bip.2.no1acat.csv", row.names = FALSE, na = "")

# BIP No Psych
df <- read.csv("original_csv/bip.nopsych.csv", header = TRUE, stringsAsFactors = FALSE)
colnames(df) <- c("gene", "description", "cases", "controls", "cases.ptv", "controls.ptv", "p.ptv", "or.ptv", "cases.d.miss", "controls.d.miss", "p.d.miss", "or.d.miss")
df$n <- rowSums(df[grepl("^(cases|controls)\\.", colnames(df), perl = TRUE)], na.rm = TRUE)
df$n[is.na(df$n)] <- 0
df_na <- df
df <- df[!(is.na(df$p.d.miss) & is.na(df$p.ptv)),]
df$p.meta <- apply(df, 1, function(x) { meta.fisher(c(as.numeric(x[7]), as.numeric(x[11]))) })
minP <- .Machine$double.xmin
maxP <- 1 - .Machine$double.eps
df$p.acat <- apply(df, 1, function(x) {
  p.ptv <- as.numeric(x[7])
  p.d.miss <- as.numeric(x[11])
  if(is.na(p.ptv) & is.na(p.d.miss)) {
    return(NA)
  } else if(is.na(p.ptv)) {
    return(p.d.miss)
  } else if(is.na(p.d.miss)) {
    return(p.ptv)
  } else {
    if(p.ptv > maxP) { p.ptv <- maxP }
    if(p.d.miss > maxP) { p.d.miss <- maxP }
    return(ACAT(c(p.ptv, p.d.miss)))
  }
})
df$gene.ensembl <- ensId2EnsGene(df$gene, convert)
df_na_2 <- df
df <- df[!is.na(df$gene) & (!is.na(df$p.meta) | !is.na(df$p.acat)),]
df_1 <- df
df_no1meta <- df
df_no1acat <- df
df_1$p.meta[df_1$p.meta==1] <- (1 - .Machine$double.eps)
df_1$p.acat[df_1$p.acat==1] <- (1 - .Machine$double.eps)
df_no1meta <- df_no1meta[df_no1meta$p.meta != 1,]
df_no1acat <- df_no1acat[df_no1acat$p.acat != 1,]
write.csv(df_1, file = "bip.nopsych.csv", row.names = FALSE, na = "")
write.csv(df_no1meta, file = "bip.nopsych.no1meta.csv", row.names = FALSE, na = "")
write.csv(df_no1acat, file = "bip.nopsych.no1acat.csv", row.names = FALSE, na = "")

# BIP Psych
df <- read.csv("original_csv/bip.psych.csv", header = TRUE, stringsAsFactors = FALSE)
colnames(df) <- c("gene", "description", "cases", "controls", "cases.ptv", "controls.ptv", "p.ptv", "or.ptv", "cases.d.miss", "controls.d.miss", "p.d.miss", "or.d.miss")
df$n <- rowSums(df[grepl("^(cases|controls)\\.", colnames(df), perl = TRUE)], na.rm = TRUE)
df$n[is.na(df$n)] <- 0
df_na <- df
df <- df[!(is.na(df$p.d.miss) & is.na(df$p.ptv)),]
df$p.meta <- apply(df, 1, function(x) { meta.fisher(c(as.numeric(x[7]), as.numeric(x[11]))) })
minP <- .Machine$double.xmin
maxP <- 1 - .Machine$double.eps
df$p.acat <- apply(df, 1, function(x) {
  p.ptv <- as.numeric(x[7])
  p.d.miss <- as.numeric(x[11])
  if(is.na(p.ptv) & is.na(p.d.miss)) {
    return(NA)
  } else if(is.na(p.ptv)) {
    return(p.d.miss)
  } else if(is.na(p.d.miss)) {
    return(p.ptv)
  } else {
    if(p.ptv > maxP) { p.ptv <- maxP }
    if(p.d.miss > maxP) { p.d.miss <- maxP }
    return(ACAT(c(p.ptv, p.d.miss)))
  }
})
df$gene.ensembl <- ensId2EnsGene(df$gene, convert)
df_na_2 <- df
df <- df[!is.na(df$gene) & (!is.na(df$p.meta) | !is.na(df$p.acat)),]
df_1 <- df
df_no1meta <- df
df_no1acat <- df
df_1$p.meta[df_1$p.meta==1] <- (1 - .Machine$double.eps)
df_1$p.acat[df_1$p.acat==1] <- (1 - .Machine$double.eps)
df_no1meta <- df_no1meta[df_no1meta$p.meta != 1,]
df_no1acat <- df_no1acat[df_no1acat$p.acat != 1,]
write.csv(df_1, file = "bip.psych.csv", row.names = FALSE, na = "")
write.csv(df_no1meta, file = "bip.psych.no1meta.csv", row.names = FALSE, na = "")
write.csv(df_no1acat, file = "bip.psych.no1acat.csv", row.names = FALSE, na = "")

# BIP w/ SA
df <- read.csv("original_csv/bip.sa.csv", header = TRUE, stringsAsFactors = FALSE)
colnames(df) <- c("gene", "description", "cases", "controls", "cases.ptv", "controls.ptv", "p.ptv", "or.ptv", "cases.d.miss", "controls.d.miss", "p.d.miss", "or.d.miss")
df$n <- rowSums(df[grepl("^(cases|controls)\\.", colnames(df), perl = TRUE)], na.rm = TRUE)
df$n[is.na(df$n)] <- 0
df_na <- df
df <- df[!(is.na(df$p.d.miss) & is.na(df$p.ptv)),]
df$p.meta <- apply(df, 1, function(x) { meta.fisher(c(as.numeric(x[7]), as.numeric(x[11]))) })
minP <- .Machine$double.xmin
maxP <- 1 - .Machine$double.eps
df$p.acat <- apply(df, 1, function(x) {
  p.ptv <- as.numeric(x[7])
  p.d.miss <- as.numeric(x[11])
  if(is.na(p.ptv) & is.na(p.d.miss)) {
    return(NA)
  } else if(is.na(p.ptv)) {
    return(p.d.miss)
  } else if(is.na(p.d.miss)) {
    return(p.ptv)
  } else {
    if(p.ptv > maxP) { p.ptv <- maxP }
    if(p.d.miss > maxP) { p.d.miss <- maxP }
    return(ACAT(c(p.ptv, p.d.miss)))
  }
})
df$gene.ensembl <- ensId2EnsGene(df$gene, convert)
df_na_2 <- df
df <- df[!is.na(df$gene) & (!is.na(df$p.meta) | !is.na(df$p.acat)),]
df_1 <- df
df_no1meta <- df
df_no1acat <- df
df_1$p.meta[df_1$p.meta==1] <- (1 - .Machine$double.eps)
df_1$p.acat[df_1$p.acat==1] <- (1 - .Machine$double.eps)
df_no1meta <- df_no1meta[df_no1meta$p.meta != 1,]
df_no1acat <- df_no1acat[df_no1acat$p.acat != 1,]
write.csv(df_1, file = "bip.sa.csv", row.names = FALSE, na = "")
write.csv(df_no1meta, file = "bip.sa.no1meta.csv", row.names = FALSE, na = "")
write.csv(df_no1acat, file = "bip.sa.no1acat.csv", row.names = FALSE, na = "")
