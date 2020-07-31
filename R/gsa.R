getFinitePs <- function(x) {
	# convert P = 0 --> machine epsilon; P = 1 --> 1 - machine epsilon
	minP <- .Machine$double.eps
	maxP <- 1 - minP
	p.finite <- x
	p.finite[p.finite < minP] <- minP
	p.finite[p.finite > maxP] <- maxP
	return(p.finite)
}

zTransformPValue <- function(x) {
	# transform p to z
	return(qnorm(1 - x))
}

getGeneLengths <- function(gene.names, file, header = TRUE, row.names = "gene_name", stringsAsFactors = FALSE) {
	# return gene lengths for a set of gene names
	len <- read.delim(file, header = header, row.names = row.names)
	gene.len <- len[gene.names,]
	return(gene.len)
}

getGeneSets <- function(file) {
	# read in a .gmt gene set file
	ncol <- max(count.fields(file, sep = "\t"))
	geneSets <- read.delim(file, sep = "\t", header = FALSE, stringsAsFactors = FALSE, col.names = paste0("V", seq_len(ncol)), row.names = 1)
	geneSets <- geneSets[2:dim(geneSets)[2]]
	return(geneSets)
}

addGeneSets <- function(df, geneSets, geneHeader = "gene.ensembl") {
  # Add gene set information to data frame
  if (TRUE %in% (rownames(geneSets) %in% colnames(df))) {
    warning("Columns present in supplied data.frame matching gene sets! These will be overwritten!")
  }
  newDf <- df
  geneSetNames <- rownames(geneSets)
  for (s in geneSetNames) {
    newDf[[s]] <- as.numeric(newDf[[geneHeader]] %in% as.character(geneSets[s,]))
  }
  return(newDf)
}

runGeneSetLogit <- function(df, geneSets, gene = "gene", stat = "z", covars) {
	# given a data frame with a column of gene names, summary statistic 'stat' and covariates supplied in vector 'covars', perform gene set analysis using a logistic regression model
	if (FALSE %in% (c(gene, stat, covars) %in% colnames(df))) {
		stop("Dataframe is missing covariate columns.")
	}
	newDf <- df[c(gene, stat, covars)]
	geneSetNames <- rownames(geneSets)
	logMods <- list()
	for (s in geneSetNames) {
		newDf[[s]] <- as.numeric(newDf[[gene]] %in% as.character(geneSets[s,]))
		tmp <- newDf[c(stat, s, covars)]
		colnames(tmp) <- c("stat", "set", covars)
		fmla <- as.formula(paste("set ~ stat +", paste(covars, collapse = " + ")))
		logMods[[s]] <- glm(fmla, data = tmp, family = "binomial")
	}
	return(logMods)
}

runGeneSetLM <- function(df, geneSets, gene = "gene", stat = "z", covars) {
	# given a data frame with a column of gene names, summary statistic 'stat' and covariates supplied in vector 'covars', perform gene set analysis using a linear regression model
	if (FALSE %in% (c(gene, stat, covars) %in% colnames(df))) {
		stop("Dataframe is missing covariate columns.")
	}
	newDf <- df[c(gene, stat, covars)]
	geneSetNames <- rownames(geneSets)
	linMods <- list()
	for (s in geneSetNames) {
		newDf[[s]] <- as.numeric(newDf[[gene]] %in% as.character(geneSets[s,]))
		tmp <- newDf[c(stat, s, covars)]
		colnames(tmp) <- c("stat", "set", covars)
		fmla <- as.formula(paste("stat ~ set +", paste(covars, collapse = " + ")))
		linMods[[s]] <- lm(fmla, data = tmp)
	}
	return(linMods)
}

getLogitStats <- function(logMods) {
	# generate a data frame of stats for all the linear models in the input list "linMods"
	logModStats <- list()
	for (s in names(logMods)) {
		coefs <- summary(logMods[[s]])$coefficients
		if ("stat" %in% rownames(coefs)) {
			logModStats[[s]] <- data.frame(set = s, p = summary(logMods[[s]])$coefficients["stat", "Pr(>|z|)"])
		} else {
			logModStats[[s]] <- data.frame(set = s, p = NA)
		}
	}
	logModStats <- as.data.frame(do.call(rbind, logModStats))
	colnames(logModStats) <- c("set", "p")
	logModStats$p.adj <- p.adjust(logModStats$p, method = "BH")
	logModStats$n <- sapply(rownames(logModStats), function(x) {
		return(table(logMods[[x]]$model$set)[2])
	})
	return(logModStats)
}

getLMStats <- function(linMods) {
	# generate a data frame of stats for all the linear models in the input list "linMods"
	linModStats <- list()
	for (s in names(linMods)) {
		coefs <- summary(linMods[[s]])$coefficients
		if ("set" %in% rownames(coefs)) {
			linModStats[[s]] <- data.frame(set = s, p = summary(linMods[[s]])$coefficients["set", "Pr(>|t|)"])
		} else {
			linModStats[[s]] <- data.frame(set = s, p = NA)
		}
	}
	linModStats <- as.data.frame(do.call(rbind, linModStats))
	colnames(linModStats) <- c("set", "p")
	linModStats$p.adj <- p.adjust(linModStats$p, method = "BH")
	linModStats$n <- sapply(rownames(linModStats), function(x) {
		return(table(linMods[[x]]$model$set)[2])
	})
	return(linModStats)
}

get_bm_convert <- function() {
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
  convert <- merge(bm2, bm, by = "gene.ensembl", all.x = TRUE)
  colnames(convert) <- c("gene.ensembl", "gene.name", "gene.name.ensembl")
  return(unique(convert))
}

# Function to return a vector of Ensembl gene names for a given vector of gene symbols
gene2EnsGene <- function(x, convertDf) {
  id <- sapply(x, function(y) {
    ids <- unique(convertDf$gene.name.ensembl[convertDf$gene.name == y])
    if (length(ids) == 1) {
      return(ids)
    } else {
      if (y %in% ids) {
        return(y)
      } else {
        return(NA)
      }
    }
  })
  return(id)
}

# Function to return a vector of Ensembl gene IDs for a given vector of gene symbols
gene2EnsId <- function(x, convertDf) {
  id <- sapply(x, function(y) {
    ids <- unique(convertDf$gene.ensembl[convertDf$gene.name == y])
    if (length(ids) == 1) {
      return(ids)
    } else {
      ids <- unique(convertDf$gene.ensembl[convertDf$gene.name.ensembl == y])
      if (length(ids) == 1) {
        return(ids)
      } else {
        return(NA)
      }
    }
  })
  return(id)
}

# Function to return a vector of Ensembl gene names for a given vector of Ensembl gene IDs
ensId2EnsGene <- function(x, convertDf) {
  id <- sapply(x, function(y) {
    ids <- unique(convertDf$gene.name.ensembl[convertDf$gene.ensembl == y])
    if (length(ids) == 1) {
      return(ids)
    } else {
      return(NA)
    }
  })
  return(id)
}

# Returns a more efficient lookup list for converting gene names to Ensembl gene names
get_convert_gene2EnsGene <- function(convertDf) {
  newDf <- convertDf[order(convertDf$gene.name),]
  geneList <- unique(newDf$gene.name)
  newList <- lapply(geneList, function(x) {
    df <- unique(newDf[newDf$gene.name == x,])
    if (dim(df)[1] == 1) {
      return(c(x, df$gene.name.ensembl))
    } else {
      df <- unique(df[df$gene.name == df$gene.name.ensembl,])
      if (dim(df)[1] == 1) {
        return(c(x, df$gene.name.ensembl))
      } else {
        return(c(x, NA))
      }
    }
  })
  newDf <- as.data.frame(do.call(rbind, newList))
  colnames(newDf) <- c("gene.name", "ensemble.ID")
  rownames(newDf) <- newDf$gene.name
  return(newDf)
}

# Returns a more efficient lookup list for converting gene names to Ensembl IDs
get_convert_gene2EnsId <- function(convertDf) {
  newDf <- convertDf[order(convertDf$gene.name),]
  geneList <- unique(newDf$gene.name)
  newList <- lapply(geneList, function(x) {
    df <- unique(newDf[newDf$gene.name == x,])
    if (dim(df)[1] == 1) {
      return(c(x, df$gene.ensembl))
    } else {
      df <- unique(df[df$gene.name == df$gene.name.ensembl,])
      if (dim(df)[1] == 1) {
        return(c(x, df$gene.ensembl))
      } else {
        return(c(x, NA))
      }
    }
  })
  newDf <- as.data.frame(do.call(rbind, newList))
  colnames(newDf) <- c("gene.name", "ensemble.ID")
  rownames(newDf) <- newDf$gene.name
  return(newDf)
}

# Convert gene sets to Ensembl gene names or IDs supplied in convertDf
geneset2EnsGene <- function(x, convertDf) {
  cols <- colnames(x)
  newDf <- apply(x, c(1, 2), function(y) {
    newId <- convertDf[y,2]
    if (is.na(newId)) {
      return("")
    } else {
      return(newId)
    }
  })
  colnames(newDf) <- cols
  return(newDf)
}

runDefaultGSA <- function(inputFile, geneHeader, pHeader, outPrefix) {
  # # read in gene-level results from study
  # df <- read.csv(inputFile, header = TRUE)
  # rownames(df) <- df[[geneHeader]]
  # # remove genes with n (minor allele count) == 0
  # df <- df[df$n > 0,]
  # 
  # df$p.finite <- getFinitePs(df[[pHeader]])
  # df$z <- zTransformPValue(df$p.finite)
  # df$len <- getGeneLengths(rownames(df), "../db/len.txt")
  # df <- df[!is.na(df$len),]
  # df$log.n <- log(df$n)
  # df$log.len <- log(df$len)
  # geneSetMsigdb <- getGeneSets("../db/gene_sets/MsigDB_canonical_and_hallmark_Tclin_v_6.1.0.gmt")
  # linModsMsigdb <- runGeneSetLM(df, geneSetMsigdb, gene = "gene.ensembl", covars = c("log.n", "log.len"))
  # linModStatsMsigdb <- getLMStats(linModsMsigdb)
  # geneSetTclin <- getGeneSets("../db/gene_sets/Tclin_enriched_MsigDB_hallmark_canonical_TCRD_v.6.1.0.gmt")
  # linModsTclin <- runGeneSetLM(df, geneSetTclin, gene = "gene.ensembl", covars = c("log.n", "log.len"))
  # linModStatsTclin <- getLMStats(linModsTclin)
  # 
  # msigdbPrefix <- paste(outPrefix, "gsa.msigdb", sep = ".")
  # tclinPrefix <- paste(outPrefix, "gsa.msigdb", sep = ".")
  # 
  # write.table(linModStatsMsigdb, paste(msigdbPrefix, ".txt", sep = ""), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  # write.table(linModStatsTclin, paste(tclinPrefix, ".txt", sep = ""), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  # write.table(linModStatsMsigdb$set[linModStatsMsigdb$p < 0.05], paste(msigdbPrefix, ".sigsets.txt", sep = ""), quote = FALSE, row.names = FALSE, col.names = FALSE)
  # write.table(linModStatsTclin$set[linModStatsTclin$p < 0.05], paste(tclinPrefix, ".sigsets.txt", sep = ""), quote = FALSE, row.names = FALSE, col.names = FALSE)
  # 
  # dfMsigdb <- addGeneSets(df, geneSetMsigdb)
  # dfTclin <- addGeneSets(df, geneSetTclin)
  # 
  # write.csv(dfMsigdb, paste(msigdbPrefix, ".genes.txt", sep = ""), quote = FALSE, row.names = FALSE)
  # write.csv(dfTclin, paste(tclinPrefix, ".genes.txt", sep = ""), quote = FALSE, row.names = FALSE)
  # 
  # # save(df, linModsMsigdb, linModStatsTclin, linModStatsMsigdb, linModStatsTclin, dfMsigdb, dfTclin, file = paste(outPrefix, ".gsa.RData", sep = ""))
  # 
}