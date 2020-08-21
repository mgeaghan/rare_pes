getFinitePs <- function(x) {
	# convert P = 0 --> machine minimum positive double; P = 1 --> 1 - machine epsilon
	minP <- .Machine$double.xmin
	maxP <- 1 - .Machine$double.eps
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

runGeneSetLogit <- function(df, geneSets, gene = "gene", stat = "z", covars = NULL) {
  if (length(covars) == 1 && is.na(covars)) {
    stop("Argument to covars cannot be NA; for no covariates, supply covars = NULL")
  }
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
		if (is.null(covars)) {
		  fmla <- as.formula(paste("stat ~ set"))
		} else {
		  fmla <- as.formula(paste("stat ~ set +", paste(covars, collapse = " + ")))
		}
		logMods[[s]] <- glm(fmla, data = tmp, family = "binomial")
	}
	return(logMods)
}

runGeneSetLM <- function(df, geneSets, gene = "gene", stat = "z", covars = NULL) {
  if (length(covars) == 1 && is.na(covars)) {
    stop("Argument to covars cannot be NA; for no covariates, supply covars = NULL")
  }
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
		if (is.null(covars)) {
		  fmla <- as.formula(paste("stat ~ set"))
		} else {
		  fmla <- as.formula(paste("stat ~ set +", paste(covars, collapse = " + ")))
		}
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

getLMStats <- function(linMods, tail = "upper") {
  # generate a data frame of stats for all the linear models in the input list "linMods"
  linModStats <- list()
  for (s in names(linMods)) {
    coefs <- summary(linMods[[s]])$coefficients
    if ("set" %in% rownames(coefs)) {
      set.p <- coefs["set", "Pr(>|t|)"]
      set.beta <- coefs["set", "Estimate"]
      set.se <- coefs["set", "Std. Error"]
      set.t <- coefs["set", "t value"]
    } else {
      set.p <- NA
      set.beta <- NA
      set.se <- NA
      set.t <- NA
    }
    log.n.p <- coefs["log.n", "Pr(>|t|)"]
    log.n.beta <- coefs["log.n", "Estimate"]
    log.n.se <- coefs["log.n", "Std. Error"]
    log.n.t <- coefs["log.n", "t value"]
    log.len.p <- coefs["log.len", "Pr(>|t|)"]
    log.len.beta <- coefs["log.len", "Estimate"]
    log.len.se <- coefs["log.len", "Std. Error"]
    log.len.t <- coefs["log.len", "t value"]
    if (tail == "upper") {
      set.p <- pt(set.t, linMods[[s]]$df, lower.tail = FALSE)
      log.n.p <- pt(log.n.t, linMods[[s]]$df, lower.tail = FALSE)
      log.len.p <- pt(log.len.t, linMods[[s]]$df, lower.tail = FALSE)
    } else if (tail == "lower") {
      set.p <- pt(set.t, linMods[[s]]$df, lower.tail = TRUE)
      log.n.p <- pt(log.n.t, linMods[[s]]$df, lower.tail = TRUE)
      log.len.p <- pt(log.len.t, linMods[[s]]$df, lower.tail = TRUE)
    } else if (tail != "both") {
      stop("getLMStats: argument to 'tail' must be either 'both', 'upper', or 'lower'.")
    }
    r2 <- summary(linMods[[s]])$adj.r.squared
    linModStats[[s]] <- data.frame(set = s, set.p = set.p, set.beta = set.beta, set.se = set.se, set.t = set.t,
                                   log.n.p = log.n.p, log.n.beta = log.n.beta, log.n.se = log.n.se, log.n.t = log.n.t,
                                   log.len.p = log.len.p, log.len.beta = log.len.beta, log.len.se = log.len.se, log.len.t = log.len.t)
  }
  linModStats <- as.data.frame(do.call(rbind, linModStats))
  # colnames(linModStats) <- c("set", "p")
  linModStats$set.p.adj <- p.adjust(linModStats$set.p, method = "BH")
  linModStats$set.n <- sapply(rownames(linModStats), function(x) {
    return(table(linMods[[x]]$model$set)[2])
  })
  return(linModStats)
}

get_convert_ensId2Hgnc <- function() {
  # Download BioMart table to convert gene names
  library(reshape2)
  library(biomaRt)
  #ensembl <- useMart('ENSEMBL_MART_ENSEMBL', 'hsapiens_gene_ensembl', host = 'http://grch37.ensembl.org')
  ensembl <- useMart('ENSEMBL_MART_ENSEMBL', 'hsapiens_gene_ensembl')
  bm <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'), mart = ensembl)
  bm <- unique(bm)
  colnames(bm) <- c("gene.ensembl", "gene.name")
  bm <- bm[bm$gene.name != "",]
  bm <- bm[(!is.na(bm$gene.ensembl)) & (!is.na(bm$gene.name)),]
  bm <- bm[!(bm$gene.ensembl %in% bm$gene.ensembl[duplicated(bm$gene.ensembl)]),]
  return(unique(bm))
}

loadDf <- function(inputFile, rm.n0 = TRUE) {
  df <- read.csv(inputFile, header = TRUE)
  rownames(df) <- df[[geneHeader]]
  # remove genes with n (minor allele count) == 0
  if (rm.n0) {
    df <- df[df[[nHeader]] > 0,]
  }
  return(df)
}

addLogNLogLen <- function(df, lengthFile, nHeader) {
  df$len <- getGeneLengths(rownames(df), lengthFile)
  df <- df[!is.na(df$len),]
  df <- df[!is.na(df[[nHeader]]),]
  df$log.n <- log(df[[nHeader]])
  df$log.len <- log(df$len)
  return(df)
}

runGSA <- function(df, geneSet, geneHeader, pHeader, covars, outPrefix, pThreshold = NA, alpha = 0.05, tail = "upper") {
  # apply pvalue threshold if present
  if (!is.na(pThreshold)) {
    df <- df[df[[pHeader]] < pThreshold,]
  }
  
  df$p.finite <- getFinitePs(df[[pHeader]])
  df$z <- zTransformPValue(df$p.finite)
  
  linMods <- runGeneSetLM(df, geneSet, gene = geneHeader, covars = covars)
  linModStats <- getLMStats(linMods, tail)
  
  gsaPrefix <- paste(outPrefix, "gsa", sep = ".")
  
  write.table(linModStats, paste(gsaPrefix, ".txt", sep = ""), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  write.table(linModStats$set[!is.na(linModStats$set.p) & linModStats$set.p < alpha], paste(gsaPrefix, ".nomsigsets.txt", sep = ""), quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(linModStats$set[!is.na(linModStats$set.p.adj) & linModStats$set.p.adj < alpha], paste(gsaPrefix, ".sigsets.txt", sep = ""), quote = FALSE, row.names = FALSE, col.names = FALSE)

  dfSets <- addGeneSets(df, geneSet, geneHeader)
  
  write.csv(dfSets, paste(gsaPrefix, ".genes.txt", sep = ""), row.names = FALSE)
}

runBfGSA <- function(df, geneSet, geneHeader, bfHeader, covars, outPrefix, alpha = 0.05, tail = "upper") {
  df$log.bf <- log(df[[bfHeader]])
  
  linMods <- runGeneSetLM(df, geneSet, gene = geneHeader, stat = "log.bf", covars = covars)
  linModStats <- getLMStats(linMods, tail)
  
  gsaPrefix <- paste(outPrefix, "gsa", sep = ".")

  write.table(linModStats, paste(gsaPrefix, ".txt", sep = ""), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  write.table(linModStats$set[!is.na(linModStats$set.p) & linModStats$set.p < alpha], paste(gsaPrefix, ".nomsigsets.txt", sep = ""), quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(linModStats$set[!is.na(linModStats$set.p.adj) & linModStats$set.p.adj < alpha], paste(gsaPrefix, ".sigsets.txt", sep = ""), quote = FALSE, row.names = FALSE, col.names = FALSE)

  dfSets <- addGeneSets(df, geneSet, geneHeader)

  write.csv(dfSets, paste(gsaPrefix, ".genes.txt", sep = ""), row.names = FALSE)
}
