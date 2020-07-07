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

runGeneSetLogit <- function(df, geneSets, gene = "gene", stat = "z", covars) {
	# given a data frame with a column of gene names, summary statistic 'stat' and covariates supplied in vector 'covars', perform gene set analysis using a logistic regression model
	if (FALSE %in% (c(gene, stat, covars) %in% colnames(df))) {
		stop("Dataframe is missing covariate columns.")
	}
	newDf <- df[c(gene, stat, covars)]
	geneSetNames <- rownames(geneSets)
	logMods <- list()
	for (s in geneSetNames) {
		newDf[[s]] <- as.numeric(newDf$gene %in% as.character(geneSets[s,]))
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
		newDf[[s]] <- as.numeric(newDf$gene %in% as.character(geneSets[s,]))
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
			logModStats[[s]] <- c(s, summary(logMods[[s]])$coefficients["stat", "Pr(>|z|)"])
		} else {
			logModStats[[s]] <- c(s, NA)
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
			linModStats[[s]] <- c(s, summary(linMods[[s]])$coefficients["set", "Pr(>|t|)"])
		} else {
			linModStats[[s]] <- c(s, NA)
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
