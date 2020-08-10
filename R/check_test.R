load("../R/convert.sets.RData")
geneSetMsigdbEns <- geneSetMsigdbEnsGene
pes <- read.csv("test_pes_out.pes.txt", row.names = 1, stringsAsFactors = FALSE)
dosages <- read.delim("../db/test_plink/test_plink.pes_score.subset.dosages.raw", row.names = 2, header = TRUE, stringsAsFactors = FALSE)
colnames(dosages) <- gsub("(.*)_[^_]+$", "\\1", colnames(dosages))
var_genes <- read.delim("../db/test_plink/test_plink.set.table", header = FALSE, stringsAsFactors = FALSE)
colnames(var_genes) <- c("variant", "gene")
genes <- read.csv("../out/_old_R_4_out/asc.gsa.msigdb.genes.txt", row.names = 1, header = TRUE, stringsAsFactors = FALSE)
annot <- read.delim("../db/annot/wgs.eur.snv.indel.g_indel.annot.rare_0.001.txt", header = TRUE, stringsAsFactors = FALSE)
bim <- read.delim("../db/test_plink/test_plink.bim", row.names = 2, header = FALSE, stringsAsFactors = FALSE)
colnames(bim) <- c("CHR", "CM", "BP", "ALT", "REF")

for (i in rownames(pes)) { if (max(pes, na.rm = TRUE) %in% pes[i,]) { x <- i } }
pathway <- "PID_IL2_PI3K_PATHWAY"

pathway_genes <- geneSetMsigdb[pathway,]
pathway_genes <- pathway_genes[pathway_genes != ""]

pathway_genes_info <- genes[pathway_genes,]
pathway_genes_info <- pathway_genes_info[!is.na(pathway_genes_info$gene),]

pathway_vars <- var_genes[var_genes$gene %in% pathway_genes,]

pathway_dosages <- dosages[x, unique(pathway_vars$variant)]

pathway_rare_vars <- colnames(pathway_dosages)[pathway_dosages != 0]
pathway_rare_var_dosages <- pathway_dosages[pathway_rare_vars]
pathway_rare_var_genes <- var_genes$gene[var_genes$variant %in% pathway_rare_vars]
pathway_rare_vars_info <- bim[pathway_rare_vars,]
pathway_rare_vars_annot <- list()
for (v in pathway_rare_vars) {
  info <- bim[v,]
  chr <- info$CHR
  pos <- info$BP
  ref <- info$REF
  alt <- info$ALT
  var_annot <- annot[annot$Chr == chr & annot$Start == pos & annot$Ref == ref & annot$Alt == alt,]
  pathway_rare_vars_annot[[v]] <- var_annot
}

pathway_rare_var_genes_info <- pathway_genes_info[pathway_rare_var_genes,]
pathway_rare_var_genes_z <- pathway_rare_var_genes_info$z

Ng <- list()
for (g in pathway_rare_var_genes) {
  Ng[[g]] <- length(var_genes$gene[var_genes$gene == g])
}
N <- list()
for (g in var_genes$gene) {
  N[[g]] <- length(var_genes$gene[var_genes$gene == g])
}
Nmin <- min(do.call(c, N))
Nmax <- max(do.call(c, N))

Wg <- list()
for (g in names(Ng)) {
  Wg[[g]] <- 1 - ((Ng[[g]] - Nmin)/(1 + Nmax - Nmin))
}

w_hsp <- pathway_rare_var_dosages$var21194854 * (1 - 10^(pathway_rare_vars_annot$var21194854[1,9]/(-10))) * (dbeta(as.numeric(pathway_rare_vars_annot$var21194854[1,6]),1,25)/25) * Wg$HSP90AA1
w_pik <- pathway_rare_var_dosages$var08700340 * (1 - 10^(pathway_rare_vars_annot$var08700340[1,9]/(-10))) * (dbeta(as.numeric(pathway_rare_vars_annot$var08700340[1,6]),1,25)/25) * Wg$PIK3R1
z_hsp <- pathway_rare_var_genes_info["HSP90AA1",]$z
z_pik <- pathway_rare_var_genes_info["PIK3R1",]$z
pes_theoretical <- ((w_hsp * z_hsp) + (w_pik * z_pik)) / sqrt(w_hsp^2 + w_pik^2)
pes_actual <- max(pes, na.rm = TRUE)
pes_theoretical - pes_actual
