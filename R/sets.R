source("../../../rare_pes/R/gsa.R")
geneSetMsigdb <- getGeneSets("./MsigDB_canonical_and_hallmark_Tclin_v_6.1.0.gmt")
geneSetTclin <- getGeneSets("./Tclin_enriched_MsigDB_hallmark_canonical_TCRD_v.6.1.0.gmt")
save(geneSetMsigdb, geneSetTclin, file = "sets.RData")
