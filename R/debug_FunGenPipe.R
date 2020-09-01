if (FALSE) {
    outPath <- file.path(resultDirOut, "3.PGSEA.KEGGREST.limma-comparisons")
    gscPGSEA<-gscKeggrestRno
    esetsPGSEA<-esetsKeggrest
    fitsPGSEA<-fitsKeggrest
    comparisons <- comparisons
    setIDCol="KeggID"
    useNameCol="keggPathNameRno2"
    fitsProbes<-fits
    annsProbes<-annPrb
    pVals <- pVals
    targetsOrder <- targetsOrder
    pValDE=0.05
    maxHeatMapProbes=50
    compName <- "dH_S.tN.gM~dtgFstat"
}

    #' REMOVED: Input/Output: eset with Biobase::fData with 3 columns: PROBEID SYMBOL ENTREZID
    #' subset eset to ENTREZIDs (they also have Symbols)
    #' Collapse multiple annotations into ;-separated values based on PROBEID
    # esetAnnSEcollapsedSym <- function(eset) {
    #     require(BiocGenerics::annotation(eset), character.only=TRUE)
    #     stopifnot(dim(Biobase::fData(esetAnnSE))[[2]] == 3) # required Biobase::fData with 3 columns: PROBEID SYMBOL ENTREZID
    #     annSE <- Biobase::fData(eset)
    #     ## Annotations with ENTREZIDs (they also have Symbols)
    #     annSym <- annSE[!is.na(annSE$ENTREZID),]
    #     ## Collapse multiple annotations into ;-separated values based on PROBEID
    #     ## 1-to-1 relation with eset
    #     annSEcollapsed <- base::data.frame("PROBEID"=annSE[!duplicated(annSE$PROBEID),1])
    #     annSEcollapsed$PROBEID <- as.character(annSEcollapsed$PROBEID)
    #     rownames(annSEcollapsed) <- annSEcollapsed$PROBEID
    #     annSEcollapsed$"SYMBOLs" <- base::sapply(annSEcollapsed$PROBEID, FUN=function(pid, myAnnSE) {paste(sort(unique(as.character(myAnnSE[myAnnSE$PROBEID==pid,"SYMBOL"]))), collapse=";")}, annSE)
    #     annSEcollapsed$"ENTREZIDs" <- base::sapply(annSEcollapsed$PROBEID, FUN=function(pid, myAnnSE) {paste(sort(unique(as.character(myAnnSE[myAnnSE$PROBEID==pid,"ENTREZID"]))), collapse=";")}, annSE)
    # 
    #     ## SUBSET collapsed annotations to SYMBOLs: dim 16822,3
    #     annSEcollapsedSym <- subset(annSEcollapsed, SYMBOLs!="")
    #     esetSym <- subset(eset, Biobase::featureNames(eset) %in% as.character(annSEcollapsedSym$PROBEID)) #  16822 features, 30 samples
    # 
    #     ## MERGE data & annotations
    #     Biobase::fData(esetSym) <- annSEcollapsedSym[Biobase::featureNames(esetSym),]
    #     esetSym
    # }
