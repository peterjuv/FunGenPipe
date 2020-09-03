if (FALSE) {

## getWriteHeatmap_PgseaTTcomparisons
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

## getColPats
aaa<-getColPats(eset)
for (aa in aaa) print(names(aa))

## getColPats2
ccc<-getColPats2(targets, HybName)
ccc
for (cc in ccc) print(names(cc))
bbb<-getColPats2(targets)
bbb
for (bb in bbb) print(names(bb))

## getColPats2 extended parameters
# default parameters
targets <- tibble(color_a=c(1,2,3), color_b=c("4"=4,"5"=5,"6"=6), a=c(11,22,33), b=c(44,55,66), c=c(77,88,99))
names(targets$color_a)
names(targets$color_c)
getColPats2(targets)
names(getColPats2(targets)$color_a)
names(getColPats2(targets)$color_b)
# using namesFrom
getColPats2(targets, namesFrom=c)
names(getColPats2(targets, namesFrom=c)$color_a)
names(getColPats2(targets, namesFrom=c)$color_b)
# using collapse=TRUE
targets
names(targets$color_b)
colorsFrom="color_"
namesFrom=NULL
collapse=FALSE
getColPats2(targets, namesFrom=a)
map(getColPats2(targets), ~pluck(.))#, unique(names(.))))
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
