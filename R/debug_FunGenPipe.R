if (FALSE) {

## gr device
dev.off()
dev.list()
dev.new()

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
#' (targets <- tibble(color_aa=c(1,1,2), color_bb=c("3"=3,"4"=4,"4"=4), aa=c(11,11,22), bb=c(33,33,44), cc=c(55,66,55)))
#' names(targets$color_aa)
#' names(targets$color_bb)
#' getColPats2(targets)
#' # using namesFrom
#' getColPats2(targets, namesFrom=NULL)
#' getColPats2(targets, namesFrom=cc)
#' cc <- NULL
#' getColPats2(targets, namesFrom=cc)
#' d <- NULL
#' getColPats2(targets, namesFrom=d)
#' # using unique
#' getColPats2(targets, unique=TRUE)
#' getColPats2(targets, namesFrom=c, unique=TRUE)

# parameters
targets
names(targets$color_b)
colorsFrom="color_"
namesFrom=NULL
collapse=FALSE
getColPats2(targets, namesFrom=a)
map(getColPats2(targets), ~pluck(.)) #, unique(names(.))))
map(getColPats2(targets), ~pluck(unique(.)))
map(getColPats2(targets), ~pluck(unique(names(.))))
map(getColPats2(targets), ~keep(., !duplicated(names(.))))#, unique(names(.))))

## pheatmapTargets
xxx<-pheatmapTargets(Biobase::exprs(esetAnn), targets=pData(esetAnn) %>% select(-color_RIN), methods=c("manhattan", "euclidean"),
    colorsFrom="color_", namesFrom=NULL, filePath=pdf(file.path(resultDir, "_debug_pheatmap_norm.pdf")),
    width=7, height=7, treeheight_row = 0, legend=FALSE)
## parameters
expLog2=Biobase::exprs(esetAnn) 
targets <- pData(esetAnn) %>% select(-color_RIN)
methods=c("manhattan", "euclidean")
method <- methods[[1]]
colorsFrom="color_"
namesFrom=NULL
filePath=pdf(file.path(resultDir, "_debug_pheatmap_norm.pdf"))
treeheight_row = 0
annotation_row = as.data.frame(getColPats2(targets, colorsFrom=colorsFrom, namesFrom=!!namesFrom))
names(annotation_row)
annotation_colors = getColPats2(targets, unique=TRUE, colorsFrom=colorsFrom, namesFrom=!!namesFrom)
names(annotation_colors)

#         annotation_for_heatmap <- data.frame(Sex = pData(esetAnn)$Sex, Fs=pData(esetAnn)$Fs)
#         row.names(annotation_for_heatmap) <- row.names(pData(esetAnn))
#         distsM <- as.matrix(dist(t(Biobase::exprs(esetAnn)), method = "manhattan"))
#         distsE <- as.matrix(dist(t(Biobase::exprs(esetAnn)), method = "euclidean"))
#         rownames(distsM) <- row.names(pData(esetAnn))
#         rownames(distsE) <- row.names(pData(esetAnn))
#         hmcol <- colorRampPalette(RColorBrewer::brewer.pal(9, "YlOrRd"))(255)
#         colnames(distsM) <- NULL
#         colnames(distsE) <- NULL
#         diag(distsM) <- NA
#         diag(distsE) <- NA
#         pdf(file.path(curDir, "pheatmap_norm.pdf"))
#         pheatmap(distsM, col = (hmcol), 
#            annotation_row = annotation_for_heatmap,
#            annotation_colors = getColPats2(pData(esetAnn), unique=TRUE),
#            legend = TRUE, 
#            treeheight_row = 0,
#            legend_breaks = c(min(distsM, na.rm = TRUE), 
#                max(distsM, na.rm = TRUE)), 
#            legend_labels = (c("similar", "diverse")),
#            main = "Manhattan heatmap for the calibrated samples")
#         pheatmap(distsE, col = (hmcol), 
#            annotation_row = annotation_for_heatmap,
#            annotation_colors = getColPats2(pData(esetAnn), unique=TRUE),
#            legend = TRUE, 
#            treeheight_row = 0,
#            legend_breaks = c(min(distsE, na.rm = TRUE), 
#                max(distsE, na.rm = TRUE)), 
#            legend_labels = (c("similar", "diverse")),
#            main = "Euclidean heatmap for the calibrated samples")
#         dev.off() 
# }

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
}
