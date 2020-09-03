## 2020-07-17 Function Genomics pipeline by peter.juvan@gmail.com
##

#' Get color patterns with names from another column from ExpressionSet slot phenoData.
#' 
#' Returns columns from phenoData that start with "color_" and have their suffix in common with another column;
#' Names of colors are set from the associated column, e.g., list(color_Sex = setNames(color_Sex, Sex), ...)
#' @param eset
#' @return List of colors, each of length equal to the number of rows in phenoData
#' @example
#' /donotrun{
#' getColPats(Biobase::pData(dataRaw/eset))
#' getColPats(Biobase::pData(eset))
#' }
#' @section Old implementation:
#' getColPats <- function(eset) {
#'     require(assertthat)
#'     require(Biobase)
#'     require(tibble)
#'     assertthat::has_attr(eset, "phenoData")
#'     colPats <- list()
#'     aNames <- colnames(Biobase::pData(eset))
#'     cNames <- aNames[grep("^color_", aNames, perl=TRUE)]
#'     fNames <- sub("^color_", "", cNames)
#'     for (i in seq_len(length(fNames))) {
#'         colPats[[fNames[i]]] <- setNames(Biobase::pData(eset)[[cNames[i]]], Biobase::pData(eset)[[fNames[i]]])
#'     }
#'     colPats
#' }
#' @export
getColPats <- function(eset) {
    require(dplyr)
    require(assertthat)
    require(Biobase)
    require(rlang)
    require(purrr)
    assertthat::has_attr(eset, "phenoData")
    cols <- Biobase::pData(eset) %>% 
        as_tibble %>% 
        select(starts_with("color_"))
    names <- Biobase::pData(eset) %>% 
        as_tibble %>%
        select(all_of(sub("color_", "", colnames(cols))))
    map2(cols, names, setNames)
}


#' Get color patterns with names from another column from targets
#' 
#' New implementation allowing for alternative naming using parameter namesFrom
#' Returns columns from phenoData that start with "color_" and have their suffix in common with another column;
#' Names of colors are set from the associated column, e.g., tibble(color_Sex = setNames(color_Sex, Sex), ...)
#' @param targets Table with columns names with a prefix colorsFrom='color_'
#' @param colorsFrom String prefix of names of columns with colors, default 'color_'
#' @param namesFrom Either NULL for using variables that share suffixes with colors
#'                  or NA for using preset names of colors
#'                  or a variable from targets that is used for setting names to colors
#' @param unique If TRUE, returns unique names and associated colors;
#'               Note that colors may be dropped, see the last example
#' @return List of colors of length equal to the number of rows;
#'         If unique=TRUE, each list is truncated to unique names
#' @examples
#' \donotrun{
#' # default parameters
#' targets <- tibble(color_a=c(1,1,2), color_b=c("3"=3,"4"=4,"4"=4), a=c(11,11,22), b=c(33,33,44), c=c(55,66,55))
#' names(targets$color_a)
#' names(targets$color_b)
#' getColPats2(targets)
#' # using namesFrom
#' getColPats2(targets, namesFrom=c)
#' # using unique
#' getColPats2(targets, unique=TRUE)
#' getColPats2(targets, namesFrom=c, unique=TRUE)
#' }
#' @export
getColPats2 <- function(targets, colorsFrom="color_", namesFrom=NULL, unique=FALSE) {
    require(magrittr)
    require(dplyr)
    require(assertthat)
    require(rlang)
    require(purrr)
    targets %<>% as_tibble
    namesFrom <- enquo(namesFrom)
    # TOFIX assertthat::assert_that(quo_is_null(namesFrom) | assertthat::has_name(targets, !!namesFrom))
    cols <- targets %>% 
        select(starts_with(colorsFrom))
    # case_when(
    #     quo_is_null(namesFrom) ~ map2(cols, targets %>% select(all_of(sub(colorsFrom, "", colnames(cols)))), setNames),
    #     #quo_is_na(namesFrom) ~ map(cols, identity),
    #     !quo_is_null(namesFrom) ~ map2(cols, targets %>% select(!!namesFrom), setNames)
    # ) %>% as_tibble
    if (quo_is_null(namesFrom))
        ncols <- map2(cols, targets %>% select(all_of(sub(colorsFrom, "", colnames(cols)))), setNames)
    else
        ncols <- map2(cols, targets %>% select(!!namesFrom), setNames)
    if (unique)
        ncols %<>% map(~keep(., !duplicated(names(.))))
    ncols
}


#' Boxplot RLE (Relative Log Expression) and correlate to RIN
#' 
#' Correlate RLE to median and IQR of RIN.
#' If RIN is not giver, correlate it to straight line, i.e. test if RLEs are equal
#' 
#' @param expLog2 Expression matrix with genes in rows and samples in columns; 
#'  columns must be named
#' @param filePath If given, output a PDF
#' @param RIN If given, correlate it to median and IQR of RLE; 
#'  otherwise correlate RLE to a straight line (RLE==1 for all samples)
#' @param width PDF width
#' @param height PDF height
#' @param ... Passed to geom_boxplot
#' @return ggplot2 object and PDF if filePath is given
#' @export
boxplotRLE <- function(expLog2, filePath=NULL, RIN=NULL, width=7, height=7, ...) {
    require(Biobase)
    require(tidyverse)
    require(magrittr)
    require(RColorBrewer)
    require(assertthat)
    require(myHelpers)
    assertthat::assert_that(is.null(RIN) | assertthat::are_equal(length(RIN), dim(expLog2)[[2]]))
    RLE_row_medians <- Biobase::rowMedians(as.matrix(expLog2))
    RLE_data <- as.data.frame(sweep(expLog2, 1, RLE_row_medians))
    RLE_data_long <- RLE_data %>% 
        tidyr::pivot_longer(cols=everything(), names_to="Sample", values_to="log2_expression_deviation")
    if (!is.null(RIN)) {
        cor_RLEmed_RIN <- cor(sapply(RLE_data, median, na.rm=TRUE), RIN, use="na.or.complete")
        cor_RLEiqr_RIN <- cor(sapply(RLE_data, IQR, na.rm=TRUE), RIN, use="na.or.complete")
        myCols <- myHelpers::brewPalCont(RIN, n=9, name="OrRd", digits=2)
    } else {
        myCols <- myHelpers::brewPalCont(rep(5,dim(expLog2)[[2]]), n=9, name="OrRd", digits=0)
    }
    prle <- ggplot2::ggplot(RLE_data_long, aes(Sample, log2_expression_deviation)) + 
        geom_boxplot(outlier.shape = NA, alpha=0.3, fill=myCols, ...) +
        scale_fill_gradient() +
        ylim(c(-2, 2)) +
        theme(axis.text.x = element_text(colour = "aquamarine4", angle = 60, size = 6.5, hjust = 1, face = "bold"),
              plot.caption = element_text(hjust=0.5))
    if (!is.null(RIN)) prle <- prle + 
        labs(caption=paste("Correlation between RLE (median, IQR) and RIN:", format(cor_RLEmed_RIN), ",", format(cor_RLEiqr_RIN)))
    if (!is.null(filePath)) {
        pdf(filePath, width=width, height=height)
        print(prle)
        dev.off()
    }
    prle
}

#' Plot PCA for gene expression and phenotype data
#' 
#' @param expLog2 Gene expression matrix in log2 scale
#' @param targets Phenotype data with columns shape, color, fill, size, 
#'  as well as color_color and color_fill
#' @param shape Variable from targets for ggplot2::aes
#' @param color Variable from targets for ggplot2::aes
#' @param fill Variable from targets for ggplot2::aes
#' @param size Variable from targets for ggplot2::aes
#' @param scale_color Variable from targets, named vector for scale_color_manual(value=scale_color)
#' @param scale_fill  Variable from targets, named vector for scale_fill_manual(value=scale_fill)
#' @param filePath PDF
#' @param width PDF width
#' @param height PDF height
#' @param stroke Line thickness, passed to geom_point
#' @param guides_fill Default "none"; use "legend" to display it
#' @param ... Passed to geom_point
#' @return ggplot2 object and PDF if filePath is given
#' @section TODO:
#'  FIX: add return value
#'  FIX: Discrete scale_color supplied for continuous color variable; scale_color not used
#'  Fix: case that scale_color and scale_fill are not named vectors or NULL
#'  Fix scale_size: Error: Discrete value supplied to continuous scale
#'  Use scale_colour_viridis_c
#' @section Implementation:
#' Put ggplot() inside {} to be able to access data using dot (.)
#' @export
plotPCAtargets <- function(expLog2, targets, shape, color, fill, size,
    scale_color = NULL, scale_fill = NULL, filePath=NULL, width=7, height=7, stroke=1, 
    guides_fill = "none", ...) {
    require(rlang)
    require(ggplot2)
    require(assertthat)
    require(magrittr)
    assertthat::are_equal(dim(expLog2)[[2]], dim(targets)[[1]])
    assertthat::are_equal(colnames(expLog2), rownames(targets))
    ## enquo arguments
    shape <- enquo(shape); color <- enquo(color); fill <- enquo(fill); size <- enquo(size)
    scale_color <- enquo(scale_color); scale_fill <- enquo(scale_fill)
    ## PCA
    PCA <- prcomp(t(expLog2), scale = FALSE)
    percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
    sd_ratio <- sqrt(percentVar[2] / percentVar[1])
    ## plot
    p <- targets %>% 
        as_tibble %>% 
        mutate(PC1 = PCA$x[,1], PC2 = PCA$x[,2]) %>% 
    {
        ggplot(., aes(PC2, PC1)) +
        geom_point(aes(shape=!!shape, color=!!color, fill=!!fill, size=!!size), stroke=stroke) +
        ggtitle("PCA plot of log2 expression data") +
        ylab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
        xlab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
        theme(plot.title = element_text(hjust = 0.5))+
        coord_fixed(ratio = sd_ratio) +
        scale_shape_manual(values = c(21:25, 0:14)[1:length(unique(pull(., !!shape)))]) +
        guides(fill = guides_fill)
        # scale_size(range = c(2,5)) 
    }
    if (!quo_is_null(scale_color))
        if (!is.numeric(pull(p$data, !!color)))
            p %<>% { . + scale_color_manual(values = getColPats2(pull(.$data, !!scale_color), unique=TRUE))}
        else
            warning("Discrete scale_color supplied for continuous color variable; scale_color not used")
    if (!quo_is_null(scale_fill))  p %<>% { . + scale_fill_manual(values = getColPats2(pull(.$data, !!scale_fill), unique=TRUE))}
    #if (!quo_is_null(scale_fill))  p %<>% { . + scale_fill_manual(pull(.$data, !!scale_fill))}
    ## PDF
    if (!is.null(filePath)) {
        pdf(filePath, width=width, height=height)
        print(p)
        dev.off()
    }
    ## return
    #p
}


##############################
#### From Osijek_HFHSD.R #####
#### For KBlag_doxy_Cla.R ####
##############################

#' Add collapsed annotations to eset using columns SYMBOLSs, ENTREZIDs
#' 
#' Keep probes w/o annotations
#' 
#' @param eset
#' @return eset with collapsed annotations in slot annotations(eset)
#' @export
esetAnnSEs <- function(eset, annDb=clariomsmousetranscriptcluster.db) {
    require(annDb$packageName, character.only=TRUE)
    annSE <- AnnotationDbi::select(annDb, keys=Biobase::featureNames(eset), columns=c("SYMBOL","ENTREZID"), keytype="PROBEID")
    ## MERGE data & annotations
    ## The number of rows in featureData must match the number of rows in assayData. 
    ## Row names of featureData must match row names of the matrix / matricies in assayData
    #Biobase::fData(eset) <- annSE
    annSEcollapsed <- base::data.frame(PROBEID=annSE[!duplicated(annSE$PROBEID),1])
    annSEcollapsed$"PROBEID" <- as.character(annSEcollapsed$PROBEID) # PROBEID: factor -> char
    rownames(annSEcollapsed) <- annSEcollapsed$"PROBEID"
    annSEcollapsed$"SYMBOLs" <- base::sapply(annSEcollapsed$"PROBEID", FUN=function(pid, myAnnSE) {
        paste(sort(unique(as.character(myAnnSE[myAnnSE$PROBEID==pid,"SYMBOL"]))), collapse=";")
    }, annSE)
    annSEcollapsed$"ENTREZIDs" <- base::sapply(annSEcollapsed$"PROBEID", FUN=function(pid, myAnnSE) {
        paste(sort(unique(as.character(myAnnSE[myAnnSE$PROBEID==pid,"ENTREZID"]))), collapse=";")
    }, annSE)
    Biobase::fData(eset) <- annSEcollapsed[Biobase::featureNames(eset),] # ensure ordered (not really needed)
    BiocGenerics::annotation(eset) <- annDb$packageName
    eset
}


#' Write DE genes from contrasts & return DE genes
#' 
#' @param outPath Path to write TTS
#' @param esets esets[[fitName]]: a list of esets for BiocGenerics::annotation of probes with columns PROBEID, SYMBOL, ENTREZID, HSA_ENTREZID
#' @param fits fits[[fitName]]: a list of PGSEA fits
#' @param contMatrices contMatrices[[fitName]]: a list of contrast matrices with Levels & Contrasts
#' @returns A list TTs of enriched sets: annTT[[fitName]]
#' @rdname getWriteHeatmap_probesTT
#' @export
getWriteHeatmap_probesTTcontrasts <- function(outPath, esets, fits, contMatrices, pVals,
                                              targetsOrder=NA, pValDE=0.05, maxHeatMapProbes=50) {
    require(gplots)
    require(Biobase)
    require(limma)
    if(!file.exists(outPath)) dir.create(outPath)
    annTT <- list()
    colPats <- list()
    for (fitName in names(fits)) {
        colPats[[fitName]] <- getColPats(esets[[fitName]])
        if (is.na(targetsOrder)) targetsOrder <- 1:dim(Biobase::pData(esets[[fitName]]))[[1]]
        ## WRITE DE genes for contrasts + draw HeatMaps
        for(contr in colnames(contMatrices[[fitName]])) {
            ## DE: write all genes
            tt <- limma::topTable(fits[[fitName]], coef=contr, number=Inf)
            write.table(format(tt,digits=3), file.path(outPath, paste(fitName, contr, "txt", sep=".")), sep="\t", quote=FALSE, row.names=FALSE)
            ## Venn: store DE genes + HSA_ENTREZID annotations
            ttVenn <- limma::topTable(fits[[fitName]], coef=contr, number=Inf, p.value=pValDE)
            #if ("HSA_ENTREZIDs" %in% colnames(ttVenn)) ttVenn <- dplyr::select(ttVenn, -HSA_ENTREZIDs)
            ttVenn$PROBEID <- as.character(ttVenn$PROBEID)
            annTT[[fitName]][[contr]] <- dplyr::inner_join(Biobase::fData(esets[[fitName]]), ttVenn, by = "PROBEID")
            ## HeatMap
            heatMapWithMaxProbesDrawn <- F
            for (pVal in pVals) {
                if (!heatMapWithMaxProbesDrawn) {
                    ## select genes, order samples
                    eSubSetAdj <- esets[[fitName]][rownames(tt[tt$adj.P.Val < pVal, ]), targetsOrder]
                    if (dim(eSubSetAdj)[[1]] > maxHeatMapProbes) {
                        eSubSetAdj <- eSubSetAdj[1:maxHeatMapProbes,] # Note: tt needs to be SORTED by adj.P.Vals !!!
                        heatMapWithMaxProbesDrawn <- T
                    }
                    if (dim(eSubSetAdj)[[1]] == 1) eSubSetAdj <- eSubSetAdj[c(1,1),]    # Note: in case of one DE probe, duplicate that probe in order to draw heatmap
                    if (dim(eSubSetAdj)[[1]] > 1) {
                        if (heatMapWithMaxProbesDrawn) {pdfFileName <- paste(fitName, contr, "fdr", pVal, "top", maxHeatMapProbes, "pdf", sep=".")} else {pdfFileName <- paste(fitName, contr, "fdr", pVal, "pdf", sep=".")}
                            pdf(file.path(outPath, pdfFileName))
                            ## order samples $ colPats using dendrogram by setting constraints == weights / no reordering
                            gplots::heatmap.2(Biobase::exprs(eSubSetAdj), col=topo.colors(100), ColSideColors=colPats[[fitName]][targetsOrder], labRow=unlist(mget(Biobase::featureNames(eSubSetAdj), Biobase::fData(eSubSetAdj)$SYMBOLs)), margins=c(8,5), Colv=targetsOrder)
                            gplots::heatmap.2(Biobase::exprs(eSubSetAdj), col=topo.colors(100), ColSideColors=colPats[[fitName]][targetsOrder], labRow=unlist(mget(Biobase::featureNames(eSubSetAdj), Biobase::fData(eSubSetAdj)$SYMBOLs)), margins=c(8,5), Colv=NULL, dendrogram='row')
                            heatmap(Biobase::exprs(eSubSetAdj), col=topo.colors(100), ColSideColors=colPats[[fitName]][targetsOrder], labRow=unlist(mget(Biobase::featureNames(eSubSetAdj), Biobase::fData(eSubSetAdj)$SYMBOLs)), margins=c(8,5), Colv=targetsOrder)
                            heatmap(Biobase::exprs(eSubSetAdj), col=topo.colors(100), ColSideColors=colPats[[fitName]][targetsOrder], labRow=unlist(mget(Biobase::featureNames(eSubSetAdj), Biobase::fData(eSubSetAdj)$SYMBOLs)), margins=c(8,5), Colv=NA)
                            dev.off()
                        } else print(paste("WARNING: no DE probes found. Fit", fitName, ", contr", contr, ", pVal", pVal))
                } # end if heatMapWithMaxProbesDrawn
            } # end pVals
            } #end contr
    } #end fits
    ## RETURN
    annTT
} #end function
        

#' Write DE genes from comparisons & return DE genes
#' 
#' @param comparisons comparisons[[fitName]]: a list of comparisons made from names of contrasts
#' @rdname getWriteHeatmap_probesTT
#' @export
getWriteHeatmap_probesTTcomparisons <- function(outPath, esets, fits, comparisons, pVals,
                                                targetsOrder=NA, pValDE=0.05, maxHeatMapProbes=50) {
    require(gplots)
    require(Biobase)
    require(limma)
    if(!file.exists(outPath)) dir.create(outPath)
    annTT <- list()
    colPats <- list()
    for (fitName in names(fits)) {
        colPats[[fitName]] <- getColPats(esets[[fitName]])
        if (is.na(targetsOrder)) targetsOrder <- 1:dim(Biobase::pData(esets[[fitName]]))[[1]]
        ## WRITE DE genes for comparisons + draw HeatMaps
        for(compName in names(comparisons[[fitName]])) {
            ## make intersection of DE genes from contrasts within coefList
            coefList <- comparisons[[fitName]][[compName]]  # [[1]][1] "sM_F.id1_0"; $sM_F.id1234_0 [1] "sM_F.id1_0" "sM_F.id2_0" "sM_F.id3_0" "sM_F.id4_0"
            ## make names for coefList
            names(coefList) <- ifelse(grepl("^$", names(coefList), perl=TRUE), make.names(coefList), names(coefList))
            ## ttComp: DE genes for a comparison compName: store intersection of DE genes from tt
            ttComp <- NULL
            ## intersertion of DE genes from a list of coeficients (coefList)
            for (coefs in coefList) {
                ## DE
                ttDE <- limma::topTable(fits[[fitName]], coef=coefs, number=Inf, p.value=pValDE)
                if (is.null(ttComp)) {ttComp <- ttDE} else {ttComp <- ttComp[ttComp$PROBEID %in% ttDE$PROBEID, ]}
            }
            #if ("HSA_ENTREZIDs" %in% colnames(ttComp)) ttComp <- dplyr::select(ttComp, -HSA_ENTREZIDs)
            write.table(format(ttComp,digits=3), file.path(outPath, paste(fitName, compName, "txt", sep=".")), sep="\t", quote=FALSE, row.names=FALSE)
            ## Store DE genes
            ttComp$PROBEID <- as.character(ttComp$PROBEID)
            annTT[[fitName]][[compName]] <- dplyr::inner_join(Biobase::fData(esets[[fitName]]), ttComp, by = "PROBEID")
            ## HeatMap
            heatMapWithMaxProbesDrawn <- F
            for (pVal in pVals) {
                if (!heatMapWithMaxProbesDrawn) {
                    ## select genes, order samples
                    eSubSetAdj <- esets[[fitName]][rownames(ttComp), targetsOrder]
                    if (dim(eSubSetAdj)[[1]] > maxHeatMapProbes) {
                        eSubSetAdj <- eSubSetAdj[1:maxHeatMapProbes,] # Note: ttComp needs to be SORTED by adj.P.Vals !!!
                        heatMapWithMaxProbesDrawn <- T
                    }
                    if (dim(eSubSetAdj)[[1]] == 1) eSubSetAdj <- eSubSetAdj[c(1,1),]    # Note: in case of one DE probe, duplicate that probe in order to draw heatmap
                    if (dim(eSubSetAdj)[[1]] > 1) {
                        if (heatMapWithMaxProbesDrawn) {pdfFileName <- paste(fitName, compName, "fdr", pVal, "top", maxHeatMapProbes, "pdf", sep=".")} else {pdfFileName <- paste(fitName, compName, "fdr", pVal, "pdf", sep=".")}
                        pdf(file.path(outPath, pdfFileName))
                        ## order samples $ colPats using dendrogram by setting constraints == weights / no reordering
                        gplots::heatmap.2(Biobase::exprs(eSubSetAdj), col=topo.colors(100), ColSideColors=colPats[[fitName]][targetsOrder], labRow=unlist(mget(Biobase::featureNames(eSubSetAdj), Biobase::fData(eSubSetAdj)$SYMBOLs)), margins=c(8,5), Colv=targetsOrder)
                        gplots::heatmap.2(Biobase::exprs(eSubSetAdj), col=topo.colors(100), ColSideColors=colPats[[fitName]][targetsOrder], labRow=unlist(mget(Biobase::featureNames(eSubSetAdj), Biobase::fData(eSubSetAdj)$SYMBOLs)), margins=c(8,5), Colv=NULL, dendrogram='row')
                        heatmap(Biobase::exprs(eSubSetAdj), col=topo.colors(100), ColSideColors=colPats[[fitName]][targetsOrder], labRow=unlist(mget(Biobase::featureNames(eSubSetAdj), Biobase::fData(eSubSetAdj)$SYMBOLs)), margins=c(8,5), Colv=targetsOrder)
                        heatmap(Biobase::exprs(eSubSetAdj), col=topo.colors(100), ColSideColors=colPats[[fitName]][targetsOrder], labRow=unlist(mget(Biobase::featureNames(eSubSetAdj), Biobase::fData(eSubSetAdj)$SYMBOLs)), margins=c(8,5), Colv=NA)
                        dev.off()
                    } else print(paste("WARNING: no DE probes found. Fit", fitName, ", comparison", compName, ", pVal", pVal))
                } # end heatMapWithMaxProbesDrawn
            } # end pVals
        } #end compName (comparisons)
    } #end fits
  ## RETURN
  annTT
} #end function



#' Write enriched patways from contrasts & return patways
#' 
#' @param outPath Path to write TTS
#' @param gscPGSEA gscPGSEA[[fitName]]: a list of gene set collection
#' @param esetsPGSEA[[fitName]]: a list of PGSEA esets
#' @param fitsPGSEA[[fitName]]: a list of PGSEA fits
#' @param contMatrices contMatrices[[fitName]]: a list of contrast matrices with Levels & Contrasts
#' @param setIDCol: add a column to TT with set IDs and name it setIDCol
#' @param setName: name of the column from Biobase::fData(esetsPGSEA[[?]]) to show with heatmaps
#' @param fitsProbes[[fitName]]: a list of probe fits
#' @param esetsProbes[[fitName]]: a list of esets for BiocGenerics::annotation of probes with columns PROBEID, SYMBOL, ENTREZID, HSA_ENTREZID
#' @return a list TTs of enriched sets: annTT[[fitName]]
#' @examples
#' \dontrun{
#' getWriteHeatmap_PgseaTTcontrasts(outPath=file.path(resultDirOut, "3.PGSEA.KEGGREST.limma-contrasts"), gscPGSEA=gscKeggrestHsa, esetsPGSEA=esetsKeggrest, fitsPGSEA=fitsKeggrest,
#'                                      setIDCol="KeggID", useNameCol="keggPathNameHsa2", fitsProbes=fits, esetsProbes=esets, pVals)
#' getWriteHeatmap_PgseaTTcontrasts(file.path(resultDirOut, "3.PGSEA.TRANSFAC.2016.1.byFA.limma-contrasts"), gscTF.facFA, esetsTF.facFA, fitsTF.facFA,
#'                                          setIDCol="facFA", useNameCol="factorFA", fits, esets, pVals)
#' }
#' @rdname getWriteHeatmap_PgseaTT
#' @export
getWriteHeatmap_PgseaTTcontrasts <- function(outPath, gscPGSEA, esetsPGSEA, fitsPGSEA, contMatrices, setIDCol, useNameCol,
                                             fitsProbes, esetsProbes, pVals, targetsOrder=NA, pValDE=0.05, maxHeatMapProbes=50) {
    require(gplots)
    require(Biobase)
    require(limma)
    if(!file.exists(outPath)) dir.create(outPath)
    annTT <- list()
    colPats <- list()
    for (fitName in names(fitsPGSEA)) {
        colPats[[fitName]] <- getColPats(esetsPGSEA[[fitName]])
        if (is.na(targetsOrder)) targetsOrder <- 1:dim(Biobase::pData(esetsProbes[[fitName]]))[[1]]
        for(contr in colnames(contMatrices[[fitName]])) {
            ttContr <- limma::topTable(fitsPGSEA[[fitName]], coef=contr, number=999999)
            ## TODO: fix tables for TF, identical columns setIDCol & factorFA
            ttContr <- cbind(setIDCol=rownames(ttContr), ttContr)
            colnames(ttContr)[1] <- setIDCol
            ## add DE genes
            if (fitName %in% names(fitsPGSEA)) {
                ttDE <- limma::topTable(fitsProbes[[fitName]], coef=contr, number=999999, p.value=pValDE)
                annOut <- t(sapply(rownames(ttContr), FUN=function(setId){
                    probeIDlst <- GSEABase::geneIds(gscPGSEA[[fitName]][[setId]]);
                    probeIDlst_DE <- intersect(probeIDlst, rownames(ttDE));
                    av <- Biobase::fData(esetsProbes[[fitName]])[Biobase::fData(esetsProbes[[fitName]])$PROBEID %in% probeIDlst_DE,];
                        c("ProbeIDsDE"=     paste(probeIDlst_DE, collapse=";"),
                        "SymbolsDE"=      paste(sort(unique(av$SYMBOL)), collapse=";"),
                        "EntrezIDsDE"=    paste(sort(unique(av$ENTREZID)), collapse=";"),
                        #"HSA_EntrezIDsDE"=paste(sort(unique(av$HSA_ENTREZID)), collapse=";"),
                        "ProbeCountDE"=   length(probeIDlst_DE),
                        "ProbeCountAll"=  length(probeIDlst),
                        "ProbeRatioDE"= length(probeIDlst_DE) / length(probeIDlst))}))
                ttContr <- cbind(ttContr, annOut)
            }
            write.table(ttContr, file.path(outPath, paste(fitName, contr, "txt", sep=".")), sep="\t", quote=FALSE, row.names=FALSE)
            ## Venn: store enriched HSA pathways
            annTT[[fitName]][[contr]] <- ttContr[ttContr$adj.P.Val < pValDE & !is.na(ttContr$adj.P.Val),]
            ## HeatMap
            heatMapWithMaxProbesDrawn <- F
            for (pVal in pVals) {
                if (!heatMapWithMaxProbesDrawn) {
                    esetSubAdj <- esetsPGSEA[[fitName]][rownames(ttContr[ttContr$adj.P.Val < pVal & !is.na(ttContr$adj.P.Val), ]), targetsOrder]
                    if (dim(esetSubAdj)[[1]] > maxHeatMapProbes) {
                        esetSubAdj <- esetSubAdj[1:maxHeatMapProbes,] # Note: tt needs to be SORTED by adj.P.Vals !!!
                        heatMapWithMaxProbesDrawn <- T
                    }
                    if (dim(esetSubAdj)[[1]] == 1) esetSubAdj <- esetSubAdj[c(1,1),]    # Note: in case of one DE probe, duplicate that probe in order to draw heatmap
                    if (dim(esetSubAdj)[[1]] > 1) {
                        if (heatMapWithMaxProbesDrawn) {pdfFileName <- paste(fitName, contr, "fdr", pVal, "top", maxHeatMapProbes, "pdf", sep=".")} else {pdfFileName <- paste(fitName, contr, "fdr", pVal, "pdf", sep=".")}
                        pdf(file.path(outPath, pdfFileName))
                        gplots::heatmap.2(Biobase::exprs(esetSubAdj), col=topo.colors(100), ColSideColors=colPats[[fitName]][targetsOrder], labRow=Biobase::fData(esetSubAdj)[,useNameCol], margins=c(8,8), Colv=targetsOrder)
                        gplots::heatmap.2(Biobase::exprs(esetSubAdj), col=topo.colors(100), ColSideColors=colPats[[fitName]][targetsOrder], labRow=Biobase::fData(esetSubAdj)[,useNameCol], margins=c(8,8), Colv=NULL, dendrogram='row')
                        heatmap(Biobase::exprs(esetSubAdj), col=topo.colors(100), ColSideColors=colPats[[fitName]][targetsOrder], labRow=Biobase::fData(esetSubAdj)[,useNameCol], margins=c(8,8), Colv=targetsOrder)
                        heatmap(Biobase::exprs(esetSubAdj), col=topo.colors(100), ColSideColors=colPats[[fitName]][targetsOrder], labRow=Biobase::fData(esetSubAdj)[,useNameCol], margins=c(8,8), Colv=NA)
                        dev.off()
                    } else print(paste("WARNING: no DE sets found. Fit", fitName, ", contr", contr, ", pVal", pVal))
                } # end heatMapWithMaxProbesDrawn
            } # end pVals
        } #end contr
    } #end fits
    ## RETURN
    annTT
} #end function


#' Write enriched patways from comparisons & return patways
#' 
#' @param comparisons comparisons[[fitName]]: a list of comparisons made from names of contrasts
#' @examples
#' \dontrun{
#' writeHeatmap_PgseaTTcomparisons(outPath=file.path(resultDirOut, "3.PGSEA.KEGGREST.limma-comparisons"), gscPGSEA=gscKeggrestHsa, esetsPGSEA=esetsKeggrest, fitsPGSEA=fitsKeggrest,
#'                           setIDCol="KeggID", useNameCol="keggPathNameHsa2", fitsProbes=fits, esetsProbes=esets, pVals, targetsOrder)
#' writeHeatmap_PgseaTTcomparisons(file.path(resultDirOut, "3.PGSEA.TRANSFAC.2016.1.byFA.limma-comparisons"), gscTF.facFA, esetsTF.facFA, fitsTF.facFA,
#'                           setIDCol="facFA", useNameCol="factorFA", fits, esets, pVals, targetsOrder)
#' }
#' @rdname getWriteHeatmap_PgseaTT
#' @export
getWriteHeatmap_PgseaTTcomparisons <- function(outPath, gscPGSEA, esetsPGSEA, fitsPGSEA, comparisons, setIDCol, useNameCol,
                                               fitsProbes, esetsProbes, pVals, targetsOrder=NA, pValDE=0.05, maxHeatMapProbes=50) {
    require(gplots)
    require(Biobase)
    require(limma)
    pVals <- pVals[pVals<=pValDE]   # pVals must be <= to pValDE since we generate limma::topTable for Enriched sets with adj.p.val <= pValDE
    if(!file.exists(outPath)) dir.create(outPath)
    annTT <- list()
    colPats <- list()
    for (fitName in names(fitsPGSEA)) {
        colPats[[fitName]] <- getColPats(gscPGSEA[[fitName]])
        if (is.na(targetsOrder)) targetsOrder <- 1:dim(Biobase::pData(esetsProbes[[fitName]]))[[1]]
        for(compName in names(comparisons[[fitName]])) {
            ## make intersection of DE from contrasts within coefList
            coefList <- comparisons[[fitName]][[compName]]  # [[1]][1] "sM_F.id1_0"; $sM_F.id1234_0 [1] "sM_F.id1_0" "sM_F.id2_0" "sM_F.id3_0" "sM_F.id4_0"
            ## make names for coefList
            names(coefList) <- ifelse(grepl("^$", names(coefList), perl=TRUE), make.names(coefList), names(coefList))
            ## ttComp: Enriched sets for a comparison compName: store intersection of Enriched sets from ttC1
            ttComp <- NULL
            for (coefs in coefList) {
                ## Enriched sets
                ttC1 <- limma::topTable(fitsPGSEA[[fitName]], coef=coefs, number=Inf, p.value=pValDE)
                ttC1 <- cbind(setIDCol=rownames(ttC1), ttC1)
                colnames(ttC1)[1] <- setIDCol
                if (is.null(ttComp)) {ttComp <- ttC1} else if (dim(ttComp)[[1]] > 0) {ttComp <- ttComp[ttComp[,setIDCol] %in% ttC1[,setIDCol], ]}
            }
            ## TODO: fix tables for TF, identical columns setIDCol & factorFA
            ## add DE genes
            if (fitName %in% names(fitsPGSEA) & dim(ttComp)[[1]] > 0) {
                ttDEComp <- NULL
                for (coefs in coefList) {
                    ## DE
                    ttDE <- limma::topTable(fitsProbes[[fitName]], coef=coefs, number=Inf, p.value=pValDE)
                    if (is.null(ttDEComp)) {ttDEComp <- ttDE} else {ttDEComp <- ttDEComp[ttDEComp$PROBEID %in% ttDE$PROBEID, ]}
                }
                #if ("HSA_ENTREZIDs" %in% colnames(ttDEComp)) ttDEComp <- dplyr::select(ttDEComp, -HSA_ENTREZIDs)
                annOut <- t(base::sapply(rownames(ttComp), FUN=function(setId) {
                    probeIDlst <- GSEABase::geneIds(gscPGSEA[[fitName]][[setId]]);
                    probeIDlst_DE <- intersect(probeIDlst, rownames(ttDEComp));
                    av <- Biobase::fData(esetsProbes[[fitName]])[Biobase::fData(esetsProbes[[fitName]])$PROBEID %in% probeIDlst_DE,];
                        c("ProbeIDsDE"=     paste(probeIDlst_DE, collapse=";"),
                        "SymbolsDE"=      paste(sort(unique(av$SYMBOL)), collapse=";"),
                        "EntrezIDsDE"=    paste(sort(unique(av$ENTREZID)), collapse=";"),
                        #"HSA_EntrezIDsDE"=paste(sort(unique(av$HSA_ENTREZID)), collapse=";"),
                        "ProbeCountDE"=   length(probeIDlst_DE),
                        "ProbeCountAll"=  length(probeIDlst),
                        "ProbeRatioDE"= length(probeIDlst_DE) / length(probeIDlst))}))
                ttComp <- cbind(ttComp, annOut)
            }
            write.table(ttComp, file.path(outPath, paste(fitName, compName, "txt", sep=".")), sep="\t", quote=FALSE, row.names=FALSE)
            ## Venn: store enriched pathways
            annTT[[fitName]][[compName]] <- ttComp[ttComp$adj.P.Val < pValDE & !is.na(ttComp$adj.P.Val),]
            ## HeatMap
            heatMapWithMaxProbesDrawn <- F
            for (pVal in pVals) {
                if (!heatMapWithMaxProbesDrawn) {
                    esetSubAdj <- esetsPGSEA[[fitName]][rownames(ttComp), targetsOrder]
                    if (dim(esetSubAdj)[[1]] > maxHeatMapProbes) {
                        esetSubAdj <- esetSubAdj[1:maxHeatMapProbes,] # Note: tt needs to be SORTED by adj.P.Vals !!!
                        heatMapWithMaxProbesDrawn <- T
                    }
                    if (dim(esetSubAdj)[[1]] == 1) esetSubAdj <- esetSubAdj[c(1,1),]    # Note: in case of one DE probe, duplicate that probe in order to draw heatmap
                    if (dim(esetSubAdj)[[1]] > 1) {
                        if (heatMapWithMaxProbesDrawn) {pdfFileName <- paste(fitName, compName, "fdr", pVal, "top", maxHeatMapProbes, "pdf", sep=".")} else {pdfFileName <- paste(fitName, compName, "fdr", pVal, "pdf", sep=".")}
                        pdf(file.path(outPath, pdfFileName))
                        gplots::heatmap.2(Biobase::exprs(esetSubAdj), col=topo.colors(100), ColSideColors=colPats[[fitName]][targetsOrder], labRow=Biobase::fData(esetSubAdj)[,useNameCol], margins=c(8,8), Colv=targetsOrder)
                        gplots::heatmap.2(Biobase::exprs(esetSubAdj), col=topo.colors(100), ColSideColors=colPats[[fitName]][targetsOrder], labRow=Biobase::fData(esetSubAdj)[,useNameCol], margins=c(8,8), Colv=NULL, dendrogram='row')
                        heatmap(Biobase::exprs(esetSubAdj), col=topo.colors(100), ColSideColors=colPats[[fitName]][targetsOrder], labRow=Biobase::fData(esetSubAdj)[,useNameCol], margins=c(8,8), Colv=targetsOrder)
                        heatmap(Biobase::exprs(esetSubAdj), col=topo.colors(100), ColSideColors=colPats[[fitName]][targetsOrder], labRow=Biobase::fData(esetSubAdj)[,useNameCol], margins=c(8,8), Colv=NA)
                        dev.off()
                    } else print(paste("WARNING: no DE sets found. Fit", fitName, ", comparison", compName, ", pVal", pVal))
                } # end heatMapWithMaxProbesDrawn
            } # end pVals
        } #end compName
    } #end fits
    ## RETURN
    annTT
} #end function



#' Get annotations for PGSEA from KEGGREST using ENTREZID.
#' 
# #' @return Data.frame with annotations and collapsed gene Symbols and ENTREZIDs
#' @export
annKeggEntrez <- function(organism="mmu", annPackage="org.Mm.eg.db") {
    require(KEGGREST)
    require(annPackage, character.only=TRUE)
    ## Get KEGG pathways for organism "org"
    keggPathOrg <- KEGGREST::keggLink("pathway", organism) # list (Named chr) org:Entrez_gene_ID -> pathway
    
    ## Transform KEGG Entrez gene IDs (org:###...###)to probe IDs (##...##) and aggregate by KEGG path IDs (path:org##...##)
    ## additionally, Remove "org:" in front of KEGG Entrez IDs and "path:" in front of KEGG path IDs
    keggPathOrgIDs <- unique(keggPathOrg) # names of paths @ organism
    keggPathOrg2EntrezIDs <- sapply(keggPathOrgIDs, function(pid) {sub(paste0(organism,":"), "", names(keggPathOrg[as.logical(keggPathOrg == pid)]))}) # list path:org04144: -> chr [1:231] "100017" "101056305" "103967" "105513" ...
    # remove "path:" in front of kegg path IDs
    names(keggPathOrg2EntrezIDs) <- sub("path:", "", names(keggPathOrg2EntrezIDs))
    
    ## create annotations of KEGG pathways
    ## get names of KEGG pathways: rename "path:map##...##" to "org##...##"
    keggPathNameMap <- KEGGREST::keggList("pathway") # list 530 paths, e.g. path:map00010 -> "Glycolysis / Gluconeogenesis"
    keggPathNameOrg <- keggPathNameMap
    names(keggPathNameOrg) <- sub("path:map", organism, names(keggPathNameOrg))
    keggPathNameOrg2 <- keggPathNameOrg[names(keggPathOrg2EntrezIDs)]
    
    ## create BiocGenerics::annotation data frame: path_id (rownames), path_name (keggPathNameOrg2), url, EntrezIDs, Symbols
    annKeggrestOrg <- as.data.frame(keggPathNameOrg2)
    annKeggrestOrg$url <- paste0("http://www.genome.jp/dbget-bin/www_bget?", rownames(annKeggrestOrg))
    ## TODO
    #annKeggrestOrg$"HSA_url" <- paste0("http://www.genome.jp/dbget-bin/www_bget?", sub(organism,"hsa",rownames(annKeggrestOrg)))
    
    ## add Entrez IDs and gene symbols
    annKeggrestOrg$EntrezIDs <- sapply(keggPathOrg2EntrezIDs[rownames(annKeggrestOrg)], paste, collapse=";")
    ## old: environment approach
    annPckgLs <- ls(paste0("package:",annPackage))
    annPckgEnv <- get(annPckgLs[grep("SYMBOL$", annPckgLs)])
    annKeggrestOrg$Symbols <- sapply(sapply(keggPathOrg2EntrezIDs[rownames(annKeggrestOrg)], function(entrezIDs) {unlist(mget(entrezIDs, annPckgEnv, ifnotfound=NA))}), paste, collapse=";")
    ## new but slow: select() and/or mapIds()
    #annKeggrestOrg$SYMBOLs <- sapply(sapply(keggPathOrg2EntrezIDs[rownames(annKeggrestOrg)], function(entrezIDs) {AnnotationDbi::mapIds(org.Mm.eg.db, keys=entrezIDs, column="SYMBOL", keytype="ENTREZID", multiVals=function(x) {paste0(x, collapse = ";")})}), paste, collapse=";")

    ## return
    annKeggrestOrg
}

#' Create a list of GeneSetCollection for a list of designs and probe annotations using ENTREZID.
#' 
#' @param organism
#' @param annPrb
#' @param annKeggEntrez
#' @param keggPathOrgIDs
#' @param design
#' @return GeneSetCollection
#' @export
gscKeggEntrez <- function(organism="mmu", annPrb, annKeggEntrez, keggPathOrgIDs, design) {
    require(GSEABase)
    require(KEGGREST)
    ## Get KEGG pathways for organism "org"
    keggPathOrg <- KEGGREST::keggLink("pathway", organism) # list (Named chr) org:Entrez_gene_ID -> pathway
    ## construct GeneSetCollection for a single design
    keggPathOrg2ProbeIDs <- list()
    gscKeggrestOrg <- list()
    for (dName in names(designs)) {
        keggPathOrg2ProbeIDs[[dName]] <- sapply(keggPathOrgIDs, function(pid) {
                                           unique(annPrb[[dName]][annPrb[[dName]]$ENTREZID %in% sub(paste0(organism,":"), "", names(keggPathOrg[as.logical(keggPathOrg == pid)])), "PROBEID"])
                                         }) # list of lists of path:org04144: -> chr [1:231] "17214142" "17215576" "17218733"  ...
        # remove "path:" in front of kegg path IDs
        names(keggPathOrg2ProbeIDs[[dName]]) <- sub("path:", "", names(keggPathOrg2ProbeIDs[[dName]]))
        ## construct GeneSetCollection
        gscKeggrestOrg[[dName]] <- GSEABase::GeneSetCollection(mapply(function(pIds, keggPathOrgId) {
                                    GSEABase::GeneSet(pIds, geneIdType=GSEABase::AnnotationIdentifier(), collectionType=GSEABase::KEGGCollection(), setName=keggPathOrgId)
                                    }, keggPathOrg2ProbeIDs[[dName]], names(keggPathOrg2ProbeIDs[[dName]])))
    }
    ## return
    gscKeggrestOrg
}

        # ## run PGSEA KEGG ##
        # pgKeggrestAll <- list()
        # esetsKeggrest <- list()
        # fitsKeggrest <- list()
        # for (dName in names(designs)) {
        #     # store original class, change it to "ExpressionSet"
        #     classEset <- class(esets[[dName]])
        #     class(esets[[dName]]) <- "ExpressionSet"
        #     pgKeggrestAll[[dName]] <- PGSEA(esets[[dName]], gscKeggrestRno[[dName]], ref=NULL, range=c(minConceptLen,  max(as.numeric(lapply(GSEABase::geneIds(gscKeggrestRno[[dName]]), length)))), center=FALSE, p.value=NA, weighted=TRUE, enforceRange=TRUE)
        #     # restore original class
        #     class(esets[[dName]]) <- classEset
        #     # new PGSEA eset
        #     esetsKeggrest[[dName]] <- new("ExpressionSet", Biobase::exprs=pgKeggrestAll[[dName]], phenoData=phenoData(esets[[dName]]), BiocGenerics::annotation="KEGGREST", experimentData=experimentData(esets[[dName]]))
        #     ## MERGE data & annotations
        #     Biobase::fData(esetsKeggrest[[dName]]) <- annKeggrestRno[Biobase::featureNames(esetsKeggrest[[dName]]),]
        #     ## fit
        #     fitsKeggrest[[dName]] <- eBayes(contrasts.fit(lmFit(esetsKeggrest[[dName]], designs[[dName]]), cont.matrices[[dName]]))
        # }



        # #######################################################################
        # ## Test enrichment of TRANSFAC gene sets sharing equal TFs using ENTREZ
        # ## 1. parse TRANSFAC using transfac-parse_v2.R
        # ## 2. load table 2019-07-05_13-30-13_fac2drDF_osRNentrez.RData
        # ## 3. subset to ENTREZ genes
        # #######################################################################

        # ## fac2drDF_entrez: import and filter TRANSFAC data.frame
        # ##  factor.AC  site.AC  factor.BSq  gene.AC  gene.DRac2  gene.DRdbName2  factorURL  siteURL  geneURL
        # load("~/CFGBC/transfac/2016.1/2019-07-05_13-30-13_fac2drDF_osRNentrez.RData")
        # dim(fac2drDF_osRNentrez)      # 3895    10

        # ## collect all ENTREZ IDs from esets and remove duplicates: from 29122, 32783, 29835, 24133 in annPrb to 22309 unique in entrezRnoEsets
        # entrezRnoEsets <- c()
        # for (nm in names(annPrb)) entrezRnoEsets <- union(entrezRnoEsets, annPrb[[nm]]$ENTREZID)
        # length(entrezRnoEsets)  # 16818

        # ## limit TF table (fac) to ENTREZ IDs from esets
        # geneDRac_entrez <- intersect(entrezRnoEsets, fac2drDF_osRNentrez$gene.DRac2)
        # length(geneDRac_entrez)   # 639
        # fac2drDF_entrez <- fac2drDF_osRNentrez[fac2drDF_osRNentrez$gene.DRac2 %in% geneDRac_entrez,]
        # dim(fac2drDF_entrez)      # 3839    10

        # ## remove factor.FA == "hsa.*" and "rno.*" from fac2drDF_entrez, retain thoese with Symbols
        # ## removed e.g.: rno-miR-23a-3p, hsa-let-7a-5p, rno-let-7a-5p, rno-miR-193a-3p, hsa-miR-127-3p
        # fac2drDF_entrez_sym <- fac2drDF_entrez[grep("^(?!(hsa|rno))", fac2drDF_entrez$factor.FA, perl=TRUE),]
        # length(unique(fac2drDF_osRNentrez$factor.FA))   # 785
        # length(unique(fac2drDF_entrez$factor.FA))       # 779
        # length(unique(fac2drDF_entrez_sym$factor.FA))   # 670


        # #######################################################################
        # ## Test enrichment of TRANSFAC gene sets sharing equal TF factor names (factor.FA)
        # ## Added 2015-12-01
        # ##  - create new factors by joining records by factor.FA
        # ##  - TODO: consider creating new factors for fac2drDF_entrez_sym or even more general, e.g. using fac2drDF or fac2drDF_osRNentrez
        # #######################################################################

        # ## TF factor names used as IDs
        # facFAs <- as.character(sort(unique(fac2drDF_entrez_sym$factor.FA))) # 1825 unique TF names

        # ## create annotations: factor.FA, gene.SD, gene.DE, URL
        # annTF.facFA <- data.frame("factorFA"=facFAs, row.names=facFAs)
        # annTF.facFA$factorACs <- base::sapply(facFAs, FUN=function(fac, myFac2drDF) {paste(sort(unique(as.character(myFac2drDF[myFac2drDF$factor.FA==fac,"factor.AC"]))), collapse=";")}, fac2drDF_entrez_sym)
        # annTF.facFA$siteACs <-   base::sapply(facFAs, FUN=function(fac, myFac2drDF) {paste(sort(unique(as.character(myFac2drDF[myFac2drDF$factor.FA==fac,"site.AC"]))),   collapse=";")}, fac2drDF_entrez_sym)
        # annTF.facFA$geneACs <-   base::sapply(facFAs, FUN=function(fac, myFac2drDF) {paste(sort(unique(as.character(myFac2drDF[myFac2drDF$factor.FA==fac,"gene.AC"]))),   collapse=";")}, fac2drDF_entrez_sym)
        # annTF.facFA$geneSDs <-   base::sapply(facFAs, FUN=function(fac, myFac2drDF) {paste(sort(unique(as.character(myFac2drDF[myFac2drDF$factor.FA==fac,"gene.SD"]))),   collapse=";")}, fac2drDF_entrez_sym)
        # #annTF.facFA$geneDEs <-   base::sapply(facFAs, FUN=function(fac, myFac2drDF) {paste(sort(unique(as.character(myFac2drDF[myFac2drDF$factor.FA==fac,"gene.DE"]))),   collapse=";")}, fac2drDF_entrez_sym)

        # ### construct gene set collection by merging factor.FA (factor name), not factor.AC (accession number)
        # facFA2ProbeIDs <- list()
        # gscTF.facFA <- list()
        # for (dName in names(designs)) {
        #     facFA2ProbeIDs[[dName]] <- sapply(facFAs, FUN=function(tffa) {
        #                                   unique(annPrb[[dName]][annPrb[[dName]]$ENTREZID %in% fac2drDF_entrez_sym[fac2drDF_entrez_sym$factor.FA==tffa, "gene.DRac2"], "PROBEID"])
        #                                 }) # list of lists of 1825 `(c-Jun)2` -> "17231461" "17257444" "17305034" "17353699" "17406586" "17427312" "17453819" "17463718"
        #     ## construct GeneSetCollection
        #     gscTF.facFA[[dName]] <- GeneSetCollection(mapply(function(pIds, tffa) {
        #                                 GeneSet(pIds, geneIdType=AnnotationIdentifier(), collectionType=NullCollection(), setName=tffa)
        #                                 }, facFA2ProbeIDs[[dName]], names(facFA2ProbeIDs[[dName]])))
        # }

        # ## run PGSEA TF.facFA ##
        # pgTFAll.facFA <- list()
        # esetsTF.facFA <- list()
        # fitsTF.facFA <- list()
        # for (dName in names(designs)) {
        #     # store original class, change it to "ExpressionSet"
        #     classEset <- class(esets[[dName]])
        #     class(esets[[dName]]) <- "ExpressionSet"
        #     pgTFAll.facFA[[dName]] <- PGSEA(esets[[dName]], gscTF.facFA[[dName]], ref=NULL, range=c(minConceptLen,  max(as.numeric(lapply(GSEABase::geneIds(gscTF.facFA[[dName]]), length)))), center=FALSE, p.value=NA, weighted=TRUE, enforceRange=TRUE)
        #     # restore original class
        #     class(esets[[dName]]) <- classEset
        #     # new PGSEA eset
        #     esetsTF.facFA[[dName]] <- new("ExpressionSet", Biobase::exprs=pgTFAll.facFA[[dName]], phenoData=phenoData(esets[[dName]]), BiocGenerics::annotation="", experimentData=experimentData(esets[[dName]]))
        #     ## MERGE data & annotations
        #     Biobase::fData(esetsTF.facFA[[dName]]) <- annTF.facFA[Biobase::featureNames(esetsTF.facFA[[dName]]),]
        #     ## fit
        #     fitsTF.facFA[[dName]] <- eBayes(contrasts.fit(lmFit(esetsTF.facFA[[dName]], designs[[dName]]), cont.matrices[[dName]]))
        # }




        # ###############################
        # #### From Winnie_meta_v3.R ####
        # ###############################

        # ### DEBUG
        # #tabList=annKeggTT
        # #tabNames=contr4venn
        # #id="DB_ID"
        # #logFC="logFC"
        # #nameCol="KeggName"
        # #addCol="SymbolsDE"
        # ### writeAnnTTlogFCcounts(file.path(resultDir, "2.limma_logFCgenesDE_p0.05.tab"), annTT, contr4venn, id="ENTREZID", nameCol="SYMBOL")
        # #tabList=annTT
        # #tabNames=contr4vennSex
        # #id="ENTREZID"
        # #logFC="logFC"
        # #nameCol="SYMBOL"
        # #addCol=NA
        # ### DEBUG END

        # writeAnnTTlogFCcounts <- function(filepath, tabList, tabNames=NULL, id="DB_ID", logFC="logFC", nameCol="SYMBOL", addCol=NA) {
        #     #### Count genes/pathways that are DE/enriched (in 3+ datasets) and output DE/enriched logFC values, else NA 
        #     #### addCol useful for adding a list of DE gene symbols to enriched sets
        #     require(dplyr)
        #     ## test parameter tabNames
        #     if (is.null(tabNames)) {
        #         tabNames <- names(tabList)
        #     } else if (!all(tabNames %in% names(tabList))) {
        #         stop("Table names (tabNames) do not match names of tables (tabList).")
        #     }
        #     ## names for venn diagrams stored as names(tabNames)
        #     if (is.null(names(tabNames))) names(tabNames) <- tabNames
        #     ## table to return
        #     logFCs <- data.frame(id = character(0), nameCol=character(0))
        #     colnames(logFCs) <- c(id, nameCol)
        #     for (tabName in tabNames) {
        #         if (!is.na(addCol)) {
        #             logFCmeans <- aggregate(tabList[[tabName]][,logFC], by=list(id=tabList[[tabName]][,id], nameCol=tabList[[tabName]][,nameCol], addCol=tabList[[tabName]][,addCol]), mean)
        #             colnames(logFCmeans) <- c(id, nameCol, addCol, logFC)
        #         } else {
        #             logFCmeans <- aggregate(tabList[[tabName]][,logFC], by=list(id=tabList[[tabName]][,id], nameCol=tabList[[tabName]][,nameCol]), mean)
        #             colnames(logFCmeans) <- c(id, nameCol, logFC)
        #         }
        #         logFCs <- dplyr::dplyr::full_join(logFCs, logFCmeans, by=c(id, nameCol))
        #     }
        #     if (is.na(addCol)) colnames(logFCs) <- c(id, nameCol, paste(names(tabNames), logFC)) else colnames(logFCs) <- c(id, nameCol, paste(rep(names(tabNames),each=2), c(addCol, logFC)))
        #     logFCs$CountBoth <- apply(logFCs[,paste(names(tabNames), logFC)], MARGIN=1, FUN=function(x){sum(!is.na(x))})
        #     logFCs$CountPos <-  apply(logFCs[,paste(names(tabNames), logFC)], MARGIN=1, FUN=function(x){sum(x>0 & !is.na(x))})
        #     logFCs$CountNeg <-  apply(logFCs[,paste(names(tabNames), logFC)], MARGIN=1, FUN=function(x){sum(x<0 & !is.na(x))})
        #     write.table(logFCs, filepath, sep="\t", quote=FALSE, row.names=FALSE)
        # }


        # ### DEBUG
        # ## plot_venneuler(file.path(curDir,"Kb.gK_W.sF~St.dSHep_NonSHep.gF~Gl.tcP_cM.dA.gF"), annTT, c("Kb.gK_W.sF", "St.dSHep_NonSHep.gF", "Gl.tcP_cM.dA.gF"))
        # # tabList <- annTT
        # # tabNames <- c("Kb.gK_W.sF", "St.dSHep_NonSHep.gF", "Gl.tcP_cM.dA.gF")
        # # id="HSA_ENTREZID"
        # # directional=TRUE
        # # logFC="logFC"
        # # vennPrint.mode=c("raw","percent")
        # # vennHeight=900
        # # vennWidth=900
        # # vennRes=150
        # # vennUnits="px"
        # # vennMargin=0.06
        # # vennImagetype="png"
        # # eulerHeight=9
        # # eulerWidth=9
        # ### DEBUG END ###

        # plot_venneuler <- function(filepath, tabList, tabNames=NULL, id="HSA_ENTREZID", directional=TRUE, logFC="logFC", 
        #                            vennPrint.mode=c("raw","percent"), vennHeight=900, vennWidth=900, vennRes=150, vennUnits="px", vennMargin=0.06, vennImagetype="png",
        #                            eulerHeight=9, eulerWidth=9) {
        #     ## TODO: unique(tab$id)
        #     ## TODO: use logFC
        #     ## TODO: write partitions
        #     require(ggplot2)
        #     require(svglite)
        #     require(venneuler)
        #     require(VennDiagram)
        #     require(RColorBrewer)
        #     require(futile.logger)
        #     ## functions
        #     col.dist <- function(inp, comp) sum( abs(inp - col2rgb(comp) ) )
        #     ## settings: no log output
        #     futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
        #     ## test parameter tabNames
        #     if (is.null(tabNames)) {
        #         tabNames <- names(tabList)
        #     } else if (!all(tabNames %in% names(tabList))) {
        #         stop("Table names (tabNames) do not match names of tables (tabList).")
        #     }
        #     ## names for venn diagrams stored as names(tabNames)
        #     if (is.null(names(tabNames))) names(tabNames) <- tabNames
        #     ## arrange data
        #     tabListPN <- list()
        #     if (directional) {
        #         sgns <- c("neg","pos")
        #         ## split tables to neg and pos DE genes, remove genes with logFC NA
        #         for (tName in tabNames) {
        #             tabListPN[["neg"]][[tName]] <- tabList[[tName]][tabList[[tName]][,logFC] < 0  & !is.na(tabList[[tName]][,logFC]), c(id,logFC)]
        #             tabListPN[["pos"]][[tName]] <- tabList[[tName]][tabList[[tName]][,logFC] >= 0 & !is.na(tabList[[tName]][,logFC]), c(id,logFC)]
        #         }
        #     } else {
        #         sgns <- c("both")
        #         ## remove genes with logFC NA
        #         for (tName in tabNames) tabListPN[["both"]][[tName]] <- tabList[[tName]][!is.na(tabList[[tName]][,logFC]), c(id,logFC)]
        #     }
        #     for (sgn in sgns) {
        #         dfP <- data.frame()  # for venneuler
        #         lnsP <- c()          # for venneuler
        #         vdlP <- list()       # for VennDiagram
        #         for (tNameLab in names(tabNames)) {
        #             tName <- tabNames[[tNameLab]]
        #             tabP <- tabListPN[[sgn]][[tName]][!duplicated(tabListPN[[sgn]][[tName]][,id]),]
        #             lnP <- dim(tabP)[[1]]
        #             dfP <- rbind(dfP, cbind(tabP, sets=rep(tNameLab, lnP), elements=tabP[,id]))
        #             lnsP <- c(lnsP, lnP)
        #             vdlP[[tName]] <- tabP[,id]
        #         }
        #         ## venneuler
        #         titP <- paste(paste(names(tabNames), collapse=" ~ "), "(logFC", sgn, ")")
        #         if (dim(dfP)[[1]] > 0) {
        #             veP <- venneuler::venneuler(dfP[,c("elements","sets")])
        #             veP$labels <- paste(veP$labels, lnsP)
        #             pdf(paste(filepath,sgn,"pdf",sep="."), height=eulerHeight, width=eulerWidth)
        #             plot(veP)
        #             title(titP)
        #             dev.off()
        #         }
        #         else warning(paste("venneuler",titP,"empty set",sep=":"))
        #         ## VennDiagram
        #         if (length(tabNames) %in% c(2,3)) euler.d <- TRUE else euler.d <- FALSE
        #         brewVD <- RColorBrewer::brewer.pal(length(tabNames), "Set2")[1:length(tabNames)]
        #         vdPlot <- VennDiagram::venn.diagram(filename=NULL, height=vennHeight, width=vennWidth, resolution=vennRes,
        #                                   x=vdlP, main=paste(paste(names(tabNames), collapse=" ~ "),"(logFC",sgn,")"), 
        #                                   force.unique=TRUE, print.mode=vennPrint.mode, euler.d=euler.d,
        #                                   fill=colors()[apply(col2rgb(brewVD), 2, function(z) which.min(sapply(colors(), function(x) col.dist(inp=z, comp=x))))],
        #                                   category.names=names(tabNames), margin=vennMargin, imagetype=vennImagetype, units=vennUnits)
        #         ggsave(vdPlot, file=paste(filepath,sgn,vennImagetype,sep="."), device=vennImagetype)
        #     }#end for signs
        # }


        # ### DEBUG
        # #filepath=file.path(resultDir, "2.vennParts", "2.CGdebug")
        # #filepath=file.path(resultDir, "2.vennParts", "3.RCGdebug")
        # ##filepath=file.path(resultDir, "2.vennParts", "4.NRCGdebug")
        # #tabList=annTT
        # #tabNames=c("C.tP_N.aw19.sM", "G.tKO_WT")
        # #tabNames=c("R.Rbpj_WT", "C.tP_N.aw19.sM", "G.tKO_WT")
        # ##tabNames=c("N.NEMO_WT", "R.Rbpj_WT", "C.tP_N.aw19.sM", "G.tKO_WT")
        # #id="ENTREZID"
        # #directional=TRUE
        # #logFC="logFC"
        # ### DEBUG END

        # vennPartition <- function(filepath, tabList, tabNames=NULL, id="HSA_ENTREZID", directional=TRUE, logFC="logFC") {
        #     ## directional==TRUE: split sets to "pos" and "neg" expressed using column logFC
        #     ##              FALSE: ignore direction of expression
        #     ## ADDED 2020-03-09: output "*.ratioArr.tab" with pairwise ratios between intersection and the 1st set size; diagonals equal set sizes
        #     require(dplyr)
        #     ## test parameter tabNames
        #     if (is.null(tabNames)) {
        #         tabNames <- names(tabList)
        #     } else if (!all(tabNames %in% names(tabList))) {
        #         stop("Table names (tabNames) do not match names of tables (tabList).")
        #     }
        #     ## names for venn diagrams stored as names(tabNames)
        #     if (is.null(names(tabNames))) names(tabNames) <- tabNames
        #     ## split tables to pos/neg logFC, remove logFC NA
        #     tabListPN <- list()
        #     if (directional) {
        #         sgns <- c("neg","pos")
        #         ## split tables to neg and pos DE genes, remove genes with logFC NA
        #         for (tName in tabNames) {
        #             tabListPN[["neg"]][[tName]] <- tabList[[tName]][tabList[[tName]][,logFC] < 0  & !is.na(tabList[[tName]][,logFC]), ]
        #             tabListPN[["pos"]][[tName]] <- tabList[[tName]][tabList[[tName]][,logFC] >= 0 & !is.na(tabList[[tName]][,logFC]), ]
        #         }
        #     } else {
        #         sgns <- c("both")
        #         ## remove genes with logFC NA
        #         for (tName in tabNames) tabListPN[["both"]][[tName]] <- tabList[[tName]][!is.na(tabList[[tName]][,logFC]), ]
        #     }
        #     for (sgn in sgns) {
        #         tls <- tabListPN[[sgn]]
        #         ## pairwise ratio between intersection and the 1st set size (e.g., 1st set corresponds to the first index); diagonals equal set sizes
        #         ratioArr <- array(dim=c(length(tabNames),length(tabNames)))
        #         colnames(ratioArr) <- rownames(ratioArr) <- names(tabNames)
        #         diag(ratioArr) <- sapply(tls, function(x)length(unique(x[,id])))
        #         for (i in 1:(dim(ratioArr)[[1]]-1)) {
        #             for (j in (i+1):dim(ratioArr)[[1]]) {
        #                 lInt <- length(intersect(tls[[i]][,id], tls[[j]][,id]))
        #                 ratioArr[i,j] <- lInt/ratioArr[i,i]
        #                 ratioArr[j,i] <- lInt/ratioArr[j,j]
        #             }
        #         }
        #         write.table(cbind("tabName"=rownames(ratioArr),format(ratioArr,digits=3)), paste(filepath, sgn, "ratioArr", "tab", sep="."), sep="\t", quote=FALSE, row.names=FALSE)
        #         ## write partitions
        #         if (length(tabNames) == 1) {
        #             if (dim(tls[[1]])[[1]] > 0) write.table(format(tls[[1]],digits=3), paste(filepath, sgn, "txt", sep="."), sep="\t", quote=FALSE, row.names=FALSE)
        #         } else if (length(tabNames) == 2) {
        #             ids1 <-  setdiff(tls[[1]][,id], tls[[2]][,id])
        #             ids2 <-  setdiff(tls[[2]][,id], tls[[1]][,id])
        #             ids12 <- intersect(tls[[1]][,id], tls[[2]][,id])

        #             tls1 <- tls[[1]][tls[[1]][,id] %in% ids1,]
        #             if (dim(tls1)[[1]] > 0) write.table(format(tls1,digits=3), paste(filepath, sgn, "A", "txt", sep="."), sep="\t", quote=FALSE, row.names=FALSE)
        #             tls2 <- tls[[2]][tls[[2]][,id] %in% ids2,]
        #             if (dim(tls2)[[1]] > 0) write.table(format(tls2,digits=3), paste(filepath, sgn, "B", "txt", sep="."), sep="\t", quote=FALSE, row.names=FALSE)
        #             t12 <- dplyr::full_join(tls[[1]][tls[[1]][,id] %in% ids12,], tls[[2]][tls[[2]][,id] %in% ids12,], by=id)
        #             if (dim(t12)[[1]] > 0) write.table(format(t12,digits=3), paste(filepath, sgn, "A&B", "txt", sep="."), sep="\t", quote=FALSE, row.names=FALSE)
        #         } else if (length(tabNames) == 3) {
        #             ids1 <-  setdiff(setdiff(tls[[1]][,id], tls[[2]][,id]), tls[[3]][,id])
        #             ids2 <-  setdiff(setdiff(tls[[2]][,id], tls[[1]][,id]), tls[[3]][,id])
        #             ids3 <-  setdiff(setdiff(tls[[3]][,id], tls[[1]][,id]), tls[[2]][,id])
        #             ids12 <-  setdiff(intersect(tls[[1]][,id], tls[[2]][,id]), tls[[3]][,id])
        #             ids13 <-  setdiff(intersect(tls[[1]][,id], tls[[3]][,id]), tls[[2]][,id])
        #             ids23 <-  setdiff(intersect(tls[[2]][,id], tls[[3]][,id]), tls[[1]][,id])
        #             ids123 <-  intersect(intersect(tls[[1]][,id], tls[[2]][,id]), tls[[3]][,id])

        #             tls1 <- tls[[1]][tls[[1]][,id] %in% ids1,]
        #             if (dim(tls1)[[1]] > 0) write.table(format(tls1,digits=3), paste(filepath, sgn, "A", "txt", sep="."), sep="\t", quote=FALSE, row.names=FALSE)
        #             tls2 <- tls[[2]][tls[[2]][,id] %in% ids2,]
        #             if (dim(tls2)[[1]] > 0) write.table(format(tls2,digits=3), paste(filepath, sgn, "B", "txt", sep="."), sep="\t", quote=FALSE, row.names=FALSE)
        #             tls3 <- tls[[3]][tls[[3]][,id] %in% ids3,]
        #             if (dim(tls3)[[1]] > 0) write.table(format(tls3,digits=3), paste(filepath, sgn, "C", "txt", sep="."), sep="\t", quote=FALSE, row.names=FALSE)

        #             tls12 <- tls[[1]][tls[[1]][,id] %in% ids12,]
        #             tls21 <- tls[[2]][tls[[2]][,id] %in% ids12,]
        #             t12 <- dplyr::full_join(tls12,tls21,by=id)
        #             if (dim(t12)[[1]] > 0) write.table(format(t12,digits=3), paste(filepath, sgn, "A&B", "txt", sep="."), sep="\t", quote=FALSE, row.names=FALSE)

        #             tls13 <- tls[[1]][tls[[1]][,id] %in% ids13,]
        #             tls31 <- tls[[3]][tls[[3]][,id] %in% ids13,]
        #             t13 <- dplyr::full_join(tls13,tls31,by=id)
        #             if (dim(t13)[[1]] > 0) write.table(format(t13,digits=3), paste(filepath, sgn, "A&C", "txt", sep="."), sep="\t", quote=FALSE, row.names=FALSE)

        #             tls23 <- tls[[2]][tls[[2]][,id] %in% ids23,]
        #             tls32 <- tls[[3]][tls[[3]][,id] %in% ids23,]
        #             t23 <- dplyr::full_join(tls23,tls32,by=id)
        #             if (dim(t23)[[1]] > 0) write.table(format(t23,digits=3), paste(filepath, sgn, "B&C", "txt", sep="."), sep="\t", quote=FALSE, row.names=FALSE)

        #             tls123 <- tls[[1]][tls[[1]][,id] %in% ids123,]
        #             tls213 <- tls[[2]][tls[[2]][,id] %in% ids123,]
        #             tls312 <- tls[[3]][tls[[3]][,id] %in% ids123,]
        #             t123 <- dplyr::full_join(dplyr::full_join(tls123,tls213,by=id),tls312,by=id)
        #             if (dim(t123)[[1]] > 0) write.table(format(t123,digits=3), paste(filepath, sgn, "A&B&C", "txt", sep="."), sep="\t", quote=FALSE, row.names=FALSE)

        #         } else if (length(tabNames) == 4) {
        #             ids1 <-  setdiff(setdiff(setdiff(tls[[1]][,id], tls[[2]][,id]), tls[[3]][,id]), tls[[4]][,id])
        #             ids2 <-  setdiff(setdiff(setdiff(tls[[2]][,id], tls[[1]][,id]), tls[[3]][,id]), tls[[4]][,id])
        #             ids3 <-  setdiff(setdiff(setdiff(tls[[3]][,id], tls[[1]][,id]), tls[[2]][,id]), tls[[4]][,id])
        #             ids4 <-  setdiff(setdiff(setdiff(tls[[4]][,id], tls[[1]][,id]), tls[[2]][,id]), tls[[3]][,id])

        #             t1 <- tls[[1]][tls[[1]][,id] %in% ids1,]
        #             if (dim(t1)[[1]] > 0) write.table(format(t1,digits=3), paste(filepath, sgn, "A", "txt", sep="."), sep="\t", quote=FALSE, row.names=FALSE)
        #             t2 <- tls[[2]][tls[[2]][,id] %in% ids2,]
        #             if (dim(t2)[[1]] > 0) write.table(format(t2,digits=3), paste(filepath, sgn, "B", "txt", sep="."), sep="\t", quote=FALSE, row.names=FALSE)
        #             t3 <- tls[[3]][tls[[3]][,id] %in% ids3,]
        #             if (dim(t3)[[1]] > 0) write.table(format(t3,digits=3), paste(filepath, sgn, "C", "txt", sep="."), sep="\t", quote=FALSE, row.names=FALSE)
        #             t4 <- tls[[4]][tls[[4]][,id] %in% ids4,]
        #             if (dim(t4)[[1]] > 0) write.table(format(t4,digits=3), paste(filepath, sgn, "D", "txt", sep="."), sep="\t", quote=FALSE, row.names=FALSE)

        #             ids12 <-  setdiff(setdiff(intersect(tls[[1]][,id], tls[[2]][,id]), tls[[3]][,id]), tls[[4]][,id])
        #             ids13 <-  setdiff(setdiff(intersect(tls[[1]][,id], tls[[3]][,id]), tls[[2]][,id]), tls[[4]][,id])
        #             ids14 <-  setdiff(setdiff(intersect(tls[[1]][,id], tls[[4]][,id]), tls[[2]][,id]), tls[[3]][,id])
        #             ids23 <-  setdiff(setdiff(intersect(tls[[2]][,id], tls[[3]][,id]), tls[[1]][,id]), tls[[4]][,id])
        #             ids24 <-  setdiff(setdiff(intersect(tls[[2]][,id], tls[[4]][,id]), tls[[1]][,id]), tls[[3]][,id])
        #             ids34 <-  setdiff(setdiff(intersect(tls[[3]][,id], tls[[4]][,id]), tls[[1]][,id]), tls[[2]][,id])

        #             t12 <- dplyr::full_join(tls[[1]][tls[[1]][,id] %in% ids12,], tls[[2]][tls[[2]][,id] %in% ids12,], by=id)
        #             if (dim(t12)[[1]] > 0) write.table(format(t12,digits=3), paste(filepath, sgn, "A&B", "txt", sep="."), sep="\t", quote=FALSE, row.names=FALSE)
        #             t13 <- dplyr::full_join(tls[[1]][tls[[1]][,id] %in% ids13,], tls[[3]][tls[[3]][,id] %in% ids13,], by=id)
        #             if (dim(t13)[[1]] > 0) write.table(format(t13,digits=3), paste(filepath, sgn, "A&C", "txt", sep="."), sep="\t", quote=FALSE, row.names=FALSE)
        #             t14 <- dplyr::full_join(tls[[1]][tls[[1]][,id] %in% ids14,], tls[[4]][tls[[4]][,id] %in% ids14,], by=id)
        #             if (dim(t14)[[1]] > 0) write.table(format(t14,digits=3), paste(filepath, sgn, "A&D", "txt", sep="."), sep="\t", quote=FALSE, row.names=FALSE)
        #             t23 <- dplyr::full_join(tls[[2]][tls[[2]][,id] %in% ids23,], tls[[3]][tls[[3]][,id] %in% ids23,], by=id)
        #             if (dim(t23)[[1]] > 0) write.table(format(t23,digits=3), paste(filepath, sgn, "B&C", "txt", sep="."), sep="\t", quote=FALSE, row.names=FALSE)
        #             t24 <- dplyr::full_join(tls[[2]][tls[[2]][,id] %in% ids24,], tls[[4]][tls[[4]][,id] %in% ids24,], by=id)
        #             if (dim(t24)[[1]] > 0) write.table(format(t24,digits=3), paste(filepath, sgn, "B&D", "txt", sep="."), sep="\t", quote=FALSE, row.names=FALSE)
        #             t34 <- dplyr::full_join(tls[[3]][tls[[3]][,id] %in% ids34,], tls[[4]][tls[[4]][,id] %in% ids34,], by=id)
        #             if (dim(t34)[[1]] > 0) write.table(format(t34,digits=3), paste(filepath, sgn, "C&D", "txt", sep="."), sep="\t", quote=FALSE, row.names=FALSE)

        #             ids123 <-  setdiff(intersect(intersect(tls[[1]][,id], tls[[2]][,id]), tls[[3]][,id]), tls[[4]][,id])
        #             ids124 <-  setdiff(intersect(intersect(tls[[1]][,id], tls[[2]][,id]), tls[[4]][,id]), tls[[3]][,id])
        #             ids134 <-  setdiff(intersect(intersect(tls[[1]][,id], tls[[3]][,id]), tls[[4]][,id]), tls[[2]][,id])
        #             ids234 <-  setdiff(intersect(intersect(tls[[2]][,id], tls[[3]][,id]), tls[[4]][,id]), tls[[1]][,id])

        #             t123 <- dplyr::full_join(dplyr::full_join(tls[[1]][tls[[1]][,id] %in% ids123,], tls[[2]][tls[[2]][,id] %in% ids123,], by=id), tls[[3]][tls[[3]][,id] %in% ids123,], by=id)
        #             if (dim(t123)[[1]] > 0) write.table(format(t123,digits=3), paste(filepath, sgn, "A&B&C", "txt", sep="."), sep="\t", quote=FALSE, row.names=FALSE)
        #             t124 <- dplyr::full_join(dplyr::full_join(tls[[1]][tls[[1]][,id] %in% ids124,], tls[[2]][tls[[2]][,id] %in% ids124,], by=id), tls[[4]][tls[[4]][,id] %in% ids124,], by=id)
        #             if (dim(t124)[[1]] > 0) write.table(format(t124,digits=3), paste(filepath, sgn, "A&B&D", "txt", sep="."), sep="\t", quote=FALSE, row.names=FALSE)
        #             t134 <- dplyr::full_join(dplyr::full_join(tls[[1]][tls[[1]][,id] %in% ids134,], tls[[3]][tls[[3]][,id] %in% ids134,], by=id), tls[[4]][tls[[4]][,id] %in% ids134,], by=id)
        #             if (dim(t134)[[1]] > 0) write.table(format(t134,digits=3), paste(filepath, sgn, "A&C&D", "txt", sep="."), sep="\t", quote=FALSE, row.names=FALSE)
        #             t234 <- dplyr::full_join(dplyr::full_join(tls[[2]][tls[[2]][,id] %in% ids234,], tls[[3]][tls[[3]][,id] %in% ids234,], by=id), tls[[4]][tls[[4]][,id] %in% ids234,], by=id)
        #             if (dim(t234)[[1]] > 0) write.table(format(t234,digits=3), paste(filepath, sgn, "B&C&D", "txt", sep="."), sep="\t", quote=FALSE, row.names=FALSE)

        #             ids1234 <-  intersect(intersect(intersect(tls[[1]][,id], tls[[2]][,id]), tls[[3]][,id]), tls[[4]][,id])

        #             t1234 <- dplyr::full_join(dplyr::full_join(dplyr::full_join(tls[[1]][tls[[1]][,id] %in% ids1234,], tls[[2]][tls[[2]][,id] %in% ids1234,], by=id), tls[[3]][tls[[3]][,id] %in% ids1234,], by=id), tls[[4]][tls[[4]][,id] %in% ids1234,], by=id)
        #             if (dim(t1234)[[1]] > 0) write.table(format(t1234,digits=3), paste(filepath, sgn, "A&B&C&D", "txt", sep="."), sep="\t", quote=FALSE, row.names=FALSE)
        #         } else {
        #             warning("Not fully implemented for > 4 sets")
        #         }
        #     }
        # }
