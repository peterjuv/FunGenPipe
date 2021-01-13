#' @section Genomics
#' Functions that operate with ranges and sequences
#' dechiperMeltDNAtime
#' @section Tagets meta-data
#' Functions for manipulating colors in targets.
#' getColPats - depricated
#' getColPats2
#' @section Plot distributions, heatmap, PCA & MDS using targets meta-data
#' Note: the order od data columns should (always) be kept in the order of targets.
#' plotDistrTargets2
#' pheatmapTargets2
#' plotPCAtargets2
#' plotPCAtargets - depricated
#' plotMDStargets2
#' plotMDStargets2_MASS - depricated
#' @section Plot distributions using limma-style data
#' plotDistr_RGList
#' boxplotRLE
#' @section From Osijek_HFHSD
#' For KBlag_doxy_Cla.R
NULL

###################
#### Genomics #####
###################

#' Calculate melting temperatures from ranges using a specified genome.
#' 
#' Tm is estimated from a range of temperatures (Trange) by finding a max of ThetaDerivative over that range.
#' Tm is set to NA in case of multiple Tm values.
#' @param myGRanges GRanges, ranges of sequences with a specified genome
#' @param myBSgenome BSgenome object
#' @param Trange Numeric vector of melting temperatures for calculation of Theta derivatives, forwarded to \code{DECIPHER::MeltDNA}; optional, default 50:100
#' @param ions Numeric molar sodium equivalent ionic concentration, forwarded to \code{DECIPHER::MeltDNA}; optional, default 0.2
#' @param dropMcols Intermediate columns to drop, default c("seq","width",Trange","ThetaDerivative","maxThetaDerivative","TmList"), NULL to include all
#' @return GRanges with melting temperatures added to mcols column "Tm" and optional intermediate columns
#' sequences and widths, (max)ThetaDerivative & list of Tms in case max Tm not unique; message execution time
#' @importFrom magrittr %>%
#' @importFrom assertthat assert_that
#' @importFrom GenomeInfoDb genome seqnames
#' @importFrom Biostrings getSeq
#' @importFrom tibble tibble
#' @importFrom dplyr mutate
#' @importFrom purrr map2 map_dbl pmap map_if
#' @importFrom DECIPHER MeltDNA
#' @export
dechiperMeltDNAtime <- function(myGRanges, myBSgenome, Trange=seq(50,100,1), ions=0.2,
    dropMcols=c("seq","width","Trange","ThetaDerivative","maxThetaDerivative","TmList")) {
    require(myBSgenome)
    assertthat::assert_that(all(GenomeInfoDb::genome(myGRanges) %in% GenomeInfoDb::genome(myBSgenome))) # test genome names
    assertthat::assert_that(all(names(GenomeInfoDb::genome(myGRanges)) %in% GenomeInfoDb::seqnames(myBSgenome))) # test seqnames
    timeDECHIPER = list()
    timeDECHIPER[["start"]] <- Sys.time()
    dss <- Biostrings::getSeq(myBSgenome, myGRanges)
    wdss <- which(width(dss) > 2)
    # dim(thDer) == Trange x myGRanges[wdss]
    thDer.w <- DECIPHER::MeltDNA(dss[wdss], type="derivative", temps=Trange, ions=ions)
    thDerList.w <- lapply(seq_len(ncol(thDer.w)), function(i) thDer.w[,i])
    tbSeqTm.w <- dss[wdss] %>% {
        tibble::tibble(name = names(.), seq = as.character(.), width = width(.), Trange = list(Trange))
        } %>%
        dplyr::mutate(
            ThetaDerivative = tibble::enframe(thDerList.w)$value,
            maxThetaDerivative = purrr::map_dbl(ThetaDerivative, max),
            TmList = purrr::pmap(list(Trange, ThetaDerivative, maxThetaDerivative), 
                ~ tryCatch(..1[which(..2 == ..3)], error=function(e) NA)),
            Tm = unlist(purrr::map_if(TmList, ~ length(.) != 1, ~ NA)),
            tmp_idx = wdss
        ) %>% 
        dplyr::relocate(Tm)
    timeDECHIPER[["finish"]] <- Sys.time()
    message(paste("DECHIPER::MeldDNA start:", timeDECHIPER[["start"]], "finish: ", timeDECHIPER[["finish"]]))
    S4Vectors::mcols(myGRanges) <- S4Vectors::mcols(myGRanges) %>%
        tibble::as_tibble() %>% 
        dplyr::mutate(tmp_idx = seq(n())) %>% 
        dplyr::left_join(tbSeqTm.w, by="tmp_idx") %>% 
        dplyr::select(-tmp_idx, -name) %>% 
        dplyr::select(-tidyselect::any_of(dropMcols))
    return(myGRanges)
}

############################
#### Targets meta-data #####
############################
 
#' Get color patterns with names from another column from ExpressionSet slot phenoData.
#' 
#' Returns columns from phenoData that start with "color_" and have their suffix in common with another column;
#' Names of colors are set from the associated column, e.g., list(color_Sex = setNames(color_Sex, Sex), ...)
#' @param eset ExpressionSet
#' @return List of colors, each of length equal to the number of rows in phenoData
#' @examples
#' /donotrun{
#' getColPats(Biobase::pData(dataRaw/eset))
#' getColPats(Biobase::pData(eset))
#' }
#' @section Old implementation:
#' \dontrun{
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
        tibble::as_tibble() %>% 
        select(starts_with("color_"))
    names <- Biobase::pData(eset) %>% 
        tibble::as_tibble() %>%
        select(all_of(sub("color_", "", colnames(cols))))
    map2(cols, names, setNames)
}


#' Get named color patterns with names preset or from other column(s) from targets
#' 
#' Returns columns from targets that start with "color_" and have their suffix in common with another column;
#' column names are dropped prefix defined by parameter colorsFrom, e.g., "color_"
#' New implementation allowing for alternative naming using parameter namesFrom.
#' Names of colors may be:
#' - preset (if rename = FALSE) or 
#' - set from the associated column, e.g., from Var for color_Var
#' - set from a common coulmn for all colors set by namesFrom parameter, e.g., namesFrom=HybName 
#' Minimum 2 columns are required, e.g. Var and color_Var
#' 
#' @param targets Table with columns names with a prefix colorsFrom='color_'
#' @param colorsFrom Character, non-empty prefix of names of columns with colors, default 'color_'
#' @param rename Logical, default TRUE: rename colors using associated variables or variable namesFrom (default) 
#'               FALSE: use existing names from 'color_' variables, or, if missing, use values from associated variables
#' @param namesFrom Either NULL for using variables that share suffixes with colors
#'                  or a variable from targets that is used for setting names to colors, e.g. HybName
#' @param unique If TRUE, returns unique names and associated colors;
#'               Note that colors may be dropped, see the last example
#' @param pullVar A variable from targets, e.g. Sex
#' @return Named list of named colors for each variable from targets with an associated color.
#'         Names of variables correspond to color variables w/o colorsFrom prefix, e.g. "VarName" for "color_VarName"
#'         Names of colors correspond to values of an associated variable by default; 
#'         alternatively, names are set to values form namesFrom variable.
#'         The lenght of each list is equal to the number of rows;
#'         alternatively, lists are truncated to unique names if unique=TRUE.
#'         If pullVar given, returns a single list of colors for that variable from targets 
#' @examples
#' \dontrun{
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
#' getColPats2(targets, namesFrom=cc, unique=TRUE)
#' # using preset names, and if missing, substiting from associated variables
#' getColPats2(targets, rename=FALSE, namesFrom=NULL)
#' # using preset names, and if missing, substiting from associated variables, not using namesFrom=cc, issueing a warning
#' getColPats2(targets, rename=FALSE, namesFrom=cc)
#' }
#' @importFrom magrittr %>% %<>%
#' @importFrom assertthat assert_that are_equal has_name
#' @importFrom tibble as_tibble
#' @importFrom rlang enquo quo_is_null quo_name
#' @importFrom dplyr select rename_with
#' @importFrom tidyselect starts_with any_of
#' @importFrom purrr map2 map keep
#' @export
getColPats2 <- function(targets, colorsFrom="color_", rename=TRUE, namesFrom=NULL, unique=FALSE, pullVar=NULL) {
    require(magrittr)
    assertthat::assert_that(is.character(colorsFrom) & colorsFrom != "", 
        msg="Parameter colorsFrom must be a non-empty character prefix of names of variables from targets")
    targets %<>% tibble::as_tibble()
    namesFrom <- rlang::enquo(namesFrom)
    pullVar <- rlang::enquo(pullVar)
    if (rename==FALSE & !rlang::quo_is_null(namesFrom)) 
        warning("Parameter namesFrom=", rlang::quo_name(namesFrom), " not used if rename=FALSE")
    assertthat::assert_that(rlang::quo_is_null(namesFrom) | assertthat::has_name(targets, rlang::quo_name(namesFrom)))
    cols <- targets %>% 
        dplyr::select(tidyselect::starts_with(colorsFrom)) %>% 
        dplyr::rename_with(~sub(colorsFrom, "", .))
    if (rename) {
        if (rlang::quo_is_null(namesFrom))
            cnames <- targets %>% dplyr::select(tidyselect::any_of(sub(colorsFrom, "", colnames(cols))))
        else
            cnames <- targets %>% dplyr::select(!!namesFrom)
        cnames[is.na(cnames)] <- "<NA>"
        assertthat::assert_that(assertthat::are_equal(dim(cols)[[1]], dim(cnames)[[1]]) & 
            (dim(cnames)[[2]] %in% c(1,dim(cols)[[2]])),
            msg = paste("The number of colors and variables is not the same. Not all variables with a prefix", colorsFrom, "have an associated variable. Consider using namesFrom parameter."))
        ncols <- purrr::map2(cols, cnames, setNames)
    } else {
        ncols <- as.list(cols)
        # add names (only if missing) from associated variables
        for (nc in names(ncols)) if (is.null(names(ncols[[nc]]))) ncols[[nc]] <- setNames(ncols[[nc]], targets[[nc]])
    }
    if (!rlang::quo_is_null(pullVar)) {
        ncol <- ncols %>% tibble::as_tibble() %>% pull(!!pullVar)
        if (unique)
            ncol <- ncol[!duplicated(names(ncol))]
        return(ncol)
    } else {
        if (unique)
            ncols %<>% purrr::map(~purrr::keep(., !duplicated(names(.))))
        return(ncols)
    }
}

########################################################################
#### Plot distributions, heatmap, PCA, MDS using targets meta-data #####
########################################################################

#' Plot distributions of signals from limma RGList, MAList, EList or EListRaw objects 
#' using multiple colors for groups of samples
#' 
#' Plots distributions for each combination of FUNS, channels and sampleColors.
#' By default, samples are ordered according to the sampleColors (if provided).
#' Colors should be provided as factors with a preset order of levels; the order of colors is determined by the order of levels
#' 
#' @param expLog2 Expression matrix (preferably in log2 scale) with genes in rows and samples in columns
#' @inheritParams getColPats2
#' @param FUNS Vector, geom_functions from ggplot2, default \code{c(geom_density, geom_boxplot)}, alternatives \code{c(geom_violin, geom_histogram, ggridges::geom_density_ridges)}
#' @param probeTypeVec Vector of probe types
#' @param probeTypeValue Value from probeTypeVec to use for plotting
#' @param numProbes Number of randomly select probes
# @param sampleColors Tibble of named factors of colors with names matching sample names
#' @param orderByColors Logical for ordering samples by colors (coded as a factor with levels in order); default TRUE
#' @param scale_x_limits NULL for auto-scale; use c(0,16) or less for log2(intensities)
#' @param filePath NULL or character; if given, output a PDF; default NULL
#' @param width PDF width, default 16/9*7=12.44
#' @param height PDF height, default 7
#' @param ... Passed to ggplot2::FUN, e.g.: bins, binwidth, show.legend
#' @return Invisibly a list of ggplot2 objects, one per a combination of FUNS and colors;
#' PDF if filePath is given
#' @examples
#' \dontrun{
#' expLog2 <- log2(dataRG$R)[1:100,]
#' ## rename=TRUE (default)
#' ## namesFrom=NULL,     orderByColors=TRUE
#' plotDistrTargets2(expLog2, targets, colorsFrom="color_", namesFrom=NULL, FUNS = c(geom_density, geom_boxplot),
#'     probeTypeVec = probeTypeVec, probeTypeValue = 0, numProbes = NULL, 
#'     scale_x_limits = NULL, filePath = file.path(resultDir, "_debug3_plotDistrTargets2-NULL.pdf"))
#' ## namesFrom=HybName
#' plotDistrTargets2(expLog2, targets, colorsFrom="color_", namesFrom=HybName, FUNS = c(geom_density, geom_boxplot),
#'     probeTypeVec = probeTypeVec, probeTypeValue = 0, numProbes = NULL, 
#'     scale_x_limits = NULL, filePath = file.path(resultDir, "_debug3_plotDistrTargets2-HybName.pdf"))
#' ## namesFrom=Birth3WHO
#' plotDistrTargets2(expLog2, targets, colorsFrom="color_", namesFrom=Birth3WHO, FUNS = c(geom_density, geom_boxplot),
#'     probeTypeVec = probeTypeVec, probeTypeValue = 0, numProbes = NULL, 
#'     scale_x_limits = NULL, filePath = file.path(resultDir, "_debug3_plotDistrTargets2-Birth3WHO.pdf"))
#' ## namesFrom=HybName,  orderByColors=FALSE
#' plotDistrTargets2(expLog2, targets, colorsFrom="color_", namesFrom=HybName, FUNS = c(geom_density, geom_boxplot),
#'     probeTypeVec = probeTypeVec, probeTypeValue = 0, numProbes = NULL, orderByColors = FALSE,
#'     scale_x_limits = NULL, filePath = file.path(resultDir, "_debug3_plotDistrTargets2-HybName-noReorder.pdf"))
#' ## rename=FALSE () (implies namesFrom=NULL)
#' ## orderByColors=TRUE
#' plotDistrTargets2(expLog2, targets, colorsFrom="color_", rename=FALSE, namesFrom=NULL, FUNS = c(geom_density, geom_boxplot),
#'     probeTypeVec = probeTypeVec, probeTypeValue = 0, numProbes = NULL, 
#'     scale_x_limits = NULL, filePath = file.path(resultDir, "_debug3_plotDistrTargets2-noRename.pdf"))
#' ## orderByColors=FALSE
#' plotDistrTargets2(expLog2, targets, colorsFrom="color_", rename=FALSE, namesFrom=NULL, FUNS = c(geom_density, geom_boxplot),
#'     probeTypeVec = probeTypeVec, probeTypeValue = 0, numProbes = NULL, orderByColors = FALSE,
#'     scale_x_limits = NULL, filePath = file.path(resultDir, "_debug3_plotDistrTargets2-noRename-noReorder.pdf"))
#' ## warning: rename=FALSE $ namesFrom!=NULL
#' plotDistrTargets2(expLog2, targets, colorsFrom="color_", rename=FALSE, namesFrom=Birth2, FUNS = c(geom_density, geom_boxplot),
#'     probeTypeVec = probeTypeVec, probeTypeValue = 0, numProbes = NULL, orderByColors = FALSE,
#'     scale_x_limits = NULL, filePath = file.path(resultDir, "_debug3_plotDistrTargets2-warning.pdf"))
#' ## w/o targets
#' getColPats2(NULL, unique=TRUE)
#' plotDistrTargets2(expLog2, NULL, colorsFrom="color_", rename=TRUE, namesFrom=NULL, FUNS = c(geom_density, geom_boxplot),
#'     probeTypeVec = probeTypeVec, probeTypeValue = 0, numProbes = NULL, orderByColors = FALSE,
#'     scale_x_limits = NULL, filePath = file.path(resultDir, "_debug3_plotDistrTargets2-noTargets.pdf"))
#' ## non-unique column names
#' exxx <- expLog2
#' colnames(exxx) <- names(getColPats2(targets, namesFrom=Birth3WHO, pullVar=Birth3WHO))
#' plotDistrTargets2(exxx, NULL, colorsFrom="color_", rename=TRUE, namesFrom=NULL, FUNS = c(geom_density, geom_boxplot),
#'     probeTypeVec = probeTypeVec, probeTypeValue = 0, numProbes = NULL, orderByColors = FALSE,
#'     scale_x_limits = NULL, filePath = file.path(resultDir, "_debug3_plotDistrTargets2-euqalSampleNames-noTargets-.pdf"))
#' plotDistrTargets2(exxx, targets, colorsFrom="color_", rename=TRUE, namesFrom=NULL, FUNS = c(geom_density, geom_boxplot),
#'     probeTypeVec = probeTypeVec, probeTypeValue = 0, numProbes = NULL, orderByColors = FALSE,
#'     scale_x_limits = NULL, filePath = file.path(resultDir, "_debug3_plotDistrTargets2-euqalSampleNames.pdf"))
#' }
#' @section Implementation:
#' See \url{https://www.tidyverse.org/blog/2020/02/glue-strings-and-tidy-eval/} for bracing variables.
#' Evaluation of expression: \url{https://adv-r.hadley.nz/evaluation.html}.
#' For \code{bquote(a +. (b))}: \url{http://adv-r.had.co.nz/Expressions.html}.
#' Errors: - dplyr::rename(!!purrr::set_names(oldColNames, newColNames))
#'         - colnames(expLog2) <- newColNames
#' TODO: -  titlePfx <- paste("Distribution of", as.character(rlang::fn_fmls()$'expLog2'), "on", numProbes, "probes", titleSfx)
#' @import ggplot2
#' @importFrom assertthat assert_that are_equal
#' @importFrom magrittr %>% %<>%
#' @importFrom tibble add_column tibble as_tibble
#' @importFrom dplyr filter sample_n right_join mutate
#' @importFrom tidyr pivot_longer
#' @importFrom tidyselect everything
#' @importFrom rlang enquo eval_tidy quo_name quo_is_null
#' @importFrom forcats fct_reorder
#' @importFrom purrr set_names
#' @importFrom vctrs vec_as_names
#' @importFrom myHelpers defactorChr
#' @export
plotDistrTargets2 <- function(expLog2, targets, colorsFrom = "color_", rename = TRUE, namesFrom = NULL, 
    FUNS = c(geom_density, geom_boxplot), probeTypeVec = NULL, probeTypeValue = 0, numProbes = NULL, 
    orderByColors = TRUE, scale_x_limits = NULL, filePath = NULL, width=16/9*7, height=7, ...) {
    # enquos
    nExpLog2 <- rlang::quo_name(rlang::enquo(expLog2))
    namesFrom <- rlang::enquo(namesFrom)
    # asserts
    assertthat::assert_that(is.matrix(expLog2))
    assertthat::are_equal(dim(expLog2)[[1]], length(probeTypeVec))
    assertthat::assert_that(is.list(FUNS))
    assertthat::assert_that(is.logical(orderByColors))
    # parameters and globals
    colnames(expLog2) <- vctrs::vec_as_names(colnames(expLog2), repair="unique")  # needed because of buggy pivot_longer(names_repair="unique")
    nNamesFrom <- rlang::quo_name(namesFrom)
    FUNS <- setNames(FUNS, as.character(1:length(FUNS)))
    if (is.null(probeTypeVec)) { 
        probeTypeVec <- rep(probeTypeValue, dim(expLog2)[[1]])
        titleSfx <- "" 
    } else titleSfx <- paste("of type", probeTypeValue, collapse=" ")
    numProbes <- min(numProbes, dim(expLog2)[[1]], sum(probeTypeVec == probeTypeValue))
    # globals
    colPats <- getColPats2(targets, colorsFrom=colorsFrom, rename=rename, namesFrom={{namesFrom}})
    colPatsUn <- getColPats2(targets, colorsFrom=colorsFrom, rename=rename, namesFrom={{namesFrom}}, unique=TRUE)
    plots <- list()
    d1 <- expLog2 %>%
        tibble::as_tibble() %>%
        tibble::add_column(probeType = probeTypeVec, .before=1) %>%
        dplyr::filter(probeType == probeTypeValue) %>%
        dplyr::select(-probeType) %>% 
        dplyr::sample_n(numProbes) %>% 
        tidyr::pivot_longer(tidyselect::everything(), names_to="Sample", values_to="Value", names_repair="unique")
    titlePfx <- paste("Distribution of", nExpLog2, "on", numProbes, "probes", titleSfx)
    for (nFUN in names(FUNS)) {
        if (length(colPats) > 0) {
            for (nColPat in names(colPats)) {
                d2 <- colPats[[nColPat]] %>%
                    tibble::tibble(nColPat = names(.), colors=.)
                # recode d1
                d1r <- d1 %>% 
                    dplyr::mutate(across(matches("Sample"), 
                        ~dplyr::recode(.x, !!!setNames(pluck(d2, "nColPat"), vctrs::vec_as_names(colnames(expLog2), repair="unique")))))
                d2 %<>%
                    dplyr::right_join(d1r, by=c("nColPat" = "Sample")) %>% 
                    dplyr::mutate(across(matches(nColPat), ~factor(.x, levels=unique(.x)))) %>% 
                    dplyr::rename(!!purrr::set_names(c("nColPat"), c(nColPat)))
                if (orderByColors) d2 %<>% 
                    dplyr::mutate(across(matches(nColPat), ~forcats::fct_reorder(.x, as.numeric(d2$colors))))
                plots[[paste0(nFUN, nColPat)]] <- d2 %>% 
                    ggplot(aes_string(x="Value", color=nColPat)) + 
                        FUNS[[nFUN]](...) +  
                        ggtitle(paste0(titlePfx, ", colors ", nColPat)) +
                        ggplot2::scale_color_manual(values = myHelpers::defactorChr(colPatsUn[[nColPat]]))
                # re-set legend title if using namesFrom 
                if (rename & !rlang::quo_is_null(namesFrom))
                    plots[[paste0(nFUN, nColPat)]] <- plots[[paste0(nFUN, nColPat)]] + 
                        guides(color=guide_legend(title=nNamesFrom))
            }
        } else {
            plots[[nFUN]] <- d1 %>% 
                ggplot(aes_string(x="Value", color="Sample")) + 
                    FUNS[[nFUN]](...) + 
                    ggtitle(titlePfx)
        }
    }
    if (!is.null(scale_x_limits))
        for (pn in names(plots)) plots[[pn]] <- plots[[pn]] + scale_x_continuous(limits=scale_x_limits)
    if (!is.null(filePath)) {
        pdf(filePath, width=width, height=height)
        for (p1 in plots) print(p1)
        dev.off()
    }
    invisible(plots)
}


#' Plot heatmap(s) of gene expression using different distance calculation methods.
#' 
#' @param expLog2 Expression matrix (preferably in log2 scale) with genes in rows and samples in named columns
#' @param targets Phenotype data with columns with colors, e.g., color_Sex
#' @param methods List of methods for distance calulation, individually passed to stats::dist(..., method)
#' @param colorsFrom String passed to getColPats2
#' @param namesFrom Variable passed to getColPats2
#' @param filePath NULL or character; if given, output a PDF; default NULL
#' @param width PDF width
#' @param height PDF height
#' @param treeheight_row Passed to pheatmap; default 0; 50 for pheatmap
#' @param ... Passed to pheatmap
#' @return ggplot2 object and PDF if filePath is given
#' @importFrom assertthat assert_that
#' @importFrom magrittr %>%
#' @importFrom rlang enquo
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @importFrom pheatmap pheatmap
#' @importFrom grid grid.newpage grid.draw
#' @export
pheatmapTargets2 <- function(expLog2, targets, methods=c("manhattan", "euclidean"),
    colorsFrom="color_", namesFrom=NULL, filePath=NULL, width=7, height=7, treeheight_row = 0, ...) {
    require(pheatmap)
    require(stringr)
    require(rlang)
    require(gtable)
    assertthat::assert_that(is.matrix(expLog2))
    namesFrom <- rlang::enquo(namesFrom)
    p <- list()
    for (method in methods) {
        dists <- as.matrix(dist(t(expLog2), method=method))
        rownames(dists) <- colnames(expLog2)
        colnames(dists) <- NULL
        diag(dists) <- NA
        hmcol <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "YlOrRd"))(255)
        annRows <- getColPats2(targets, colorsFrom=colorsFrom, namesFrom={{namesFrom}})
        annRowsN <- as.data.frame(lapply(annRows, FUN=names))
        row.names(annRowsN) <- colnames(expLog2)
        p[[method]] <- pheatmap::pheatmap(
           dists,
           col = (hmcol), 
           annotation_row = annRowsN,
           annotation_colors = getColPats2(targets, unique=TRUE, colorsFrom=colorsFrom, namesFrom={{namesFrom}}),
           treeheight_row = treeheight_row,
           legend_breaks = c(min(dists, na.rm = TRUE), max(dists, na.rm = TRUE)), 
           legend_labels = (c("similar", "diverse")),
           main = paste(str_to_title(method), "heatmap for", dim(expLog2)[[1]], "probes"), 
           silent=TRUE,
           ...)$gtable
    }
    ## PDF
    if (!is.null(filePath)) {
        pdf(filePath, width=width, height=height)
        for (p1 in p) {
            grid::grid.newpage()
            grid::grid.draw(p1)
        }
        dev.off()
    }
    ## return
    invisible(p)
}


#' Plot PCA for gene expression (expLog2) and phenotype data (targets)
#' 
#' Colnames of expLog2 should be equal as rownames from targets.
#' Up to 4 variables from targets can be used for color, fill, shape and size.
#' Colors and fill of points is collected from targets using variable names with a prefix colorsFrom (default: "color_").
#' Continuous variables may be used; colors from targets are nor used thus.
#' Shape for a continuous variable must not be used together with fill.
#' 
#' @param expLog2 Expression matrix (preferably in log2 scale) with genes in rows and samples in columns
#' @param targets Phenotype data with columns shape, color, fill, size, 
#'                as well as color_color and color_fill, e.g., Sex and color_Sex
#' @param shape Variable from targets for ggplot2::aes
#' @param color Variable from targets for ggplot2::aes
#' @param fill Variable from targets for ggplot2::aes
#' @param size Variable from targets for ggplot2::aes
#' @param colorsFrom Character prefix of variable names from targets with colors
#' @param filePath NULL or character; if given, output a PDF; default NULL
#' @param width PDF width
#' @param height PDF height
#' @param stroke Line thickness, passed to geom_point
#' @param alpha Transparency of points, passed to geom_point
#' @param sizeUniform Integer point size, used when size parameter is NULL; default 4
#' @param sizeRange Integer c(min,max) for scaling point size; used when size parameter is given
#' @param ... Passed to geom_point
#' @return ggplot2 object and PDF if filePath is given
#' @section TODO:
#'  FIX: add guides_*
#'  FIX: Discrete color_targets supplied for continuous color variable; color_targets not used
#'  Fix scale_size: Error: Discrete value supplied to continuous scale
#'  Consider using scale_colour_viridis_c
#' @section Implementation:
#' Put ggplot() inside {} to be able to access data using dot (.), e.g. mutata(...) %>% { ggplot(., ...) }
#' Access data of an existing ggplot using .$data inside {}, e.g. { . + scale_color_manual(..., .$data) }
#' @import ggplot2
#' @importFrom assertthat assert_that are_equal
#' @importFrom magrittr %>% %<>%
#' @importFrom tibble as_tibble
#' @importFrom dplyr mutate pull
#' @importFrom rlang enquo eval_tidy quo_is_null
#' @importFrom purrr map_lgl
#' @importFrom myHelpers defactorChr
#' @export
plotPCAtargets2 <- function(expLog2, targets, shape = NULL, color = NULL, fill = NULL, size = NULL,
    colorsFrom = "color_", filePath = NULL, width = 7, height = 7, 
    stroke = 1, alpha = 0.75, sizeUniform = 4, sizeRange = c(3,5), ...) {
    assertthat::are_equal(dim(expLog2)[[2]], dim(targets)[[1]])
    assertthat::are_equal(colnames(expLog2), rownames(targets))
    ## enquo vars for aes
    varsAes <- rlang::enquos(shape = shape, color = color, fill = fill, size = size, .ignore_empty = "all")
    varsAesNotNull <- varsAes[purrr::map_lgl(varsAes, ~!rlang::quo_is_null(.x))]
    callAes <- rlang::call2("aes", !!!varsAesNotNull)
    # PCA
    PCA <- prcomp(t(expLog2), scale = FALSE)
    percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
    sd_ratio <- sqrt(percentVar[2] / percentVar[1])
    ## plot
    p <- targets %>% 
        tibble::as_tibble() %>% 
        dplyr::mutate(PC1 = PCA$x[,1], PC2 = PCA$x[,2]) %>% 
    {
        ggplot(., aes(PC2, PC1)) +
        geom_point(mapping=rlang::eval_tidy(callAes), stroke=stroke, alpha=alpha, ...) +
        ggtitle("PCA plot") +
        ylab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
        xlab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
        theme(plot.title = element_text(hjust = 0.5))+
        coord_fixed(ratio = sd_ratio)
    }
    # shape
    if (!rlang::quo_is_null(rlang::enquo(shape))) {
        if (!is.numeric(pull(p$data, !!rlang::enquo(shape))))
            p %<>%  { . + 
                scale_shape_manual(values = c(21:25, 0:14)[1:length(unique(pull(.$data, {{shape}})))])
            }
        else {
            assertthat::assert_that(rlang::quo_is_null(rlang::enquo(fill)), msg = "Using fill together with shape for a continuous variable is not implememnted.")
            warn("Warning: shape for a continuous variable must not be used together with fill.")
            p <- p + scale_shape_binned(solid = FALSE)
        }
        p <- p + guides(shape = guide_legend(override.aes = list(alpha = 1, size = sizeUniform)))
    }
    else 
        p$layers[[1]]$aes_params$shape = 21
    # colors
    if (!rlang::quo_is_null(rlang::enquo(color))) {
        if (!is.numeric(dplyr::pull(p$data, !!rlang::enquo(color))))
            p %<>% { . + 
                scale_color_manual(values = myHelpers::defactorChr(getColPats2(.$data, colorsFrom=colorsFrom, unique=TRUE, pullVar={{color}})))
            }
        else 
            warning("Cannot use discrete colors for continuous color variable")
        p <- p + guides(color = guide_legend(override.aes = list(alpha = 1, size = sizeUniform, shape=21)))
    }
    # fill
    if (!rlang::quo_is_null(rlang::enquo(fill))) {
        if (!is.numeric(dplyr::pull(p$data, !!rlang::enquo(fill)))) 
            p %<>% { . + 
                scale_fill_manual(values = myHelpers::defactorChr(getColPats2(.$data, colorsFrom=colorsFrom, unique=TRUE, pullVar={{fill}})))
            }
        else 
            warning("Cannot use discrete colors for continuous fill variable")
        p <- p + guides(fill = guide_legend(override.aes = list(alpha = 1, size = sizeUniform, shape=21)))
    }
    # size
    if (rlang::quo_is_null(rlang::enquo(size))) 
        p$layers[[1]]$aes_params$size = sizeUniform
    else  {
        if (is.numeric(dplyr::pull(p$data, !!rlang::enquo(size)))) 
            p <- p + scale_size(range = sizeRange)
        else  
            p <- p + scale_size_discrete(range = sizeRange)
        p <- p + guides(size = guide_legend(override.aes = list(alpha = 1, shape=21)))
    }
    # PDF
    if (!is.null(filePath)) {
        pdf(filePath, width=width, height=height)
        print(p)
        dev.off()
    }
    invisible(p)
}


#' Depricated: Plot PCA for gene expression and phenotype data
#' 
#' @rdname plotPCAtargets2
#' @param scale_color Variable from targets, named vector for scale_color_manual(value=scale_color)
#' @param scale_fill  Variable from targets, named vector for scale_fill_manual(value=scale_fill)
#' @param guides_fill Default "none"; use "legend" to display it
#' @return ggplot2 object and PDF if filePath is given
#' @section TODO:
#'  FIX: add return value
#'  FIX: Discrete scale_color supplied for continuous color variable; scale_color not used
#'  Fix: case that scale_color and scale_fill are not named vectors or NULL
#'  Fix scale_size: Error: Discrete value supplied to continuous scale
#'  Use scale_colour_viridis_c
#' @section Implementation:
#' Put ggplot() inside {} to be able to access data using dot (.)
#' @importFrom assertthat are_equal
#' @export
plotPCAtargets <- function(expLog2, targets, shape, color, fill, size,
    scale_color = NULL, scale_fill = NULL, filePath=NULL, width=7, height=7, stroke=1, 
    guides_fill = "none", alpha=0.5, ...) {
    require(rlang)
    require(ggplot2)
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
        tibble::as_tibble() %>% 
        mutate(PC1 = PCA$x[,1], PC2 = PCA$x[,2]) %>% 
    {
        ggplot(., aes(PC2, PC1)) +
        geom_point(aes(shape=!!shape, color=!!color, fill=!!fill, size=!!size), stroke=stroke, alpha=alpha, ...) +
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
    invisible(p)
}


#' Plot MDS(es) of gene expression using different distance calculation methods and alternative colors.
#' 
#' @param expLog2 Expression matrix (preferably in log2 scale) with genes in rows and samples in columns
#' @param targets Phenotype data with columns with colors, e.g., color_Sex
#' @param colorsFrom String variable prefix(es) from targets, passed to \code{getColPats2()}; default "color_"
#' @param namesFrom Variable from targets, passed to \code{getColPats2()}; default NULL, alternative suggested HybName
#' @inheritParams myHelpers::MDScols
#' @param FUNS Character vector MDS functions, passed individually to \code{MDScols()}
#' @param tops Integer vector number of rows (genes), passed individually to \code{MDScols()}
#' @param size Size of ggplot labels, passed to \code{geom_text()}
#' @param filePath NULL or character; if given, output a PDF; default NULL
#' @param width PDF width
#' @param height PDF height
#' @param ... Passed to ggplot2::geom_text
#' @return Invisibly a list of ggplot2 objects and PDF if filePath is given
#' @import dplyr
#' @importFrom tibble as_tibble
#' @importFrom ggplot2 ggplot aes geom_text labs theme element_text element_blank scale_color_manual
#' @importFrom myHelpers MDScols defactorChr
#' @export
#' @section TODO: add parameter dim.plot as in limma::plotMDS
plotMDStargets2 <- function(expLog2, targets, colorsFrom="color_", rename = TRUE, namesFrom=NULL, 
    scale=FALSE, center=FALSE, FUNS = c("cmdscale", "isoMDS", "sammon"), p = 2, selection = "pairwise", tops = 500,
    maxit = 50, trace = FALSE, tol = 1e-3, size=3, filePath=NULL, width=7, height=7, ...) {
    # require(myHelpers)
    colPats <- getColPats2(targets, colorsFrom=colorsFrom, rename=rename, namesFrom={{namesFrom}})
    colPatsUn <- getColPats2(targets, colorsFrom=colorsFrom, rename=rename, namesFrom={{namesFrom}}, unique=TRUE)
    if (is.null(selection)) tops <- dim(expLog2)[[1]]
    plots <- list()
    fits <- list()
    for (FUN in FUNS)
    for (top in tops) {
        fits[[paste0(FUN, top)]] <- MDScols(expLog2, scale=scale, center=center, FUN=FUN, p=p, 
            selection=selection, top=top, k=2, maxit=maxit, trace=trace, tol=tol, plot=FALSE)
        for (nColPat in names(colPats)) {
            plots[[paste0(FUN, top, nColPat)]] <- fits[[paste0(FUN, top)]] %>% 
                tibble::as_tibble() %>% 
                mutate(color = colPats[[nColPat]],
                       name = names(colPats[[nColPat]])) %>% 
                ggplot(aes(x=V1, y=V2, color=name, label=name)) + 
                geom_text(show.legend=FALSE, size=size, ...) +
                labs(title = paste(FUN, dim(expLog2)[[2]], "objects,", min(top, dim(expLog2)[[1]]), selection, "parameters,", nColPat, "colors")) + 
                theme(plot.title = element_text(hjust = 0.5),
                      axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
                      axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
                scale_color_manual(values = defactorChr(colPatsUn[[nColPat]]))
        }
    }
    ## PDF
    if (!is.null(filePath)) {
        pdf(filePath, width=width, height=height)
        for (pt in plots) print(pt)
        dev.off()
    }
    invisible(plots)
}


#' Depricated: Plot MDS(es) of gene expression using different distance calculation methods and alternative colors.
#' 
#' @rdname plotMDStargets2
#' @param expLog2 Expression matrix (preferably in log2 scale) with genes in rows and samples in columns
#' @param targets Phenotype data with columns with colors, e.g., color_Sex
#' @param colorsFrom String variable prefix(es) from targets, passed to \code{getColPats2()}; default "color_"
#' @param namesFrom Variable from targets, passed to \code{getColPats2()}; suggested HybName
#' @param scale Logical scale data, passed to \code{MASS_MDScols()}
#' @param center Logical center data, passed to \code{MASS_MDScols()}
#' @param methods List of methods for distance calulation, passed individually to \code{stats::dist(method)} through \code{MASS_MDScols()}
#' @param FUNS List of MDS function from pacakge MASS, passed individually to \code{MASS_MDScols()}
#' @param p Power of the Minkowski distance, passed to \code{dist()} through \code{MASS_MDScols()}
#' @param maxit Passed to \code{MASS_MDScols()}, 
#' @param trace Trace progress, passed to \code{MASS_MDScols()}, default FALSE
#' @param tol Tolerance, passed to \code{MASS_MDScols()}
#' @param size Size of ggplot labels, passed to \code{geom_text()}
#' @param filePath NULL or character; if given, output a PDF; default NULL
#' @param width PDF width
#' @param height PDF height
#' @param ... Passed to ggplot2::geom_text
#' @return Invisibly a list of ggplot2 objects and PDF if filePath is given
#' @import dplyr
#' @importFrom ggplot2 aes labs theme element_text element_blank
#' @importFrom myHelpers MASS_MDScols defactorChr
#' @export
plotMDStargets2_MASS <- function(expLog2, targets, colorsFrom="color_", namesFrom=NULL, 
    scale=FALSE, center=FALSE, methods=c("euclidean", "manhattan"), FUNS = c("isoMDS", "sammon"),
    p = 2, maxit = 50, trace = FALSE, tol = 1e-3, size=3, filePath=NULL, width=7, height=7, ...) {
    # require(myHelpers)
    colPats <- getColPats2(targets, colorsFrom=colorsFrom, namesFrom={{namesFrom}})
    colPatsUn <- getColPats2(targets, colorsFrom=colorsFrom, namesFrom={{namesFrom}}, unique=TRUE)
    plots <- list()
    fits <- list()
    for (method in methods) {
        for (FUN in FUNS) {
            fits[[paste0(method, FUN)]] <- MASS_MDScols(expLog2, scale=scale, center=center, method=method, 
                FUN=FUN, p=p, k=2, maxit=maxit, trace=trace, tol=tol, plot=FALSE)
            for (nColPat in names(colPats)) {
                plots[[paste0(method, FUN, nColPat)]] <- fits[[paste0(method, FUN)]] %>% 
                    tibble::as_tibble() %>% 
                    mutate(color = colPats[[nColPat]],
                           name = names(colPats[[nColPat]])) %>% 
                    ggplot(aes(x=V1, y=V2, color=name, label=name)) + 
                    geom_text(show.legend=FALSE, size=size, ...) +
                    labs(title = paste(FUN, method, dim(expLog2)[[2]], "objects,", dim(expLog2)[[1]], "parameters", nColPat, "colors")) + 
                    theme(plot.title = element_text(hjust = 0.5),
                          axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
                          axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
                    scale_color_manual(values = defactorChr(colPatsUn[[nColPat]]))
            }
        }
    }
    ## PDF
    if (!is.null(filePath)) {
        pdf(filePath, width=width, height=height)
        for (pt in plots) print(pt)
        dev.off()
    }
    invisible(plots)
}


###########################################################
#### Plot distributions from limma::RGList style data #####
###########################################################

#' Plot distributions of signals in a RGList object using multiple colors for groups of samples
#' 
#' Plots distributions for each combination of FUNS, channels and sampleColors.
#' By default, samples are ordered according to the sampleColors (if provided).
#' Colors should be provided as factors with a preset order of levels; the order of colors is determined by the order of levels
#' 
#' @param RGList A list of matrices of signal intensities per channel
#' @param filePath NULL or character; if given, output a PDF; default NULL
#' @param channels List of data.frames from RGList, e.g. list("log2(R)", "log2(G)")
#' @param FUNS geom_functions from ggplot2, default \code{c(geom_density, geom_boxplot)}, others: \code{c(geom_violin, geom_histogram, ggridges::geom_density_ridges)}
#' @param probeTypeVec Vector of probe types
#' @param probeTypeValue Value from probeTypeVec to use for plotting
#' @param numProbes Number of randomly select probes
#' @param sampleColors Tibble of named factors of colors with names matching sample names
#' @param orderByColors Logical for ordering samples by colors (coded as a factor with levels in order); default TRUE
#' @param scale_x_limits NULL for auto-scale; use c(0,16) or less for log2(intensities)
#' @param width PDF width, default 16/9*7=12.44
#' @param height PDF height, default 7
#' @param ... Passed to ggplot2::FUN, e.g.: bins, binwidth, show.legend
#' @return
#' A list of ggplots, one per a combination of FUNS, channels and sampleColors
#' PDF if filePath is not NULL
#' @examples
#' \dontrun{
#' plotDistr_RGList(dataRG.bgc, file.path(curDir, "density-boxplot_dataRG.bgc_genes.pdf"), 
#' channels = c("log2R","log2G","beta"), probeTypeVec = dataRG.bgc$genes$ControlType, probeTypeValue = 0,
#' numProbes = 10000, sampleColors = targets %>% select(starts_with("color_")))
#' ## Removed parameters/code:
#' # transform = log2: transformation function applied to each channel (default: log2; identity for none)
#' mutate("{{transform}}({ch})" := transform(!!parse_expr(ch)))
#' ggplot(aes(x=transform(!!parse_expr(ch)), color=HybName)) +
#' ggplot(d2, aes(x=!!parse_expr(ch), color=fct_reorder(HybName, order(colors)))) + #, order=colors)) + 
#' }
#' @section Implementation:
#' See \url{https://www.tidyverse.org/blog/2020/02/glue-strings-and-tidy-eval/} for bracing variables.
#' Evaluation of expression: \url{https://adv-r.hadley.nz/evaluation.html}.
#' For \code{bquote(a +. (b))}: \url{http://adv-r.had.co.nz/Expressions.html}.
#' @import ggplot2
#' @importFrom assertthat assert_that are_equal
#' @importFrom tibble add_column tibble
#' @importFrom dplyr filter sample_n right_join mutate
#' @importFrom tidyr pivot_longer
#' @importFrom rlang parse_expr
#' @importFrom forcats fct_reorder
#' @export
plotDistr_RGList <- function(RGList, filePath = NULL, channels = c("R","G"), FUNS = c(geom_density, geom_boxplot),
    probeTypeVec = NULL, probeTypeValue = 0, numProbes = NULL, sampleColors = NULL, orderByColors = TRUE, 
    scale_x_limits = NULL, width=16/9*7, height=7, ...) {
    assertthat::assert_that(is(RGList, "RGList"))
    assertthat::assert_that(all(channels %in% names(RGList)))
    FUNS <- setNames(FUNS, as.character(1:length(FUNS)))
    if (is.null(probeTypeVec)) { probeTypeVec <- rep(probeTypeValue, dim(RGList[[1]])[[1]]); titleSfx <- "" } else { titleSfx <- paste("of type", probeTypeValue, collapse=" ") }
    assertthat::are_equal(dim(RGList[[1]])[[1]], length(probeTypeVec))
    numProbes <- min(numProbes, dim(RGList[[1]])[[1]], sum(probeTypeVec == probeTypeValue))
    assertthat::assert_that(is.null(sampleColors) |
                            assertthat::are_equal(colnames(RGList[[channels[[1]]]]), names(sampleColors[[1]])))
    assertthat::assert_that(is.logical(orderByColors))
    plots <- list()
    for (ch in channels) {
        d1 <- RGList[[ch]] %>%
                as_tibble %>%
                add_column(probeType = probeTypeVec, .before=1) %>%
                filter(probeType == probeTypeValue) %>%
                sample_n(numProbes) %>%
                pivot_longer(-probeType, names_to="HybName", values_to=ch)
        titlePfx <- paste("Distribution of ", ch, "on", numProbes, "probes", titleSfx)
        for (nFUN in names(FUNS)) {
            if (!is.null(sampleColors))
                for (cn in colnames(sampleColors)) {
                    d2 <- sampleColors[[cn]] %>%
                        tibble(HybName = names(.), colors=.) %>% 
                        right_join(d1, by="HybName") %>% 
                        mutate(HybName = factor(HybName, levels=unique(HybName)))
                    if (orderByColors) d2 %<>% 
                        mutate(HybName = fct_reorder(HybName, as.numeric(colors)))
                    plots[[paste(ch,nFUN,cn, sep="_")]] <- d2 %>% 
                        ggplot(aes(x=!!parse_expr(ch), color=HybName)) + 
                            FUNS[[nFUN]](...) +  
                            ggplot2::scale_color_manual(values=setNames(as.character(sampleColors[[cn]]), names(sampleColors[[cn]]))) +
                            ggtitle(paste(titlePfx, cn))
                }
            else 
                plots[[paste(ch,nFUN, sep="_")]] <- d1 %>% 
                    mutate(HybName = as_factor(HybName)) %>% 
                    ggplot(aes(x=!!parse_expr(ch), color=HybName)) + 
                        FUNS[[nFUN]](...) + 
                        ggtitle(titlePfx)
        }
    }
    if (!is.null(scale_x_limits))
        for (pn in names(plots)) plots[[pn]] <- plots[[pn]] + scale_x_continuous(limits=scale_x_limits)
    if (!is.null(filePath)) {
        pdf(filePath, width=width, height=height)
        for (p1 in plots) print(p1)
        dev.off()
    }
    invisible(plots)
}
 

#' Boxplot RLE (Relative Log Expression) and correlate to RIN
#' 
#' Correlate RLE to median and IQR of RIN.
#' If RIN is not giver, correlate it to straight line, i.e. test if RLEs are equal
#' 
#' @param expLog2 Expression matrix (preferably in log2 scale) with genes in rows and samples in columns
#' @param filePath NULL or character; if given, output a PDF; default NULL
#' @param RIN If given, correlate it to median and IQR of RLE; 
#'  otherwise correlate RLE to a straight line (RLE==1 for all samples)
#' @param width PDF width
#' @param height PDF height
#' @param alpha Transparency of points, passed to geom_point
#' @param ... Passed to geom_boxplot
#' @return ggplot2 object and PDF if filePath is given
#' @importFrom assertthat assert_that are_equal has_name
#' @export
boxplotRLE <- function(expLog2, filePath=NULL, RIN=NULL, width=7, height=7, alpha=0.5, ...) {
    require(Biobase)
    require(magrittr)
    require(RColorBrewer)
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
        geom_boxplot(outlier.shape = NA, alpha=alpha, fill=myCols, ...) +
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
    invisible(prle)
}


##############################
#### From Osijek_HFHSD.R #####
#### For KBlag_doxy_Cla.R ####
##############################

#' Add collapsed annotations to eset using columns SYMBOLSs, ENTREZIDs
#' 
#' Keep probes w/o annotations
#' 
#' @param eset ExpressionSet
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
#' @rdname getWriteHeatmap_probesTT
#' @param outPath Path to write TTs
#' @param esets Named list of esets for BiocGenerics::annotation of probes with columns PROBEID, SYMBOL, ENTREZID, HSA_ENTREZID
#' @param fits Named list of PGSEA fits
#' @param contMatrices Named list of contrast matrices with Levels & Contrasts
#' @returns Named list of \code{limma::topTable()} results (TTs) of DE genes; names corresponds to names of fits
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
#' @rdname getWriteHeatmap_probesTT
#' @param comparisons Named list of comparisons made from names of contrasts
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



#' Write enriched pathways from contrasts & return patways
#' 
#' @rdname getWriteHeatmap_PgseaTT
#' @param outPath Path to write TTS
#' @param gscPGSEA Named list of gene set collection
#' @param esetsPGSEA Named list of PGSEA esets
#' @param fitsPGSEA Named list of PGSEA fits
#' @param contMatrices Named list of contrast matrices with Levels & Contrasts
#' @param setIDCol Add a column to TT with set IDs and name it setIDCol
#' @param setName Name of the column from Biobase::fData(esetsPGSEA[[?]]) to show with heatmaps
#' @param fitsProbes Named list of probe fits
#' @param esetsProbes Named list of esets for BiocGenerics::annotation of probes with columns PROBEID, SYMBOL, ENTREZID, HSA_ENTREZID
#' @returns Named list of \code{limma::topTable()} results (TTs) of enriched sets; names corresponds to names of fits
#' @examples
#' \dontrun{
#' getWriteHeatmap_PgseaTTcontrasts(outPath=file.path(resultDirOut, "3.PGSEA.KEGGREST.limma-contrasts"), gscPGSEA=gscKeggrestHsa, esetsPGSEA=esetsKeggrest, fitsPGSEA=fitsKeggrest,
#'                                      setIDCol="KeggID", useNameCol="keggPathNameHsa2", fitsProbes=fits, esetsProbes=esets, pVals)
#' getWriteHeatmap_PgseaTTcontrasts(file.path(resultDirOut, "3.PGSEA.TRANSFAC.2016.1.byFA.limma-contrasts"), gscTF.facFA, esetsTF.facFA, fitsTF.facFA,
#'                                          setIDCol="facFA", useNameCol="factorFA", fits, esets, pVals)
#' }
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
#' @rdname getWriteHeatmap_PgseaTT
#' @param comparisons Named list of comparisons made from names of contrasts
#' @examples
#' \dontrun{
#' writeHeatmap_PgseaTTcomparisons(outPath=file.path(resultDirOut, "3.PGSEA.KEGGREST.limma-comparisons"), gscPGSEA=gscKeggrestHsa, esetsPGSEA=esetsKeggrest, fitsPGSEA=fitsKeggrest,
#'                           setIDCol="KeggID", useNameCol="keggPathNameHsa2", fitsProbes=fits, esetsProbes=esets, pVals, targetsOrder)
#' writeHeatmap_PgseaTTcomparisons(file.path(resultDirOut, "3.PGSEA.TRANSFAC.2016.1.byFA.limma-comparisons"), gscTF.facFA, esetsTF.facFA, fitsTF.facFA,
#'                           setIDCol="facFA", useNameCol="factorFA", fits, esets, pVals, targetsOrder)
#' }
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
#' @param organism Character organism 3-letter code, e.g. "mmu" or "hsa"
#' @param annPrb Data frame probe annotations
#' @param annKeggEntrez Data frame KEGG annotations
#' @param keggPathOrgIDs TOWRITE
#' @param design TOWRITE
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

        # getWriteAnnTTlogFCcounts <- function(filepath, tabList, tabNames=NULL, id="DB_ID", logFC="logFC", nameCol="SYMBOL", addCol=NA) {
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
        #         logFCs <- dplyr::full_join(logFCs, logFCmeans, by=c(id, nameCol))
        #     }
        #     if (is.na(addCol)) colnames(logFCs) <- c(id, nameCol, paste(names(tabNames), logFC)) else colnames(logFCs) <- c(id, nameCol, paste(rep(names(tabNames),each=2), c(addCol, logFC)))
        #     logFCs$CountBoth <- apply(logFCs[,paste(names(tabNames), logFC)], MARGIN=1, FUN=function(x){sum(!is.na(x))})
        #     logFCs$CountPos <-  apply(logFCs[,paste(names(tabNames), logFC)], MARGIN=1, FUN=function(x){sum(x>0 & !is.na(x))})
        #     logFCs$CountNeg <-  apply(logFCs[,paste(names(tabNames), logFC)], MARGIN=1, FUN=function(x){sum(x<0 & !is.na(x))})
        #     write.table(logFCs, filepath, sep="\t", quote=FALSE, row.names=FALSE)
        #     return(logFCs)
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
