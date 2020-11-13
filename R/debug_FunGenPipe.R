#### plotDistrTargets2
if (FALSE) {
## limmaList
ll <- list(R=c(1,2), G=c(2,4))
library(tidyverse)
ch <- rlang::quo(log2(R))
ch <- rlang::quo(R)
rlang::eval_tidy(ch, ll %>% as_tibble) 
rlang::expr_text(ch)
rlang::quo_name(ch)
#colPats <- getColPats2(targets)
colorsFrom="color_"
(colPats <- getColPats2(targets, colorsFrom=colorsFrom, namesFrom=HybName))
(colPats <- getColPats2(targets, colorsFrom=colorsFrom, namesFrom=Birth3WHO))
(colPats <- getColPats2(targets, colorsFrom=colorsFrom, rename=FALSE))
nColPat <- names(colPats)[[1]]
expLog2 <- log2(dataRG$R)[1:100,]
colnames(expLog2)
te <- expLog2 %>% as_tibble %>% sample_n(33)
colnames(te)
colnames(te) <- myHelpers::defactorChr(pluck(targets, "Birth3"))
probeTypeVec <- rep(0,100)
probeTypeValue = 0
numProbes = NULL
orderByColors = TRUE
FUNS = c(geom_density, geom_boxplot)

le <- te %>% tidyr::pivot_longer(tidyselect::everything(), names_to="HybName", values_to="Values")
print(myHelpers::defactorChr(pluck(le, "HybName")))
str(pluck(le, "HybName"))
oldColNames <- colnames(te)
newColNames <- sub("(.*)_mibi...", "\\1", oldColNames)
newColNames <- as.factor(newColNames)
#' ter <- te %>% dplyr::rename(!!!purrr::set_names(oldColNames, newColNames))
#' ter <- te %>% dplyr::rename_with(~ sub("(.*)_mibi...", "\\1", .x))
ter <- te
colnames(ter) <- newColNames
ter
(ter %>% tidyr::pivot_longer(tidyselect::everything(), names_to="HybName", values_to="Values"))$HybName[300:1000]
        dplyr::mutate("{{namesFrom}}" := dplyr::recode_factor({{namesFrom}}, setNames(colnames(.x), pluck(targets, nNamesFrom))))

tec <- te %>% tidyr::pivot_longer(tidyselect::everything(), names_to="HybName", values_to="Values")
dplyr::recode_factor(tec$HybName, !!!setNames(pluck(targets, "Birth3WHO"), colnames(te)))
dplyr::recode(tec$HybName, !!!setNames(pluck(targets, "Birth3WHO"), colnames(te)))
tec %>% dplyr::mutate(HybName = dplyr::recode(HybName, !!!setNames(pluck(targets, "Birth3WHO"), colnames(te))))

for (nColPat in names(colPats)) {
    print(colPats[[nColPat]] %>% tibble::tibble("{nColPat}" := names(.), colors=.))
}

#' @examples
expLog2 <- log2(dataRG$R)[1:100,]
## rename=TRUE (default)
## namesFrom=NULL,     orderByColors=TRUE
plotDistrTargets2(expLog2, targets, colorsFrom="color_", namesFrom=NULL, FUNS = c(geom_density, geom_boxplot),
    probeTypeVec = probeTypeVec, probeTypeValue = 0, numProbes = NULL, 
    scale_x_limits = NULL, filePath = file.path(resultDir, "_debug3_plotDistrTargets2-NULL.pdf"))
## namesFrom=HybName
plotDistrTargets2(expLog2, targets, colorsFrom="color_", namesFrom=HybName, FUNS = c(geom_density, geom_boxplot),
    probeTypeVec = probeTypeVec, probeTypeValue = 0, numProbes = NULL, 
    scale_x_limits = NULL, filePath = file.path(resultDir, "_debug3_plotDistrTargets2-HybName.pdf"))
## namesFrom=Birth3WHO
plotDistrTargets2(expLog2, targets, colorsFrom="color_", namesFrom=Birth3WHO, FUNS = c(geom_density, geom_boxplot),
    probeTypeVec = probeTypeVec, probeTypeValue = 0, numProbes = NULL, 
    scale_x_limits = NULL, filePath = file.path(resultDir, "_debug3_plotDistrTargets2-Birth3WHO.pdf"))
## namesFrom=HybName,  orderByColors=FALSE
plotDistrTargets2(expLog2, targets, colorsFrom="color_", namesFrom=HybName, FUNS = c(geom_density, geom_boxplot),
    probeTypeVec = probeTypeVec, probeTypeValue = 0, numProbes = NULL, orderByColors = FALSE,
    scale_x_limits = NULL, filePath = file.path(resultDir, "_debug3_plotDistrTargets2-HybName-noReorder.pdf"))
## rename=FALSE () (implies namesFrom=NULL)
## orderByColors=TRUE
plotDistrTargets2(expLog2, targets, colorsFrom="color_", rename=FALSE, namesFrom=NULL, FUNS = c(geom_density, geom_boxplot),
    probeTypeVec = probeTypeVec, probeTypeValue = 0, numProbes = NULL, 
    scale_x_limits = NULL, filePath = file.path(resultDir, "_debug3_plotDistrTargets2-noRename.pdf"))
## orderByColors=FALSE
plotDistrTargets2(expLog2, targets, colorsFrom="color_", rename=FALSE, namesFrom=NULL, FUNS = c(geom_density, geom_boxplot),
    probeTypeVec = probeTypeVec, probeTypeValue = 0, numProbes = NULL, orderByColors = FALSE,
    scale_x_limits = NULL, filePath = file.path(resultDir, "_debug3_plotDistrTargets2-noRename-noReorder.pdf"))
## warning: rename=FALSE $ namesFrom!=NULL
plotDistrTargets2(expLog2, targets, colorsFrom="color_", rename=FALSE, namesFrom=Birth2, FUNS = c(geom_density, geom_boxplot),
    probeTypeVec = probeTypeVec, probeTypeValue = 0, numProbes = NULL, orderByColors = FALSE,
    scale_x_limits = NULL, filePath = file.path(resultDir, "_debug3_plotDistrTargets2-warning.pdf"))
## w/o targets
getColPats2(NULL, unique=TRUE)
plotDistrTargets2(expLog2, NULL, colorsFrom="color_", rename=TRUE, namesFrom=NULL, FUNS = c(geom_density, geom_boxplot),
    probeTypeVec = probeTypeVec, probeTypeValue = 0, numProbes = NULL, orderByColors = FALSE,
    scale_x_limits = NULL, filePath = file.path(resultDir, "_debug3_plotDistrTargets2-noTargets.pdf"))
## non-unique column names
exxx <- expLog2
colnames(exxx) <- names(getColPats2(targets, namesFrom=Birth3WHO, pullVar=Birth3WHO))
plotDistrTargets2(exxx, NULL, colorsFrom="color_", rename=TRUE, namesFrom=NULL, FUNS = c(geom_density, geom_boxplot),
    probeTypeVec = probeTypeVec, probeTypeValue = 0, numProbes = NULL, orderByColors = FALSE,
    scale_x_limits = NULL, filePath = file.path(resultDir, "_debug4_plotDistrTargets2-euqalSampleNames-noTargets-.pdf"))
plotDistrTargets2(exxx, targets, colorsFrom="color_", rename=TRUE, namesFrom=NULL, FUNS = c(geom_density, geom_boxplot),
    probeTypeVec = probeTypeVec, probeTypeValue = 0, numProbes = NULL, orderByColors = FALSE,
    scale_x_limits = NULL, filePath = file.path(resultDir, "_debug4_plotDistrTargets2-euqalSampleNames.pdf"))
#' -------------

ff <- function(expLog2) {
    print(as.character(rlang::fn_fmls()))
    print(names(formals()))
    print(methods::formalArgs(sys.function(sys.parent())))
    print(rlang::quo_name(rlang::enquo(expLog2)))
    print(rlang::quo_name(rlang::quo(expLog2)))
    print(substitute(expLog2))
}
ff(c(1,2))
aa <- c(1,2)
ff(aa)
ff(exxx)



}
## WORKING w/o namesFrom=NULL
# #' Plot distributions of signals from limma RGList, MAList, EList or EListRaw objects 
# #' using multiple colors for groups of samples
# #' 
# #' Plots distributions for each combination of FUNS, channels and sampleColors.
# #' By default, samples are ordered according to the sampleColors (if provided).
# #' Colors should be provided as factors with a preset order of levels; the order of colors is determined by the order of levels
# #' 
# #' @param expLog2 Expression matrix (preferably in log2 scale) with genes in rows and samples in columns
# #' @param targets Phenotype data with columns with colors, e.g., color_Sex
# #' @param colorsFrom String variable prefix(es) from targets, passed to \code{getColPats2()}; default "color_"
# #' @param namesFrom Variable from targets for (re)setting names to colors; its values are used for guides;
# #' passed to \code{getColPats2()}; default NULL, alternative suggested HybName;
# # @param limmaList Object of class RGList, MAList, EList or EListRaw; a list of matrices of signal intensities or processed signals
# # @param channels Vector, components from limmaList object, possibly with transformations; default \code{c(log2(R), log2(G))}
# #' @param FUNS Vector, geom_functions from ggplot2, default \code{c(geom_density, geom_boxplot)}, alternatives \code{c(geom_violin, geom_histogram, ggridges::geom_density_ridges)}
# #' @param probeTypeVec Vector of probe types
# #' @param probeTypeValue Value from probeTypeVec to use for plotting
# #' @param numProbes Number of randomly select probes
# # @param sampleColors Tibble of named factors of colors with names matching sample names
# #' @param orderByColors Logical for ordering samples by colors (coded as a factor with levels in order); default TRUE
# #' @param scale_x_limits NULL for auto-scale; use c(0,16) or less for log2(intensities)
# #' @param filePath NULL or character; if given, output a PDF; default NULL
# #' @param width PDF width, default 16/9*7=12.44
# #' @param height PDF height, default 7
# #' @param ... Passed to ggplot2::FUN, e.g.: bins, binwidth, show.legend
# #' @return Invisibly a list of ggplot2 objects, one per a combination of FUNS and colors;
# #' PDF if filePath is given
# #' @examples
# #' 
# #' \dontrun{
# #' plotDistr_limmaList(dataRG.bgc, file.path(curDir, "density-boxplot_dataRG.bgc_genes.pdf"), 
# #' channels = c("log2R","log2G","beta"), probeTypeVec = dataRG.bgc$genes$ControlType, probeTypeValue = 0,
# #' numProbes = 10000, sampleColors = targets %>% select(starts_with("color_")))
# #' ## Removed parameters/code:
# #' # transform = log2: transformation function applied to each channel (default: log2; identity for none)
# #' mutate("{{transform}}({ch})" := transform(!!parse_expr(ch)))
# #' ggplot(aes(x=transform(!!parse_expr(ch)), color=HybName)) +
# #' ggplot(d2, aes(x=!!parse_expr(ch), color=fct_reorder(HybName, order(colors)))) + #, order=colors)) + 
# #' }
# #' @section Implementation:
# #' See \url{https://www.tidyverse.org/blog/2020/02/glue-strings-and-tidy-eval/} for bracing variables.
# #' Evaluation of expression: \url{https://adv-r.hadley.nz/evaluation.html}.
# #' For \code{bquote(a +. (b))}: \url{http://adv-r.had.co.nz/Expressions.html}.
# #' @import ggplot2
# #' @importFrom tibble add_column tibble as_tibble
# #' @importFrom dplyr filter sample_n right_join mutate
# #' @importFrom tidyr pivot_longer
# #' @importFrom tidyselect everything
# #' @importFrom rlang enquos eval_tidy quo_name quo_is_null
# #' @importFrom forcats fct_reorder
# #' @importFrom myHelpers defactorChr
# #' @export
# plotDistrTargets2 <- function(expLog2, targets, colorsFrom = "color_", namesFrom = HybName, FUNS = c(geom_density, geom_boxplot),
#     probeTypeVec = NULL, probeTypeValue = 0, numProbes = NULL, orderByColors = TRUE, 
#     scale_x_limits = NULL, filePath = NULL, width=16/9*7, height=7, ...) {
#     FUNS <- setNames(FUNS, as.character(1:length(FUNS)))
#     if (is.null(probeTypeVec)) { probeTypeVec <- rep(probeTypeValue, dim(expLog2)[[1]]); titleSfx <- "" } else { titleSfx <- paste("of type", probeTypeValue, collapse=" ") }
#     assertthat::are_equal(dim(expLog2)[[1]], length(probeTypeVec))
#     numProbes <- min(numProbes, dim(expLog2)[[1]], sum(probeTypeVec == probeTypeValue))
# #    assertthat::assert_that(is.logical(orderByColors))

#     namesFrom <- enquo(namesFrom)
#     nNamesFrom <- rlang::quo_name(namesFrom)
#     colPats <- getColPats2(targets, colorsFrom=colorsFrom, namesFrom={{namesFrom}})
#     colPatsUn <- getColPats2(targets, colorsFrom=colorsFrom, namesFrom={{namesFrom}}, unique=TRUE)
#     print(colPats)

#     plots <- list()
#     message("=============namesFrom==============================")
#     print(rlang::quo_name(namesFrom))
#     print(nNamesFrom)
#     print(pluck(targets, nNamesFrom))
#     print(myHelpers::defactorChr(pluck(targets, nNamesFrom)))
#     # # if (!rlang::quo_is_null(namesFrom))
#     # #     oldColNames <- colnames(expLog2)
#     # #     newColNames <- myHelpers::defactorChr(pluck(targets, rlang::quo_name(namesFrom)))
#     # if (!rlang::quo_is_null(namesFrom))
#     #     newColNames <- pluck(targets, nNamesFrom)
#     #     colnames(expLog2) <- newColNames
#     message("=============expLog2==============================")
#     print(head(expLog2))
#     d1 <- expLog2 %>%
#         tibble::as_tibble() %>%
#         tibble::add_column(probeType = probeTypeVec, .before=1) %>%
#         dplyr::filter(probeType == probeTypeValue) %>%
#         dplyr::select(-probeType) %>% 
#         dplyr::sample_n(numProbes) %>% 
#         tidyr::pivot_longer(tidyselect::everything(), names_to=nNamesFrom, values_to="Values", names_repair="minimal") %>% 
#         dplyr::mutate({{namesFrom}} := dplyr::recode({{namesFrom}}, !!!setNames(pluck(targets, nNamesFrom), colnames(expLog2))))
#         # dplyr::rename(!!purrr::set_names(oldColNames, newColNames)) %>% 
#     message("=============d1==============================")
#     print(tail(d1,33))
#     print(myHelpers::defactorChr(pluck(d1, nNamesFrom)))
#     titlePfx <- paste("Distribution on", numProbes, "probes", titleSfx)
#     message("=============d2s==============================")
#     for (nFUN in names(FUNS)) {
#         if (length(colPats) > 0) {
#             for (nColPat in names(colPats)) {
#                 d2 <- colPats[[nColPat]] %>%
#                     tibble::tibble("{{namesFrom}}" := names(.), colors=.)
#                 d2 %<>% 
#                     dplyr::right_join(d1, by=nNamesFrom) %>% 
#                     dplyr::mutate("{{namesFrom}}" := factor({{namesFrom}}, levels=unique({{namesFrom}})))
#                 if (orderByColors) d2 %<>% 
#                     dplyr::mutate("{{namesFrom}}" := forcats::fct_reorder({{namesFrom}}, as.numeric(colors)))
#                 print(tail(d2,33))
#                 plots[[paste(nFUN,nColPat, sep="_")]] <- d2 %>% 
#                     ggplot(aes(x=Values, color={{namesFrom}})) + 
#                         FUNS[[nFUN]](...) +  
#                         ggtitle(paste(titlePfx, nColPat)) +
#                         ggplot2::scale_color_manual(values = defactorChr(colPatsUn[[nColPat]]))
#             }
#         } else {
#             plots[[nFUN]] <- d1 %>% 
#                 ggplot(aes(x=Values, color={{namesFrom}})) + 
#                     FUNS[[nFUN]](...) + 
#                     ggtitle(titlePfx)
#         }
#     }
#     if (!is.null(scale_x_limits))
#         for (pn in names(plots)) plots[[pn]] <- plots[[pn]] + scale_x_continuous(limits=scale_x_limits)
#     if (!is.null(filePath)) {
#         pdf(filePath, width=width, height=height)
#         for (p1 in plots) print(p1)
#         dev.off()
#     }
#     invisible(plots)
# }
#         # if (!is.null(sampleColors))
#         #     for (cn in colnames(sampleColors)) {
#         #         d2 <- sampleColors[[cn]] %>%
#         #             tibble::tibble(HybName = names(.), colors=.) %>% 
#         #             dplyr::right_join(d1, by="HybName") %>% 
#         #             dplyr::mutate(HybName = factor(HybName, levels=unique(HybName)))
#         #         if (orderByColors) d2 %<>% 
#         #             dplyr::mutate(HybName = forcats::fct_reorder(HybName, as.numeric(colors)))
#         #         plots[[paste(ch,nFUN,cn, sep="_")]] <- d2 %>% 
#         #             ggplot(aes(x=!!parse_expr(ch), color=HybName)) + 
#         #                 FUNS[[nFUN]](...) +  
#         #                 ggplot2::scale_color_manual(values=setNames(as.character(sampleColors[[cn]]), names(sampleColors[[cn]]))) +
#         #                 ggtitle(paste(titlePfx, cn))
#         #     }
#         # else 
#         #     plots[[paste(ch,nFUN, sep="_")]] <- d1 %>% 
#         #         mutate(HybName = as_factor(HybName)) %>% 
#         #         ggplot(aes(x=!!parse_expr(ch), color=HybName)) + 
#         #             FUNS[[nFUN]](...) + 
#         #             ggtitle(titlePfx)
## WORKING with namesFrom=NULL
# plotDistrTargets2 <- function(expLog2, targets, colorsFrom = "color_", namesFrom = HybName, FUNS = c(geom_density, geom_boxplot),
#     probeTypeVec = NULL, probeTypeValue = 0, numProbes = NULL, orderByColors = TRUE, 
#     scale_x_limits = NULL, filePath = NULL, width=16/9*7, height=7, ...) {
#     FUNS <- setNames(FUNS, as.character(1:length(FUNS)))
#     if (is.null(probeTypeVec)) { probeTypeVec <- rep(probeTypeValue, dim(expLog2)[[1]]); titleSfx <- "" } else { titleSfx <- paste("of type", probeTypeValue, collapse=" ") }
#     assertthat::are_equal(dim(expLog2)[[1]], length(probeTypeVec))
#     numProbes <- min(numProbes, dim(expLog2)[[1]], sum(probeTypeVec == probeTypeValue))
# #    assertthat::assert_that(is.logical(orderByColors))

#     namesFrom <- enquo(namesFrom)
#     nNamesFrom <- rlang::quo_name(namesFrom)
#     colPats <- getColPats2(targets, colorsFrom=colorsFrom, namesFrom={{namesFrom}})
#     colPatsUn <- getColPats2(targets, colorsFrom=colorsFrom, namesFrom={{namesFrom}}, unique=TRUE)
#     print(colPats)

#     plots <- list()
#     message("=============namesFrom==============================")
#     print(rlang::quo_name(namesFrom))
#     print(nNamesFrom)
#     print(pluck(targets, nNamesFrom))
#     print(myHelpers::defactorChr(pluck(targets, nNamesFrom)))
#     message("=============expLog2==============================")
#     print(head(expLog2))
#     d1 <- expLog2 %>%
#         tibble::as_tibble() %>%
#         tibble::add_column(probeType = probeTypeVec, .before=1) %>%
#         dplyr::filter(probeType == probeTypeValue) %>%
#         dplyr::select(-probeType) %>% 
#         dplyr::sample_n(numProbes) %>% 
#         tidyr::pivot_longer(tidyselect::everything(), names_to=nNamesFrom, values_to="Values", names_repair="minimal")
#     if (!rlang::quo_is_null(namesFrom)) d1 %<>%  
#         dplyr::mutate({{namesFrom}} := dplyr::recode({{namesFrom}}, !!!setNames(pluck(targets, nNamesFrom), colnames(expLog2))))
#     message("=============d1==============================")
#     print(tail(d1,33))
#     print(myHelpers::defactorChr(pluck(d1, nNamesFrom)))
#     titlePfx <- paste("Distribution on", numProbes, "probes", titleSfx)
#     message("=============d2s==============================")
#     for (nFUN in names(FUNS)) {
#         if (length(colPats) > 0) {
#             for (nColPat in names(colPats)) {
#                 d2 <- colPats[[nColPat]] %>%
#                     tibble::tibble("{{namesFrom}}" := names(.), colors=.)
#                 d2 %<>% 
#                     dplyr::right_join(d1, by=nNamesFrom) %>% 
#                     dplyr::mutate(across(matches(nNamesFrom), ~factor(.x, levels=unique(.x))))
#                     # dplyr::mutate(across(matches(nNamesFrom), factor))
#                     # dplyr::mutate("{{namesFrom}}" := factor({{namesFrom}}, levels=unique({{namesFrom}})))
#                 print(tail(d2,33))
#                 if (orderByColors) d2 %<>% 
#                     dplyr::mutate(across(matches(nNamesFrom), ~forcats::fct_reorder(.x, as.numeric(d2$colors))))
#                     # dplyr::mutate("{{namesFrom}}" := forcats::fct_reorder({{namesFrom}}, as.numeric(colors)))
#                 print(tail(d2,33))
#                 plots[[paste(nFUN, nColPat, sep="_")]] <- d2 %>% 
#                     ggplot(aes(x=Values, color={{namesFrom}})) + 
#                         FUNS[[nFUN]](...) +  
#                         ggtitle(paste(titlePfx, nColPat)) +
#                         ggplot2::scale_color_manual(values = defactorChr(colPatsUn[[nColPat]]))
#             }
#         } else {
#             plots[[nFUN]] <- d1 %>% 
#                 ggplot(aes(x=Values, color={{namesFrom}})) + 
#                     FUNS[[nFUN]](...) + 
#                     ggtitle(titlePfx)
#         }
#     }
#     if (!is.null(scale_x_limits))
#         for (pn in names(plots)) plots[[pn]] <- plots[[pn]] + scale_x_continuous(limits=scale_x_limits)
#     if (!is.null(filePath)) {
#         pdf(filePath, width=width, height=height)
#         for (p1 in plots) print(p1)
#         dev.off()
#     }
#     invisible(plots)
# }
## WORKING
# plotDistrTargets2 <- function(expLog2, targets, colorsFrom = "color_", rename = TRUE, namesFrom = HybName, 
#     FUNS = c(geom_density, geom_boxplot),
#     probeTypeVec = NULL, probeTypeValue = 0, numProbes = NULL, orderByColors = TRUE, 
#     scale_x_limits = NULL, filePath = NULL, width=16/9*7, height=7, ...) {
#     FUNS <- setNames(FUNS, as.character(1:length(FUNS)))
#     if (is.null(probeTypeVec)) { probeTypeVec <- rep(probeTypeValue, dim(expLog2)[[1]]); titleSfx <- "" } else { titleSfx <- paste("of type", probeTypeValue, collapse=" ") }
#     assertthat::are_equal(dim(expLog2)[[1]], length(probeTypeVec))
#     numProbes <- min(numProbes, dim(expLog2)[[1]], sum(probeTypeVec == probeTypeValue))
# #    assertthat::assert_that(is.logical(orderByColors))

#     namesFrom <- enquo(namesFrom)
#     nNamesFrom <- rlang::quo_name(namesFrom)
#     colPats <- getColPats2(targets, colorsFrom=colorsFrom, rename=rename, namesFrom={{namesFrom}})
#     colPatsUn <- getColPats2(targets, colorsFrom=colorsFrom, rename=rename, namesFrom={{namesFrom}}, unique=TRUE)
#     print(colPats)

#     plots <- list()
#     message("=============namesFrom==============================")
#     print(rlang::quo_name(namesFrom))
#     print(myHelpers::defactorChr(pluck(targets, nNamesFrom)))
#     message("=============expLog2==============================")
#     print(head(expLog2))
#     d1 <- expLog2 %>%
#         tibble::as_tibble() %>%
#         tibble::add_column(probeType = probeTypeVec, .before=1) %>%
#         dplyr::filter(probeType == probeTypeValue) %>%
#         dplyr::select(-probeType) %>% 
#         dplyr::sample_n(numProbes) %>% 
#         tidyr::pivot_longer(tidyselect::everything(), names_to="Sample", values_to="Value", names_repair="minimal")
#     ## not needed; always rename using names from colors inside a loop
#     # if (rename & !rlang::quo_is_null(namesFrom)) d1 %<>%  
#     #     dplyr::mutate({{namesFrom}} := dplyr::recode({{namesFrom}}, !!!setNames(pluck(targets, nNamesFrom), colnames(expLog2))))
#     message("=============d1==============================")
#     print(tail(d1,33))
#     print(myHelpers::defactorChr(pluck(d1, "Sample")))
#     titlePfx <- paste("Distribution on", numProbes, "probes", titleSfx)
#     message("=============d2s==============================")
#     for (nFUN in names(FUNS)) {
#         if (length(colPats) > 0) {
#             for (nColPat in names(colPats)) {
#                 message(paste("=============", nColPat, "d2s=============================="))
#                 d2 <- colPats[[nColPat]] %>%
#                     tibble::tibble(nColPat = names(.), colors=.)
#                 print(tail(d2,33))
#                 d1r <- d1 %>% 
#                     dplyr::mutate(across(matches("Sample"), ~dplyr::recode(.x, !!!setNames(pluck(d2, "nColPat"), colnames(expLog2)))))
#                 print(tail(d1r,33))
#                 d2 %<>%
#                     dplyr::right_join(d1r, by=c("nColPat" = "Sample")) %>% 
#                     dplyr::mutate(across(matches(nColPat), ~factor(.x, levels=unique(.x)))) %>% 
#                     dplyr::rename(!!purrr::set_names(c("nColPat"), c(nColPat)))

#                     # dplyr::mutate(across(matches(nNamesFrom), factor))
#                     # dplyr::mutate("{{namesFrom}}" := factor({{namesFrom}}, levels=unique({{namesFrom}})))
#                 if (orderByColors) d2 %<>% 
#                     dplyr::mutate(across(matches(nColPat), ~forcats::fct_reorder(.x, as.numeric(d2$colors))))
#                     # dplyr::mutate("{{namesFrom}}" := forcats::fct_reorder({{namesFrom}}, as.numeric(colors)))
#                 ## moved to guides
#                 # if (rename & !rlang::quo_is_null(namesFrom)) d2 %<>%  
#                 #     dplyr::rename("{nColPats}" := nNamesFrom)
#                 print(tail(d2,33))
#                 plots[[paste0(nFUN, nColPat)]] <- d2 %>% 
#                     ggplot(aes_string(x="Value", color=nColPat)) + 
#                         FUNS[[nFUN]](...) +  
#                         ggtitle(paste0(titlePfx, ", colors ", nColPat)) +
#                         ggplot2::scale_color_manual(values = myHelpers::defactorChr(colPatsUn[[nColPat]]))
#                 # re-set legend title if using namesFrom 
#                 if (rename & !rlang::quo_is_null(namesFrom))
#                     plots[[paste0(nFUN, nColPat)]] <- plots[[paste0(nFUN, nColPat)]] + 
#                         guides(color=guide_legend(title=nNamesFrom))
#             }
#         } else {
#             plots[[nFUN]] <- d1 %>% 
#                 ggplot(aes_string(x="Value", color="Sample")) + 
#                     FUNS[[nFUN]](...) + 
#                     ggtitle(titlePfx)
#         }
#     }
#     if (!is.null(scale_x_limits))
#         for (pn in names(plots)) plots[[pn]] <- plots[[pn]] + scale_x_continuous(limits=scale_x_limits)
#     if (!is.null(filePath)) {
#         pdf(filePath, width=width, height=height)
#         for (p1 in plots) print(p1)
#         dev.off()
#     }
#     invisible(plots)
# }

#### getColPats2
if (FALSE) {
    ## getColPats2 <- function(targets, colorsFrom="color_", rename=TRUE, namesFrom=NULL, unique=FALSE, pullVar=NULL) {
    colorsFrom="color_"
    rename=TRUE
    namesFrom=NULL
    unique=FALSE
    pullVar=NULL
    targets %>% {select(., matches("3WHO"))}
    targets3WHO <- targets %>% {bind_cols(select(., -matches("^color_.*$")), select(., matches("color_Birth3WHO")))}
    targets3WHO %<>% tibble::as_tibble()
    namesFrom <- rlang::enquo(namesFrom)
    pullVar <- rlang::enquo(pullVar)

    # TOFIX assertthat::assert_that(quo_is_null(namesFrom) | assertthat::has_name(targets3WHO, !!namesFrom))
    cols <- targets3WHO %>% 
        dplyr::select(tidyselect::starts_with(colorsFrom)) %>% 
        dplyr::rename_with(~sub(colorsFrom, "", .))
    if (rename) {
        if (rlang::quo_is_null(namesFrom))
            names <- targets3WHO %>% dplyr::select(tidyselect::any_of(sub(colorsFrom, "", colnames(cols))))
        else
            names <- targets3WHO %>% dplyr::select(!!namesFrom)
        assertthat::assert_that(assertthat::are_equal(dim(cols)[[1]], dim(names)[[1]]) & 
            (dim(names)[[2]] %in% c(1,dim(cols)[[2]])),
            msg = paste("The number of colors and variables is not the same. Not all variables with a prefix", colorsFrom, "have an associated variable. Consider using namesFrom parameter."))
        ncols <- purrr::map2(cols, names, setNames)
    } else 
        ncols <- as.list(cols)
    if (!rlang::quo_is_null(pullVar)) {
        ncol <- ncols %>% tibble::as_tibble() %>% pull(!!pullVar)
        if (unique)
            ncol <- ncol[!duplicated(names(ncol))]
        return(ncol)
    } else {
        if (unique)
            ncols %<>% purrr::map(~keep(., !duplicated(names(.))))
        return(ncols)
    }
}

#### plotMDStargets 2
if (FALSE) {

sd <- dataRG$Mproc %>% as_tibble %>% slice_sample(n=1000)
pp2<-plotMDStargets2(sd, targets, trace=TRUE, filePath = file.path(resultDir, "_debug_MDS2.pdf"))
pp3<-plotMDStargets2(sd, targets, trace=TRUE, filePath = file.path(resultDir, "_debug_MDS3.pdf"), namesFrom=HybName)

# plotMDStargets bug p not numeric
sd <- dataRG$Mproc %>% as_tibble %>% slice_sample(n=10)
s2c(sd,"sd")
s2c(targets,"targets")
exit()
ff<-plotMDStargets2(sd, targets, methods=c("euclidean"))
cp <- getColPats2(targets, namesFrom=NULL)
cpu <- getColPats2(targets, namesFrom=NULL, unique=TRUE)
im <- isoMDScols(sd)

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

}

#### pheatmapTargets
if (FALSE) {
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
}

#### from CFGBC pipelines

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

