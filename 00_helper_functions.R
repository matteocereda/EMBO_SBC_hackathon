
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y, use = 'complete.obs')) # specify complete obs to remove NAs
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}

# Adapted from MarcElosua https://github.com/drighelli/SpatialExperiment/issues/115
seurat_to_spe <- function(seu, sample_id, img_id) {
  require(SPIAT)
  ## Convert to SCE
  sce <- Seurat::as.SingleCellExperiment(seu)
  
  ## Extract spatial coordinates
  spatialCoords =  GetTissueCoordinates(object = seu, image = 'slice')[,1:2] %>% as.matrix()
  colnames(spatialCoords) = c('Cell.X.Position', 'Cell.Y.Position')
  ## Extract and process image data
  img <- SpatialExperiment::SpatialImage(
    x = as.raster(seu@images[[img_id]]@image))
  
  imgData <- DataFrame(
    sample_id = sample_id,
    image_id = img_id,
    data = I(list(img)),
    scaleFactor = seu@images[[img_id]]@scale.factors$lowres)
  
  # Convert to SpatialExperiment
  spe <- SpatialExperiment(
    assays = assays(sce),
    rowData = rowData(sce),
    colData = colData(sce),
    metadata = metadata(sce),
    reducedDims = reducedDims(sce),
    altExps = altExps(sce),
    sample_id = sample_id,
    spatialCoords = spatialCoords,
    imgData = imgData
  )
  # indicate all spots are on the tissue
  spe$in_tissue <- 1
  spe$sample_id <- sample_id
  # Return Spatial Experiment object
  spe
}


# Function readpted from https://rdrr.io/cran/clustermole/src/R/enrichment.R
clustermole_enrichment_cosmx = function (expr_mat, species, method = "gsva") 
{
  if (!is(expr_mat, "matrix")) {
    stop("expression matrix is not a matrix")
  }
  # if (nrow(expr_mat) < 5000) { this would not work with cosmx because we have ~ 900 genes) , so we meodfied this id=f in order to use this function on cosmx data either 
  
  #   stop("expression matrix does not appear to be complete (too few rows)")
  # }
  if (nrow(expr_mat) < 500) {
    stop("expression matrix does not appear to be complete (too few rows)")
  }
  if (ncol(expr_mat) < 5) {
    stop("expression matrix does not appear to be complete (too few columns)")
  }
  if (max(expr_mat) > 100) {
    stop("expression values do not appear to be log-scaled")
  }
  expr_mat <- expr_mat[rowMeans(expr_mat) > min(expr_mat), 
  ]
  markers_tbl <- clustermole_markers(species = species)
  markers_tbl <- dplyr::filter(markers_tbl, .data$gene %in% 
                                 rownames(expr_mat))
  markers_tbl <- dplyr::add_count(markers_tbl, .data$celltype_full, 
                                  name = "n_genes")
  markers_tbl <- dplyr::filter(markers_tbl, .data$n_genes >= 
                                 5)
  markers_list <- dplyr::distinct(markers_tbl, .data$celltype_full, 
                                  .data$gene)
  markers_list <- split(x = markers_list$gene, f = markers_list$celltype_full)
  celltypes_tbl <- markers_tbl %>% dplyr::select(!dplyr::starts_with("gene")) %>% 
    dplyr::distinct()
  scores_tbl <- get_scores(expr_mat = expr_mat, markers_list = markers_list, 
                           method = method)
  scores_tbl <- scores_tbl %>% dplyr::filter(.data$score_rank <= 
                                               100) %>% dplyr::inner_join(celltypes_tbl, by = "celltype_full") %>% 
    dplyr::arrange(.data$cluster, .data$score_rank)
  scores_tbl
}


# Function taken from https://rdrr.io/cran/clustermole/src/R/enrichment.R
get_scores <- function(expr_mat, markers_list, method = c("gsva", "ssgsea", "singscore", "all")) {
  method <- match.arg(method)
  
  if (method == "gsva" || method == "all") {
    gsva_param <- GSVA::gsvaParam(exprData = expr_mat, geneSets = markers_list, kcdf = "Gaussian")
    scores_mat <- GSVA::gsva(gsva_param, verbose = FALSE)
    scores_tbl <- lengthen_scores(scores_mat)
    scores_gsva_tbl <- dplyr::select(scores_tbl, "cluster", "celltype_full", score_rank_gsva = "score_rank")
  }
  
  if (method == "ssgsea" || method == "all") {
    ssgsea_param <- GSVA::ssgseaParam(exprData = expr_mat, geneSets = markers_list)
    scores_mat <- GSVA::gsva(ssgsea_param, verbose = FALSE)
    scores_tbl <- lengthen_scores(scores_mat)
    scores_ssgsea_tbl <- dplyr::select(scores_tbl, "cluster", "celltype_full", score_rank_ssgsea = "score_rank")
  }
  
  if (method == "singscore" || method == "all") {
    markers_gsc <- Map(function(x, y) GSEABase::GeneSet(x, setName = y), markers_list, names(markers_list))
    markers_gsc <- GSEABase::GeneSetCollection(markers_gsc)
    scores_mat <- singscore::multiScore(rankData = rankGenes(expr_mat), upSetColc = markers_gsc)
    scores_mat <- scores_mat$Scores
    scores_tbl <- lengthen_scores(scores_mat)
    scores_singscore_tbl <- dplyr::select(scores_tbl, "cluster", "celltype_full", score_rank_singscore = "score_rank")
  }
  
  if (method == "all") {
    # combine all scores into a single table
    scores_tbl <- scores_gsva_tbl
    scores_tbl <- dplyr::full_join(scores_tbl, scores_ssgsea_tbl, by = c("cluster", "celltype_full"))
    scores_tbl <- dplyr::full_join(scores_tbl, scores_singscore_tbl, by = c("cluster", "celltype_full"))
    # create a score matrix for easier stats
    scores_mat <- dplyr::select(scores_tbl, dplyr::starts_with("score_rank_"))
    scores_mat <- as.matrix(scores_mat)
    # calculate stats
    scores_tbl$score_ranks_min <- apply(scores_mat, 1, min)
    scores_tbl$score_ranks_mean <- round(apply(scores_mat, 1, mean), 3)
    scores_tbl$score_ranks_median <- round(apply(scores_mat, 1, median), 3)
    # set the average rank as the default rank
    # not using the median as ssGSEA and singscore ranks tend to correlate well
    scores_tbl$score_rank <- scores_tbl$score_ranks_mean
  }
  
  scores_tbl
}

# Function taken from https://rdrr.io/cran/clustermole/src/R/enrichment.R
lengthen_scores <- function(scores_mat) {
  scores_mat %>%
    round(10) %>%
    tibble::as_tibble(rownames = "celltype_full") %>%
    tidyr::gather(key = "cluster", value = "score", -"celltype_full") %>%
    dplyr::select("cluster", "celltype_full", "score") %>%
    dplyr::group_by(.data$cluster) %>%
    dplyr::mutate(score_rank = rank(desc(.data$score), ties.method = "first")) %>%
    dplyr::ungroup()
}
