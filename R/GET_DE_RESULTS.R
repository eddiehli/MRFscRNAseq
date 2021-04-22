
#' Get sensitivity, specificity and FDR for the obtained simulation results.
#'
#' @param results A list that contains MRF results
#' @param alpha Significance level. Default = 0.05
#' @return Returns a binary matrix that corresponds to cell-type specific DE results. 
#' @export
#'

get_DE_results = function(MRF_Results, alpha = 0.05) {
  results = MRF_Results$postDE
  res = matrix(results, nrow = dim(results)[2], ncol = dim(results)[3])
  posterior_threshold_vec = apply(res, 2, get_posterior_threshold, alpha = alpha)
  DE_states = matrix(NA, nrow = nrow(res), ncol = ncol(res))
  for (i in 1:ncol(DE_states)) { DE_states[,i] = as.numeric(res[,i] >= posterior_threshold_vec[i]) }
  rownames(DE_states) = rownames(MRF_Results$gene_gene)
  colnames(DE_states) = rownames(MRF_Results$cell_cell)
  return(DE_states)
}

get_posterior_threshold = function(data, alpha = alpha) {
  q_sorted = sort(1-data, decreasing = F)
  index = which.min((cumsum(q_sorted)/seq(1, length(data), by = 1) <= alpha) == TRUE)
  posterior.threshold = 1 - q_sorted[index]
}

