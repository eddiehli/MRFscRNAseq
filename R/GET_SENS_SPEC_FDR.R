
#' Get sensitivity, specificity and FDR for the obtained simulation results.
#'
#' @param data A list that contains
#' A 4-dimensional expression array: (2 groups)*(num of replicates)*(num of genes)*(num of cell types)
#' and other information such as number of genes, number of DE genes, etc.
#' @param test_type Choose from "ttest", "hurdle" and "MRF"
#' @param results Simulation results. For "ttest" and "hurdle", a matrix of p_values. For "MRF". See getMRFDE().
#' @param alpha Significance level. Default=0.01
#' @return Returns sensitivity, specificity, and fdr
#' @import caret
#' @export
#'

get_sens_spec_fdr = function(data, test_type, results, alpha = 0.05) {
  if (test_type == "ttest") {
    results = get_noMRF_sens_spec_fdr(data, results, alpha = alpha)
  }
  else if (test_type == "MRF") {
    results = get_MRF_sens_spec_fdr(data, results, alpha = alpha)
  }
  return(results)
}

get_noMRF_sens_spec_fdr = function(data, results, alpha) {

  expression_mat = data$expression_mat
  n_group = data$n_group
  n_sample = data$n_sample
  n_gene = data$n_gene
  n_cell = data$n_cell
  n_de_gene = data$n_de_gene
  n_de_cell = data$n_de_cell
  true_states = data$states
  p_vals = results

  p_vals_adj = matrix(p.adjust(p_vals, method = "BH"), nrow = n_gene, ncol = n_cell)
  total_DE_genes_mat = array(as.numeric(p_vals_adj < alpha), dim = c(1, n_gene, n_cell))
  
  dat = ifelse(total_DE_genes_mat == 1, "A", "B")
  ref = ifelse(true_states == 1, "A", "B")
  
  # Sensitivity
  sensitivity = sensitivity(factor(dat), factor(ref))
  # Specificity
  specificity = specificity(factor(dat), factor(ref))
  # FDR
  fdr = 1- posPredValue(factor(dat), factor(ref))

  noMRF_results = c(sensitivity = sensitivity, specificity = specificity, fdr = fdr)
  return(noMRF_results)
}

get_MRF_sens_spec_fdr = function(data, MRFDE_results, alpha) {

  expression_mat = data$expression_mat
  n_group = data$n_group
  n_sample = data$n_sample
  n_gene = data$n_gene
  n_cell = data$n_cell
  n_de_gene = data$n_de_gene
  n_de_cell = data$n_de_cell
  g_g = data$g_g
  c_c = data$c_c
  MRFDE_results = results

  q_sorted = sort(as.vector(1-MRFDE_results$postDE), decreasing = F)
  index = which.min((cumsum(q_sorted)/seq(1, n_gene*n_cell, by = 1) <= alpha) == TRUE)
  posterior.threshold = 1 - q_sorted[index]
  posterior.threshold
  results_postDE = matrix(MRFDE_results$postDE, nrow = n_gene, ncol = n_cell)
  total_DE_genes_mat = array(as.numeric(results_postDE > posterior.threshold), dim = c(1, n_gene, n_cell))
  
  dat = ifelse(total_DE_genes_mat == 1, "A", "B")
  ref = ifelse(true_states == 1, "A", "B")
  
  # Sensitivity
  sensitivity = sensitivity(factor(dat), factor(ref))
  # Specificity
  specificity = specificity(factor(dat), factor(ref))
  # FDR
  fdr = 1- posPredValue(factor(dat), factor(ref))
  
  MRF_results = c(sensitivity = sensitivity, specificity = specificity, fdr = fdr)
  return(MRF_results)
}


