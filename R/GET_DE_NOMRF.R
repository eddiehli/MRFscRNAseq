
#' Get simulation results for three types of tests
#'
#' @param data A list that contains
#' A 4-dimensional expression array: (2 groups)*(num of replicates)*(num of genes)*(num of cell types)
#' and other information such as number of genes, number of DE genes, etc.
#' @param test_type Choose from "ttest" and "hurdle"
#' @param alpha significance level. Default=0.01
#' @return Returns simulation results
#' @import pscl
#' @export
#'

get_DE_noMRF = function(data, test_type, alpha = 0.01) {
  if (test_type == "ttest") {
    results = get_noMRF_results(data, alpha = alpha)
  }
  else if (test_type == "hurdle") {
    results = get_noMRF_MAST_results(data, alpha = alpha)
  }
  return(results)
}

get_noMRF_results = function(data, alpha) {

  expression_mat = data$expression_mat
  n_group = data$n_group
  n_sample = data$n_sample
  n_gene = data$n_gene
  n_cell = data$n_cell
  n_de_gene = data$n_de_gene
  n_de_cell = data$n_de_cell

  p_vals = matrix(NA, nrow = n_gene, ncol = n_cell)
  for (g in 1:n_gene) for (c in 1:n_cell) p_vals[g,c] = t.test(expression_mat[1,,g,c], expression_mat[2,,g,c])$p.value

  return(p_vals)
}

get_noMRF_MAST_results = function(data, alpha) {

  expression_mat = data$expression_mat
  n_group = data$n_group
  n_sample = data$n_sample
  n_gene = data$n_gene
  n_cell = data$n_cell
  n_de_gene = data$n_de_gene
  n_de_cell = data$n_de_cell

  p_vals = matrix(NA, nrow = n_gene, ncol = n_cell)
  for(c in 1:n_cell){
    for (g in 1:n_gene) {
      reg.data = data.frame(y = c(expression_mat[1,,g,c], expression_mat[2,,g,c]), x = c(rep(0,100), rep(1,100)))
      reg.model = hurdle(y ~ x, data = reg.data, dist = "negbin", zero.dist = "poisson")
      reg.summary = summary(reg.model)
      x = reg.summary$coefficients$count[2,3]
      y = reg.summary$coefficients$zero[2,3]
      if (is.na(x) == TRUE) { x = 0 }
      if (is.na(y) == TRUE) { y = 0 }
      test_stat = x^2 + y^2
      p_val = 1 - pchisq(test_stat, 2)
      p_vals[g,c] = p_val
    }
  }

  return(p_vals)
}

