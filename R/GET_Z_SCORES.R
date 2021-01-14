
#' Get z scores
#'
#' @param data A list that contains
#' A 4-dimensional expression array: (2 groups)*(num of replicates)*(num of genes)*(num of cell types)
#' and other information such as number of genes, number of DE genes, etc.
#' @return Returns z scores
#' @export
#'

get_z_scores = function(data) {
  
  expression_mat = data$expression_mat
  n_group = data$n_group
  n_sample = data$n_sample
  n_gene = data$n_gene
  n_cell = data$n_cell
  n_de_gene = data$n_de_gene
  n_de_cell = data$n_de_cell
  g_g = data$g_g
  c_c = data$c_c
  
  zz = array(NA, dim = c(1, n_gene, n_cell))
  for (g in 1:n_gene) for (c in 1:n_cell) zz[1,g,c] = qnorm(pt(t.test(expression_mat[1,,g,c], expression_mat[2,,g,c])$statistic,
                                                               df = sum(length(expression_mat[1,,g,c]), length(expression_mat[2,,g,c]))))
  return(zz)
}
