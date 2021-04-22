
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
  gene_gene = data$gene_gene
  cell_cell = data$cell_cell
  n_gene = dim(gene_gene)[1]
  n_cell = dim(cell_cell)[1]
  
  zz = array(NA, dim = c(1, n_gene, n_cell))
  for (g in 1:n_gene) for (c in 1:n_cell) zz[1,g,c] = qnorm(pt(t.test(expression_mat[1,,g,c], expression_mat[2,,g,c])$statistic,
                                                               df = sum(length(expression_mat[1,,g,c]), length(expression_mat[2,,g,c]))))
  return(zz)
}
