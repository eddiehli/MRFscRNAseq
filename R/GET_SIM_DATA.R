
#' Get simulated data
#'
#' @param gg_info A list that contains gene to gene network, DE genes and non-DE genes
#' @param random.seed A random seed
#' @param gamma_alpha Shape parameter of Gamma distribution
#' @param gamma_beta Scale parameter of Gamma distribution
#' @param tau DE effect
#' @return A simulated data set
#' @export
#' @examples
#' data = get_sim_data(gg_info, gamma_alpha = 2, gamma_beta = 4, tau = 5, random.seed = 2020)

get_sim_data = function(gg_info, gamma_alpha, gamma_beta, tau, random.seed) {

  # Set Random Seed
  set.seed(random.seed)

  # G_G Matrix
  g_g = gg_info$g_g
  gene_names_DE = gg_info$gene_names_DE
  gene_names_EE = gg_info$gene_names_EE
  gene_names_all = gg_info$gene_names_all
  n_gene = dim(g_g)[1]

  de_gene_index = which(rownames(g_g) %in% gene_names_DE)

  # DE Genes and Cell Types
  n_de_gene = length(gene_names_DE)
  n_de_cell = 10

  # C_C Matrix
  n_cell = 20
  c_c = matrix(0, nrow = n_cell, ncol = n_cell)
  c_c[1:n_de_cell, 1:n_de_cell] = 1

  de_cell_index = 1:n_de_cell

  # Latent States
  statesI = array(0, dim = c(1, n_gene, n_cell))
  statesI[1, de_gene_index, de_cell_index] = 1
  sum(statesI)
  sum(statesI[1,1:length(gene_names_DE),1:10])

  paraMRF = c(194, 0.2, 0.4)
  gamma = paraMRF[1]; beta_gene = paraMRF[2]; beta_cell = paraMRF[3]
  iter = 10

  states = array(NA, dim = c(1, n_gene, n_cell))
  for (it in 1:iter) {
    for (g in 1:n_gene) {
      for (c in 1:n_cell) {
        a = gamma + beta_gene*sum(statesI[1,,c]*g_g[g,]) + beta_cell*sum(statesI[1,g,]*c_c[c,])
        b = beta_gene*sum(1-statesI[1,,c]*g_g[g,]) + beta_cell*sum(1-statesI[1,g,]*c_c[c,])
        prob = exp(a)/(exp(a)+ exp(b))
        states[1,g,c] <- (prob >= runif(length(prob))) + 0
      }
    }
  }

  sum(states)
  sum(states[1,1:length(gene_names_DE),1:10])

  # Generate Expression Matrix
  n_group = 2
  n_sample = 100

  expression_mat = array(NA, dim = c(n_group, n_sample, n_gene, n_cell))

  for (g in 1:n_gene) {
    for (c in 1:n_cell) {
      mu = rgamma(1, gamma_alpha, gamma_beta)
      disp = 0.1
      expression_mat[1,,g,c] = rnbinom(n_sample, mu = mu, size = disp)
      if (states[1,g,c] == 1) {
        lambda = sample(c(1/tau, tau), 1, prob = c(0.5, 0.5), replace = T)
        expression_mat[2,,g,c] = rnbinom(n_sample, mu = lambda*mu, size = lambda*disp)
      } else {
        expression_mat[2,,g,c] = rnbinom(n_sample, mu = mu, size = disp)
      }
    }
  }

  sim_data = list(expression_mat = expression_mat,
                  n_group = n_group, n_sample = n_sample,
                  n_gene = n_gene, n_cell = n_cell, g_g = g_g, c_c = c_c,
                  n_de_gene = n_de_gene, n_de_cell = n_de_cell,
                  states = states)
  return(sim_data)
}











