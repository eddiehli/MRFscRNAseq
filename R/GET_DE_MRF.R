
#' Run the MRF model to estimate posterior probabilities of differential expression for each gene across each cell type
#'
#' @param data Summary statistics matrix, rows are genes, columns are cell types
#' @param g_g Gene to gene network matrix
#' @param c_c Cell to cell dependency matrix
#' @param nulltype Type of null hypothesis assumed in estimating f0, see locfdr package.Default is the MLE (nulltype=1)
#' @param df Degrees of freedom for fitting the estimated density, see locfdr package. Default df=15
#' @param iterEM Max number of iterations for the EM algorithm. Default=200
#' @param iterGibbsPost Number of Gibbs posterior samples. Default=20,000
#' @param brPost Number of burn-in for the posterior samples. Default=10,000
#' @return The estimated model parameters and the posterior probabilities of differential expression
#' \item{postDE}{Posterior probabilities of differential expression. A 2-dimensional array: (num of genes)*(num of cell types)}
#' \item{paraMRF}{Estimated model parameters}
#' \item{paraMRFTrace}{Trace of the estimated model parameters in the EM algorithm}
#' @import locfdr
#' @export
#'

get_DE_MRF = function(data, g_g, c_c,
                      nulltype = 1, df = 15,
                      iterEM = 200, iterGibbsPost = 20000, brPost = 10000) {

  n_gene = dim(g_g)[1]
  n_cell = dim(c_c)[1]

  # Check non-finite summary statistics
  z_scores = array(data, dim = c(1, n_gene, n_cell))
  z_scores[z_scores == Inf] = max(z_scores[is.finite(z_scores)])
  z_scores[z_scores == -Inf] = min(z_scores[is.finite(z_scores)])
  z_scores[is.na(z_scores)] = 0

  # Use locfdr to estimate f0 and f1
  results_locfdr = suppressWarnings(locfdr(z_scores, nulltype = nulltype, df = df, plot = 0))

  p0 = (results_locfdr$fp0)[3, 3]
  p1 = 1 - p0
  fdr = array(results_locfdr$fdr, dim = dim(z_scores))
  fy0 = fdr/p0
  fy1 = (1 - fdr)/p1

  ## Initialization
  mf1 = ((p1*fy1/(p1*fy1 + (1-p1)*fy0)) > runif(length(fy0))) + 0
  if ((sum(mf1)/prod(dim(mf1))) < 0.0005) return("The total number of DE states is too small for MRF.")
  paraMRF = rep(0,6)
  paraMRFTrace = c()
  converge = c()
  cat("\n")

  # EM Algorithm
  for (j in 1:iterEM) {
    cat("\r", "Estimating model parameters,", floor(j/iterEM*100), "%", "completed")
    mf1 = calmf1(meanf = mf1, fy1 = fy1, fy0 = fy0, paraMRF = paraMRF, g_g = g_g, c_c = c_c, n_gene = n_gene, n_cell = n_cell)
    w1 = calw(fy1, fy0, paraMRF, mf1, g_g = g_g, c_c = c_c, n_gene = n_gene, n_cell = n_cell)
    tmp = optimNR(paraIni = rep(0, 3), meanf = mf1, w1 = w1, alpha = 10^(-6), maxiter = 100, g_g = g_g, c_c = c_c, n_gene = n_gene, n_cell = n_cell)
    if ((tmp=="NA")[1]){
      return("NA")
    }
    paraMRF = as.numeric(tmp$para)
    paraMRFTrace = rbind(paraMRFTrace, paraMRF)
    converge = c(converge, tmp$converged)
  }

  statesI = ((p1*fy1/(p1*fy1 + (1-p1)*fy0)) > runif(length(fy0))) + 0
  pfdr1 = get_states_nopara(fy1, fy0, paraMRF, statesI, iterGibbsPost, brPost, skip = 5, g_g = g_g, c_c = c_c, n_gene = n_gene, n_cell = n_cell)

  # Save the results in a list
  names(paraMRF) = c("Gamma", "Beta_Gene", "Beta_Cell")
  colnames(paraMRFTrace) = names(paraMRF)
  rownames(paraMRFTrace) = 1:iterEM
  results = list(postDE = pfdr1, paraMRF = paraMRF, paraMRFTrace = paraMRFTrace)
  return(results)
}

calmf1 = function(meanf, fy1, fy0, paraMRF, g_g, c_c, n_gene, n_cell) {
  gamma = paraMRF[1]; beta_gene = paraMRF[2]; beta_cell = paraMRF[3]
  for (g in 1:n_gene) {
    for (c in 1:n_cell) {
      a = gamma + beta_gene*sum(meanf[1,,c]*g_g[g,]) + beta_cell*sum(meanf[1,g,]*c_c[c,])
      b = beta_gene*sum(1-meanf[1,,c]*g_g[g,]) + beta_cell*sum(1-meanf[1,g,]*c_c[c,])
      prob = exp(a)*fy1[1,g,c]/(exp(a)*fy1[1,g,c]+ exp(b)*fy0[1,g,c])
      meanf[1,g,c] <- (prob >= runif(length(prob))) + 0
    }
  }
  return(meanf)
}

calw = function(fy1, fy0, paraMRF, meanf, g_g, c_c, n_gene, n_cell) {
  gamma = paraMRF[1]; beta_gene = paraMRF[2]; beta_cell = paraMRF[3]
  w1 = meanf*0
  for (g in 1:n_gene) {
    for (c in 1:n_cell) {
      a = gamma + beta_gene*sum(meanf[1,,c]*g_g[g,]) + beta_cell*sum(meanf[1,g,]*c_c[c,])
      b = beta_gene*sum(1-meanf[1,,c]*g_g[g,]) + beta_cell*sum(1-meanf[1,g,]*c_c[c,])
      w1[1,g,c] = exp(a)*fy1[1,g,c]/(exp(a)*fy1[1,g,c]+ exp(b)*fy0[1,g,c])
    }
  }
  return(w1)
}

optimNR = function(paraIni = rep(0, 3), meanf, w1, alpha = 10^(-6), maxiter = 100, g_g, c_c, n_gene, n_cell) {
  w0 = 1 - w1
  iter = 1
  para = rep(1, 3)
  parapre = paraIni

  while (max(abs(para - parapre)) >= alpha & iter <= maxiter) {
    print(iter)
    if (iter != 1) {
      parapre = para
    }
    gamma = parapre[1]; beta_gene = parapre[2]; beta_cell = parapre[3]
    hes = matrix(0, nrow = 3, ncol = 3)
    del = rep(0, 3)
    for (g in 1:n_gene) {
      for (c in 1:n_cell) {
        a = gamma + beta_gene*sum(meanf[1,,c]*g_g[g,]) + beta_cell*sum(meanf[1,g,]*c_c[c,])
        b = beta_gene*sum(1-meanf[1,,c]*g_g[g,]) + beta_cell*sum(1-meanf[1,g,]*c_c[c,])
        g_a = sum(meanf[1,,c]*g_g[g,])
        g_b = sum(1-meanf[1,,c]*g_g[g,])
        c_a = sum(meanf[1,g,]*c_c[c,])
        c_b = sum(1-meanf[1,g,]*c_c[c,])
        denom = exp(a) + exp(b)

        ## Gamma: Intercept
        numer_gamma = exp(a)
        del[1] = del[1] + sum(w1[1,g,c]) - sum(numer_gamma/denom)
        ## Beta: Gene
        numer_gene = g_a*exp(a) + g_b*exp(b)
        del[2] = del[2] + sum(w1[1,g,c]*g_a) + sum(w0[1,g,c]*g_b) - sum(numer_gene/denom)
        ## Beta: Cell
        numer_cell = c_a*exp(a) + c_b*exp(b)
        del[3] = del[3] + sum(w1[1,g,c]*c_a) + sum(w0[1,g,c]*c_b) - sum(numer_cell/denom)

        ## Diagonal of the Hessian Matrix
        # gamma, gamma
        numer = exp(a) * exp(b)
        hes[1,1] = hes[1,1] + sum(numer/denom^2)
        # gene, gene
        term1 = ( g_a^2*exp(a) + g_b^2*exp(b) ) / denom
        term2 = - numer_gene^2 /denom^2
        hes[2,2] = hes[2,2] + sum(term1 + term2)
        # cell, cell
        term1 = ( c_a^2*exp(a) + c_b^2*exp(b) ) / denom
        term2 = - numer_cell^2 /denom^2
        hes[3,3] = hes[3,3] + sum(term1 + term2)

        ## Off-diagonal of the Hessian Matrix
        # gamma, gene
        term1 = g_a*exp(a) / denom
        term2 = - numer_gene * numer_gamma / denom^2
        hes[1,2] = hes[1,2] + sum(term1 + term2)
        # gamma, cell
        term1 = c_a*exp(a) / denom
        term2 = - numer_cell * numer_gamma / denom^2
        hes[1,3] = hes[1,3] + sum(term1 + term2)
        # gene, cell
        term1 = ( g_a*c_a*exp(a) + g_b*c_b*exp(b) ) / denom
        term2 = - numer_gene * numer_cell / denom^2
        hes[2,3] = hes[2,3] + sum(term1 + term2)
      }
    }
    hes = - hes
    hes = t(hes) + hes
    diag(hes) = diag(hes)/2
    if (is.na(det(hes))) return("NA")
    if (abs(det(hes)) <= 1/10^6) return("NA")
    para = parapre - solve(hes) %*% del
    iter = iter + 1
    if (sum(is.na(para))) return("NA")
  }
  converge = 1
  if (iter == (maxiter + 1)){
    print("Not converged, please increase alpha or maxiter")
    converge = 0
  }
  return(list(para = para, converge = converge))
}

get_states_nopara = function(fy1, fy0, paraMRF, statesI, iterGibbs, br, skip, g_g, c_c, n_gene, n_cell){
  gamma = paraMRF[1]; beta_gene = paraMRF[2]; beta_cell = paraMRF[3]
  count = 0
  statessum = statesI*0
  cat("\n")
  for (iter in 1:iterGibbs){
    if (iter%%200 == 0){
      cat("\r", "Posterior sampling,", floor(iter/iterGibbs*100), "%", "completed")
    }
    for (g in 1:n_gene) {
      for (c in 1:n_cell) {
        a = gamma + beta_gene*sum(statesI[1,,c]*g_g[g,]) + beta_cell*sum(statesI[1,g,]*c_c[c,])
        b = beta_gene*sum(1-statesI[1,,c]*g_g[g,]) + beta_cell*sum(1-statesI[1,g,]*c_c[c,])
        prob = exp(a)*fy1[1,g,c]/(exp(a)*fy1[1,g,c]+ exp(b)*fy0[1,g,c])
        statesI[1,g,c] <- (prob >= runif(length(prob))) + 0
      }
    }
    if (iter >= 2 && (iter%%skip)==0)
    {
      count = count+1
      statessum = statessum + statesI
    }
  }
  pfdr1 <- statessum/count
  return(pfdr1)
}


