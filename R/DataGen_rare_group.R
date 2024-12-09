#' Function used to generate input data (used only for Simulations)
#' Generate SPPMIs, dummy matrices based on prior group structures, and code-code pairs for tuning and evaluation
#'
#' @importFrom MASS mvrnorm
#' @importFrom mvtnorm rmvnorm
#' @importFrom utils combn
#' @importFrom stats rnorm
#' @importFrom fastDummies dummy_cols
#' @importFrom methods signature
#'
#' @param seed for reproducibility
#' @param p the length of an embedding
#' @param n1 the number of codes in site 1
#' @param n2 the number of codes in site 2
#' @param n.common common: the number of overlapping codes
#' @param n.group the number of groups
#' @param sigma.eps.1 the sd of error in site 1
#' @param sigma.eps.2 the sd of error in site 2
#' @param ratio.delta the proportion of codes in each site that have site-specific effects applied to them
#' @param network.k the number of distinct blocks within each site for which unique inter-code correlations are modeled
#' @param rho.beta AR parameter for the group effects covariance matrix
#' @param rho.U0 AR parameter for the code effects covariance matrix
#' @param rho.delta AR parameter for the code-site effects covariance matrix
#' @param sigma.rare the sd of error for rare codes (usually larger than sigma.eps.1 and sigma.eps.2)
#' @param n.rare The number of rare codes
#' @param group.size the size of each group
#'
#' @return Returns input data, SPPMIs, dummy matrices based on prior group structures and code-code pairs for tuning and evaluation
#' @export
#'


DataGen_rare_group <- function(seed, p, n1, n2, n.common, n.group, sigma.eps.1, sigma.eps.2, ratio.delta, network.k, rho.beta, rho.U0, rho.delta, sigma.rare, n.rare, group.size){
  required_packages <- c("Matrix", "MASS", "fastDummies", "rsvd", "Rcpp", "Rcpp", "RcppArmadillo", "inline")

  # Check for missing packages and stop if any are not installed
  missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
  if (length(missing_packages) > 0) {
    stop("Missing packages need to be installed: ", paste(missing_packages, collapse = ", "), call. = FALSE)
  }

  # Load required packages
  lapply(required_packages, library, character.only = TRUE)

  code <- '
  using namespace Rcpp;
  int n = as<int>(n_);
  arma::vec mu = as<arma::vec>(mu_);
  arma::mat sigma = as<arma::mat>(sigma_);
  int ncols = sigma.n_cols; // Corrected syntax
  arma::mat Y = arma::randn(n, ncols);
  return wrap(arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma));
  '
  rmvnorm.rcpp <- cxxfunction(signature(n_="integer", mu_="numeric", sigma_="matrix"), code, plugin="RcppArmadillo", verbose=TRUE)

  set.seed(seed)
  N = n1 + n2 - n.common
  n1.no = n1 - n.common
  n2.no = n2 - n.common
  #### Group effect
  ar1_cor <- function(n, rho) {
    exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) -
                      (1:n - 1))
    rho^exponent
  }
  M.beta.1 <- ar1_cor(3,rho.beta)
  M.beta.2 <- Matrix::bdiag(M.beta.1, diag(7))
  Sigma.beta <- Matrix::bdiag(replicate(40, M.beta.2 , simplify=FALSE))
  beta <- mvrnorm(p, rep(0,n.group), Sigma.beta)
  beta <- t(beta)
  #### Code effect
  M.U0.1 <- ar1_cor(6,rho.U0)
  Sigma.U0.1 <- Matrix::bdiag(replicate(N/6, M.U0.1 , simplify=FALSE))
  M.U0.2 <- matrix(0.05, 6, 6)
  diag(M.U0.2) <- 0
  Sigma.U0.2.temp <-  Matrix::bdiag(replicate(300/6, M.U0.2 , simplify=FALSE))
  Sigma.U0.2 <- matrix(0, N, N)
  Sigma.U0.2[501:800, 501:800] <- as.matrix(Sigma.U0.2.temp)
  Sigma.U0 <- Sigma.U0.1 + Sigma.U0.2
  u0 <- rmvnorm.rcpp(p, rep(0,N), as.matrix(Sigma.U0))
  u0 <- t(u0)
  temp <- seq(501, 900, 2)
  group.only.codes <- rep(0, length(temp)*group.size)
  for (i in 1: length(temp)){
    group.only.codes[(group.size*(i-1)+1):(group.size*i) ] <- seq(temp[i], 2500, n.group)
  }
  u0[group.only.codes,] =0
  #### code-site effect
  Sigma.delta.0 <- ar1_cor(network.k, rho.delta)
  Sigma.delta.1 <- Matrix::bdiag(replicate(n1*ratio.delta/network.k, Sigma.delta.0, simplify=FALSE))
  Sigma.delta.2 <- Matrix::bdiag(replicate(n2*ratio.delta/network.k, Sigma.delta.0, simplify=FALSE))
  delta1.temp <- mvrnorm(p, rep(0,n1*ratio.delta), Sigma.delta.1)
  delta2.temp <- mvrnorm(p, rep(0,n2*ratio.delta), Sigma.delta.2)
  set.seed(1)
  idx1 <- sample(n.common, ncol(delta1.temp)) + 1000
  idx2 <- sample(n.common, ncol(delta2.temp))
  delta1 <- matrix(0, n1, p)
  delta2 <- matrix(0, n2, p)
  delta1[as.vector(idx1),] <-  t(delta1.temp)
  delta2[as.vector(idx2),] <-  t(delta2.temp)
  #### Code-Code pairs for tuning and evaluation
  ## similar pairs
  name.beta.full.1 <- matrix(501:2500, 400, 5)
  name.beta.full <-c(t(name.beta.full.1 ))
  supp.cov.beta <- Matrix::bdiag(replicate(40, Matrix::bdiag(matrix(1, 3*5, 3*5), Matrix::bdiag(replicate(7, matrix(1, 1*5, 1*5), simplify=FALSE))), simplify=FALSE))
  supp.cov.beta.full <- matrix(0, N, N)
  supp.cov.beta.full[name.beta.full, name.beta.full] <- as.matrix(supp.cov.beta)
  supp.cov.beta.full.upper <- supp.cov.beta.full
  supp.cov.beta.full.upper[lower.tri(supp.cov.beta.full.upper, diag = TRUE)] <- 0
  pairs.sim <- which(supp.cov.beta.full.upper!=0, arr.ind = T)
  pairs.sim <- as.data.frame(pairs.sim)
  pairs.sim$type <- rep('similarity', length(pairs.sim$row))
  pairs.sim$row <- as.character(pairs.sim$row)
  pairs.sim$col <- as.character(pairs.sim$col)
  ## related pairs
  supp.cov.u0 <- as.matrix(Sigma.U0)
  supp.cov.u0[as.matrix(Sigma.U0)!=0] <- 1
  supp.cov.u0.upper <- supp.cov.u0
  supp.cov.u0.upper[lower.tri(supp.cov.u0.upper, diag = TRUE)] <- 0
  pairs.rel.new <- matrix(0,n2*ratio.delta/network.k*(network.k*(network.k-1)/2),2)
  for (j in 1: (n2*ratio.delta/network.k)){
    pairs.rel.new[((network.k*(network.k-1)/2)*(j-1) + 1) :(j*network.k*(network.k-1)/2),] <- t(combn(idx2[((j-1)*network.k+1):(j*network.k)],2))
  }
  pairs.rel.new <- as.data.frame(pairs.rel.new+ (n1 - (n1-n.common)))
  supp.cov.delta.2 <- matrix(0, n2, n2)
  supp.cov.delta.2[cbind(as.vector((pairs.rel.new-(n1 - (n1-n.common)))[,1]), as.vector((pairs.rel.new-(n1 - (n1-n.common)))[,2]))] <- 1
  supp.cov.delta.2 <- supp.cov.delta.2 + t(supp.cov.delta.2 )
  supp.cov.delta.2.upper <- supp.cov.delta.2
  supp.cov.delta.2.upper[lower.tri(supp.cov.delta.2.upper , diag = TRUE)] <- 0
  n.rel <- sum(supp.cov.u0.upper!=0) + n2*ratio.delta/network.k*(network.k*(network.k-1)/2)
  pairs.rel.shared <- which(supp.cov.u0.upper!=0, arr.ind = T)
  pairs.rel.2 <- which(supp.cov.delta.2.upper!=0, arr.ind = T) + 1000
  ## Merge similar and related pairs
  pairs.rel <- rbind(pairs.rel.shared, pairs.rel.2)
  pairs.rel <- as.data.frame(pairs.rel)
  pairs.rel$type <- rep('related', dim(pairs.rel)[1])
  pairs.rel$row <- as.character(pairs.rel$row)
  pairs.rel$col <- as.character(pairs.rel$col)
  pairs.rel <- pairs.rel[!(duplicated(pairs.rel)),]
  pairs.rel.full <- rbind(pairs.sim, pairs.rel)
  pairs.rel.full$type <- 'related'
  set.seed(1)
  idx.rel <- sample(1:dim( pairs.rel.full)[1])
  pairs.rel.CV <-  pairs.rel.full[idx.rel[1:floor(dim( pairs.rel.full)[1]/2) ], ]
  pairs.rel.EV <-  pairs.rel.full[idx.rel[(floor(dim(pairs.rel)[1]/2)+ 1):dim( pairs.rel.full)[1]], ]
  idx.rare <- as.numeric(names(sort(table(c(pairs.rel.EV[,1], pairs.rel.EV[,2])), decreasing = T)))
  idx.rare <- (intersect(idx.rare, 1001:3000))[(1+100):(n.rare+100)]
  #### Embeddings
  beta.full.1 = do.call(rbind, replicate(5, beta, simplify=FALSE))
  dim(beta.full.1)
  beta.full <- rbind(rbind(matrix(0, 500,p), beta.full.1), matrix(0,500,p))
  # site 1
  u.1 = beta.full[1:n1,] + u0[1:n1,] + delta1
  ## site 2
  u.2 = beta.full[(n1 - n1.no +1):N,] + u0[(n1 - n1.no +1):N,] +  delta2
  ## True S matrices
  S.1.0 <- u.1%*%t(u.1)
  S.2.0 <- u.2%*%t(u.2)
  S.miss.0 <- u.1[1:n1.no, ]%*%t(u.2[(n.common+1):n2, ])
  ## Add noise
  set.seed(seed)
  err.1 <- matrix(rnorm(n1^2, 0, sigma.eps.1), n1, n1)
  S.1 <- S.1.0 + err.1
  #### Add noises for rare codes in target
  idx.freq <- setdiff(1001:3000, idx.rare)
  n.rare <- length(idx.rare)
  n.freq <- length(idx.freq)
  err.rare.1 <- matrix(rnorm(n.rare^2, 0, sigma.rare), n.rare, n.rare)
  err.rare.2 <- matrix(rnorm(n.rare*n.freq, 0, sqrt(sigma.rare*sigma.eps.2)), n.rare, n.freq)
  err.freq <- matrix(rnorm(n.freq^2, 0, sigma.eps.2), n.freq, n.freq)
  err.2.temp <- matrix(0, n2, n2)
  err.2.temp[1:n.rare, 1:n.rare] <- err.rare.1
  err.2.temp[1:n.rare, (1+n.rare):n2] <- err.rare.2
  err.2.temp[(1+n.rare):n2, 1:n.rare] <- t(err.rare.2)
  err.2.temp[(1+n.rare):n2, (1+n.rare):n2] <- err.freq
  err.2 <- err.2.temp[match(1001:3000, c(idx.rare, idx.freq)), match(1001:3000, c(idx.rare, idx.freq))]
  S.2 <- S.2.0 + err.2
  S.1[lower.tri(S.1)] = t(S.1)[lower.tri(S.1)]
  S.2[lower.tri(S.2)] = t(S.2)[lower.tri(S.2)]
  rownames(S.1) <- as.character(1:n1)
  colnames(S.1) <- as.character(1:n1)
  rownames(S.2) <- as.character((n1.no +1):N)
  colnames(S.2) <- as.character((n1.no +1):N)
  ######## Dummy matrices based on prior group structures ########
  group.names <- as.character(1:n.group)
  group.names.full <- rep(group.names, group.size)
  X.group.0 <- dummy_cols(group.names.full)
  X.group.1 <- X.group.0[,-1]
  X.group <- rbind(rbind(matrix(0,500,n.group), as.matrix(X.group.1)),matrix(0,500,n.group))
  X.group.source <- X.group[1:n1, ]
  X.group.target <- X.group[1001:N, ]
  rownames(X.group.source) <-  as.character(1:n1)
  rownames(X.group.target) <- as.character(1001:N)
  output <- list('delta1' = delta1, 'delta2' = delta2, 'u.1' = u.1, 'u.2' = u.2, 'S.1' = S.1, 'S.2' = S.2,
                 'S.1.0' = S.1.0, 'S.2.0' = S.2.0, 'X.group.source' = X.group.source, 'X.group.target' = X.group.target,
                 'pairs.rel.CV' = pairs.rel.CV, 'pairs.rel.EV' = pairs.rel.EV)
  return(output)
}
