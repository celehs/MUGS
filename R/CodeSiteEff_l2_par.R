#' Function Used To Estimate Code-Site Effects Parallelly
#'
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doSNOW registerDoSNOW
#' @importFrom foreach %dopar% foreach
#' @importFrom grpreg grpreg
#' @importFrom glmnet glmnet
#' @importFrom stats na.omit
#' @importFrom dplyr anti_join
#' @importFrom mvtnorm rmvnorm
#'
#' @param S.1 SPPMI from the source site
#' @param S.2 SPPMI from the target site
#' @param n1 the number of codes from the source site
#' @param n2 the number of codes from the target site
#' @param U.1 the left embeddings (left singular vectors times the square root of the singular values) from the source site
#' @param U.2 the left embeddings (left singular vectors times the square root of the singular values) from the target site
#' @param V.1 the right embeddings (right singular vectors times the square root of the singular values) from the source site
#' @param V.2  the right embeddings (right singular vectors times the square root of the singular values) from the target site
#' @param delta.int the initial estimator for the code-site effect
#' @param lambda.delta the tuning parameter controls the intensity of penalization on the code-site effects
#' @param p the length of an embedding
#' @param common_codes the list of overlapping codes
#' @param n.common the number of overlapping codes
#' @param n.core the number of cored used for parallel computation
#'
#' @return The output for the estimation of code-site effects
#' @export
#'
#' @examples
#' S1_data <- system.file("extdata", "S.1.Rdata", package = "MUGS")
#' load(S1_data)
#' S2_data <- system.file("extdata", "S.2.Rdata", package = "MUGS")
#' load(S2_data)
#' U1_data <- system.file("extdata", "U.1.Rdata", package = "MUGS")
#' load(U1_data)
#' U2_data <- system.file("extdata", "U.2.Rdata", package = "MUGS")
#' load(U2_data)
#' X_group_source_data <- system.file("extdata", "X.group.source.Rdata", package = "MUGS")
#' load(X_group_source_data)
#' X_group_target_data <- system.file("extdata", "X.group.target.Rdata", package = "MUGS")
#' load(X_group_target_data)
# 'names.list.1 <- rownames(S.1)
# 'names.list.2 <- rownames(S.2)
# 'common_codes <- intersect(names.list.1, names.list.2)
# 'full.name.list <- c(names.list.1, names.list.2)
# 'delta.int <- matrix(0, 4000, 10)
# 'rownames(delta.int) <- full.name.list
# 'CodeSiteEff_l2_par.out <-  CodeSiteEff_l2_par(S.1 = S.1, S.2 = S.2, n1 = 2000, n2 =2000, U.1 = U.1, U.2 = U.2, V.1= U.1, V.2 = U.2,
# '                                              delta.int = delta.int, lambda.delta = 1000, p=10, common_codes = common_codes, n.common = 1000, n.core=2)


CodeSiteEff_l2_par <- function(S.1, S.2, n1, n2, U.1, U.2, V.1, V.2, delta.int, lambda.delta, p, common_codes, n.common, n.core){
  utils::globalVariables("i")
  delta.int.1 <- delta.int[1:n1,]
  Y.1 <- S.1 - U.1 %*% t(V.1 - delta.int.1)
  delta.1.common <- matrix(0, n.common, p)
  Y.1.common <- Y.1[,colnames(Y.1)%in%common_codes]
  cl <- makeCluster(n.core, type="SOCK")
  registerDoSNOW (cl)
  DELTA.1 = foreach(i = 1:n.common, .packages = c("grplasso"), .combine = 'cbind') %dopar%{
    fit <- grplasso(U.1, y =  Y.1.common[,i], index = rep(1, p), lambda = lambda.delta*sqrt(log(p)/n1)/sqrt(p), model = LinReg(),
                    penscale = sqrt, center = FALSE, standardize = FALSE)
    return(fit$coefficients)
  }
  stopCluster(cl)
  delta.1.common <- t(as.matrix(DELTA.1))
  rownames(delta.1.common) <- rownames(delta.int.1)[rownames(delta.int.1)%in%common_codes]
  delta.1.temp <- rbind(delta.1.common, delta.int.1[!(rownames(delta.int.1)%in%common_codes),])
  delta.1 <- delta.1.temp[match(rownames(delta.int.1), rownames(delta.1.temp)),]
  delta.int.2 <- delta.int[(1+n1):(n1 + n2),]
  Y.2 <- S.2 - U.2%*% t(V.2 -delta.int.2)
  delta.2.common <- matrix(0, n.common, p)
  Y.2.common <- Y.2[,colnames(Y.2)%in%common_codes]
  cl <- makeCluster(n.core, type="SOCK")
  registerDoSNOW (cl)
  DELTA.2= foreach(i = 1:n.common, .packages = c("grplasso"), .combine = 'cbind') %dopar%{
    fit <- grplasso(U.2, y =  Y.2.common[,i], index = rep(1, p), lambda = lambda.delta*sqrt(log(p)/n2)/sqrt(p), model = LinReg(),
                    penscale = sqrt, center = FALSE, standardize = FALSE)
    return(fit$coefficients)
  }
  stopCluster(cl)
  delta.2.common <- t(as.matrix(DELTA.2))
  rownames(delta.2.common) <- rownames(delta.int.2)[rownames(delta.int.2)%in%common_codes]
  delta.2.temp <- rbind(delta.2.common, delta.int.2[!(rownames(delta.int.2)%in%common_codes),])
  delta.2 <- delta.2.temp[match(rownames(delta.int.2), rownames(delta.2.temp)),]

  delta <- rbind(delta.1, delta.2)
  V.1.new <- V.1 - delta.int.1 + delta.1
  V.2.new <- V.2 - delta.int.2 + delta.2
  output <- list("delta" = delta, 'V.1.new' = V.1.new, 'V.2.new' = V.2.new)
  return(output)
}
