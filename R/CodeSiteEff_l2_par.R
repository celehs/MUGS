utils::globalVariables(c(
  "S.1", "S.2", "U.1", "U.2",
  "X.group.source", "X.group.target",
  "pairs.rel.CV", "pairs.rel.EV","i","LinReg"
))
#' Function Used To Estimate Code-Site Effects Parallelly
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


CodeSiteEff_l2_par <- function(S.1, S.2, n1, n2, U.1, U.2, V.1, V.2, delta.int, lambda.delta, p, common_codes, n.common, n.core){

  required_packages <- c("parallel", "doSNOW")

  # Check for missing packages and stop if any are not installed
  missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
  if (length(missing_packages) > 0) {
    stop("Missing packages need to be installed: ", paste(missing_packages, collapse = ", "), call. = FALSE)
  }

  # Load required packages
  lapply(required_packages, library, character.only = TRUE)

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
