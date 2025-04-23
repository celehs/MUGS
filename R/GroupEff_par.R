utils::globalVariables(c(
  "S.1", "S.2", "U.1", "U.2",
  "X.group.source", "X.group.target",
  "pairs.rel.CV", "pairs.rel.EV"
))
#' Function Used To Estimate Group Effects Parallelly
#'
#' @param S.MGB SPPMI from the source site
#' @param S.BCH SPPMI from the target site
#' @param n.MGB the number of codes from the source site
#' @param n.BCH the number of codes from the target site
#' @param U.MGB the left embeddings (left singular vectors times the square root of the singular values) from the source site
#' @param U.BCH the left embeddings (left singular vectors times the square root of the singular values) from the target site
#' @param V.MGB the right embeddings (right singular vectors times the square root of the singular values) from the source site
#' @param V.BCH the right embeddings (right singular vectors times the square root of the singular values) from the target site
#' @param X.MGB.group the dummy matrix based on prior group structures at the source site
#' @param X.BCH.group the dummy matrix based on prior group structures at the target site
#' @param n.group the number of groups
#' @param name.list the full list of code names from the source site and the target site with repeated names of overlapping codes
#' @param beta.int the initial estimator for the group effects
#' @param lambda the tuning parameter controls the intensity of penalization on the group effect; by default we set it to 0
#' @param p the length of an embedding
#' @param n.core the number of cored used for parallel computation
#'
#' @return The output of estimating group effects parallelly
#' @export


GroupEff_par <- function(S.MGB, S.BCH, n.MGB, n.BCH, U.MGB, U.BCH, V.MGB, V.BCH, X.MGB.group, X.BCH.group, n.group, name.list, beta.int, lambda=0, p, n.core){
  required_packages <- c("parallel", "doSNOW", "Matrix")

  # Check for missing packages and stop if any are not installed
  missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
  if (length(missing_packages) > 0) {
    stop("Missing packages need to be installed: ", paste(missing_packages, collapse = ", "), call. = FALSE)
  }

  # Load required packages
  lapply(required_packages, library, character.only = TRUE)


  n.j.MGB <- apply(X.MGB.group>0, 2, sum)
  n.j.BCH <- apply(X.BCH.group>0, 2, sum)
  BCH.only <- names(which(n.j.MGB == 0))
  MGB.only <- names(which(n.j.BCH == 0))
  n.BCH.only <- length(BCH.only)
  n.MGB.only <- length(MGB.only)
  BETA.BCHonly <- matrix(0, p, n.BCH.only)
  colnames(BETA.BCHonly) <- BCH.only
  BETA.MGBonly <- matrix(0, p, n.MGB.only )
  colnames(BETA.MGBonly) <- MGB.only
  if (n.BCH.only >0){
    for (j in 1:n.BCH.only){
      nn <- n.j.BCH[names(n.j.BCH)==BCH.only[j]]
      temp.BCH <- matrix(0,nn,p)
      for (jj in 1: nn){
        temp.BCH[jj,] <- as.vector((V.BCH[X.BCH.group[, colnames(X.BCH.group)==BCH.only[j] ]>0,])[jj,] - beta.int[rownames(beta.int)==BCH.only[j],] )
      }
      name.BCH <- (name.list[(1+n.MGB):(n.MGB + n.BCH)])[X.BCH.group[,colnames(X.BCH.group)==BCH.only[j]]>0]
      Y.BCH.1 <-as.vector(S.BCH[,rownames(S.BCH)%in%name.BCH])
      Y.BCH <- Y.BCH.1 - as.vector(as.matrix(U.BCH)%*%t(temp.BCH))
      V.j.BCH <- do.call(rbind, replicate(nn, U.BCH, simplify=FALSE))
      m <- glmnet(V.j.BCH, Y.BCH, alpha = 0, family = 'gaussian', lambda = lambda*2/length(Y.BCH), intercept=F)
      BETA.BCHonly[,j] <- as.vector(m$beta)
    }
  }
  if (n.MGB.only >0){
    for (j in 1:n.MGB.only){
      nn <- n.j.MGB[names(n.j.MGB)==MGB.only[j]]
      temp.MGB <- matrix(0,nn,p)
      for (jj in 1: nn){
        temp.MGB[jj,] <- as.vector((V.MGB[X.MGB.group[, colnames(X.MGB.group)==MGB.only[j] ]>0,])[jj,] - beta.int[rownames(beta.int)==MGB.only[j],] )
      }
      name.MGB <- (name.list[1:n.MGB])[X.MGB.group[,colnames(X.MGB.group)==MGB.only[j]]>0]
      Y.MGB.1 <- as.vector(S.MGB[,rownames(S.MGB)%in%name.MGB])
      Y.MGB <- Y.MGB.1 - as.vector(as.matrix(U.MGB)%*%t(temp.MGB))
      V.j.MGB <- do.call(rbind, replicate(nn, U.MGB, simplify=FALSE))
      m <- glmnet(V.j.MGB, Y.MGB, alpha = 0, family = 'gaussian', lambda = lambda*2/length(Y.MGB), intercept=F)
      BETA.MGBonly[,j] <- as.vector(m$beta)
    }
  }
  cl <- makeCluster(n.core, type="SOCK")
  registerDoSNOW (cl)
  n.both <- n.group - n.BCH.only - n.MGB.only
  name.both <- setdiff(colnames(X.MGB.group),c(BCH.only, MGB.only))
  BETA = foreach(j = 1:n.both, .packages = c("glmnet"), .combine = 'cbind') %dopar%{
    name.MGB <- (name.list[1:n.MGB])[X.MGB.group[, colnames(X.MGB.group)==name.both[j]]>0]
    name.BCH <- (name.list[(1+n.MGB):(n.MGB + n.BCH)])[X.BCH.group[,colnames(X.BCH.group)==name.both[j]]>0]
    temp.MGB <- V.MGB[X.MGB.group[,colnames(X.MGB.group)==name.both[j]]>0,] - do.call(rbind, replicate(length(name.MGB),  beta.int[rownames(beta.int)==name.both[j],], simplify=FALSE))
    Y.MGB <- as.vector(S.MGB[,rownames(S.MGB)%in%name.MGB]) - as.vector(as.matrix(U.MGB)%*%t(temp.MGB))
    temp.BCH<- V.BCH[X.BCH.group[,colnames(X.BCH.group)==name.both[j]]>0,] - do.call(rbind, replicate(length(name.BCH),  beta.int[rownames(beta.int)==name.both[j],], simplify=FALSE))
    Y.BCH <- as.vector(S.BCH[,rownames(S.BCH)%in%name.BCH]) - as.vector(as.matrix(U.BCH)%*%t(temp.BCH))
    Y <- c(Y.MGB, Y.BCH)
    V.j <- rbind(do.call(rbind, replicate(length(name.MGB), U.MGB, simplify=FALSE)), do.call(rbind, replicate(length(name.BCH), U.BCH, simplify=FALSE)))
    m <- glmnet(V.j, Y, alpha = 0, family = 'gaussian', lambda = lambda*2/length(Y), intercept=F)
    return(m$beta)
  }
  stopCluster(cl)
  colnames(BETA) <- name.both
  BETA.1 <- t(cbind(BETA.BCHonly, BETA.MGBonly, as.matrix(BETA)))
  beta <- BETA.1[match(rownames(beta.int),rownames(BETA.1)),]
  dif_F <- norm(beta- beta.int, type = "F")^2/(p*n.group)
  X.group <- rbind(X.MGB.group,X.BCH.group)
  V.new <- rbind(V.MGB, V.BCH) - as.matrix(X.group)%*%(beta.int- beta)
  V.MGB.new  <- V.new[1:n.MGB,]
  V.BCH.new  <- V.new[(n.MGB+1):(n.MGB + n.BCH),]
  output <- list("beta" = beta, "dif_F" = dif_F,'V.MGB.new' = V.MGB.new, 'V.BCH.new' = V.BCH.new)
  return(output)
}
