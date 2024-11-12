#' Function Used To Estimate Code Effects
#'
#' @param S.1 SPPMI from the source site
#' @param S.2 SPPMI from the target site
#' @param n1 the number of codes from the source site
#' @param n2 the number of codes from the target site
#' @param U.1 the left embeddings (left singular vectors times the square root of the singular values) from the source site
#' @param U.2 the left embeddings (left singular vectors times the square root of the singular values) from the target site
#' @param V.1 the right embeddings (right singular vectors times the square root of the singular values) from the source site
#' @param V.2 the right embeddings (right singular vectors times the square root of the singular values) from the target site
#' @param common_codes the list of overlapping codes
#' @param zeta.int the initial estimator for the code effects
#' @param lambda the tuning parameter controls the intensity of penalization on the code effect
#' @param p the length of an embedding
#'
#' @return The output for the estimation of code effects
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
#‘ names.list.1 <- rownames(S.1)
#‘ names.list.2 <- rownames(S.2)
#‘ common_codes <- intersect(names.list.1, names.list.2)
#‘ full.name.list <- c(names.list.1, names.list.2)

#‘ source('CodeEff_Matrix.R')
# Any 4000*10 matrix with row names being the code names
#‘ zeta.int <- matrix(0, 4000, 10)
#‘ rownames(zeta.int) <- full.name.list
#‘ CodeEff_Matrix.out <- CodeEff_Matrix(S.1=S.1, S.2=S.2, n1=2000, n2=2000, U.1=U.1 , U.2=U.2, V.1=U.1, V.2=U.2,
#‘                                      common_codes = common_codes, zeta.int = zeta.int, lambda=10, p=10)

CodeEff_Matrix<- function(S.1, S.2, n1, n2, U.1, U.2, V.1, V.2, common_codes, zeta.int, lambda, p){
  zeta.int.1 <- zeta.int[1:n1,]
  zeta.int.2 <- zeta.int[(1+n1):(n1+n2),]
  Y.1 <- S.1 - U.1%*%t(V.1 -  zeta.int.1)
  Y.2 <- S.2 - U.2%*%t(V.2 -  zeta.int.2)
  zeta.1.only <- t(solve(t(U.1)%*%U.1 + lambda*sqrt(log(p)/n1)*diag(p))%*%t(U.1)%*%Y.1[,!(colnames(Y.1)%in%common_codes)])
  zeta.2.only <- t(solve(t(U.2)%*%U.2 + lambda*sqrt(log(p)/n2)*diag(p))%*%t(U.2)%*%Y.2[,!(colnames(Y.2)%in%common_codes)])
  X <- rbind(U.1, U.2)
  Y <- rbind(Y.1[,colnames(Y.1)%in%common_codes], Y.2[,colnames(Y.2)%in%common_codes])
  zeta.common <- t(solve(t(X)%*%X + lambda*sqrt(log(p)/(n1+n2))*diag(p))%*%t(X)%*%Y)
  zeta.1.temp <- rbind(zeta.1.only, zeta.common)
  zeta.1 <- zeta.1.temp[match(rownames(U.1),rownames(zeta.1.temp)),]
  zeta.2.temp <- rbind(zeta.2.only, zeta.common)
  zeta.2 <- zeta.2.temp[match(rownames(U.2),rownames(zeta.2.temp)),]
  zeta <- rbind(zeta.1, zeta.2)

  dif_F <- norm(zeta- zeta.int, type = "F")^2/(p*(n1+n2))
  V.new <- rbind(V.1, V.2) - zeta.int + zeta
  V.1.new  <- V.new[1:n1,]
  V.2.new  <- V.new[(n1+1):(n1 + n2),]
  output <- list("zeta" = zeta, "dif_F" = dif_F,'V.1.new' = V.1.new, 'V.2.new' = V.2.new)
  return(output)
}
