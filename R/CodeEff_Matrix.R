utils::globalVariables(c(
  "S.1", "S.2", "U.1", "U.2",
  "X.group.source", "X.group.target",
  "pairs.rel.CV", "pairs.rel.EV"
))
#' Function Used To Estimate Code Effects
#'
#' This function estimates code effects using left and right embeddings from source and target sites.
#'
#' @param S.1 SPPMI from the source site.
#' @param S.2 SPPMI from the target site.
#' @param n1 The number of codes from the source site.
#' @param n2 The number of codes from the target site.
#' @param U.1 The left embeddings left singular vectors times the square root of the singular values from the source site.
#' @param U.2 The left embeddings left singular vectors times the square root of the singular values from the target site.
#' @param V.1 The right embeddings right singular vectors times the square root of the singular values from the source site.
#' @param V.2 The right embeddings right singular vectors times the square root of the singular values from the target site.
#' @param common_codes The list of overlapping codes.
#' @param zeta.int The initial estimator for the code effects.
#' @param lambda The tuning parameter controls the intensity of penalization on the code effect.
#' @param p The length of an embedding.
#'
#' @return A list with the following elements:
#' \item{zeta}{The estimated code effects.}
#' \item{dif_F}{The Frobenius norm difference between the updated and initial estimators.}
#' \item{V.1.new}{Updated right embeddings for the source site.}
#' \item{V.2.new}{Updated right embeddings for the target site.}
#'
#' @export



CodeEff_Matrix <- function(S.1, S.2, n1, n2, U.1, U.2, V.1, V.2, common_codes, zeta.int, lambda, p) {
  zeta.int.1 <- zeta.int[1:n1, ]
  zeta.int.2 <- zeta.int[(1 + n1):(n1 + n2), ]

  Y.1 <- S.1 - U.1 %*% t(V.1 - zeta.int.1)
  Y.2 <- S.2 - U.2 %*% t(V.2 - zeta.int.2)

  zeta.1.only <- t(solve(t(U.1) %*% U.1 + lambda * sqrt(log(p) / n1) * diag(p)) %*% t(U.1) %*% Y.1[, !(colnames(Y.1) %in% common_codes)])
  zeta.2.only <- t(solve(t(U.2) %*% U.2 + lambda * sqrt(log(p) / n2) * diag(p)) %*% t(U.2) %*% Y.2[, !(colnames(Y.2) %in% common_codes)])

  X <- rbind(U.1, U.2)
  Y <- rbind(Y.1[, colnames(Y.1) %in% common_codes], Y.2[, colnames(Y.2) %in% common_codes])
  zeta.common <- t(solve(t(X) %*% X + lambda * sqrt(log(p) / (n1 + n2)) * diag(p)) %*% t(X) %*% Y)

  zeta.1.temp <- rbind(zeta.1.only, zeta.common)
  rownames(zeta.1.temp) <- rownames(U.1)
  zeta.1 <- zeta.1.temp[rownames(U.1), , drop = FALSE]

  zeta.2.temp <- rbind(zeta.2.only, zeta.common)
  rownames(zeta.2.temp) <- rownames(U.2)
  zeta.2 <- zeta.2.temp[rownames(U.2), , drop = FALSE]

  zeta <- rbind(zeta.1, zeta.2)
  dif_F <- norm(zeta - zeta.int, type = "F")^2 / (p * (n1 + n2))
  V.new <- rbind(V.1, V.2) - zeta.int + zeta
  V.1.new <- V.new[1:n1, ]
  V.2.new <- V.new[(n1 + 1):(n1 + n2), ]

  list(zeta = zeta, dif_F = dif_F, V.1.new = V.1.new, V.2.new = V.2.new)
}
