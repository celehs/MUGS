utils::globalVariables(c(
  "S.1", "S.2", "U.1", "U.2",
  "X.group.source", "X.group.target",
  "pairs.rel.CV", "pairs.rel.EV"
))
#' Function Used For Tuning And Evaluation
#'
#' @param pairs.rel the known code-code pairs
#' @param U the code embedding matrix
#' @param seed Optional integer for reproducibility of sampling.
#'
#' @return The output of tuning and evaluation
#' @export


evaluation.sim <- function(pairs.rel, U, seed = NULL) {
  required_packages <- c("dplyr", "pROC")
  missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
  if (length(missing_packages) > 0) {
    stop("Missing packages need to be installed: ", paste(missing_packages, collapse = ", "), call. = FALSE)
  }
  lapply(required_packages, library, character.only = TRUE)

  if (!is.null(seed)) set.seed(seed)

  names <- rownames(U)
  pairs <- data.frame(matrix(0, length(pairs.rel$row), 2))
  colnames(pairs) <- c('id1', 'id2')
  pairs$id1 <- match(pairs.rel$row, names)
  pairs$id2 <- match(pairs.rel$col, names)
  pairs <- na.omit(pairs)
  n.rel.neg <- nrow(pairs)

  pairs$id1 <- pmin(pairs$id1, pairs$id2)
  pairs$id2 <- pmax(pairs$id1, pairs$id2)

  names.rel.1 <- pairs.rel[, 1][(pairs.rel[, 1] %in% names) & (pairs.rel[, 2] %in% names)]
  names.rel.2 <- pairs.rel[, 2][(pairs.rel[, 1] %in% names) & (pairs.rel[, 2] %in% names)]
  names.rel <- union(names.rel.1, names.rel.2)

  pairs.rel.neg <- data.frame(
    id1 = sample(which(names %in% names.rel), 3 * n.rel.neg, replace = TRUE),
    id2 = sample(which(names %in% names.rel), 3 * n.rel.neg, replace = TRUE)
  )
  pairs.rel.neg <- subset(pairs.rel.neg, id1 != id2)
  pairs.rel.neg$id1 <- pmin(pairs.rel.neg$id1, pairs.rel.neg$id2)
  pairs.rel.neg$id2 <- pmax(pairs.rel.neg$id1, pairs.rel.neg$id2)
  pairs.rel.neg <- pairs.rel.neg[!duplicated(pairs.rel.neg), ]
  pairs.rel.neg <- dplyr::anti_join(pairs.rel.neg, pairs, by = c("id1", "id2"))

  if (nrow(pairs.rel.neg) < n.rel.neg) {
    message("Too few negative samples generated. Returning NA.")
    return(list(n.rel = 0, AUC.rel = NA))
  }

  pairs.rel.neg <- pairs.rel.neg[sample(nrow(pairs.rel.neg), n.rel.neg, replace = FALSE), ]
  y <- c(rep(1, n.rel.neg), rep(0, n.rel.neg))

  id1 <- c(pairs$id1, pairs.rel.neg$id1)
  id2 <- c(pairs$id2, pairs.rel.neg$id2)

  p <- mapply(function(i, j) sum(U[i, ] * U[j, ]), id1, id2)
  AUC.rel <- pROC::auc(y, p)

  list(n.rel = n.rel.neg, AUC.rel = AUC.rel)
}
