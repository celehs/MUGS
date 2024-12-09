#' Function Used For Tuning And Evaluation
#'
#' @importFrom stats na.omit
#' @importFrom dplyr anti_join
#' @importFrom pROC auc
#'
#' @param pairs.rel the known code-code pairs
#' @param U the code embedding matrix
#'
#' @return The output of tuning and evaluation
#' @export
#'

evaluation.sim <- function(pairs.rel, U){
  required_packages <- c("dplyr", "pROC")

  # Check for missing packages and stop if any are not installed
  missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
  if (length(missing_packages) > 0) {
    stop("Missing packages need to be installed: ", paste(missing_packages, collapse = ", "), call. = FALSE)
  }

  # Load required packages
  lapply(required_packages, library, character.only = TRUE)

  names <- rownames(U)
  pairs <- data.frame(matrix(0, length(pairs.rel$row), 2) )
  colnames(pairs) <- c('id1', 'id2')
  pairs$id1 <- match(pairs.rel$row,names)
  pairs$id2 <- match(pairs.rel$col, names)
  pairs <- na.omit(pairs)
  n.rel.neg <- nrow(pairs)
  pairs.temp <- pairs
  pairs$id1 <- apply(pairs.temp, 1, min)
  pairs$id2 <- apply(pairs.temp, 1, max)

  names.rel.1 <- (pairs.rel[,1])[(pairs.rel[,1]%in%names)&(pairs.rel[,2]%in%names)]
  names.rel.2 <-(pairs.rel[,2])[(pairs.rel[,1]%in%names)&(pairs.rel[,2]%in%names)]
  names.rel <- union(names.rel.1, names.rel.2)
  set.seed(1)
  pairs.rel.neg = data.frame(id1 = sample(which(names%in%names.rel), 3*n.rel.neg, replace = TRUE),
                             id2 = sample(which(names%in%names.rel), 3*n.rel.neg, replace = TRUE))
  pairs.rel.neg = subset(pairs.rel.neg, id1!=id2)
  pairs.rel.neg_temp = pairs.rel.neg
  pairs.rel.neg$id1 = apply(pairs.rel.neg_temp, 1, min)
  pairs.rel.neg$id2 = apply(pairs.rel.neg_temp, 1, max)
  pairs.rel.neg = pairs.rel.neg[!(duplicated(pairs.rel.neg)), ]
  pairs.rel.neg = anti_join(pairs.rel.neg, pairs, by = c("id1", "id2")) #19273     2

  if(nrow(pairs.rel.neg)<n.rel.neg){
    cat("Warning! Something Wrong Here! Too little Null sample!\n")
    return(0)
  }

  pairs.rel.pos <- pairs
  pairs.rel.neg <- pairs.rel.neg[sample(nrow(pairs.rel.neg),n.rel.neg,replace = FALSE), c('id1','id2')]
  y = c(rep(1, nrow(pairs.rel.pos)), rep(0, nrow(pairs.rel.neg)))
  id1 = pairs.rel.pos$id1; id2 = pairs.rel.pos$id2; num = length(id1)
  p1 = unlist(lapply(1:num, function(i){
    return(sum(U[id1[i],]*U[id2[i],]))
  }))
  id1 = pairs.rel.neg$id1; id2 = pairs.rel.neg$id2; num = length(id1)
  p2 = unlist(lapply(1:num, function(i){
    return(sum(U[id1[i],]*U[id2[i],]))
  }))
  p = c(p1,p2)
  AUC.rel = auc(y,p)
  output <- list("n.rel" = n.rel.neg, "AUC.rel" = AUC.rel)
  return(output)
}
