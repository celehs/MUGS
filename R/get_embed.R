utils::globalVariables(c(
  "S.1", "S.2", "U.1", "U.2",
  "X.group.source", "X.group.target",
  "pairs.rel.CV", "pairs.rel.EV"
))
#'Function For Getting Embedding From SVD
#'
#' @param mysvd the (managed) svd result (adding an element with 'names')
#' @param d dim of the final embedding
#' @param normalize if the output embeddings have l2 norm equal to 1
#'
#' @return The embedding from SVD
#' @export


get_embed = function(mysvd, d=2000, normalize=TRUE){
  id = which(sign(mysvd$u[1,])==sign(mysvd$v[1,]))
  id = id[1:min(d,length(id))]
  embed = mysvd$u[,id]%*%diag(sqrt(mysvd$d[id]))
  if(normalize){
    embed = embed/apply(embed,1,norm,'2')
  }
  rownames(embed) = mysvd$names
  return(embed)
}
