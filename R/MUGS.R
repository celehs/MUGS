utils::globalVariables(c(
  "S.1", "S.2", "U.1", "U.2",
  "X.group.source", "X.group.target",
  "pairs.rel.CV", "pairs.rel.EV"
))

#' Main function for MUGS algorithm
#'
#' @param TUNE Logical value indicating whether the function should tune parameters TRUE or use predefined parameters FALSE.
#' @param Eva Logical value indicating whether to perform evaluation (TRUE) or skip it (FALSE).
#' @param Lambda The candidate values for the tuning parameter controlling the intensity of penalization on the code effects.
#' @param Lambda.delta The candidate values for the tuning parameter controlling the intensity of penalization on the code-site effects.
#' @param n.core Integer specifying the number of cores to use for parallel processing.
#' @param tol Numeric value representing the tolerance level for convergence in the algorithm.
#' @param seed Integer used to set the seed for random number generation, ensuring reproducibility. Set to NULL to disable.
#' @param S.1 The SPPMI matrix from site 1.
#' @param S.2 The SPPMI matrix from site 2.
#' @param X.group.source The dummy matrix representing the group structure of codes at site 1.
#' @param X.group.target The dummy matrix representing the group structure of codes at site 2.
#' @param pairs.rel.CV Code-code pairs used for tuning via cross-validation.
#' @param pairs.rel.EV Code-code pairs used for evaluation.
#' @param p Integer indicating the length of embeddings.
#' @param n.group The number of groups.
#' @param outdir Optional directory to write output files. Defaults to a temporary directory.
#'
#' @return A list or saved files containing the embedding matrices, similarity matrices, and site-heterogeneous code analysis.
#' @export
MUGS <- function(TUNE = FALSE, Eva = TRUE,
                 Lambda = c(10), Lambda.delta = c(1000), n.core = 4, tol = 1, seed = NULL,
                 S.1 = NULL, S.2 = NULL, X.group.source = NULL, X.group.target = NULL,
                 pairs.rel.CV = NULL, pairs.rel.EV = NULL, p = 100, n.group = 400, outdir = NULL) {

  required_packages <- c("rsvd", "dplyr", "pROC")
  missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
  if (length(missing_packages) > 0) {
    stop("Missing packages need to be installed: ", paste(missing_packages, collapse = ", "), call. = FALSE)
  }
  lapply(required_packages, library, character.only = TRUE)

  if (!is.null(seed)) set.seed(seed)

  ################  Data Preperation ################
  # Ensure data is provided for non-simulated cases
  stopifnot(!is.null(S.1), !is.null(S.2), !is.null(pairs.rel.CV), !is.null(X.group.source), !is.null(X.group.target))

  # Names for embeddings, common code handling etc.
  names.list.1 <- rownames(S.1)
  names.list.2 <- rownames(S.2)
  n1 <- length(names.list.1)
  n2 <- length(names.list.2)
  common_codes <- intersect(names.list.1, names.list.2)
  n.common <- length(common_codes)
  K1 <-length(Lambda)
  K2 <- length(Lambda.delta)
  N = n1 + n2 - n.common
  n1.no = n1 - n.common
  n2.no = n2 - n.common

  ################  Initialization ################
  # Single-site svd to obtain initial embeddings at each site
  svd.1 <- rsvd(S.1, 3*p)
  svd.2 <- rsvd(S.2, 3*p)
  emb.1 <- get_embed(svd.1,d=p)
  emb.2 <- get_embed(svd.2,d=p)
  rownames(emb.1) <- names.list.1
  rownames(emb.2) <- names.list.2
  # Align two sets of initial embeddings
  comm.embed.1 <- emb.1[match( common_codes, rownames(emb.1)) ,]
  comm.embed.2 <- emb.2[match( common_codes, rownames(emb.2)) ,]
  # orthogonal Procrustes problem based on overlapped parts
  A <- t(comm.embed.1)%*%comm.embed.2
  eigenA<-eigen(t(A)%*%A + 10^-5*diag(p), symmetric = TRUE)
  W <- A%*%eigenA$vectors%*%diag(eigenA$values^(-1/2), p, p)%*%t(eigenA$vectors)
  embed.1.init <- emb.1%*%W
  embed.2.init <- emb.2
  #  Calculate initial group effects
  beta.names.1 <- unique(c(colnames(X.group.source), colnames(X.group.target)))
  beta.int <- matrix(0, n.group, p)
  rownames(beta.int) <- beta.names.1
  X.group <- rbind(X.group.source, X.group.target)
  embed.init <- rbind(embed.1.init,embed.2.init)
  name.list <- rownames(embed.init)
  for (i in 1:n.group){
    beta.int[i, ] <- colMeans(embed.init[X.group[,i]==TRUE ,])
  }
  # Calculate initial code effect
  name.list.unique <- unique(c(names.list.1, names.list.2))
  zeta.int.temp  <- matrix(0, N, p)
  rownames(zeta.int.temp ) <- name.list.unique
  zeta.int.temp[match(common_codes, rownames(zeta.int.temp )) , ] <- (embed.1.init[match(common_codes, rownames(embed.1.init)),] + embed.2.init[match(common_codes, rownames(embed.2.init)),] )/2
  zeta.int.temp[match(setdiff(names.list.1, common_codes),  rownames(zeta.int.temp ))  ,] <- embed.1.init[match(setdiff(names.list.1, common_codes), rownames(embed.1.init)),]
  zeta.int.temp[match(setdiff(names.list.2, common_codes),  rownames(zeta.int.temp ))  ,] <- embed.2.init[match(setdiff(names.list.2, common_codes), rownames(embed.2.init)),]
  unique_rows <- !duplicated(rownames(X.group))
  X.group.unique <- X.group[unique_rows, ]
  beta.int.full <- as.matrix(X.group.unique)%*%beta.int
  zeta.int.1 <- zeta.int.temp[match(names.list.1, rownames(zeta.int.temp)), ] -beta.int.full[match(names.list.1, rownames(beta.int.full)), ]
  zeta.int.2 <- zeta.int.temp[match(names.list.2, rownames(zeta.int.temp)), ] -beta.int.full[match(names.list.2, rownames(beta.int.full)), ]
  zeta.int <- rbind(zeta.int.1, zeta.int.2)
  # Calculate initial code-site effect
  delta.int <- embed.init - as.matrix(X.group)%*%beta.int - zeta.int

  ################  Main Loop for MUGS Algorithm ################
  # Create arrays to store results
  delta.spst.ovl <- matrix(0,K1, K2)
  CV.res <- matrix(0, K1, K2)

  for (k1 in 1:K1) {
    message("k1 = ", k1)
    lambda = Lambda[k1]
    for (k2 in 1:K2) {
      message("k2 = ", k2)
      lambda.delta = Lambda.delta[k2]
      # initialization
      U.1 <- embed.1.init
      U.2 <- embed.2.init
      V.1 <- U.1
      V.2 <- U.2
      delta.int.V <- delta.int
      zeta.int.V <- zeta.int
      beta.int.V <-beta.int
      delta.int.U <- delta.int
      zeta.int.U <-zeta.int
      beta.int.U <-beta.int
      dif_loss <- 10
      while (abs(dif_loss) > tol){
        loss_main <- (norm(S.1 - U.1%*%t(V.1), type = 'F')^2 + norm(S.2 - U.2%*%t(V.2), type = 'F')^2)/(n1^2 + n2^2)

        loss_zeta_V_ovl <- lambda*sqrt(log(p)/(n1+n2))*norm((zeta.int.V[1:n1, ])[rownames(zeta.int.V[1:n1, ])%in%common_codes, ],type = 'F')^2/(n1^2 + n2^2)
        loss_zeta_V_1 <- lambda*sqrt(log(p)/n1)*norm((zeta.int.V[1:n1, ])[!(rownames(zeta.int.V[1:n1, ])%in%common_codes), ],type = 'F')^2/(n1^2 + n2^2)
        loss_zeta_V_2 <- lambda*sqrt(log(p)/n2)*norm((zeta.int.V[(1+n1):(n1+n2), ])[!(rownames(zeta.int.V[(1+n1):(n1+n2), ])%in%common_codes), ],type = 'F')^2/(n1^2 + n2^2)
        loss_zeta_V <- loss_zeta_V_ovl + loss_zeta_V_1 + loss_zeta_V_2
        loss_zeta_U_ovl <- lambda*sqrt(log(p)/(n1+n2))*norm((zeta.int.U[1:n1, ])[rownames(zeta.int.U[1:n1, ])%in%common_codes, ],type = 'F')^2/(n1^2 + n2^2)
        loss_zeta_U_1 <- lambda*sqrt(log(p)/n1)*norm((zeta.int.U[1:n1, ])[!(rownames(zeta.int.U[1:n1, ])%in%common_codes), ],type = 'F')^2/(n1^2 + n2^2)
        loss_zeta_U_2 <- lambda*sqrt(log(p)/n2)*norm((zeta.int.U[(1+n1):(n1+n2), ])[!(rownames(zeta.int.U[(1+n1):(n1+n2), ])%in%common_codes), ],type = 'F')^2/(n1^2 + n2^2)
        loss_zeta_U <- loss_zeta_U_ovl + loss_zeta_U_1 + loss_zeta_U_2

        loss_delta_V.1 <- lambda.delta*sqrt(log(p)/n1)*sum(apply(delta.int.V[1:n1,], 1, norm, "2"))/(n1^2 + n2^2)
        loss_delta_V.2 <- lambda.delta*sqrt(log(p)/n2)*sum(apply(delta.int.V[(1+n1):(n1+n2),], 1, norm, "2"))/(n1^2 + n2^2)
        loss_delta_V <- loss_delta_V.1 + loss_delta_V.2
        loss_delta_U.1 <- lambda.delta*sqrt(log(p)/n1)*sum(apply(delta.int.U[1:n1,], 1, norm, "2"))/(n1^2 + n2^2)
        loss_delta_U.2 <- lambda.delta*sqrt(log(p)/n2)*sum(apply(delta.int.U[(1+n1):(n1+n2),], 1, norm, "2"))/(n1^2 + n2^2)
        loss_delta_U <- loss_delta_U.1 + loss_delta_U.2

        loss_p_V <- loss_zeta_V +  loss_delta_V
        loss_p_U <- loss_zeta_U +  loss_delta_U
        loss <- loss_main +  loss_p_V + loss_p_U
        message('\n loss:',loss)
        loss.V  <- loss
        dif_loss_V <-10
        # Update V
        while (abs(dif_loss_V) > tol)
        {
          # Update group effects
          GroupEff.out <- GroupEff_par(S.1, S.2, n1, n2, U.1, U.2, V.1, V.2, X.group.source, X.group.target, n.group, name.list, beta.int.V, 0, p, n.core)
          beta.int.V <- as.matrix(GroupEff.out$beta)
          V.1 <- as.matrix(GroupEff.out$V.MGB.new)
          V.2 <- as.matrix(GroupEff.out$V.BCH.new)
          temp.group <- (norm(S.1 - U.1%*%t(V.1), type = 'F')^2 + norm(S.2 - U.2%*%t(V.2), type = 'F')^2)/(n1^2 + n2^2)
          message('\n diff_groupEff=', temp.group - loss_main)
          # Update code effects
          CodeEff.out <-CodeEff_Matrix(S.1, S.2, n1, n2, U.1, U.2, V.1, V.2, common_codes, zeta.int.V, lambda, p)
          zeta.int.V <- as.matrix(CodeEff.out$zeta)
          V.1 <- as.matrix(CodeEff.out$V.1.new)
          V.2 <-as.matrix(CodeEff.out$V.2.new)
          temp.code.main <- (norm(S.1 - U.1%*%t(V.1), type = 'F')^2 + norm(S.2 - U.2%*%t(V.2), type = 'F')^2)/(n1^2 + n2^2)
          temp.code.p.ovl <- lambda*sqrt(log(p)/(n1+n2))*norm((zeta.int.V[1:n1, ])[rownames(zeta.int.V[1:n1, ])%in%common_codes, ],type = 'F')^2/(n1^2 + n2^2)
          temp.code.p.1 <-lambda*sqrt(log(p)/n1)*norm((zeta.int.V[1:n1, ])[!(rownames(zeta.int.V[1:n1, ])%in%common_codes), ],type = 'F')^2/(n1^2 + n2^2)
          temp.code.p.2 <-lambda*sqrt(log(p)/n2)*norm((zeta.int.V[(1+n1):(n1+n2), ])[!(rownames(zeta.int.V[(1+n1):(n1+n2), ])%in%common_codes), ],type = 'F')^2/(n1^2 + n2^2)
          message('\n diff_CodeEff=', temp.code.main-temp.group+temp.code.p.ovl+temp.code.p.1+temp.code.p.2-loss_zeta_V)
          # Update code-site effects
          CodeSiteEff.out <- CodeSiteEff_l2_par(S.1, S.2, n1, n2, U.1, U.2, V.1, V.2, delta.int.V, lambda.delta, p, common_codes, n.common, n.core)
          delta.int.V <- as.matrix(CodeSiteEff.out$delta)
          V.1 <- as.matrix(CodeSiteEff.out$V.1.new)
          V.2 <- as.matrix(CodeSiteEff.out$V.2.new)
          temp.codesite.main <- (norm(S.1 - U.1%*%t(V.1), type = 'F')^2 + norm(S.2 - U.2%*%t(V.2), type = 'F')^2)/(n1^2 + n2^2)
          temp.codesite.p <- lambda.delta*sqrt(log(p)/n1)*sum(apply(delta.int.V[1:n1,], 1, norm, "2"))/(n1^2 + n2^2) +lambda.delta*sqrt(log(p)/n2)*sum(apply(delta.int.V[(1+n1):(n1+n2),], 1, norm, "2"))/(n1^2 + n2^2)
          message('\n diff_CodeSiteEff=', temp.codesite.main - temp.code.main + temp.codesite.p -loss_delta_V)

          loss_main <- (norm(S.1 - U.1%*%t(V.1), type = 'F')^2 + norm(S.2 - U.2%*%t(V.2), type = 'F')^2)/(n1^2 + n2^2)
          loss_zeta_V_ovl <- lambda*sqrt(log(p)/(n1+n2))*norm((zeta.int.V[1:n1, ])[rownames(zeta.int.V[1:n1, ])%in%common_codes, ],type = 'F')^2/(n1^2 + n2^2)
          loss_zeta_V_1 <- lambda*sqrt(log(p)/n1)*norm((zeta.int.V[1:n1, ])[!(rownames(zeta.int.V[1:n1, ])%in%common_codes), ],type = 'F')^2/(n1^2 + n2^2)
          loss_zeta_V_2 <- lambda*sqrt(log(p)/n2)*norm((zeta.int.V[(1+n1):(n1+n2), ])[!(rownames(zeta.int.V[(1+n1):(n1+n2), ])%in%common_codes), ],type = 'F')^2/(n1^2 + n2^2)
          loss_zeta_V <- loss_zeta_V_ovl + loss_zeta_V_1 + loss_zeta_V_2

          loss_delta_V.1 <- lambda.delta*sqrt(log(p)/n1)*sum(apply(delta.int.V[1:n1,], 1, norm, "2"))/(n1^2 + n2^2)
          loss_delta_V.2 <- lambda.delta*sqrt(log(p)/n2)*sum(apply(delta.int.V[(1+n1):(n1+n2),], 1, norm, "2"))/(n1^2 + n2^2)

          loss_delta_V <- loss_delta_V.1 + loss_delta_V.2
          loss_p_V <- loss_zeta_V +  loss_delta_V
          loss.new.V <- loss_main + loss_p_V + loss_p_U
          dif_loss_V <- (loss.new.V - loss.V)/loss.V
          message('\n dif_loss_V',dif_loss_V)
          loss.V <- loss.new.V
        }
        # Update U
        loss.U <- loss.new.V
        dif_loss_U <- 10
        while (abs(dif_loss_U) > 1)
        {
          # Update group effect
          GroupEff.out <- GroupEff_par(t(S.1), t(S.2), n1, n2, V.1, V.2, U.1, U.2, X.group.source, X.group.target, n.group, name.list, beta.int.U, 0, p, n.core)
          beta.int.U <- as.matrix(GroupEff.out$beta)
          U.1 <- as.matrix(GroupEff.out$V.MGB.new)
          U.2 <- as.matrix(GroupEff.out$V.BCH.new)
          temp.group <- (norm(S.1 - U.1%*%t(V.1), type = 'F')^2 + norm(S.2 - U.2%*%t(V.2), type = 'F')^2)/(n1^2 + n2^2)
          message('\n diff_groupEff=', temp.group - loss_main)

          # Update code effect
          CodeEff.out <-CodeEff_Matrix(t(S.1), t(S.2), n1, n2, V.1, V.2, U.1, U.2, common_codes, zeta.int.U, lambda, p)
          zeta.int.U <- as.matrix(CodeEff.out$zeta)
          U.1 <- as.matrix(CodeEff.out$V.1.new)
          U.2 <-as.matrix(CodeEff.out$V.2.new)
          temp.code.main <- (norm(S.1 - U.1%*%t(V.1), type = 'F')^2 + norm(S.2 - U.2%*%t(V.2), type = 'F')^2)/(n1^2 + n2^2)
          temp.code.p.ovl <- lambda*sqrt(log(p)/(n1+n2))*norm((zeta.int.U[1:n1, ])[rownames(zeta.int.U[1:n1, ])%in%common_codes, ],type = 'F')^2/(n1^2 + n2^2)
          temp.code.p.1 <-lambda*sqrt(log(p)/n1)*norm((zeta.int.U[1:n1, ])[!(rownames(zeta.int.U[1:n1, ])%in%common_codes), ],type = 'F')^2/(n1^2 + n2^2)
          temp.code.p.2 <-lambda*sqrt(log(p)/n2)*norm((zeta.int.U[(1+n1):(n1+n2), ])[!(rownames(zeta.int.U[(1+n1):(n1+n2), ])%in%common_codes), ],type = 'F')^2/(n1^2 + n2^2)
          message('\n diff_codeEff=', temp.code.main-temp.group+temp.code.p.ovl+temp.code.p.1+temp.code.p.2-loss_zeta_U)

          # Update code-site effect
          CodeSiteEff.out <- CodeSiteEff_l2_par(t(S.1), t(S.2), n1, n2, V.1, V.2, U.1, U.2, delta.int.U, lambda.delta, p, common_codes, n.common, n.core)
          delta.int.U <- as.matrix(CodeSiteEff.out$delta)
          U.1 <- as.matrix(CodeSiteEff.out$V.1.new)
          U.2 <- as.matrix(CodeSiteEff.out$V.2.new)
          temp.codesite.main <- (norm(S.1 - U.1%*%t(V.1), type = 'F')^2 + norm(S.2 - U.2%*%t(V.2), type = 'F')^2)/(n1^2 + n2^2)
          temp.codesite.p <- lambda.delta*sqrt(log(p)/n1)*sum(apply(delta.int.U[1:n1,], 1, norm, "2"))/(n1^2 + n2^2) +lambda.delta*sqrt(log(p)/n2)*sum(apply(delta.int.U[(1+n1):(n1+n2),], 1, norm, "2"))/(n1^2 + n2^2)
          message('\n diff_CodeSiteEff=', temp.codesite.main - temp.code.main + temp.codesite.p - loss_delta_U)

          loss_main <- (norm(S.1 - U.1%*%t(V.1), type = 'F')^2 + norm(S.2 - U.2%*%t(V.2), type = 'F')^2)/(n1^2 + n2^2)
          loss_zeta_U_ovl <- lambda*sqrt(log(p)/(n1+n2))*norm((zeta.int.U[1:n1, ])[rownames(zeta.int.U[1:n1, ])%in%common_codes, ],type = 'F')^2/(n1^2 + n2^2)
          loss_zeta_U_1 <- lambda*sqrt(log(p)/n1)*norm((zeta.int.U[1:n1, ])[!(rownames(zeta.int.U[1:n1, ])%in%common_codes), ],type = 'F')^2/(n1^2 + n2^2)
          loss_zeta_U_2 <- lambda*sqrt(log(p)/n2)*norm((zeta.int.U[(1+n1):(n1+n2), ])[!(rownames(zeta.int.U[(1+n1):(n1+n2), ])%in%common_codes), ],type = 'F')^2/(n1^2 + n2^2)
          loss_zeta_U <- loss_zeta_U_ovl + loss_zeta_U_1 + loss_zeta_U_2

          loss_delta_U.1 <- lambda.delta*sqrt(log(p)/n1)*sum(apply(delta.int.U[1:n1,], 1, norm, "2"))/(n1^2 + n2^2)
          loss_delta_U.2 <- lambda.delta*sqrt(log(p)/n2)*sum(apply(delta.int.U[(1+n1):(n1+n2),], 1, norm, "2"))/(n1^2 + n2^2)

          loss_delta_U <- loss_delta_U.1 + loss_delta_U.2

          loss_p_U <- loss_zeta_U + loss_delta_U
          loss.new.U <- loss_main + loss_p_V + loss_p_U
          dif_loss_U <- (loss.new.U - loss.U)/loss.U
          message('\n dif_loss_U',dif_loss_U)
          loss.U <- loss.new.U
        }
        dif_loss <-  (loss.new.U - loss)/loss
        message('\n dif_loss=', dif_loss)
      }
      delta.int.U.2 <- delta.int.U[(n1.no+1):(n1+n2) ,]
      delta.spst.ovl[k1, k2]  <- sum(rowSums(abs(delta.int.U.2[rownames(delta.int.U.2)%in%common_codes,]))!=0)
      # validation
      if (TUNE==TRUE){
        ans <- evaluation.sim(pairs.rel.CV,  U.2)
        CV.res[k1, k2] <- ans$AUC.rel
      }
    }
  }

  # Tuning if TUNE==T and output the selected tuning parameters and the number of site-dissimilar codes
  if (TUNE==TRUE){
    idx <- which(CV.res == max(CV.res), arr.ind = TRUE)
    lambda.opt <- Lambda[idx[1,1]]
    lambda.delta.opt <- Lambda.delta[idx[1,2]]
    message("lambda1 = ", lambda.opt)
    message("lambda2 = ", lambda.delta.opt)
    message("The number of site-heterogeneous codes = ", delta.spst.ovl)
  } else {
    # Output the embedding matrices, cosine similarity matrices, and the names of similar codes and dissimilar codes if TUNE==F
    if (Eva ==TRUE ){
      ans <- evaluation.sim(pairs.rel.EV,  U.2)
      message('\n evaluation AUC =', ans$AUC.rel )
    }
    if (is.null(outdir)) outdir <- tempdir()
    save(U.1, file = file.path(outdir, "U.1.Rdata"))
    save(U.2, file = file.path(outdir, "U.2.Rdata"))
    CS.1 <- (U.1/apply(U.1,1,norm,'2'))%*%t(U.1/apply(U.1,1,norm,'2'))
    CS.2 <- (U.2/apply(U.2,1,norm,'2'))%*%t(U.2/apply(U.2,1,norm,'2'))
    save(CS.1, file = file.path(outdir, "CS.1.Rdata"))
    save(CS.2, file = file.path(outdir, "CS.2.Rdata"))
    save(beta.int.U, file = file.path(outdir, "beta.Rdata"))
    # Similar codes between sites vs dissimilar codes between sites
    similar.codes <-intersect (rownames(delta.int.U)[rowSums(abs(delta.int.U))==0], common_codes )
    save(similar.codes, file = file.path(outdir, "similar.codes.Rdata"))
    dissimilar.codes <- rownames(delta.int.U)[rowSums(abs(delta.int.U))!=0]
    save(dissimilar.codes, file = file.path(outdir, "dissimilar.codes.Rdata"))
    #Sparsity of code-site effect
    message("The number of site-heterogeneous codes = ", delta.spst.ovl)
  }
}
