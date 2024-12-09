#' S.1 Dataset
#'
#' @description A matrix containing SPPMI data from the source site.
#' This dataset is used as input for analysis in the package.
#'
#' @format A matrix with 2000 rows and 10 columns:
#' \describe{
#'   \item{Row Names}{Unique identifiers for each row.}
#'   \item{Columns}{Numeric values representing SPPMI data.}
#' }
"S.1"

#' S.2 Dataset
#'
#' @description A matrix containing SPPMI data from the target site.
#' This dataset is used as input for analysis in the package.
#'
#' @format A matrix with 2000 rows and 10 columns:
#' \describe{
#'   \item{Row Names}{Unique identifiers for each row.}
#'   \item{Columns}{Numeric values representing SPPMI data.}
#' }
"S.2"

#' U.1 Dataset
#'
#' @description A matrix containing left embeddings from the source site.
#' These embeddings are used for embedding-based computations.
#'
#' @format A matrix with 2000 rows and 10 columns:
#' \describe{
#'   \item{Row Names}{Unique identifiers for each row.}
#'   \item{Columns}{Numeric values representing embeddings.}
#' }
"U.1"

#' U.2 Dataset
#'
#' @description A matrix containing left embeddings from the target site.
#' These embeddings are used for embedding-based computations.
#'
#' @format A matrix with 2000 rows and 10 columns:
#' \describe{
#'   \item{Row Names}{Unique identifiers for each row.}
#'   \item{Columns}{Numeric values representing embeddings.}
#' }
"U.2"

#' X.group.source Dataset
#'
#' @description A matrix containing group structures at the source site.
#' It represents binary group membership of entities at the source.
#'
#' @format A matrix with 2000 rows and 50 columns:
#' \describe{
#'   \item{Rows}{Entities at the source site.}
#'   \item{Columns}{Binary values (0 or 1) indicating group membership.}
#' }
"X.group.source"

#' X.group.target Dataset
#'
#' @description A matrix containing group structures at the target site.
#' It represents binary group membership of entities at the target.
#'
#' @format A matrix with 2000 rows and 50 columns:
#' \describe{
#'   \item{Rows}{Entities at the target site.}
#'   \item{Columns}{Binary values (0 or 1) indicating group membership.}
#' }
"X.group.target"

#' pairs.rel.CV Dataset
#'
#' @description A data frame containing cross-validation pairs for relative comparisons.
#'
#' @format A data frame with multiple columns:
#' \describe{
#'   \item{col}{Integer representing the column index of a pair.}
#'   \item{row}{Integer representing the row index of a pair.}
#'   \item{type}{Character string indicating the type of data (e.g., "train", "test").}
#' }
"pairs.rel.CV"

#' pairs.rel.EV Dataset
#'
#' @description A data frame containing evaluation pairs for relative comparisons.
#'
#' @format A data frame with multiple columns:
#' \describe{
#'   \item{col}{Integer representing the column index of a pair.}
#'   \item{row}{Integer representing the row index of a pair.}
#'   \item{type}{Character string indicating the type of data (e.g., "validation").}
#' }
"pairs.rel.EV"
