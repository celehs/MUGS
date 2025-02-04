# R/load_data.R

#' Download Datasets from Zenodo
#'
#' This function downloads the datasets from Zenodo and saves them to a specified directory.
#'
#' @param destdir (optional) The directory where the files should be saved.
#'                If `NULL`, the files will be saved in a temporary directory.
#' @return A list of file paths for the downloaded datasets.
#' @export
#'
#' @examples
#' \dontrun{
#' files <- download_zenodo_data()
#' }
download_zenodo_data <- function(destdir = NULL) {
  # Zenodo record URL (replace with your Zenodo record URL)
  zenodo_base_url <- "https://zenodo.org/record/14776064/files/"

  # List of dataset files
  dataset_files <- c(
    "pairs.rel.CV.Rdata",
    "pairs.rel.EV.Rdata",
    "S.1.Rdata",
    "S.2.Rdata",
    "U.1.Rdata",
    "U.2.Rdata",
    "X.group.source.Rdata",
    "X.group.target.Rdata"
  )

  # Set the destination directory
  if (is.null(destdir)) {
    destdir <- tempdir()
  }

  # Download each file
  file_paths <- sapply(dataset_files, function(file) {
    url <- paste0(zenodo_base_url, file)
    destfile <- file.path(destdir, file)
    utils::download.file(url, destfile, mode = "wb")
    return(destfile)
  })

  message("Datasets downloaded successfully to: ", destdir)
  return(file_paths)
}

#' Load pairs.rel.CV Dataset
#'
#' This function loads the `pairs.rel.CV` dataset into R.
#'
#' @param filepath (optional) Path to the downloaded `.Rdata` file.
#'                 If `NULL`, the file will be downloaded automatically.
#' @return The loaded dataset.
#' @export
#'
#' @examples
#' \dontrun{
#' pairs.rel.CV <- load_pairs_rel_CV()
#' }
load_pairs_rel_CV <- function(filepath = NULL) {
  if (is.null(filepath)) {
    filepath <- download_zenodo_data()[["pairs.rel.CV.Rdata"]]
  }
  load(filepath)
  return(pairs.rel.CV)
}

#' Load pairs.rel.EV Dataset
#'
#' This function loads the `pairs.rel.EV` dataset into R.
#'
#' @param filepath (optional) Path to the downloaded `.Rdata` file.
#'                 If `NULL`, the file will be downloaded automatically.
#' @return The loaded dataset.
#' @export
#'
#' @examples
#' \dontrun{
#' pairs.rel.EV <- load_pairs_rel_EV()
#' }
load_pairs_rel_EV <- function(filepath = NULL) {
  if (is.null(filepath)) {
    filepath <- download_zenodo_data()[["pairs.rel.EV.Rdata"]]
  }
  load(filepath)
  return(pairs.rel.EV)
}

#' Load S.1 Dataset
#'
#' This function loads the `S.1` dataset into R.
#'
#' @param filepath (optional) Path to the downloaded `.Rdata` file.
#'                 If `NULL`, the file will be downloaded automatically.
#' @return The loaded dataset.
#' @export
#'
#' @examples
#' \dontrun{
#' S.1 <- load_S_1()
#' }
load_S_1 <- function(filepath = NULL) {
  if (is.null(filepath)) {
    filepath <- download_zenodo_data()[["S.1.Rdata"]]
  }
  load(filepath)
  return(S.1)
}

#' Load S.2 Dataset
#'
#' This function loads the `S.2` dataset into R.
#'
#' @param filepath (optional) Path to the downloaded `.Rdata` file.
#'                 If `NULL`, the file will be downloaded automatically.
#' @return The loaded dataset.
#' @export
#'
#' @examples
#' \dontrun{
#' S.2 <- load_S_2()
#' }
load_S_2 <- function(filepath = NULL) {
  if (is.null(filepath)) {
    filepath <- download_zenodo_data()[["S.2.Rdata"]]
  }
  load(filepath)
  return(S.2)
}

#' Load U.1 Dataset
#'
#' This function loads the `U.1` dataset into R.
#'
#' @param filepath (optional) Path to the downloaded `.Rdata` file.
#'                 If `NULL`, the file will be downloaded automatically.
#' @return The loaded dataset.
#' @export
#'
#' @examples
#' \dontrun{
#' U.1 <- load_U_1()
#' }
load_U_1 <- function(filepath = NULL) {
  if (is.null(filepath)) {
    filepath <- download_zenodo_data()[["U.1.Rdata"]]
  }
  load(filepath)
  return(U.1)
}

#' Load U.2 Dataset
#'
#' This function loads the `U.2` dataset into R.
#'
#' @param filepath (optional) Path to the downloaded `.Rdata` file.
#'                 If `NULL`, the file will be downloaded automatically.
#' @return The loaded dataset.
#' @export
#'
#' @examples
#' \dontrun{
#' U.2 <- load_U_2()
#' }
load_U_2 <- function(filepath = NULL) {
  if (is.null(filepath)) {
    filepath <- download_zenodo_data()[["U.2.Rdata"]]
  }
  load(filepath)
  return(U.2)
}

#' Load X.group.source Dataset
#'
#' This function loads the `X.group.source` dataset into R.
#'
#' @param filepath (optional) Path to the downloaded `.Rdata` file.
#'                 If `NULL`, the file will be downloaded automatically.
#' @return The loaded dataset.
#' @export
#'
#' @examples
#' \dontrun{
#' X.group.source <- load_X_group_source()
#' }
load_X_group_source <- function(filepath = NULL) {
  if (is.null(filepath)) {
    filepath <- download_zenodo_data()[["X.group.source.Rdata"]]
  }
  load(filepath)
  return(X.group.source)
}

#' Load X.group.target Dataset
#'
#' This function loads the `X.group.target` dataset into R.
#'
#' @param filepath (optional) Path to the downloaded `.Rdata` file.
#'                 If `NULL`, the file will be downloaded automatically.
#' @return The loaded dataset.
#' @export
#'
#' @examples
#' \dontrun{
#' X.group.target <- load_X_group_target()
#' }
load_X_group_target <- function(filepath = NULL) {
  if (is.null(filepath)) {
    filepath <- download_zenodo_data()[["X.group.target.Rdata"]]
  }
  load(filepath)
  return(X.group.target)
}
