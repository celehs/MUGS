#' Download Datasets from Zenodo
#'
#' This function downloads large datasets hosted on Zenodo for use with the MUGS package.
#'
#' @param dataset The name of the dataset to download (e.g., "pairs.rel.CV").
#' @param destfile The destination file path where the dataset will be saved. Defaults to the current working directory.
#' @return A message confirming the dataset has been downloaded.
#' @examples
#' \dontrun{
#'   download_dataset("pairs.rel.CV", "data/pairs.rel.CV.rda")
#' }
#' @export
download_dataset <- function(dataset, destfile = paste0(dataset, ".rda")) {
  base_url <- "https://zenodo.org/record/<Zenodo_Record_ID>/files/"
  url <- paste0(base_url, dataset, ".rda")
  download.file(url, destfile, mode = "wb")
  message("Dataset ", dataset, " has been downloaded and saved to ", destfile)
}
