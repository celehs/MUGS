#' Download and Load Example Data from Zenodo
#'
#' @param file Name of the .Rdata file to download (e.g., "S.1.Rdata").
#' @param destdir Directory to store the downloaded data. Defaults to a temporary directory.
#'
#' @return A list containing the loaded dataset.
#' @export
#'

download_example_data <- function(file, destdir = tempdir()) {
  base_url <- "https://zenodo.org/records/14776064/files/"
  dest <- file.path(destdir, file)

  if (!file.exists(dest)) {
    url <- paste0(base_url, file, "?download=1")
    message("Downloading ", file, " from Zenodo...")
    utils::download.file(url, dest, mode = "wb")
  }

  env <- new.env()
  load(dest, envir = env)
  return(as.list(env))
}
