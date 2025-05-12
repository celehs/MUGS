## Test environments
- Local: Windows 11, R 4.4.0
- R-hub:
  - Windows (R-devel)
  - macOS (R-devel)
  - Linux (R-devel)
- Win-builder: Windows (R-devel)

## R CMD check results
All checks passed:
- 0 errors ✔
- 0 warnings ✔
- 0 notes ✔

## Comments
This is the initial CRAN submission of the package.
* All acronyms in the DESCRIPTION (e.g., EHR) have been explained.
* The DESCRIPTION file includes a properly formatted reference to the associated publication:
  Li et al. (2024) <doi:10.1038/s41746-024-01320-4>.
* All examples that download data are wrapped in `\dontrun{}` to avoid timeouts.
* `TRUE` and `FALSE` are used instead of `T` and `F`.
* No functions modify `.GlobalEnv`, write to user filespace by default, or set random seeds internally.
* Messages are displayed using `message()` or `warning()` rather than `cat()` or `print()`.
