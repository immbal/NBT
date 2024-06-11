#' Check if required columns are present in the data frame
#'
#' This function checks whether the specified required columns are present in the given data frame.
#' If any of the required columns are missing, it stops the execution and returns an error message.
#'
#' @param df A data frame in which to check for required columns.
#' @param required_cols A character vector of required column names.
#' @return None. The function stops execution with an error message if any required columns are missing.
#' @examples
#' df <- data.frame(ID = 1:10, Value = rnorm(10))
#' check_columns(df, c("ID", "Value")) # No error
#' check_columns(df, c("ID", "Value", "MissingColumn")) # Error: The following columns are missing: MissingColumn
#' @export
check_columns <- function(df, required_cols) {
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop(paste("The following columns are missing:", paste(missing_cols, collapse = ", ")))
  }
} 


#' Check if required packages are installed
#'
#' This function checks if the specified packages are installed. If any of the packages are not installed,
#' it stops the execution and returns an error message.
#'
#' @param packages A character vector of package names to check.
#' @return None. The function stops execution with an error message if any packages are not installed.
#' @examples
#' check_required_packages(c("Hmisc", "progress")) # No error if both packages are installed
#' check_required_packages(c("Hmisc", "nonexistent_package")) # Error: The following packages are required but not installed: nonexistent_package
#' @export
check_required_packages <- function(packages) {
  missing_packages <- packages[!sapply(packages, requireNamespace, quietly = TRUE)]
  if (length(missing_packages) > 0) {
    stop("The following packages are required but not installed: ", paste(missing_packages, collapse = ", "))
  }
}