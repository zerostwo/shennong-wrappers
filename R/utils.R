#' Ensure a GitHub Package is Installed
#'
#' This internal utility checks whether a specified package is installed.
#' If not, it installs the package from GitHub using the `pak` package.
#'
#' @param pkg A string. The name of the package to check.
#' @param repo A string of the form `"owner/repo"`, specifying the GitHub repository.
#'
#' @details
#' This function is primarily used to ensure that packages not available on CRAN
#' (e.g., `ROGUE`) are available at runtime. It leverages `rlang::check_installed()`
#' with a custom install `action` that calls `pak::pak("github::owner/repo")`.
#'
#' @return Called for side effect. Installs the package if necessary.
#'
#' @seealso [rlang::check_installed()], [pak::pak()]
#'
#' @keywords internal
#' @noRd
check_installed_github <- function(pkg, repo) {
  check_installed(pkg = "pak")
  check_installed(pkg, action = function(pkg, ...) pak::pak(repo))
}
