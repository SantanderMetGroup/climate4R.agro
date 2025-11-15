#' @importFrom utils packageDescription

.onAttach <- function(...) {
    pkgname <- "climate4R.agro"
    pkg_desc <- utils::packageDescription(pkgname)
    ver_raw <- pkg_desc$Version
    ver <- tryCatch(package_version(ver_raw), error = function(...) NA)
    ver_label <- if (is.null(ver_raw) || is.na(ver_raw)) "<unknown>" else ver_raw
    builddate <- pkg_desc$Date
    if (is.null(builddate) || is.na(builddate)) {
        builddate <- pkg_desc$`Date/Publication`
    }
    if (is.null(builddate) || is.na(builddate)) {
        builddate <- "unknown date"
    }
    packageStartupMessage(sprintf("%s version %s (%s) is loaded", pkgname, ver_label, builddate))
    url <- paste0("https://raw.githubusercontent.com/SantanderMetGroup/", pkgname, "/master/DESCRIPTION")
    b <- tryCatch(suppressWarnings(readLines(url)), error = function(er) {
        er <- NULL
        return(er)
    })
    if (!is.null(b)) {
        latest_str <- gsub("Version: ", "", b[grep("^Version", b)])
        latest.ver <- tryCatch(package_version(latest_str), error = function(...) NA)
        if (!is.na(ver) && !is.na(latest.ver)) {
            if (ver < latest.ver) {
                ver.mess1 <- paste0("WARNING: Your current version of ", pkgname, " (v", ver_label, ") is not up-to-date")
                ver.mess <- paste0("Get the latest stable version (", latest.ver,
                                   ") using <devtools::install_github('SantanderMetGroup/", pkgname, "')>")
                packageStartupMessage(ver.mess1)
                packageStartupMessage(ver.mess)
            } else if (ver > latest.ver) {
                ver.mess1 <- paste0("WARNING: Your current version of ", pkgname, " (v", ver_label,
                                    ") is ahead of the master branch version (", latest.ver, ")")
                ver.mess <- "Development version may have unexpected behaviour"
                packageStartupMessage(ver.mess1)
                packageStartupMessage(ver.mess)
            }
        }
    }
    packageStartupMessage("Use 'agroindexShow()' for an overview of the available indices")
}
# End
