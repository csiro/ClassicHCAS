#' The HCAS colors
#'
#' @param n Integer. Number of color codes to return.
#'
#' @return a charterer vector of color codes
#' @export
hcas_color <- function(n = 10) {
    hcas_pal <- grDevices::colorRampPalette(c("#430E59", "#CCCC66", "#184F0F"))

    return(
        hcas_pal(n)
    )
}


#' The HCAS reference density colors
#'
#' @param n Integer. Number of color codes to return.
#'
#' @return a charterer vector of color codes
#' @export
ref_density_color <- function(n = 10) {
    ref_cols <- grDevices::colorRampPalette(
        c("#ffffff", "#ffffd9", "#edf8b1", "#c7e9b4", "#7fcdbb",
          "#41b6c4", "#1d91c0", "#225ea8", "#253494", "#081d58"
        )
    )

    return(
        ref_cols(n)
    )
}

