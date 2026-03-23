#' ClassicHCAS palettes
#'
#' @param n Integer. Number of color codes to return.
#' @param name Character. Palette name. One of `"hcas"` or `"ref_density"`.
#'
#' @return A character vector of color codes.
#' @export
palettes <- function(n = 10, name = c("hcas", "ref_density")) {
    name <- match.arg(name)

    palette_fn <- switch(
        name,
        hcas = grDevices::colorRampPalette(c("#430E59", "#CCCC66", "#184F0F")),
        ref_density = grDevices::colorRampPalette(
            c(
                "#ffffff", "#ffffd9", "#edf8b1", "#c7e9b4", "#7fcdbb",
                "#41b6c4", "#1d91c0", "#225ea8", "#253494", "#081d58"
            )
        )
    )

    palette_fn(n)
}

