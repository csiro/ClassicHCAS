#' ClassicHCAS colour palettes
#'
#' Returns colour palettes used by ClassicHCAS plots.
#'
#' @details
#' The \code{"hcas"} palette is intended for habitat condition maps, with low
#' condition shown in purple, intermediate values in yellow, and high condition
#' in green. The \code{"ref_density"} palette is intended for raw or normalised
#' reference density surfaces.
#'
#' @param n Integer. Number of colour codes to return.
#' @param name Character. Palette name. One of \code{"hcas"} or
#' \code{"ref_density"}.
#'
#' @return A character vector of hexadecimal colour codes.
#' @export
#'
#' @examples
#' palettes(5)
#' palettes(5, "ref_density")
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
