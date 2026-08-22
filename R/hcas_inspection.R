#' Launch an interactive HCAS point-inspection tool
#'
#' Builds a Shiny application for inspecting how individual target locations
#' are benchmarked against reference samples. The application shows the raw
#' condition estimate, the reference samples inside the geographic search
#' radius, the samples retained by the two-stage HCAS selection, and their
#' predicted and observed distances on the reference-density surface.
#'
#' @details
#' The inspection tool uses the current \code{\link{benchmark}} and
#' \code{\link{reference_use}} implementations. Consequently, its controls map
#' directly to the current benchmarking arguments: \code{k1} is the first-stage
#' predicted-distance filter, \code{k2} is the reference-density filter, and
#' \code{xy_penalty} applies the optional scaled-coordinate penalty. Feature
#' counts, geographic distance handling, density dimensions, and density
#' metadata are inferred by ClassicHCAS rather than supplied as legacy fixed
#' values.
#'
#' The user interface intentionally retains the standalone inspection tool's
#' default Shiny theme, \pkg{shinyWidgets} controls, layout, labels, map styling,
#' and reference-density plot styling.
#'
#' Matrix and data.frame inputs must be ordered as \code{x}, \code{y}, predicted
#' RS variables, then observed RS variables. Raster inputs must contain
#' predicted RS layers followed by observed RS layers in the same variable
#' order. For raster inputs, the target selector uses \code{x} and \code{y}
#' coordinates in the raster coordinate reference system. Map clicks are
#' projected back to that coordinate system and snapped to the containing raster
#' cell centre before values are extracted. Temporal reference-sample lists are
#' not supported because \code{\link{reference_use}} currently supports only a
#' single reference year.
#'
#' The map requires longitude/latitude coordinates. If the input coordinates
#' are projected, supply their coordinate reference system through \code{crs};
#' the display coordinates are then transformed to EPSG:4326. The CRS affects
#' map display only. Benchmarking uses the coordinate handling implemented by
#' \code{\link{benchmark}}.
#'
#' @inheritParams benchmark
#' @param data A matrix, data.frame, or \pkg{terra} \code{SpatRaster} of target
#' locations. Matrix and data.frame inputs use a row-ID selector. Raster inputs
#' use coordinate inputs and extract the selected cell values.
#' @param samples A matrix or data.frame of reference locations in the same
#' column order as \code{data}. For raster \code{data}, this may contain only
#' \code{x} and \code{y} coordinates, in which case raster values are extracted.
#' Temporal sample lists are not supported.
#' @param crs Optional coordinate reference system understood by
#' \code{\link[terra]{vect}}, such as \code{"EPSG:3577"}. Required for the map
#' when matrix/data.frame coordinates are projected or when raster coordinates
#' are projected and the raster has no CRS.
#' @param background Optional background raster drawn beneath the sample
#' markers, given as a file path to a GeoTIFF or a \pkg{terra}
#' \code{SpatRaster}. The map renderer reads web mercator directly, so a
#' single-band Cloud-Optimized GeoTIFF (COG) already in \code{"EPSG:3857"} or
#' \code{"EPSG:4326"} is served unchanged and its internal overviews stream for
#' fast, overview-accelerated rendering. Any other input is reduced to its
#' first layer, reprojected to \code{"EPSG:3857"} when needed, and written to a
#' temporary COG with overviews. Requires the suggested package \pkg{leafem}.
#' @param background_colors Optional character vector of colours (for example
#' hexadecimal codes such as \code{c("#430E59", "#CCCC66", "#184F0F")}) used as
#' the continuous colour ramp for \code{background}. Defaults to the
#' \code{"hcas"} palette from \code{\link{palettes}}. Ignored when
#' \code{background} is \code{NULL}.
#' @param launch Logical. If \code{TRUE}, run the application with
#' \code{\link[shiny]{runApp}}. If \code{FALSE}, return the application object
#' without running it, which is useful for testing or custom deployment.
#' @param confidence Numeric between 0 and 1. Weight given to the selected
#' maximum probability component relative to the distance-weighted mean
#' probability when computing raw condition. Ignored when \code{boost} is not
#' \code{NULL} or \code{NA}; the inspection app default is \code{boost = k2}.
#' @param boost \code{NULL}, \code{NA}, or one positive finite numeric factor.
#' The inspection app default \code{k2} multiplies the kernel weight of the
#' highest-probability retained site by \code{k2} and returns the resulting
#' weighted mean instead of the LDC blend. \code{confidence} is ignored in this
#' mode. Use \code{NULL} or \code{NA} for the unboosted LDC blend.
#' @param ... Additional arguments passed to \code{\link[shiny]{runApp}} when
#' \code{launch = TRUE}, such as \code{host}, \code{port}, or
#' \code{launch.browser}.
#'
#' @return A Shiny application object when \code{launch = FALSE}. When
#' \code{launch = TRUE}, the return value from \code{\link[shiny]{runApp}} is
#' returned invisibly after the application stops.
#' @seealso \code{\link{benchmark}}, \code{\link{reference_use}}
#' @export
#'
#' @examples
#' \dontrun{
#' app <- hcas_inspection(
#'     data = reference_samples,
#'     samples = reference_samples,
#'     ref_density = normalised_density,
#'     crs = "EPSG:3577",
#'     launch = FALSE
#' )
#' shiny::runApp(app)
#' }
hcas_inspection <- function(
        data,
        samples,
        ref_density,
        xy_stats = c(0, 0, 1, 1),
        xy_penalty = 0.0,
        radius_km = 200,
        k1 = 70,
        k2 = 10,
        bin_width = NULL,
        interpolate = TRUE,
        offset = 0,
        confidence = 0.5,
        lambda = 1.0,
        exclude_slef = FALSE,
        drop_features = NULL,
        num_threads = -1,
        crs = NULL,
        background = NULL,
        background_colors = NULL,
        launch = interactive(),
        kernel = c("Gaussian", "Cauchy"),
        boost = k2,
        ...) {

    kernel <- .check_kernel(kernel)
    required <- c("shiny", "shinyWidgets", "leaflet", "ggplot2")
    available <- vapply(required, requireNamespace, logical(1), quietly = TRUE)
    if (!all(available)) {
        stop(
            "'hcas_inspection()' requires the suggested package(s): ",
            paste(required[!available], collapse = ", "),
            "."
        )
    }
    if (!is.null(background) && !requireNamespace("leafem", quietly = TRUE)) {
        stop("The 'background' map requires the suggested package 'leafem'.")
    }

    raster_data <- .is_rast(data)
    data <- if (raster_data) {
        .check_rast(data, name = "data")
    } else if (.is_mat(data)) {
        .check_mat(data, name = "data")
    } else {
        stop("'data' must be a raster, matrix, or data.frame.")
    }
    if (raster_data && !inherits(data, "SpatRaster")) {
        stop("'data' must be convertible to a terra SpatRaster object.")
    }
    if (raster_data && .inspection_has_crs(crs)) {
        terra::crs(data) <- crs
    }
    samples <- if (.is_mat(samples)) {
        .check_mat(samples, name = "samples")
    } else {
        stop("'samples' must be a matrix or data.frame, not a temporal list.")
    }
    if (!raster_data && !nrow(data)) {
        stop("'data' must contain at least one target row.")
    }
    if (raster_data && terra::ncell(data) < 1L) {
        stop("'data' raster must contain at least one cell.")
    }
    if (!nrow(samples)) {
        stop("'samples' must contain at least one reference row.")
    }
    if (raster_data) {
        if (terra::nlyr(data) %% 2L) {
            stop("'data' raster must contain matching predicted and observed RS layers.")
        }
        samples <- .inspection_raster_samples(data, samples)
    } else {
        if (ncol(data) != ncol(samples)) {
            stop("'data' and 'samples' must have matching columns.")
        }
        .num_rs_vars_mat(data, "data")
    }
    .num_rs_vars_mat(samples, "samples")

    background_file <- .inspection_background_cog(background)
    background_colors <- .inspection_background_colors(background_colors)

    if (length(k1) != 1L || !is.finite(k1) || k1 < 1) {
        stop("'k1' must be one finite number greater than or equal to one.")
    }
    if (length(k2) != 1L || !is.finite(k2) || k2 < 1 || k2 > k1) {
        stop("'k2' must be between one and 'k1'.")
    }
    boost <- .check_boost(boost)

    density <- .inspection_density(
        ref_density = ref_density,
        bin_width = bin_width,
        offset = offset,
        interpolate = interpolate
    )
    geographic <- .is_lonlat(data)
    map_crs <- if (raster_data) {
        .inspection_data_crs(data, crs)
    } else {
        crs
    }
    if (!geographic && !.inspection_has_crs(map_crs)) {
        stop("'crs' must be supplied to display projected coordinates on the map.")
    }

    condition_mode_default <- if (is.null(boost)) "ldc" else "boost"
    boost_default <- if (is.null(boost)) as.numeric(k2) else boost
    condition_value_default <- if (condition_mode_default == "ldc") {
        confidence
    } else {
        boost_default
    }
    condition_value_label <- if (condition_mode_default == "ldc") {
        "LDC confidence"
    } else {
        "Boost factor"
    }
    condition_value_max <- if (condition_mode_default == "ldc") 1 else Inf
    condition_value_step <- if (condition_mode_default == "ldc") 0.05 else 1

    radius_choices <- sort(unique(c(seq(0, 500, 50), radius_km)))
    k1_choices <- sort(unique(c(seq(0, 100, 5), as.integer(k1))))
    k2_choices <- sort(unique(c(seq(0, 50, 5), as.integer(k2))))
    map_bounds <- if (raster_data) {
        .inspection_raster_map_bounds(data, map_crs, geographic)
    } else {
        NULL
    }
    target_selector <- if (raster_data) {
        default_xy <- .inspection_raster_default_xy(data)
        shiny::fluidRow(
            shiny::column(
                width = 4,
                shiny::numericInput(
                    "target_x", "X coordinate:",
                    value = default_xy[1],
                    min = terra::xmin(data),
                    max = terra::xmax(data)
                )
            ),
            shiny::column(
                width = 4,
                shiny::numericInput(
                    "target_y", "Y coordinate:",
                    value = default_xy[2],
                    min = terra::ymin(data),
                    max = terra::ymax(data)
                )
            ),
            shiny::column(
                width = 3,
                shiny::HTML("<br/>"),
                shinyWidgets::actionBttn(
                    inputId = "submit",
                    style = "pill",
                    color = "default",
                    size = "md",
                    label = "Run",
                    icon = shiny::icon("play")
                )
            )
        )
    } else {
        shiny::fluidRow(
            shiny::column(
                width = 6,
                shiny::numericInput(
                    "point_id", "Enter sample id:",
                    value = 1, min = 1, max = nrow(data)
                )
            ),
            shiny::column(
                width = 3,
                shiny::HTML("<br/>"),
                shinyWidgets::actionBttn(
                    inputId = "submit",
                    style = "pill",
                    color = "default",
                    size = "md",
                    label = "Run",
                    icon = shiny::icon("play")
                )
            )
        )
    }

    ui <- shiny::fluidPage(
        shiny::tags$head(
            shiny::tags$style(shiny::HTML(
                "
                #condition_mode .radiobtn.btn,
                #kernel .radiobtn.btn {
                    background-color: #e6e9ed;
                    border-color: #c7ccd1;
                    color: #2f3942;
                }

                #condition_mode .radiobtn.btn.active,
                #kernel .radiobtn.btn.active {
                    background-color: #2f80c1;
                    border-color: #246aa2;
                    color: #ffffff;
                }

                #condition_mode .radiobtn.btn:not(.active):hover,
                #kernel .radiobtn.btn:not(.active):hover {
                    background-color: #d8dde3;
                    border-color: #b7bec6;
                    color: #25313b;
                }
                "
            ))
        ),
        shiny::titlePanel("HCAS inspection tool"),
        shiny::sidebarLayout(
            shiny::sidebarPanel(
                target_selector,
                shiny::h4("Benchmarking Options:"),
                shiny::fluidRow(
                    shiny::column(
                        width = 6,
                        shinyWidgets::radioGroupButtons(
                            inputId = "condition_mode",
                            label = "Condition estimator",
                            choices = c(
                                "LDC confidence" = "ldc",
                                "Boost" = "boost"
                            ),
                            selected = condition_mode_default,
                            justified = TRUE,
                            status = "default",
                            size = "sm"
                        )
                    ),
                    shiny::column(
                        width = 3,
                        shinyWidgets::numericInputIcon(
                            inputId = "condition_value",
                            value = condition_value_default,
                            step = condition_value_step,
                            min = 0,
                            max = condition_value_max,
                            label = condition_value_label
                        )
                    )
                ),
                shiny::fluidRow(
                    shiny::column(
                        width = 3,
                        shinyWidgets::numericInputIcon(
                            inputId = "xy_penalty",
                            value = xy_penalty,
                            step = 0.5,
                            min = 0,
                            max = Inf,
                            label = "Geo-penalty"
                        )
                    ),
                    shiny::column(
                        width = 9,
                        shiny::h3(" "),
                        shiny::HTML("<br/>"),
                        shinyWidgets::materialSwitch(
                            inputId = "exclude_self",
                            value = exclude_slef,
                            label = "Exclude self-assessment",
                            status = "primary",
                            inline = FALSE,
                            right = TRUE
                        )
                    )
                ),
                shiny::fluidRow(
                    shiny::column(
                        width = 6,
                        shinyWidgets::radioGroupButtons(
                            inputId = "kernel",
                            label = "Kernel",
                            choices = c(
                                "Gaussian" = "gaussian",
                                "Cauchy" = "cauchy"
                            ),
                            selected = kernel,
                            justified = TRUE,
                            status = "default",
                            size = "sm"
                        )
                    ),
                    shiny::column(
                        width = 3,
                        shinyWidgets::numericInputIcon(
                            inputId = "lambda",
                            value = lambda,
                            step = 0.1,
                            min = 0,
                            max = Inf,
                            label = "Lambda"
                        )
                    )
                ),
                shinyWidgets::sliderTextInput(
                    inputId = "radius_km",
                    label = "Search radius",
                    choices = radius_choices,
                    selected = radius_km,
                    grid = TRUE
                ),
                shinyWidgets::sliderTextInput(
                    inputId = "k1",
                    label = "Number of ENV neighbours",
                    choices = k1_choices,
                    selected = as.integer(k1),
                    grid = TRUE
                ),
                shinyWidgets::sliderTextInput(
                    inputId = "k2",
                    label = "Number of RS neighbours",
                    choices = k2_choices,
                    selected = as.integer(k2),
                    grid = TRUE
                ),
                shiny::textOutput("condition"),
                shiny::textOutput("nearby_count"),
                shiny::plotOutput("density_plot", width = "100%", height = "300px")
            ),
            shiny::mainPanel(
                leaflet::leafletOutput("map", height = 1000)
            )
        )
    )

    server <- function(input, output, session) {
        condition_values <- shiny::reactiveValues(
            ldc = confidence,
            boost = boost_default
        )
        previous_condition_mode <- shiny::reactiveVal(condition_mode_default)
        get_condition_value <- function(mode) {
            if (mode == "ldc") {
                condition_values$ldc
            } else {
                condition_values$boost
            }
        }
        set_condition_value <- function(mode, value) {
            if (mode == "ldc") {
                condition_values$ldc <- value
            } else {
                condition_values$boost <- value
            }
        }

        shiny::observeEvent(input$condition_mode, {
            old_mode <- previous_condition_mode()
            if (!is.null(input$condition_value) &&
                    old_mode %in% c("ldc", "boost")) {
                set_condition_value(old_mode, input$condition_value)
            }

            new_mode <- input$condition_mode
            if (is.null(new_mode) || !(new_mode %in% c("ldc", "boost"))) {
                return()
            }
            previous_condition_mode(new_mode)

            new_value <- get_condition_value(new_mode)
            shiny::updateNumericInput(
                session,
                inputId = "condition_value",
                label = if (new_mode == "ldc") {
                    "LDC confidence"
                } else {
                    "Boost factor"
                },
                value = new_value,
                min = 0,
                max = if (new_mode == "ldc") 1 else Inf,
                step = if (new_mode == "ldc") 0.05 else 1
            )
        }, ignoreInit = FALSE)

        shiny::observeEvent(input$condition_value, {
            mode <- previous_condition_mode()
            if (mode %in% c("ldc", "boost")) {
                set_condition_value(mode, input$condition_value)
            }
        }, ignoreInit = TRUE)

        if (raster_data) {
            shiny::observeEvent(input$map_click, {
                clicked <- input$map_click
                target_xy <- tryCatch(
                    .inspection_raster_click_xy(
                        data,
                        lng = clicked$lng,
                        lat = clicked$lat,
                        crs = map_crs
                    ),
                    error = function(cond) {
                        shiny::showNotification(
                            conditionMessage(cond),
                            type = "warning"
                        )
                        NULL
                    }
                )
                if (is.null(target_xy)) {
                    return()
                }

                shiny::updateNumericInput(session, "target_x", value = target_xy[1])
                shiny::updateNumericInput(session, "target_y", value = target_xy[2])

                target_map <- .inspection_map_coordinates(
                    data.frame(x = target_xy[1], y = target_xy[2]),
                    map_crs,
                    geographic
                )
                leaflet::leafletProxy("map", session = session) |>
                    leaflet::clearGroup("target_preview") |>
                    leaflet::addCircleMarkers(
                        data = target_map,
                        lng = ~x,
                        lat = ~y,
                        color = "red",
                        radius = 4,
                        group = "target_preview",
                        label = "Selected target"
                    )
            }, ignoreNULL = TRUE)
        }

        inspected <- shiny::eventReactive(input$submit, {
            target <- NULL
            point_id <- NULL
            target_error <- NULL
            if (raster_data) {
                target_x <- input$target_x
                target_y <- input$target_y
                shiny::validate(
                    shiny::need(
                        length(target_x) == 1L && is.finite(target_x) &&
                            length(target_y) == 1L && is.finite(target_y),
                        "Target coordinates must be finite numbers."
                    )
                )
                target <- tryCatch(
                    .inspection_raster_target(data, target_x, target_y),
                    error = function(cond) {
                        target_error <<- conditionMessage(cond)
                        NULL
                    }
                )
                shiny::validate(shiny::need(is.null(target_error), target_error))
            } else {
                point_id <- as.integer(input$point_id)
            }
            radius_value <- as.numeric(input$radius_km)
            k1_value <- as.integer(input$k1)
            k2_value <- as.integer(input$k2)
            kernel_value <- if (is.null(input$kernel)) {
                NA_character_
            } else {
                input$kernel
            }
            lambda_value <- input$lambda
            condition_mode <- if (is.null(input$condition_mode)) {
                NA_character_
            } else {
                input$condition_mode
            }
            condition_value <- input$condition_value
            shiny::validate(
                shiny::need(
                    raster_data ||
                        (is.finite(point_id) &&
                            point_id >= 1L &&
                            point_id <= nrow(data)),
                    "Sample id is outside the available target rows."
                ),
                shiny::need(
                    condition_mode %in% c("ldc", "boost"),
                    "Choose either LDC confidence or Boost."
                ),
                shiny::need(
                    length(condition_value) == 1L && is.finite(condition_value),
                    "The estimator value must be finite."
                ),
                shiny::need(
                    condition_mode != "ldc" ||
                        (condition_value >= 0 && condition_value <= 1),
                    "LDC confidence must be between 0 and 1."
                ),
                shiny::need(
                    condition_mode != "boost" || condition_value > 0,
                    "Boost must be greater than zero."
                ),
                shiny::need(
                    kernel_value %in% c("gaussian", "cauchy"),
                    "Choose either Gaussian or Cauchy kernel."
                ),
                shiny::need(
                    length(lambda_value) == 1L &&
                        is.finite(lambda_value) &&
                        lambda_value > 0,
                    "Lambda must be greater than zero."
                ),
                shiny::need(k1_value >= 1L, "k1 must be greater than zero."),
                shiny::need(k2_value >= 1L, "k2 must be greater than zero."),
                shiny::need(k2_value <= k1_value, "k2 must be less than or equal to k1.")
            )
            active_confidence <- if (condition_mode == "ldc") {
                condition_value
            } else {
                confidence
            }
            active_boost <- if (condition_mode == "boost") {
                condition_value
            } else {
                NULL
            }
            if (!raster_data) {
                target <- data[point_id, , drop = FALSE]
            }

            .inspection_point(
                target = target,
                samples = samples,
                ref_density = density$values,
                xy_stats = xy_stats,
                xy_penalty = input$xy_penalty,
                radius_km = radius_value,
                k1 = k1_value,
                k2 = k2_value,
                bin_width = density$bin_width,
                offset = density$offset,
                confidence = active_confidence,
                boost = active_boost,
                lambda = lambda_value,
                exclude_slef = input$exclude_self,
                drop_features = drop_features,
                num_threads = num_threads,
                kernel = kernel_value,
                geographic = geographic,
                crs = map_crs
            )
        }, ignoreNULL = TRUE)

        output$condition <- shiny::renderText({
            result <- inspected()
            paste("Estimated condition:", round(result$condition, 5))
        })

        output$nearby_count <- shiny::renderText({
            result <- inspected()
            paste("Number of gray points:", nrow(result$nearby))
        })

        output$density_plot <- shiny::renderPlot({
            result <- inspected()
            .inspection_density_plot(
                density = density$values,
                selected = result$selected,
                bin_width = density$bin_width,
                offset = density$offset
            )
        }, bg = "grey96")

        # Render the base map once so it appears immediately, before any Run,
        # and is never torn down. Markers are updated in place via leafletProxy.
        output$map <- leaflet::renderLeaflet({
            .inspection_leaflet_map(
                result = NULL,
                bounds = map_bounds,
                background = background_file,
                background_colors = background_colors
            )
        })

        shiny::observeEvent(input$submit, {
            result <- tryCatch(inspected(), error = function(cond) NULL)
            proxy <- leaflet::leafletProxy("map", session = session) |>
                leaflet::clearGroup("target_preview") |>
                leaflet::clearGroup("nearby") |>
                leaflet::clearGroup("selected") |>
                leaflet::clearGroup("target")
            if (!is.null(result)) {
                proxy <- .inspection_leaflet_result_markers(proxy, result)
                zoom_bounds <- .inspection_result_bounds(result)
                if (!is.null(zoom_bounds)) {
                    leaflet::flyToBounds(
                        proxy,
                        lng1 = zoom_bounds[1],
                        lat1 = zoom_bounds[2],
                        lng2 = zoom_bounds[3],
                        lat2 = zoom_bounds[4]
                    )
                }
            }
        }, ignoreInit = TRUE)
    }

    app <- shiny::shinyApp(ui = ui, server = server)
    if (!isTRUE(launch)) {
        return(app)
    }

    invisible(shiny::runApp(app, ...))
}


.inspection_density <- function(ref_density, bin_width, offset, interpolate) {
    if (!.is_mat(ref_density)) {
        stop("'ref_density' must be a matrix or convertible to one.")
    }

    is_density <- methods::is(ref_density, "reference_density")
    if (is_density) {
        density_bin_width <- attr(ref_density, "bin.width")
        density_offset <- attr(ref_density, "offset")
        if (is.null(bin_width)) {
            bin_width <- density_bin_width
        } else if (!is.null(density_bin_width) && bin_width != density_bin_width) {
            warning("Provided 'bin_width' differs from reference density attribute.")
        }
        if (is.null(offset)) {
            offset <- density_offset
        } else if (!is.null(density_offset) && offset != density_offset) {
            warning("Provided 'offset' differs from reference density attribute.")
        }
    }

    if (length(bin_width) != 1L || !is.finite(bin_width) || bin_width <= 0) {
        stop("'bin_width' must be supplied and greater than zero.")
    }
    if (is.null(offset)) {
        offset <- 0
    }
    if (length(offset) != 1L || !is.finite(offset) || offset < 0) {
        stop("'offset' must be one non-negative number.")
    }

    values <- .check_mat(unclass(ref_density), name = "ref_density")
    if (isTRUE(interpolate)) {
        values <- terra::as.matrix(
            terra::disagg(
                terra::rast(values),
                fact = 2,
                method = "bilinear"
            ),
            wide = TRUE
        )
        bin_width <- bin_width / 2
        offset <- offset * 2
    }

    list(values = values, bin_width = bin_width, offset = as.integer(offset))
}


.inspection_has_crs <- function(x) {
    crs <- if (inherits(x, "SpatRaster")) {
        terra::crs(x)
    } else {
        x
    }

    length(crs) == 1L && !is.na(crs) && nzchar(crs)
}


.inspection_data_crs <- function(data, crs) {
    if (inherits(data, "SpatRaster") && .inspection_has_crs(data)) {
        return(terra::crs(data))
    }
    if (.inspection_has_crs(crs)) {
        return(crs)
    }
    NULL
}


.inspection_raster_samples <- function(data, samples) {
    if (ncol(samples) == 2L) {
        samples <- cbind(
            samples,
            as.matrix(terra::extract(data, samples))
        )
    } else if ((ncol(samples) - 2L) != terra::nlyr(data)) {
        stop("Sample feature count does not match number of raster layers.")
    }

    if (anyNA(samples)) {
        stop("'samples' contains missing values or locations outside the raster.")
    }

    samples
}


.inspection_raster_default_xy <- function(data) {
    # Start at the centre of the raster extent. The point need not fall on a
    # populated cell; the map renders there and the user can move it elsewhere.
    c(
        (terra::xmin(data) + terra::xmax(data)) / 2,
        (terra::ymin(data) + terra::ymax(data)) / 2
    )
}


.inspection_raster_cell_xy <- function(data, x, y) {
    .inspection_raster_cell_info(data, x, y)$xy
}


.inspection_raster_cell_info <- function(data, x, y) {
    if (length(x) != 1L || !is.finite(x) ||
            length(y) != 1L || !is.finite(y)) {
        stop("Target coordinates must be finite numbers.")
    }

    cell <- terra::cellFromXY(data, cbind(x, y))
    if (!length(cell) || is.na(cell[1])) {
        stop("Target coordinate is outside the raster extent.")
    }

    list(
        cell = cell[1],
        xy = unname(terra::xyFromCell(data, cell[1])[1, ])
    )
}


.inspection_raster_target <- function(data, x, y) {
    target <- .inspection_raster_cell_info(data, x, y)
    values <- as.matrix(terra::extract(data, target$cell))
    if (!nrow(values)) {
        stop("Target coordinate could not be extracted from the raster.")
    }
    if (anyNA(values)) {
        missing_layers <- names(data)[which(is.na(values[1, ]))]
        if (is.null(missing_layers) || any(!nzchar(missing_layers))) {
            missing_layers <- paste0("layer ", which(is.na(values[1, ])))
        }
        stop(
            "Target raster cell ", target$cell,
            " has missing values in: ",
            paste(missing_layers, collapse = ", "),
            "."
        )
    }

    cbind(
        matrix(target$xy, nrow = 1L, dimnames = list(NULL, c("x", "y"))),
        values
    )
}


.inspection_map_click_xy <- function(lng, lat, crs) {
    if (length(lng) != 1L || !is.finite(lng) ||
            length(lat) != 1L || !is.finite(lat)) {
        stop("Map click did not provide finite coordinates.")
    }
    if (!.inspection_has_crs(crs)) {
        return(c(lng, lat))
    }

    point <- terra::vect(
        data.frame(x = lng, y = lat),
        geom = c("x", "y"),
        crs = "EPSG:4326"
    )
    point <- terra::project(point, crs)
    unname(terra::crds(point)[1, ])
}


.inspection_raster_click_xy <- function(data, lng, lat, crs) {
    xy <- .inspection_map_click_xy(lng, lat, crs)
    .inspection_raster_cell_info(data, xy[1], xy[2])$xy
}


.inspection_raster_map_bounds <- function(data, crs, geographic) {
    extent <- terra::ext(data)
    corners <- data.frame(
        x = c(extent[1], extent[1], extent[2], extent[2]),
        y = c(extent[3], extent[4], extent[3], extent[4])
    )
    corners <- .inspection_map_coordinates(corners, crs, geographic)
    bounds <- c(
        lng1 = min(corners$x, na.rm = TRUE),
        lat1 = min(corners$y, na.rm = TRUE),
        lng2 = max(corners$x, na.rm = TRUE),
        lat2 = max(corners$y, na.rm = TRUE)
    )

    if (any(!is.finite(bounds))) {
        return(NULL)
    }
    if (bounds["lng1"] == bounds["lng2"]) {
        bounds[c("lng1", "lng2")] <- bounds[c("lng1", "lng2")] + c(-0.01, 0.01)
    }
    if (bounds["lat1"] == bounds["lat2"]) {
        bounds[c("lat1", "lat2")] <- bounds[c("lat1", "lat2")] + c(-0.01, 0.01)
    }

    bounds
}


# lon/lat bounding box (lng1, lat1, lng2, lat2) of the inspected target and its
# nearby/selected reference samples, used to zoom the map in on Run.
.inspection_result_bounds <- function(result) {
    xs <- c(result$target$x, result$nearby$x, result$selected$x)
    ys <- c(result$target$y, result$nearby$y, result$selected$y)
    xs <- xs[is.finite(xs)]
    ys <- ys[is.finite(ys)]
    if (!length(xs) || !length(ys)) {
        return(NULL)
    }

    bounds <- c(min(xs), min(ys), max(xs), max(ys))
    # Pad a single-point (or single-line) extent so the map does not over-zoom.
    if (bounds[1] == bounds[3]) {
        bounds[c(1, 3)] <- bounds[c(1, 3)] + c(-0.05, 0.05)
    }
    if (bounds[2] == bounds[4]) {
        bounds[c(2, 4)] <- bounds[c(2, 4)] + c(-0.05, 0.05)
    }

    bounds
}


# georaster-layer-for-leaflet renders EPSG:3857 and EPSG:4326 natively.
.inspection_is_webmercator <- function(r) {
    isTRUE(terra::same.crs(r, "EPSG:3857")) ||
        isTRUE(terra::same.crs(r, "EPSG:4326"))
}


# Normalise the 'background' argument to a single-band Cloud-Optimized GeoTIFF
# in web mercator that is served to leafem::addGeotiff() via url = (not file =).
# The url = path performs no gdal_translate/gdalwarp, so internal overviews are
# preserved and streamed. A single-band EPSG:3857/4326 file on disk is served
# as-is; anything else is reduced to its first layer, reprojected to EPSG:3857
# when needed, and written to a COG with overviews under the session tempdir.
.inspection_background_cog <- function(background) {
    if (is.null(background)) {
        return(NULL)
    }
    if (is.character(background)) {
        if (length(background) != 1L || !nzchar(background)) {
            stop("'background' must be a single file path or a terra SpatRaster.")
        }
        if (!file.exists(background)) {
            stop("'background' file does not exist: ", background)
        }
        r <- terra::rast(background)
        if (terra::nlyr(r) == 1L && .inspection_is_webmercator(r)) {
            return(background)
        }
    } else if (.is_rast(background)) {
        r <- .check_rast(background, name = "background")
    } else {
        stop("'background' must be a file path or a terra SpatRaster.")
    }
    if (terra::nlyr(r) > 1L) {
        r <- r[[1L]]
    }
    if (!.inspection_is_webmercator(r)) {
        r <- terra::project(r, "EPSG:3857")
    }
    dir <- file.path(tempdir(), "hcas_background")
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    path <- file.path(dir, "background.tif")
    terra::writeRaster(
        r, path, filetype = "COG",
        gdal = c("OVERVIEWS=AUTO", "COMPRESS=DEFLATE"), overwrite = TRUE
    )
    path
}


# Validate the optional colour ramp for the background raster.
.inspection_background_colors <- function(colors) {
    if (is.null(colors)) {
        return(NULL)
    }
    if (!is.character(colors) || !length(colors)) {
        stop(
            "'background_colors' must be a character vector of colours, ",
            "e.g. c(\"#430E59\", \"#CCCC66\", \"#184F0F\")."
        )
    }
    valid <- tryCatch({
        grDevices::col2rgb(colors)
        TRUE
    }, error = function(cond) FALSE)
    if (!valid) {
        stop("'background_colors' contains invalid colour codes.")
    }
    colors
}


.inspection_leaflet_map <- function(result = NULL, bounds = NULL,
                                    background = NULL,
                                    background_colors = NULL) {
    map <- leaflet::leaflet()
    map <- leaflet::addTiles(map)
    map <- leaflet::addProviderTiles(map, "Esri.WorldImagery")

    if (!is.null(background)) {
        palette <- if (is.null(background_colors)) palettes() else background_colors
        # Serve the COG over url = (not file =): the url = path skips leafem's
        # gdal_translate/gdalwarp, so the internal overviews survive and stream.
        # group/layerId must be set explicitly here because leafem derives them
        # from `file`, which is NULL on the url = path. autozoom = FALSE leaves
        # the view to fitBounds() below; bands = 1 renders one thematic layer.
        shiny::addResourcePath("hcas_background", dirname(background))
        map <- leafem::addGeotiff(
            map,
            url = paste0("hcas_background/", basename(background)),
            group = "background",
            layerId = "background",
            bands = 1,
            opacity = 0.8,
            autozoom = FALSE,
            colorOptions = leafem::colorOptions(
                palette = palette,
                na.color = "transparent"
            )
        )
    }

    if (!is.null(bounds)) {
        # Open at the raster extent. Strip names: leaflet serialises a named
        # scalar to a JSON object ({"lng1": 133}) instead of a bare number,
        # which breaks fitBounds.
        bounds <- unname(bounds[c("lng1", "lat1", "lng2", "lat2")])
        map <- leaflet::fitBounds(
            map,
            lng1 = bounds[1],
            lat1 = bounds[2],
            lng2 = bounds[3],
            lat2 = bounds[4]
        )
    }

    if (is.null(result)) {
        return(map)
    }

    .inspection_leaflet_result_markers(map, result)
}


.inspection_leaflet_result_markers <- function(map, result) {
    if (nrow(result$nearby)) {
        map <- leaflet::addCircleMarkers(
            map,
            data = result$nearby,
            lng = ~x,
            lat = ~y,
            color = "gray",
            radius = 2,
            group = "nearby",
            label = ~paste("Sample id:", id)
        )
    }
    if (nrow(result$selected)) {
        map <- leaflet::addCircleMarkers(
            map,
            data = result$selected,
            lng = ~x,
            lat = ~y,
            color = "blue",
            radius = 3,
            group = "selected",
            label = ~paste("Selected sample id:", id)
        )
    }
    leaflet::addCircleMarkers(
        map,
        data = result$target,
        lng = ~x,
        lat = ~y,
        color = "red",
        radius = 3,
        group = "target",
        label = "Target"
    )
}


.inspection_point <- function(
        target,
        samples,
        ref_density,
        xy_stats,
        xy_penalty,
        radius_km,
        k1,
        k2,
        bin_width,
        offset,
        confidence,
        boost,
        lambda,
        exclude_slef,
        drop_features,
        num_threads,
        kernel,
        geographic,
        crs) {

    n_vars <- .num_rs_vars_mat(samples, "samples")
    keep_features <- .keep_rs_features(drop_features, n_vars)
    samples_kept <- .subset_rs_mat(samples, keep_features)
    target_kept <- .subset_rs_mat(target, keep_features)
    kept_n_vars <- length(keep_features)
    predicted_columns <- 2L + seq_len(kept_n_vars)
    observed_columns <- predicted_columns + kept_n_vars
    bin_num <- min(dim(ref_density))

    use <- reference_use_cpp(
        target_vals = target_kept,
        sample_vals = samples_kept,
        ref_density = ref_density,
        xy_stats = xy_stats,
        xy_penalty = xy_penalty,
        geographic = geographic,
        radius_km = radius_km,
        k_env = k1,
        k_rs = k2,
        bin_width = bin_width,
        bin_num = bin_num,
        offset = offset,
        confidence = confidence,
        boost = boost,
        lambda = lambda,
        exclude_slef = exclude_slef,
        num_threads = num_threads,
        weighted_max = FALSE,
        kernel = kernel
    )
    condition <- bench_cpp(
        raster_vals = target_kept,
        sample_vals = samples_kept,
        ref_density = ref_density,
        xy_stats = xy_stats,
        xy_penalty = xy_penalty,
        geographic = geographic,
        radius_km = radius_km,
        k_env = k1,
        k_rs = k2,
        bin_width = bin_width,
        bin_num = bin_num,
        offset = offset,
        confidence = confidence,
        boost = boost,
        lambda = lambda,
        exclude_slef = exclude_slef,
        temporal_weights = NULL,
        make_su = FALSE,
        num_threads = num_threads,
        kernel = kernel
    )

    selected_ids <- which(use$density > 0)
    selected <- data.frame(
        id = selected_ids,
        x = samples_kept[selected_ids, 1],
        y = samples_kept[selected_ids, 2],
        predicted_distance = numeric(length(selected_ids)),
        observed_distance = numeric(length(selected_ids))
    )
    if (length(selected_ids)) {
        selected$predicted_distance <- rowSums(abs(
            sweep(
                samples_kept[selected_ids, predicted_columns, drop = FALSE],
                2,
                target_kept[1, predicted_columns],
                FUN = "-"
            )
        ))
        selected$observed_distance <- rowSums(abs(
            sweep(
                samples_kept[selected_ids, observed_columns, drop = FALSE],
                2,
                target_kept[1, observed_columns],
                FUN = "-"
            )
        ))
    }

    nearby_ids <- .inspection_radius_ids(
        target_xy = target[1, 1:2],
        sample_xy = samples_kept[, 1:2, drop = FALSE],
        radius_km = radius_km,
        geographic = geographic
    )
    nearby <- data.frame(
        id = nearby_ids,
        x = samples_kept[nearby_ids, 1],
        y = samples_kept[nearby_ids, 2]
    )
    target_map <- data.frame(x = target[1, 1], y = target[1, 2])

    list(
        condition = unname(condition[1, 1]),
        nearby = .inspection_map_coordinates(nearby, crs, geographic),
        selected = .inspection_map_coordinates(selected, crs, geographic),
        target = .inspection_map_coordinates(target_map, crs, geographic)
    )
}


.inspection_radius_ids <- function(target_xy, sample_xy, radius_km, geographic) {
    radius_m <- radius_km * 1000
    dx <- sample_xy[, 1] - target_xy[1]
    dy <- sample_xy[, 2] - target_xy[2]
    if (geographic) {
        dx <- dx * cos(target_xy[2] * pi / 180)
        radius <- radius_m / 111320
    } else {
        radius <- radius_m
    }
    which(dx * dx + dy * dy <= radius * radius)
}


.inspection_map_coordinates <- function(x, crs, geographic) {
    if (!nrow(x)) {
        return(x)
    }
    if (geographic && !.inspection_has_crs(crs)) {
        return(x)
    }

    source_crs <- if (.inspection_has_crs(crs)) crs else "EPSG:4326"
    points <- terra::vect(
        x[, c("x", "y"), drop = FALSE],
        geom = c("x", "y"),
        crs = source_crs
    )
    if (!terra::is.lonlat(points, perhaps = TRUE, warn = FALSE)) {
        points <- terra::project(points, "EPSG:4326")
    }
    x[, c("x", "y")] <- terra::crds(points)
    x
}


.inspection_density_plot <- function(density, selected, bin_width, offset) {
    colours <- c(
        "gray92", "#F5F2D8", "#C6E8BC", "#7ED5B8",
        "#34B8C0", "#478EC1", "#7A55AB", "#80146E"
    )
    density_x <- (seq_len(nrow(density)) - 1L + offset) * bin_width
    density_y <- (seq_len(ncol(density)) - 1L + offset) * bin_width
    density_plot <- expand.grid(x = density_x, y = density_y)
    density_plot$value <- as.vector(density)
    maximum <- if (nrow(selected)) {
        max(selected[, c("predicted_distance", "observed_distance")]) + 0.5
    } else {
        0.5
    }

    x <- y <- value <- predicted_distance <- observed_distance <- NULL
    ggplot2::ggplot(
        data = selected,
        ggplot2::aes(x = predicted_distance, y = observed_distance)
    ) +
        ggplot2::geom_tile(
            data = density_plot,
            ggplot2::aes(x = x, y = y, fill = value),
            alpha = 0.7,
            inherit.aes = FALSE
        ) +
        ggplot2::geom_point(alpha = 0.5, colour = "blue", size = 3) +
        ggplot2::scale_x_continuous(limits = c(-0.1, maximum)) +
        ggplot2::scale_y_continuous(limits = c(-0.1, maximum)) +
        ggplot2::geom_abline(intercept = 0, slope = 1, alpha = 0.5, linetype = 3) +
        ggplot2::scale_fill_gradientn(colours = colours) +
        ggplot2::theme(
            plot.background = ggplot2::element_blank(),
            panel.border = ggplot2::element_blank(),
            legend.background = ggplot2::element_blank()
        ) +
        ggplot2::coord_equal() +
        ggplot2::labs(
            x = "Predicted distance",
            y = "Observed distance",
            fill = "Condition"
        )
}
