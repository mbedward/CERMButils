#' Derive a set of progression polygons
#'
#' This function takes a spatial data frame of fire extent polygons with one or
#' more attributes defining time points, and derives the corresponding sequence
#' of fire progression polygons.
#'
#' @param x An \code{sf} spatial data frame with one or more columns defining time
#'   and polygonal features representing the total extent of a fire at each time.
#'
#' @param time_cols A character vector of one or more column names that define
#'   the time steps for the progressions.
#'
#' @param min_geom_area Threshold minimum area (square metres) for progression
#'   polygons. If the progression for a given time step consists of two or more
#'   spatially distinct parts, the threshold is applied to each individually.
#'
#' @return An \code{sf} spatial data frame of progression polygons with
#'   attribute values copied from the input data.
#'
#' @export
#'
make_progressions <- function(x, time_cols, min_geom_area = 100) {
  checkmate::assert_class(x, "sf")

  the_crs <- sf::st_crs(x)
  if (is.na(the_crs)) {
    stop("The input sf data frame must have a coordinate reference system defined")
  }

  checkmate::assert_character(time_cols, min.len = 1, any.missing = FALSE)
  ok <- time_cols %in% colnames(x)
  if (!all(ok)) {
    msg <- paste(time_cols[!ok], collapse = ", ")
    msg <- glue::glue("Missing time column(s): {msg}")
    stop(msg)
  }

  # Check that times are unique
  xtimes <- apply(x[, time_cols], 1, paste, collapse=" ")
  ntimes <- table(xtimes)

  if (any(ntimes) > 1) {
    msg <- paste(time_cols, collapse = ", ")
    msg <- glue::glue("One or more times defined by {msg} occur more than once")
  }

  # Sort data by time
  x <- dplyr::arrange(x, dplyr::across(dplyr::all_of(time_cols)))

  # Subset the input data to records with geometrically distinct features
  gunique <- unique(sf::st_geometry(x)) %>%
    sf::st_as_sfc()

  sf::st_crs(gunique) <- sf::st_crs(x)
  ii <- match(gunique, sf::st_geometry(x))
  x <- x[ii, ]

  # Initialize the outer fire extent
  gouter <- sf::st_geometry(x[1,])

  # From the second time step onwards, identify the fire progression (if any) and
  # update the cumulative extent.
  #
  gprog <- lapply(2:nrow(x), function(i) {
    res <- NULL
    g <- sf::st_geometry(x[i,])

    # Get the spatial difference (new extent, if any) and format as
    # an `sfc` geometry list object with the parent CRS applied. This
    # ensures that st_area will provide values in m^2 (below).
    gnew <- sf::st_difference(g, gouter) %>%
      sf::st_sfc(., crs = sf::st_crs(x))

    # Update outer extent so far
    gouter <<- sf::st_union(gouter, g)

    if (length(gnew) > 0) {
      # check for any non-polygonal geometries or tiny slivers
      gnew <- sf::st_cast(gnew, "POLYGON")
      a <- as.numeric( sf::st_area(gnew) )
      res <- gnew[a >= min_geom_area]
    }

    res
  })

  # Combine the progression polygons and the relevant input attribute rows
  ilen <- lengths(gprog)
  indices <- rep(seq_along(ilen), ilen) + 1

  dat_prog <- sf::st_drop_geometry( x[indices, ] )
  gprog <- do.call(c, gprog)
  dat_prog$geom <- gprog
  dat_prog <- sf::st_as_sf(dat_prog)

  dat_prog
}

