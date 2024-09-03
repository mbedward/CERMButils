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
#'   the time steps for the progressions. Normally a time column will will hold
#'   a Date or Date-Time value (e.g. POSIX.ct) but any scalar variable type that
#'   can be sorted is acceptable. When two or more time columns are specified,
#'   sorting will be based on the provided order of column names.
#'
#' @param out_epsg Integer EPSG code specifying the map projection for the
#'   output progression polygons. Default is 8058 (NSW Lambert / GDA2020).
#'
#' @param min_geom_area Threshold minimum area (square metres) for progression
#'   polygons. If the progression for a given time step consists of two or more
#'   spatially distinct parts, the threshold is applied to each individually.
#'
#' @param replicate_times A single string that specifies what to do if there is
#'   more than one input data record for any time step. Options are:
#'   \code{'merge'} (merge feature geometries within each time step);
#'   \code{'largest'} (choose the record with the largest extent polygon area);
#'   \code{'smallest'} (choose the record with the largest extent polygon area);
#'   \code{'fail'} (stop with an error message).
#'   The default is \code{'merge'}. The argument may be abbreviated and is
#'   case-sensitive.
#'
#' @param dTolerance Distance tolerance used to simplify the input extent
#'   polygons using the \code{\link[sf]{st_simplify}} function. Default is 2 for
#'   a minimum between vertex distance of 2 metres. Set to 0 to omit the
#'   simplification step.
#'
#' @return An \code{sf} spatial data frame of progression polygons with
#'   time column(s) copied from the input fire extent data.
#'
#' @export
#'
make_progressions <- function(x,
                              time_cols,
                              out_epsg = 8058,
                              min_geom_area = 100,
                              replicate_times = c('merge', 'largest', 'smallest', 'fail'),
                              dTolerance = 2) {

  checkmate::assert_class(x, "sf")

  the_crs <- sf::st_crs(x)
  if (is.na(the_crs)) {
    stop("The input sf data frame must have a coordinate reference system defined")
  }

  checkmate::assert_integerish(out_epsg, any.missing = FALSE, len = 1)
  x <- sf::st_transform(x, out_epsg)

  checkmate::assert_number(dTolerance, lower = 0)
  if (dTolerance > 0) x <- sf::st_simplify(x, dTolerance = dTolerance)

  checkmate::assert_character(time_cols, min.len = 1, any.missing = FALSE)
  ok <- time_cols %in% colnames(x)
  if (!all(ok)) {
    msg <- paste(time_cols[!ok], collapse = ", ")
    msg <- glue::glue("Missing time column(s): {msg}")
    stop(msg)
  }

  checkmate::assert_number(min_geom_area, lower = 0)

  replicate_times <- match.arg(replicate_times)

  # Subset the input data to records with geometrically distinct features
  # gunique <- unique(sf::st_geometry(x)) %>%
  #   sf::st_as_sfc()
  #
  # sf::st_crs(gunique) <- sf::st_crs(x)
  # ii <- match(gunique, sf::st_geometry(x))
  # x <- x[ii, ]

  # Check whether times are unique
  xtimes <- apply(sf::st_drop_geometry(x[, time_cols]), 1, paste, collapse=" ")
  ntimes <- table(xtimes)

  if (any(ntimes > 1)) {
    if (replicate_times == "fail") {
      msg <- paste(time_cols, collapse = ", ")
      msg <- glue::glue("One or more times defined by {msg} occur more than once
                         and replicate_times was set to 'fail'")
      stop(msg)

    } else if (replicate_times == "merge") {
      # Make sure that the geometry column is 'geom' so the dplyr::summarize
      # call is easier to code
      gname <- attr(x, "sf_column", exact = TRUE)
      i <- which(colnames(x) == gname)
      colnames(x)[i] <- "geom"

      x <- x %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(time_cols))) %>%
        dplyr::summarize(geom = st_union(geom))

    } else {  # 'largest' or 'smallest'
      # Choose one record per time based on the largest or smallest feature area
      a <- as.numeric( sf::st_area(x) )

      x <- x %>%
        dplyr::mutate(.area_ = a) %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(time_cols)))

      if (replicate_times == "largest") {
        x <- dplyr::filter(x, .area_ == max(.area_))
      } else {
        x <- dplyr::filter(x, .area_ == min(.area_))
      }

      # If there are identical feature geometries in any time step there
      # will still be multiple records per time, so just choose the first
      x <- dplyr::slice_head(x, n = 1)

      x <- x %>%
        dplyr::ungroup() %>%
        dplyr::select(-.area_)
    }
  }

  # Sort data by time
  x <- dplyr::arrange(x, dplyr::across(dplyr::all_of(time_cols)))

  # From the second time step onwards, identify the fire progression (if any) and
  # update the cumulative extent.
  #
  gprog <- lapply(1:nrow(x), function(i) {
    if (i == 1) {
      # Initialize the outer fire extent and return it as the first progression
      # polygon
      gouter <<- sf::st_geometry(x[i,])
      res <- gouter

    } else {
      # For subsequent extent polygons, derive the difference between the
      # previous outer polygon and the current extent
      #
      res <- NULL
      g <- sf::st_geometry(x[i,])

      # Get the spatial difference (new extent, if any) and format as
      # an `sfc` geometry list object with the parent CRS applied. This
      # ensures that st_area will provide values in m^2 (below).
      gnew <- sf::st_difference(g, gouter) %>%
        sf::st_sfc(., crs = sf::st_crs(x)) %>%
        sf::st_make_valid()

      # Discard any non-polygonal geometries and any polygons with area
      # less than the minimum threshold.
      is_poly <- sf::st_is(gnew, c("POLYGON", "MULTIPOLYGON"))
      if (any(is_poly)) {
        gnew <- gnew[is_poly]

        # Convert any multipolygons to polygons
        gnew <- sf::st_cast(gnew, "POLYGON")

        # Check area
        a <- as.numeric( sf::st_area(gnew) )
        ok <- a >= min_geom_area
        if (any(ok)) {
          gnew <- gnew[ok]

          # Update outer extent so far
          gouter <<- sf::st_union(gouter, gnew) %>%
            sf::st_union() %>%  # second call to dissolve internal boundaries
            sf::st_make_valid()

          res <- sf::st_cast(gnew, "POLYGON")
        }
      }
    }

    res
  })

  # Combine the progression polygons and the relevant input attribute rows
  ilen <- lengths(gprog)
  indices <- rep(seq_along(ilen), ilen) + 1

  dat_prog <- sf::st_drop_geometry( x[indices, time_cols] )
  gprog <- do.call(c, gprog)
  dat_prog$geom <- gprog
  dat_prog <- sf::st_as_sf(dat_prog)

  dat_prog
}

