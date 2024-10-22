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
#' @param cookie_cutter An optional \code{sf} spatial data frame of polygonal
#'   features that define the outer-most supposed extent of the fire. If
#'   provided, progression polygons will be clipped and or removed so that none
#'   extend outside the cookie features.
#'
#' @param out_epsg Integer EPSG code specifying the map projection for the
#'   output progression polygons. This \emph{must} be a projected coordinate
#'   system with metres as map units. The default is 8058 (NSW Lambert /
#'   GDA2020).
#'
#' @param min_geom_area Threshold minimum area (square metres) for progression
#'   polygons. If the progression for a given time step consists of two or more
#'   spatially distinct parts, the threshold is applied to each individually.
#'
#' @param unique_geom (logical) If \code{TRUE} (the default), any group of
#'   records with identical geometries will be reduced to a single record, being
#'   that with the earliest time in the group. If \code{FALSE} this step will be
#'   skipped.
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
                              cookie_cutter = NULL,
                              out_epsg = 8058,
                              min_geom_area = 100,
                              unique_geom = TRUE,
                              replicate_times = c('merge', 'largest', 'smallest', 'fail'),
                              dTolerance = 2) {

  checkmate::assert_class(x, "sf")

  if (is.na( sf::st_crs(x) )) {
    stop("The input sf data frame must have a coordinate reference system defined")
  }

  # Attempt to fix any invalid geometries
  x <- sf::st_make_valid(x)

  # Ensure that the name of the geometry column is 'geom' to make the code for
  # later steps a bit less fiddly
  sf::st_geometry(x) <- "geom"

  # Check that the output CRS is defined and has metres as map units
  checkmate::assert_integerish(out_epsg, any.missing = FALSE, len = 1)

  units_txt <- units::deparse_unit(sf::st_crs(out_epsg)$ud_unit)
  if (!units_txt == "m") {
    msg <- glue::glue("Argument out_epsg ({out_epsg}) does not have metres as map units")
    stop(msg)
  }

  # Transform the input extent polygons into the output CRS if required
  x <- sf::st_transform(x, out_epsg)

  # Check distance tolerance used to optionally generalize polygon vertices
  checkmate::assert_number(dTolerance, lower = 0)

  # Check the column(s) specified for time
  checkmate::assert_character(time_cols, min.len = 1, any.missing = FALSE)
  ok <- time_cols %in% colnames(x)
  if (!all(ok)) {
    msg <- paste(time_cols[!ok], collapse = ", ")
    msg <- glue::glue("Missing time column(s): {msg}")
    stop(msg)
  }

  checkmate::assert_number(min_geom_area, lower = 0)

  checkmate::assert_flag(unique_geom)

  replicate_times <- match.arg(replicate_times)

  # Remove any extent records with non-polygonal geometries
  n_in <- nrow(x)
  x <- .remove_non_polygonal(x)
  n_poly <- nrow(x)

  if (n_poly == 0) {
    stop("There are no valid polygonal features in the input extent data")
  } else if (n_poly < n_in) {
    msg <- glue::glue("Discarding {n_in - n_poly} input record(s) for features with empty or non-polygonal geometries")
    warning(msg, immediate. = TRUE)
  }

  # Check cookie_cutter polygons if provided
  if (!is.null(cookie_cutter)) {
    # Check that the cookie cutter object is an sf data frame
    checkmate::assert_class(cookie_cutter, "sf")

    # Check it has a CRS defined
    if (is.na( sf::st_crs(cookie_cutter) )) {
      stop("The sf data frame for cookie_cutter must have a coordinate reference system defined")
    }

    # Re-project if required
    cookie_cutter <- sf::st_transform(cookie_cutter, out_epsg)

    # Remove any non-polygonal features
    cookie_cutter <- .remove_non_polygonal(cookie_cutter)
    if (nrow(cookie_cutter) == 0) stop("No polygonal features in the sf data frame set for argument cookie_cutter")

    # Subset cookie cutter polygons to those that intersect the bounding box of the extent polygons
    bb_x <- sf::st_bbox(x)
    cookie_cutter <- sf::st_crop(cookie_cutter, bb_x)

    if (nrow(cookie_cutter) == 0) {
      # There is no overlap between the cookie cutter and the bounding box of the fire extent
      # polygons, so there's nothing to do.
      return(NULL)
    }
  }

  # If requested, subset the input data to unique geometries by taking the
  # earliest extent polygon from each group of spatially identical polygons
  if (unique_geom) {
    x <- x %>%
      sf::st_make_valid() %>%
      dplyr::group_by(geom) %>%
      dplyr::mutate(.time_order = order(.data[[time_cols]])) %>%
      dplyr::ungroup() %>%
      dplyr::filter(.time_order == 1) %>%
      dplyr::select(-.time_order)
  }

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
      x <- x %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(time_cols))) %>%
        dplyr::summarize(geom = sf::st_union(geom))

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


  # Generalize the feature geometries if requested
  if (dTolerance > 0) {
    x <- sf::st_simplify(x, dTolerance = dTolerance)
    if (!is.null(cookie_cutter)) cookie_cutter <- sf::st_simplify(cookie_cutter, dTolerance = dTolerance)
  }


  # Sort data by time
  x <- dplyr::arrange(x, dplyr::across(dplyr::all_of(time_cols)))

  # Process the extent polygons ordered by time and derive the progression
  # polygons
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

          # Cookie cutter the new feature
          gnew <- gnew

          # Guard against the new feature not being a valid polygon

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

  # Apply cookie cutter if defined
  if (!is.null(cookie_cutter)) {
    # Do fast st_intersects check first to save time
    ii <- which( lengths( sf::st_intersects(dat_prog, cookie_cutter) ) > 0 )

    if (length(ii) > 0) {
      dat_prog <- suppressWarnings({
        # Note: just intersect with the cookie geometry to avoid picking up
        # unwanted columns if `cookie_cutter` is an sf data frame
        dat_prog <- sf::st_intersection(dat_prog, sf::st_geometry(cookie_cutter[ii, ]) )

        if (nrow(dat_prog) > 0) {
          .remove_non_polygonal(dat_prog)
        } else {
          # There were no polygonal geometries :(
          NULL
        }
      })

    } else {
      # All progression polygons are outside the cookie cutter polygons
      # so discard them by returning NULL
      dat_prog <- NULL
    }
  }

  dat_prog
}


# Non-exported helper function to remove features from an sf data frame
# with non-polygonal geometries
.remove_non_polygonal <- function(x) {
  ok <- !sf::st_is_empty(x) & sf::st_geometry_type(x) %in% c("POLYGON", "MULTIPOLYGON")

  # Return valid polygon records / geometries
  if (inherits(x, "sf")) x[ok, ]
  else if (inherits(x, "sfc")) x[ok]
  else stop("Input should be an 'sf' or 'sfc' object")
}

