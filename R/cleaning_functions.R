#' Identify groups of similar line features
#'
#' Layers of line features representing back-burning sometimes contain groups of
#' lines that appear to be duplicates, i.e. having identical or generally
#' similar geometries. This function attempts to identify such groups. Two or
#' more line features are considered to be members of a group if the closest
#' distance between them is less than a threshold distance, specified via the
#' (\code{max_nbr_dist}) argument.
#'
#' \emph{Beware}: this function is not yet smart enough to distinguish between a
#' tightly clustered set of lines on the one hand, and a more dispersed set in
#' which pairs of lines are closer to each other than the specified threshold
#' distance. For example, lines forming a ladder arrangement over a substantial
#' distance, where each member is sufficiently close to at least one other,
#' would all be assigned to a single group.
#'
#' @param lines An \code{'sf'} spatial data frame of input line features. This
#'   must have a projected (i.e. non-geographic) coordinate reference system
#'   assigned. All features should have single-part geometries; either type
#'   \code{'LINESTRING'} or type \code{'MULTILINESTRING'} with one part.
#'
#' @param max_nbr_dist (numeric; default 50) Threshold distance for grouped
#'   lines. Two lines will be considered as potential members of a single group
#'   if the closest distance between them is equal to or less than this
#'   threshold distance.
#'
#' @param use_attributes A character vector of one or more attribute names
#'   corresponding to columns in the input \code{'sf'} data frame. If provided,
#'   these attributes will be considered when grouping lines features in
#'   addition to location and topological similarity. See example below.
#'
#' @return An integer vector of group ID values with length equal to the number
#'   of rows in the input \code{`sf`} data frame.
#'
#' @examples
#' \dontrun{
#' # Default behaviour: assign group IDs to features based on similarity
#' # of location
#' ids <- assign_lines_to_groups(dat_firelines)
#'
#' # Group on both location and the attributes 'IncidentName' and 'FIRELINETYPE'
#' ids <- assign_lines_to_groups(dat_firelines,
#'                               use_attributes = c("IncidentName", "FIRELINETYPE"))
#' }
#'
#' @export
#
assign_lines_to_groups <- function(lines,
                                   max_nbr_dist = 50,
                                   use_attributes = NULL) {

  CRS <- .check_sf_and_metres(lines, "lines")

  # Ensure that the name of the geometry column is 'geom' to make the code for
  # later steps a bit less fiddly
  sf::st_geometry(lines) <- "geom"

  checkmate::assert_number(max_nbr_dist, finite = TRUE, lower = 0)

  if (!is.null(use_attributes)) {
    checkmate::assert_character(use_attributes,
                                min.chars = 1,
                                any.missing = FALSE,
                                min.len = 1,
                                unique = TRUE)

    # Make sure that 'geom' is not in the 'use_attributes' vector
    use_attributes <- setdiff(use_attributes, 'geom')

    ok <- use_attributes %in% colnames(lines)

    if (!all(ok)) {
      msg <- paste(use_attributes[!ok], collapse = ", ")
      msg <- glue::glue("Column(s) specified in 'use_attributes' not found: {msg}")
      stop(msg)
    }
  }


  # Convert all features to LINESTRINGs and check that the number of features
  # does not change, i.e. there were no multi-part lines. This step will also
  # raise an error if any of the input features are not linear (e.g. points).
  #
  nfeatures <- nrow(lines)

  lines <- tryCatch(.geom_to_linestring(lines),
                    error = function(e) {
                      msg <- glue::glue("Cannot treat all input features as 'LINESTRING'
                                         {e}")
                      stop(msg)
                    })

  if (nrow(lines) != nfeatures) {
    stop("Not all input line features have a single-part geometry")
  }

  # Create a buffer around each line
  bufs <- lines  %>%
    st_geometry() %>%
    sf::st_buffer(dist = max_nbr_dist / 2,  # Note: buffer by half the threshold distance
                  endCapStyle = 'FLAT')

  # Record intersections between buffers
  xbufs <- sf::st_intersects(bufs, sparse = TRUE)

  # Identify initial groups based on intersections between line buffers
  ids <- xbufs %>%
    igraph::graph_from_adj_list() %>%
    igraph::components() %>%
    `$`(membership)

  # If grouping columns were specified via `use_attributes`, split the initial
  # groups as required
  if (!is.null(use_attributes)) {
    dat_grps <- sf::st_drop_geometry(lines) %>%
      dplyr::select(dplyr::all_of(use_attributes))

    dat_grps$.initial_id <- ids

    dat_grps <- dat_grps %>%
      dplyr::group_by(dplyr::all_of(use_attributes), .initial_id) %>%
      dplyr::mutate(.final_id = dplyr::cur_group_id())

    ids <- dat_grps[['.final_id']]
  }

  # Return group IDs
  ids
}


#' Derive an average line for a group of line features
#'
#' The purpose of this function is to reduce groups of similar line features to
#' a single average or representative line, where \emph{similar} means that all
#' lines in a group are closely located and have approximately the same
#' orientation along their length. The function takes an \code{sf} data frame of
#' line features and a vector of integer group ID values to identify which group
#' each line feature belongs to. A buffer polygon is placed around the unioned
#' lines in each group. An average line is then generated by taking the
#' centreline of the buffer polygon.
#'
#' Group IDs will generally have been created using the
#' \code{assign_lines_to_groups()} function which presently only requires that
#' lines be with a certain distance of each other at some point. Hence, it is
#' possible for a group to contain lines that would not have been grouped
#' together by manual inspection, e.g. lines meeting at a T-junction. To help
#' identify such cases, the function calculates the furthest distance between a
#' point on any of the group lines and the derived average line and records this
#' in the \code{maxdist} column of the output \code{sf} data frame. The larger
#' this value, the higher the chance that the average line does not reliably
#' represent all lines in the group.
#'
#' The value of the \code{buffer_dist} argument can greatly influence the result
#' of this function. In general, it seems to be best to use a fairly small value
#' for the buffer distance, relative to the average length of input line
#' features, to make it more likely that the buffer polygons placed around each
#' group of lines will be elongated. Polygons that are more circular or square
#' will usually result in representative lines that are unacceptable. Hopefully,
#' most such cases will be flagged for checking in the output data frame, but
#' beware of ones that slip through.
#'
#' @param lines An \code{'sf'} spatial data frame of input line features. This
#'   must have a projected (i.e. non-geographic) coordinate reference system
#'   assigned. All features should have single-part geometries; either type
#'   \code{'LINESTRING'} or type \code{'MULTILINESTRING'} with one part.
#'
#' @param group_ids An integer vector of group ID values with length equal to
#'   the number of rows in the input \code{`sf`} data frame. The ID values do
#'   not have to be consecutive.
#'
#' @param buffer_dist The width of the buffer to place around each group of line
#'   features to form the polygon from which a representative line will be
#'   generated using the centreline algorithm.
#'
#' @return An \code{sf} data frame of \code{N} line features, where \code{N} is
#'   the number of distinct input group ID values.
#'   \describe{
#'   \item{id}{Integer group ID value}
#'   \item{maxdist}{The furthest distance between a member line in the group and
#'     the average line. This will always be zero for groups with only line.
#'     Large distance values indicate groups where the average line is possibly
#'     not representative of all lines in the group.
#'     }
#'   \item{geom}{The average line geometry}
#'   }
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' library(sf)
#' library(CERMButils)
#'
#' ids <- assign_lines_to_groups(dat_firelines)
#' dat_group_lines <- get_group_line(dat_firelines, ids, buffer_dist = 25)
#'
#' # Optionally, sort by the 'maxdist' column to identify any groups where the
#' # average line does not usefully represent the member lines.
#'
#' dat_group_lines <- dat_group_lines %>%
#'   arrange(desc(maxdist))
#' }
#'
#' @export
#
get_group_line <- function(lines, group_ids = NULL, buffer_dist = 25, progress = TRUE) {
  CRS <- .check_sf_and_metres(lines, "lines")

  nfeatures <- nrow(lines)

  # Check that all features are lines (LINESTRING or MULTILINESTRING)
  gt <- sf::st_geometry_type(lines)
  ok <- grepl(gt, pattern = "LINESTRING")
  if (!all(ok)) {
    stop("One or more input features do not have LINESTRING geometry")
  }

  # Reduce any MULTILINESTRING features to LINESTRING, checking that they all
  # were single-part geometries
  if (any(grepl(gt, pattern = "MULTILINESTRING"))) {
    lines <- .geom_to_linestring(lines)
    if (nrow(lines) != nfeatures) {
      stop("Input features include one or more MULTILINESTRINGs with more than one part")
    }
  }

  if (is.null(group_ids)) {
    # Assume that we are looking for a representative line for all features in
    # the input layer
    group_ids <- rep(1, nrow(lines))

  } else {
    # With user-provided IDs there should be as many values as features
    checkmate::assert_integerish(group_ids, any.missing = FALSE, len = nrow(lines))
  }

  # Process each group of lines
  UniqueIDs <- sort(unique(group_ids))

  if (progress) pb <- progress::progress_bar$new(total = length(UniqueIDs))
  pb$tick(0)

  average_geoms <- lapply(UniqueIDs, function(id) {
    dat <- lines[which(group_ids == id), ]
    glines <- sf::st_geometry(dat)

    res <- NULL

    ZeroDistance <- units::set_units(0.0, "m")

    if (nrow(dat) == 1) {
      # Nothing more to do...
      res <- sf::st_sf(id = id, maxdist = ZeroDistance, geom = glines)

    } else {
      ulines <- sf::st_union(glines)
      bb <- sf::st_minimum_rotated_rectangle(ulines)

      buf_poly <- sf::st_buffer(ulines, dist = buffer_dist, endCapStyle = 'ROUND')

      c_line <- centerline::cnt_path_guess(buf_poly, keep = 5)

      # Return the portion of the inferred centre line within the bounding
      # rectangle
      c_line <- sf::st_intersection(c_line, bb)

      # Check if the derived line has more than one part. This can happen if the
      # current group of lines includes one or more members that are more than 2
      # x buffer_dist apart. At the moment we don't have a good way of dealing
      # with this so we just reject the group and issue a warning message
      #
      if (length(c_line) > 1) {
        msg <- glue::glue("Unable to derive a single line for group id {id}")

        if (progress) {
          pb$message(msg)
        } else {
          warning(msg, immediate. = TRUE)
        }

      } else {
        # Calculate Hausdorff distance between the average line and the union of
        # the group lines to find the worst-case distance as a measure of how
        # well the average line represents the group
        maxdist = sf::st_distance(ulines, c_line, which = "Hausdorff")

        res <- sf::st_sf(id = id, maxdist = maxdist, geom = c_line)
      }
    }

    if (progress) pb$tick()

    res
  })

  # Combine average lines into a single 'sf' data frame and return
  do.call(rbind, average_geoms)
}


# Private helper function to convert geometries in an sf-ish object to
# type LINESTRING.
#
.geom_to_linestring <- function(x) {
  # Suppress any sf warnings about repeated attributes
  suppressWarnings(
    sf::st_cast(x, "LINESTRING")
  )
}

