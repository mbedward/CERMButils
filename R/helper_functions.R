# Private helper function to use if-else blocks in a dplyr or magrittr pipeline.
#
# Black magic code slightly adapted from a StackOverflow post by Johann-Friedrich Salzmann
# https://stackoverflow.com/a/78368837/40246
#
.pipe_ifelse <- function(data, cond, a, b){
  ce <- rlang::enexpr(cond)

  if(rlang::eval_tidy(ce, data = data)) {
    e <- rlang::enexpr(a)
  } else {
    if(missing(b)) {
      return(data)
    } else {
      e <- rlang::enexpr(b)
    }
  }

  u <- rlang::expr(`%>%`((.), rlang::`!!`(e)))
  data %>% {rlang::eval_tidy(u)}
}


# Private helper function to check that 'sf' values passed to functions
# are actually class 'sf' and have a CRS defined with metres as units.
#
# Returns the CRS of the input sf object
#
.check_sf_and_metres <- function(x, .var.name) {
  checkmate::assert_class(x, "sf", .var.name = .var.name)

  # Check that the input layer has a CRS defined with metres as map units
  CRS <- sf::st_crs(x)
  if (is.na(CRS)) stop("A cooordinate reference system must be set for ", .var.name)

  units_txt <- units::deparse_unit(sf::st_crs(CRS)$ud_unit)
  if (units_txt != "m") {
    msg <- glue::glue("The CRS of {.var.name} does not have metres as map units")
    stop(msg)
  }

  CRS
}
