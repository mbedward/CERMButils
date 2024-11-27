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
