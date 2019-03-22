# horizontal nudge position adjustment
# copied from https://github.com/tidyverse/ggplot2/issues/2733
position_hnudge <- function(x = 0) {
  ggproto(NULL, PositionHNudge, x = x)
}

PositionHNudge <- ggproto("PositionHNudge", Position,
                          x = 0,
                          required_aes = "x",
                          setup_params = function(self, data) {
                            list(x = self$x)
                          },
                          compute_layer = function(data, params, panel) {
                            transform_position(data, function(x) x + params$x)
                          }
)