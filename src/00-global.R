# Global Functions and Constants
#
# Jonas Schöley

# Init ------------------------------------------------------------

# fonts
library(showtext)
font_add_google('Roboto', 'roboto')
font_add_google('Roboto Condensed', 'robotocondensed')
showtext_auto()

# figure specifications
figspec <- list()

# Figure dimensions -----------------------------------------------

# figure dimensions in mm
figspec$fig_dims <- list(width = 124.6)

# Figure theme ----------------------------------------------------

# ggplot theme by Jonas Schöley
figspec$MyGGplotTheme <-
  function (
    size = 8,
    family = 'roboto',
    scaler = 1,
    axis = 'x',
    panel_border = FALSE,
    grid = 'y',
    minor_grid = '',
    show_legend = TRUE,
    ar = NA,
    axis_title_just = 'rt',
    axis_ticks = TRUE
  ) {

    size_med = size*scaler
    size_sml = round(size*0.7)*scaler
    base_linesize = 0.3*scaler

    # justification of axis titles
    xj <- switch(tolower(substr(axis_title_just, 1, 1)), b = 0,
                 l = 0, m = 0.5, c = 0.5, r = 1, t = 1)
    yj <- switch(tolower(substr(axis_title_just, 2, 2)), b = 0,
                 l = 0, m = 0.5, c = 0.5, r = 1, t = 1)

    list(
      theme_minimal(base_size = size_med, base_family = family),
      theme(
        # basic
        text = element_text(color = 'black'),
        line = element_line(size = base_linesize, lineend = 'square'),
        # axis
        axis.title = element_text(size = size_med, face = 'bold'),
        axis.title.x = element_text(hjust = xj),
        axis.title.y = element_text(hjust = yj),
        axis.title.y.right = element_text(hjust = yj, angle = 90),
        axis.text = element_text(size = size_med, color = 'black'),
        # strips
        strip.text = element_text(color = 'black', size = size_med),
        strip.background = element_blank(),
        # plot
        title = element_text(face = 'bold'),
        plot.subtitle = element_text(color = 'black', size = size_med, face = 'bold'),
        plot.caption = element_text(color = 'black', size = size_sml, face = 'plain'),
        plot.background = element_blank(),
        panel.background = element_blank(),
        #plot.margin = unit(c(1, 0.1, 0.5, 0.5), units = 'mm'),
        # grid
        panel.grid = element_blank()
      ),
      if (isTRUE(axis_ticks)) {
        theme(axis.ticks = element_line(size = rel(0.5), color = 'black'))
      },
      if (identical(grid, 'y')) {
        theme(panel.grid.major.y =
                element_line(size = base_linesize, linetype = 3, color = 'grey80'))
      },
      if (identical(grid, 'x')) {
        theme(panel.grid.major.x =
                element_line(size = base_linesize, linetype = 3, color = 'grey80'))
      },
      if (identical(grid, 'xy') | identical(grid, 'yx')) {
        theme(panel.grid.major.y =
                element_line(size = base_linesize, linetype = 3, color = 'grey80'),
              panel.grid.major.x =
                element_line(size = base_linesize, linetype = 3, color = 'grey80'))
      },
      if (identical(minor_grid, 'y')) {
        theme(panel.grid.minor.y =
                element_line(size = base_linesize, linetype = 3, color = 'grey80'))
      },
      if (identical(minor_grid, 'x')) {
        theme(panel.grid.minor.x =
                element_line(size = base_linesize, linetype = 3, color = 'grey80'))
      },
      if (identical(minor_grid, 'xy') | identical(grid, 'yx')) {
        theme(panel.grid.minor.y =
                element_line(size = base_linesize, linetype = 3, color = 'grey80'),
              panel.grid.minor.x =
                element_line(size = base_linesize, linetype = 3, color = 'grey80'))
      },
      if (isTRUE(panel_border)) {
        theme(
          panel.border =
            element_rect(fill = NA)
        )
      },
      if (!isTRUE(show_legend)) {
        theme(legend.position = 'none')
      },
      if (axis == 'x') {
        theme(
          axis.line.x = element_line(linetype = 1, color = 'black')
        )
      },
      if (axis == 'y') {
        theme(
          axis.line.y = element_line(linetype = 1, color = 'black')
        )
      },
      if (axis == 'xy') {
        theme(
          axis.line = element_line(linetype = 1, color = 'black')
        )
      },
      if (!is.na(ar)) {
        theme(
          aspect.ratio = ar
        )
      }
    )
  }

# Figure breaks and labels ----------------------------------------

figspec$bl <- list()
figspec$bl <- within(figspec$bl, {

  # nice log-labels for plot
  NiceLogs <- function (base = 10, min = -6, max = 0,
                        between = c(2, 4, 6, 8)) {

    powers <- min:max
    primary_labels <- paste0('10^', powers)
    secondary_labels <- paste0('scriptstyle(', between, ')')

    breaks <-
      c(vapply(
        base^powers,
        FUN = function (x) x*c(1, between),
        FUN.VALUE = vector('numeric', length = length(between)+1)
      ))
    labels <-
      c(vapply(
        primary_labels,
        FUN = function (x) c(x, secondary_labels),
        FUN.VALUE = vector('character', length = length(between)+1)
      ))

    list(breaks = breaks, labels = parse(text = labels))

  }

  # days since birth
  age = c(0, 7, 14, 21, 28)

  # hazard rate breaks and labels for plot axes
  hazard1 = 10^(seq(-6, 0, 1))
  names(hazard1) <-
    parse(
      text = c('10^-6', 'scriptstyle(-5)', 'scriptstyle(-4)',
               'scriptstyle(-3)', 'scriptstyle(-2)',
               'scriptstyle(-1)', '10^0')
    )
  hazard2 = NiceLogs(10, -7, 0)$breaks
  names(hazard2) = NiceLogs(10, -7, 0)$labels

})

# Export figures --------------------------------------------------

#' Export ggplot
#'
#' @author Jonas Schöley
figspec$ExportFigure <-
  function(figure,
           path,
           filename,
           width = figspec$fig_dims$width,
           height = 100,
           scale = 1,
           device = 'png',
           dpi = 300,
           add_date = FALSE) {
    require(ggplot2)

    if (missing(filename)) {
      filename <- tolower(gsub('\\.', '_', make.names(deparse(substitute(figure)))))
    }
    if (isTRUE(add_date)) {
      filename <- paste0(Sys.Date(), '-', filename)
    }

    arguments <-
      list(
        filename = paste0(filename, '.', device),
        plot = figure,
        path = path,
        width = width,
        height = height,
        units = "mm",
        scale = scale,
        dpi = dpi,
        device = device
      )
    if (device == 'pdf') {
      arguments$useDingbats <- FALSE
    }

    do.call(ggsave, arguments)
  }

#' Export ggplots Stored in List
#'
#' @author Jonas Schöley
figspec$ExportFiguresFromList <- function(lst, path, ...) {
  figure_names <- tolower(gsub('\\.+', '_', make.names(names(lst))))
  Fun <- function (figure, filename, ...) {
    ExportFigure(figure = figure, filename = filename, ...)
  }
  purrr::pwalk(
    list(lst, figure_names),
    Fun, path = path, ...
  )
}