# Take a list of slendr_pop population boundary objects and
# "interpolate" all of them at time points specified by others
# (unless a given population is supposed to be removed at that time)
fill_maps <- function(pops, time = NULL) {

  removal_times <- sapply(pops, attr, "remove")

  # get times of all spatial maps across all populations
  all_times <- unique(sort(c(
    0,
    time,
    removal_times,
    unlist(sapply(pops, function(i) i$time))
  ))) %>% .[. != Inf & . != -1]

  all_maps <- lapply(seq_along(pops), function(i) {

    # get times where the spatial map of the current population
    # needs to be filled in
    missing_times <- all_times[
      all_times >= removal_times[i] &
        !all_times %in% pops[[i]]$time
    ]

    # generate the missing maps
    new_maps <- lapply(missing_times, function(t) {
      # get all preceding maps
      previous_map <- pops[[i]] %>% .[.$time > t, ]
      if (!nrow(previous_map)) return(NULL)
      latest_map <- previous_map[nrow(previous_map), ]
      latest_map$time <- t
      latest_map
    }) %>%
      do.call(rbind, .)

    if (!is.null(new_maps)) {
      combined_maps <-
        rbind(pops[[i]], new_maps) %>%
        .[order(-.$time), ] %>%
        .[.$time != Inf, ]

      attributes(combined_maps) <- attributes(pops[[i]])
    } else {
      combined_maps <- pops[[i]]
    }

    combined_maps

  })

  all_maps
}


# Pick the next/previous value from a vector
get_time_point <- function(times, current_value, what) {
  current_index <- which(current_value <= times & times <= current_value)

  if (!length(current_index)) {
    if (what == "previous")
      return(times[current_value <= times][1])
    else
      return(utils::tail(times[current_value >= times], 1))
  } else {
    if (what == "previous")
      new_index <- current_index + 1
    else if (what == "next")
      new_index <- current_index - 1
    else
      stop("Invalid direction for the time point selection")

    # prevent jumping out of the allowed range
    if (new_index > length(times) | new_index <= 0) new_index <- current_index

    times[new_index]
  }
}

#' Open an interactive browser of the spatial model
#'
#' @param model Compiled \code{slendr_model} model object
#'
#' @import shiny
#' @export
explore <- function(model) {

  # generate choices for the coordinate system graticules
  if (has_crs(model$world)) {
    crs <- sf::st_crs(model$world)$epsg
    coord_choice <- c("original", "internal")
    names(coord_choice) <- c("original (longitude-latitude)",
                             sprintf("internal (EPSG:%s)", crs))
  } else {
    coord_choice <- list("original (abstract coordinates" = "original")
  }

  # generate event table for manual selection of time points
  split_events <- model$splits
  split_events$event <- with(
    split_events,
    sprintf("split of %s from %s", pop, parent)
  )
  split_events <- split_events[split_events$tsplit_orig != Inf, c("tsplit_orig", "event")]
  colnames(split_events) <- c("time", "event")

  if (!is.null(model$geneflow)) {
    geneflow_starts <- model$geneflow
    geneflow_starts$event <- with(
      geneflow_starts,
      sprintf("geneflow %s -> %s, %.2f%%", from, to, 100 * rate)
    )
    geneflow_starts <- geneflow_starts[, c("tstart_orig", "event")]
    colnames(geneflow_starts) <- c("time", "event")

    geneflow_ends <- model$geneflow
    geneflow_ends$event <- with(
      geneflow_ends,
      sprintf("geneflow %s -> %s ends", from, to)
    )
    geneflow_ends <- geneflow_ends[, c("tend_orig", "event")]
    colnames(geneflow_ends) <- c("time", "event")
  } else {
    geneflow_starts <- NULL
    geneflow_ends <- NULL
  }

  cleanup_events <- do.call(rbind, lapply(model$populations, function(pop)
    data.frame(time = attr(pop, "remove"),
               event = sprintf("%s removed", pop$pop[1]),
               stringsAsFactors = FALSE)
    ))
  cleanup_events <- cleanup_events[with(cleanup_events, time != -1), ]

  events <- do.call(rbind, list(split_events, geneflow_starts,
                                geneflow_ends, cleanup_events))
  events <- stats::aggregate(event~time, data = events, FUN = paste, collapse = ", ")
  events$label <- sprintf("time %s: %s", events$time, events$event)
  events <- events[order(events$time), ]
  event_choices <- events$time
  names(event_choices) <- events$label

  # generate time points for the slider
  time_point_snapshots <-
    as.integer(c(0, event_choices, unlist(lapply(model$populations, `[[`, "time"))) %>%
    sort %>%
    unique %>%
    .[. != Inf])

  interpolated_maps <- fill_maps(model$populations, time_point_snapshots)

  ui <- fluidPage(
    tags$style(type = "text/css", ".recalculating { opacity: 1.0; }"),

    navbarPage(
      "Model explorer",

      tabPanel(
        "Spatial dynamics",

        sidebarLayout(

          sidebarPanel(

            htmlOutput(outputId = "time_label"),

            br(),

            fluidRow(
              column(2, actionButton("previous_time", label = "",
                                     icon = icon("angle-double-left", "fa-1x"))),

              column(8, shinyWidgets::sliderTextInput(
                inputId = "time_slider",
                label = "",
                choices = rev(time_point_snapshots),
                selected = max(time_point_snapshots),
                width = "100%",
                animate = animationOptions(interval = 2000, loop = FALSE)
              )),

              column(2, actionButton("next_time", label = "",
                                     icon = icon("angle-double-right", "fa-1x")))
            ),

            selectInput(
              inputId = "time_select",
              label = "Select event:",
              choices = event_choices,
              selected = max(time_point_snapshots)
            ),

            selectInput(
              inputId = "coord_system",
              label = "Coordinate system:",
              choices = coord_choice
            ),

            fluidRow(
              column(4, checkboxInput(
                inputId = "intersect",
                label = "Intersect ranges",
                value = TRUE
              )),

              if (nrow(model$world)) {
                column(4, checkboxInput(
                  inputId = "show_map",
                  label = "Show map",
                  value = TRUE
                ))
              } else NULL,

              if (!is.null(model$geneflow)) {
                column(4, checkboxInput(
                  inputId = "show_geneflows",
                  label = "Indicate geneflows",
                  value = TRUE
                ))
              } else NULL,
            ),

            p(strong("Generation time: "), model$generation_time, " time units")

          ),

          mainPanel(

            fluidRow(
              align = "center",

              plotOutput(outputId = "slendr_maps", height = 480),

              hr(),

              tableOutput("geneflows_table")

            )

          )
        )

      ), # tabPanel

      tabPanel(
        "Population history graph",

        sidebarLayout(

          sidebarPanel(

            checkboxInput(
              inputId = "show_cleanups",
              label = "Show removal times",
              value = TRUE
            ),
            width = 3

          ),

          mainPanel(

            plotOutput(outputId = "slendr_graph")

          )

        )

      ) # tabPanel

    ) # navbarPage

  )

  server <- function(input, output, session) {

    output$time_label = renderText({
      event <- events[events$time == input$time_slider, "event"]
      if (length(event))
        label <- sprintf("<i>(%s)</i>", event)
      else
        label <- ""
      sprintf("<b>Time point:</b> %s %s", input$time_slider, label)
    })

    observeEvent(input$time_select, {
      value <- as.numeric(input$time_select)
      shinyWidgets::updateSliderTextInput(session, "time_slider", selected = value)
    }, ignoreInit = TRUE)

    observeEvent(input$previous_time, {
      value <- get_time_point(time_point_snapshots, input$time_slider, "previous")
      shinyWidgets::updateSliderTextInput(session, "time_slider", selected = value)
    })

    observeEvent(input$next_time, {
      value <- get_time_point(time_point_snapshots, input$time_slider, "next")
      shinyWidgets::updateSliderTextInput(session, "time_slider", selected = value)
    })

    output$slendr_maps <- renderPlot({

      plot(
        model,
        time = input$time_slider,
        graticules = input$coord_system,
        intersect = input$intersect,
        show_map = input$show_map,
        geneflows = if (is.null(model$geneflow)) FALSE else input$show_geneflows,
        interpolated_maps = interpolated_maps
      )

    })

    output$geneflows_table <- renderTable({
      if (!is.null(model$geneflow)) {
        migr_df <- get_geneflows(model, input$time_slider)
        table <- migr_df[, c("from", "to", "tstart_orig", "tend_orig", "rate")]
        table$rate_gen <- sprintf("%.1f%%", table$rate / model$generation_time * 100)
        table$tstart_orig <- as.integer(table$tstart_orig)
        table$tend_orig <- as.integer(table$tend_orig)
        table$rate <- sprintf("%.1f%%", table$rate * 100)
        colnames(table) <- c("source", "target", "start", "end", "rate", "rate per gen.")
        table$overlapping <- ifelse(migr_df$overlap, "yes", "no")
        if (!nrow(table)) return(NULL)
        table
      } else return(NULL)
    }, sanitize.text.function = identity)

    output$slendr_graph <- renderPlot({ graph(model, input$show_cleanups) },
                                      height = 600)

  }

  shinyApp(ui, server)
}
