#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(AMORE)
library(dplyr)
library(leaflet)
library(markdown)
library(randomForest)
library(shiny)
library(shinyjs)
library(sf)
library(viridis)

options(shiny.launch.browser = TRUE)
options(rsconnect.quarto = FALSE)
options(shiny.maxRequestSize = 500 * 1024^2)

# ===============================
# GLOBAL DATA
# ===============================

df_pred_base   <- readRDS("data/df_pred.rds")
grid_unique_sf <- readRDS("data/grid_unique_sf.rds")

static_vars <- c("HW", "CW", "WR", "Depth", "days_from_sunday")

range_scale <- function(x, min_val, max_val) {
  (x - min_val) / (max_val - min_val)
}

# ===============================
# USER TEMPLATE
# ===============================

# variables users are allowed to replace
replaceable_vars <- setdiff(
  names(df_pred_base),
  c("id", "date", static_vars)
)

# small template: structure only
df_template <- df_pred_base |>
  select(id, date, all_of(replaceable_vars)) |>
  slice(1:1000)   


# ===============================
# PALETTES
# ===============================

presence_bins <- c(
  "[0,0.3]", "(0.3,0.5]", "(0.5,0.7]",
  "(0.7,0.8]", "(0.8,0.9]", "(0.9,1]"
)

presence_pal <- colorFactor(
  palette = viridis(length(presence_bins), direction = -1),
  domain  = presence_bins,
  ordered = TRUE
)

interaction_pal <- function(x) {
  colorNumeric(
    palette = viridis(256, direction = -1),
    domain  = x,
    na.color = "transparent"
  )
}

# ===============================
# VALIDATION
# ===============================

validate_user_layer <- function(x, df_pred) {
  
  if (!all(c("id", "date") %in% names(x)))
    stop("File must contain 'id' and 'date'")
  
  vars <- setdiff(names(x), c("id", "date"))
  
  if (length(vars) == 0)
    stop("No variables to replace found")
  
  forbidden <- intersect(vars, static_vars)
  if (length(forbidden) > 0)
    stop(paste("Cannot replace:", paste(forbidden, collapse = ", ")))
  
  non_numeric <- vars[!sapply(x[vars], is.numeric)]
  if (length(non_numeric) > 0)
    stop(paste("Non-numeric variables:", paste(non_numeric, collapse = ", ")))
  
  key <- df_pred |> select(id, date) |> distinct()
  if (nrow(anti_join(x, key, by = c("id", "date"))) > 0)
    stop("Invalid id/date combinations detected")
  
  TRUE
}

# ===============================
# FAST PREDICTORS 
# ===============================

range_scale_df <- function(df, minv, maxv) {
  as.data.frame(
    Map(function(x, mn, mx) (x - mn) / (mx - mn), df, minv, maxv),
    check.names = FALSE
  )
}

predict_rf_fast <- function(df, rf, thr, minv, maxv,
                            bsize = NULL,
                            progress_cb = NULL) {
  
  sel_cols <- attr(rf$terms, "term.labels")
  n <- nrow(df)
  prediction <- integer(n)
  
  if (is.null(bsize)) {
    bsize <- max(10000L, min(100000L, floor(n / 20)))
  }
  
  nblock <- (n + bsize - 1L) %/% bsize
  t0 <- Sys.time()
  
  tryCatch({
    
    for (i in seq_len(nblock)) {
      bs <- (i - 1L) * bsize + 1L
      be <- min(i * bsize, n)
      
      Xi  <- df[bs:be, sel_cols, drop = FALSE]
      Xsi <- range_scale_df(Xi, minv, maxv)
      
      p2 <- predict(rf, Xsi, type = "prob")[, 2]
      prediction[bs:be] <- as.integer(p2 > thr)
      
      if (!is.null(progress_cb)) {
        frac <- i / nblock
        elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
        eta <- round(elapsed / frac - elapsed)
        progress_cb(frac, eta)
      }
    }
    
  }, error = function(e) {
    
    Xall  <- df[, sel_cols, drop = FALSE]
    XallS <- range_scale_df(Xall, minv, maxv)
    p2 <- predict(rf, XallS, type = "prob")[, 2]
    prediction <<- as.integer(p2 > thr)
    
    if (!is.null(progress_cb)) {
      progress_cb(1, 0)
    }
  })
  
  df$prediction <- prediction
  df
}

predict_mpn_fast <- function(df, net, thr, sel_cols, minv, maxv,
                             bsize = 50000,
                             progress_cb = NULL) {
  
  X <- range_scale_df(df |> select(all_of(sel_cols)), minv, maxv)
  n <- nrow(X)
  prediction <- integer(n)
  nblock <- (n + bsize - 1L) %/% bsize
  
  t0 <- Sys.time()
  
  for (i in seq_len(nblock)) {
    bs <- (i - 1L) * bsize + 1L
    be <- min(i * bsize, n)
    
    prediction[bs:be] <-
      as.integer(sim.MLPnet(net$net, X[bs:be, , drop = FALSE]) > thr)
    
    if (!is.null(progress_cb)) {
      frac <- i / nblock
      elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      eta <- round(elapsed / frac - elapsed)
      progress_cb(frac, eta)
    }
  }
  
  df$prediction <- prediction
  df
}

# ===============================
# PRELOAD MODELS & METADATA
# ===============================

rf_Tt  <- readRDS("models/rf_Tt.RData")
thr_Tt <- readRDS("models/RF_best_threshold_Tt.RData")
min_Tt <- readRDS("models/train_min_Tt.RData")[-1]
max_Tt <- readRDS("models/train_max_Tt.RData")[-1]

rf_Cc  <- readRDS("models/rf_Cc.RData")
thr_Cc <- readRDS("models/RF_best_threshold_Cc.RData")
min_Cc <- readRDS("models/train_min_Cc.RData")[-1]
max_Cc <- readRDS("models/train_max_Cc.RData")[-1]

net_Mm      <- readRDS("models/MPN_Mm.RData")
thr_Mm      <- readRDS("models/MPN_best_threshold_Mm.RData")
sel_cols_Mm <- readRDS("models/InputVars_Mm.Rdata")
min_Mm      <- readRDS("models/train_min_Mm.RData")[-1]
max_Mm      <- readRDS("models/train_max_Mm.RData")[-1]


# ===============================
# SCENARIO FUNCTIONS
# ===============================

scenario_presence_sf <- function(df) {
  
  df_c <- df[complete.cases(df), ]
  
  agg_pred <- df_c |>
    group_by(id) |>
    summarise(sum_pred = sum(prediction), .groups = "drop")
  
  effort_df <- df |>
    group_by(id) |>
    summarise(sampling_days = n(), .groups = "drop")
  
  agg_pred <- agg_pred |>
    left_join(effort_df, by = "id")
  
  max_effort <- max(agg_pred$sampling_days, na.rm = TRUE)
  
  agg_pred <- agg_pred |>
    mutate(
      log_weight = log(sampling_days + 1),
      log_weight = log_weight / log(max_effort + 1),
      adj_pred_prob = (sum_pred / sampling_days) * log_weight
    )
  
  agg_pred$adj_prob_bin <- cut(
    agg_pred$adj_pred_prob,
    breaks = c(0, 0.3, 0.5, 0.7, 0.8, 0.9, 1),
    include.lowest = TRUE
  )
  
  agg_pred |>
    left_join(grid_unique_sf, by = "id") |>
    st_as_sf()
}

# ===============================
# INTERACTIONS
# ===============================

interaction_Tt <- function(df) {
  df |> mutate(
    int_trawlers = trawlers * prediction,
    int_purse_seines = (other_purse_seines + purse_seines + tuna_purse_seines) * prediction,
    int_fixed = (fixed_gear + set_gillnets) * prediction,
    int_ship = SD * prediction
  ) |>
    group_by(id) |>
    summarise(across(starts_with("int_"), mean), .groups = "drop") |>
    left_join(grid_unique_sf, by = "id") |>
    st_as_sf()
}

interaction_Mm <- function(df) {
  df |> mutate(
    int_trawlers = trawlers * prediction,
    int_purse_seines = other_purse_seines * prediction,
    int_ship = SD * prediction
  ) |>
    group_by(id) |>
    summarise(across(starts_with("int_"), mean), .groups = "drop") |>
    left_join(grid_unique_sf, by = "id") |>
    st_as_sf()
}

interaction_Cc <- function(df) {
  df |> mutate(
    int_trawlers = trawlers * prediction,
    int_purse_seines = (other_purse_seines + purse_seines + tuna_purse_seines) * prediction,
    int_fixed = (set_longlines + fixed_gear + pots_and_traps) * prediction,
    int_ship = SD * prediction
  ) |>
    group_by(id) |>
    summarise(across(starts_with("int_"), mean), .groups = "drop") |>
    left_join(grid_unique_sf, by = "id") |>
    st_as_sf()
}

# ===============================
# UI
# ===============================

ui <- navbarPage(
  title = "SAMESEA Activity 1.1 Explorer",
  inverse = FALSE,

  header = tags$style(HTML("
    .navbar { background-color:#7fbf7f; }
    .navbar-default .navbar-brand,
    .navbar-default .navbar-nav > li > a { color:#ffffff; }
    .navbar-default .navbar-nav > .active > a {
      background-color:#5fa85f; color:white;
    }
    body { background-color:#f4faf4; }
    h2, h3, h4 { color:#2f6f2f; }
    .btn-primary { background-color:#7fbf7f; border-color:#7fbf7f; }
  ")),
  
  # -------- EXPLORER TAB --------
  tabPanel(
    "Explorer",
    sidebarLayout(
      sidebarPanel(
        radioButtons(
          "mode", "Mode",
          choices = c(
            "Explore default results" = "default",
            "Run custom scenario" = "scenario"
          ),
          selected = "default"
        ),
        hr(),
        selectInput(
          "species", "Species",
          choices = c(
            "Tursiops truncatus" = "Tt",
            "Monachus monachus"  = "Mm",
            "Caretta caretta"   = "Cc"
          )
        ),
        radioButtons(
          "map_type", "Map type",
          choices = c(
            "Predicted presence" = "presence",
            "Interaction probability" = "interaction"
          )
        ),
        conditionalPanel(
          condition = "input.map_type == 'interaction'",
          selectInput(
            "activity", "Human activity",
            choices = c(
              "Trawlers" = "int_trawlers",
              "Purse seines" = "int_purse_seines",
              "Fixed gear" = "int_fixed",
              "Ship density" = "int_ship"
            )
          )
        ),
        conditionalPanel(
          condition = "input.mode == 'scenario'",
          hr(),
          fileInput("layer", "Upload custom layer (.rds)", accept = ".rds"),
          helpText("Static layers cannot be replaced:", paste(static_vars, collapse = ", ")),
          helpText("See Help tab for more details on input data requirements."),
          actionButton("apply", "Apply layer"),
          actionButton(
            "run", "Run prediction", class = "btn-primary",
            title = "Run prediction for the selected species. Results are cached per species."
          ),
          div(style="margin-top:5px; color:red;",
              "Note: Tt predictions can take >10 min!"),
          uiOutput("cache_status"),
        )
      ),
      mainPanel(
        leafletOutput("map", height = "750px")
      )
    )
  ),
  
  # -------- HELP TAB --------
  
  tabPanel(
    "Help",
    fluidRow(
      column(
        10, offset = 1,
        
        h3("Custom layer template"),
        p("Download a template .rds file to prepare your custom input layer."),

        downloadButton(
          "download_template",
          "Download input template (.rds)",
          class = "btn-primary"
        ),
        
        hr(),
        includeMarkdown("www/help_custom_layers.md")
      )
    )
  )
  ,
  
  # -------- ABOUT SAMESEA TAB --------
  tabPanel(
    "About SAMESEA",
    fluidRow(
      column(
        4,
        img(src = "samesea_logo.jpg", width = "100%")
      ),
      column(
        8,
        h3("Project Summary"),
        p("SAMESEA aims to develop standardized maritime monitoring practices across the Adriatic-Ionian region (EUSAIR basin), with a main emphasis on the monitoring of marine sentinel species (dolphins, sea turtles, monk seals)."),
        p("In pursuing this main objective, SAMESEA will also improve the dialogue between socio-economic activities and the authorities responsible for marine biodiversity conservation, which is key for designing effective conservation measures."),
        p("A cooperation at the macro-regional level and among different stakeholders will ensure a more functional management of the marine environment, and a long-lasting sustainable coexistence between human activities and marine wildlife."),
        p("The project will create a strategy for the monitoring of marine sentinel species, pilot actions, methods to spread its conservation practices, and an action plan to improve the sustainable management of the EUSAIR basin."),
        p("For the first time in this region, a regional task force will establish a common approach to assess interactions between marine sentinel species and human activities at multiple levels."),
        br(),
        h3("About Activity 1.1"),
        p("WP1: Joint development of the transnational strategy for the monitoring of sentinel species"),
        p("Activity 1.1 focuses on the identification of high-risk areas for sentinel species in the EUSAIR Region, providing a synoptic overview of multi-hazard susceptibility areas."),
        p("This R-based Shiny platform allows users and policy makers to visualize areas of presence and interaction between species and human activities."),
        br(),
        h3("Project Partners"),
        tags$ul(
          tags$li("CoNISMa (Italy)"),
          tags$li("Blue World Institute of Marine Research and Conservation – BWI (Croatia)"),
          tags$li("Faculty of Veterinary Medicine, University of Zagreb – VEFUNIZG (Croatia)"),
          tags$li("Aleksander Moisiu University, Durrës – UAMD (Albania)"),
          tags$li("Centre for Economic, Technological and Environmental Development – CETOR (Bosnia and Herzegovina)"),
          tags$li("Morigenos – Slovenian Marine Mammal Society (Slovenia)"),
          tags$li("Montenegro Dolphin Research – MDR (Montenegro)"),
          tags$li("Ministry of Tourism and Environment – MTE (Albania)"),
          tags$li("Archipelagos Institute of Marine Conservation (Greece)"),
          tags$li("Municipality of Neum (Bosnia and Herzegovina)"),
          tags$li("Veneto Regional Park of the Po Delta (Italy)")
        ),
        
        br(),
        tags$a(
          href = "https://samesea.interreg-ipa-adrion.eu/",
          target = "_blank",
          class = "btn btn-default",
          "Visit the SAMESEA project website"
        ),
        
        br(),
        h3("App Authorship"),
        
        p(
          "Maria Silvia Labriola, Tommaso Russo ",
          "(University of Rome Tor Vergata / CoNISMa)"
        ),
        
        br(),
        h3("Contact Information"),
        
        p(
          "Maria Silvia Labriola (App Developer): maria.silvia.labriola@uniroma2.it"),
        p(
          "CoNISMa (Lead Partner): marinemammals.bca@unipd.it"),

      )
    )
  ),
  
  # -------- TARGET SPECIES TAB --------
  tabPanel(
    "Target Species",
    fluidRow(
      column(4, img(src = "Tt.jpg", width = "100%"),
             h4("Tursiops truncatus"),
             p("The common bottlenose dolphin is a coastal top predator that plays an important role in regulating fish and cephalopod populations and maintaining trophic balance. Due to its long lifespan, high site fidelity, and frequent overlap with human activities such as fisheries and maritime traffic, it is an effective sentinel species. Changes in its distribution and habitat use provide integrated signals of ecosystem condition and anthropogenic pressure in coastal marine environments.
.")),
      column(4, img(src = "Mm.jpg", width = "100%"),
             h4("Monachus monachus"),
             p("The Mediterranean monk seal is a rare coastal apex predator closely associated with undisturbed habitats. Its strong sensitivity to human disturbance, habitat degradation, and coastal pressures makes it an important indicator of ecosystem integrity. The presence or absence of Monachus monachus reflects the level of environmental quality and conservation status of coastal marine ecosystems.")),
      column(4, img(src = "Cc.jpg", width = "100%"),
             h4("Caretta caretta"),
             p("The loggerhead sea turtle is a long-lived marine species that connects pelagic and benthic ecosystems through its wide migratory range and feeding behavior. It integrates multiple human pressures, including fisheries interactions, shipping, pollution, and climate variability. These characteristics make the loggerhead turtle a valuable sentinel for assessing cumulative impacts across large marine areas."))
    ),
    hr(),
    h4("Additional information"),
    p("The analysis highlighted great variability in survey intensity across subareas, as well as uneven representation of target species occurrence. In terms of sightings, spatial and temporal discrepancies were observed, reflecting differences in the availability and coverage of data among the target species. The common bottlenose dolphin (T. truncatus) showed a moderate number of records distributed across a relatively broad time span, whereas reports for the loggerhead sea turtle (C. caretta) and the Mediterranean monk seal (M. monachus) were fewer and more temporally restricted. This pattern suggests that observations of these latter species are either sporadic or linked to limited monitoring efforts.")
  )
)

# ===============================
# SERVER
# ===============================

server <- function(input, output, session) {
  
  applying <- reactiveVal(FALSE)
  running  <- reactiveVal(FALSE)
  
  observe({ shinyjs::toggleState("apply", !applying()) })
  observe({ shinyjs::toggleState("run", !running()) })
  
  df_pred_user <- reactiveVal(df_pred_base)
  scenario_cache <- reactiveValues(Tt = NULL, Mm = NULL, Cc = NULL)
  
  observeEvent(input$apply, {
    req(input$layer)
    applying(TRUE)
    on.exit(applying(FALSE), add = TRUE)
    
    withProgress(message = "Applying custom layer", value = 0, {
      incProgress(0.2, detail = "Reading file")
      user_data <- readRDS(input$layer$datapath)
      
      tryCatch({
        incProgress(0.4, detail = "Validating structure")
        validate_user_layer(user_data, df_pred_base)
        
        vars <- setdiff(names(user_data), c("id", "date"))
        incProgress(0.6, detail = "Joining to prediction grid")
        
        df_new <- df_pred_user() |>
          left_join(user_data, by = c("id", "date"), suffix = c("", ".user"))
        
        incProgress(0.8, detail = "Replacing variables")
        for (v in vars) {
          idx <- !is.na(df_new[[paste0(v, ".user")]])
          df_new[[v]][idx] <- df_new[[paste0(v, ".user")]][idx]
          df_new[[paste0(v, ".user")]] <- NULL
        }
        incProgress(1)
        
        df_pred_user(df_new)
        showNotification(
          paste("Layer applied:", paste(vars, collapse = ", ")),
          type = "message"
        )
        
      }, error = function(e) {
        showNotification(e$message, type = "error", duration = NULL)
      })
    })
  })
  
  observeEvent(input$run, {
    running(TRUE)
    on.exit(running(FALSE), add = TRUE)
    
    withProgress(message = "Running prediction", value = 0, {
      
      progress_cb <- function(frac, eta) {
        incProgress(
          frac * 0.6,
          detail = paste("Running species model (~", eta, "s remaining)")
        )
      }
      
      df <- switch(input$species,
                   "Tt" = predict_rf_fast(df_pred_user(), rf_Tt, thr_Tt, min_Tt, max_Tt, progress_cb = progress_cb),
                   "Cc" = predict_rf_fast(df_pred_user(), rf_Cc, thr_Cc, min_Cc, max_Cc, progress_cb = progress_cb),
                   "Mm" = predict_mpn_fast(df_pred_user(), net_Mm, thr_Mm, sel_cols_Mm, min_Mm, max_Mm, progress_cb = progress_cb)
      )
      
      incProgress(0.8, detail = "Computing presence & interactions")
      
      scenario_cache[[input$species]] <- list(
        presence    = scenario_presence_sf(df),
        interaction = switch(input$species,
                             "Tt" = interaction_Tt(df),
                             "Cc" = interaction_Cc(df),
                             "Mm" = interaction_Mm(df)
        )
      )
      
      incProgress(1)
    })
  })
  
  output$cache_status <- renderUI({
    cached <- names(scenario_cache)[sapply(scenario_cache, Negate(is.null))]
    tags$div(style = "margin-top:5px;font-style:italic;color:darkgreen;",
             if (length(cached)) paste("Cached:", paste(cached, collapse = ", "))
             else "No species cached yet.")
  })
  
  map_data <- reactive({
    if (input$mode == "default") {
      file <- if (input$map_type == "presence")
        file.path("meta/presence", paste0(input$species, "_sf.rds"))
      else
        file.path("meta/interaction", paste0(input$species, "_sf.rds"))
      readRDS(file)
    } else {
      req(scenario_cache[[input$species]])
      if (input$map_type == "presence")
        scenario_cache[[input$species]]$presence
      else
        scenario_cache[[input$species]]$interaction
    }
  })
  
  output$map <- renderLeaflet({
    dat <- map_data()
    req(inherits(dat, "sf"))

    if (input$map_type == "presence") {
      dat$adj_prob_bin <- factor(dat$adj_prob_bin,
                                 levels = presence_bins, ordered = TRUE)
      leaflet(dat) |>
        addProviderTiles("CartoDB.Positron") |>
        addPolygons(
          fillColor = ~presence_pal(adj_prob_bin),
          color = "white",
          weight = 0.3,
          fillOpacity = 0.8
        )
    } else {
      pal <- interaction_pal(dat[[input$activity]])
      leaflet(dat) |>
        addProviderTiles("CartoDB.Positron") |>
        addPolygons(
          fillColor = ~pal(dat[[input$activity]]),
          color = "white",
          weight = 0.3,
          fillOpacity = 0.8
        )
    }
  })
  

  output$download_template <- downloadHandler(
    filename = function() {
      "samesea_custom_layer_template.rds"
    },
    content = function(file) {
      saveRDS(df_template, file)
    }
  )
  
  
}

shinyApp(ui, server)