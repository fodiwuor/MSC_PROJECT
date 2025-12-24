# =========================================================
# app.R  (shinyapps.io safe — NO install.packages() in app)
# =========================================================

options(repos = c(CRAN = "https://cloud.r-project.org"))

# ---------------------------
# 0) Packages
# ---------------------------
required_packages <- c(
  "shiny","shinydashboard","dplyr","tidyr","stringr",
  "DT","readr","jsonlite","httr"
)

missing_pkgs <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkgs) > 0) {
  stop("Missing packages: ", paste(missing_pkgs, collapse = ", "))
}

suppressPackageStartupMessages({
  library(shiny)
  library(shinydashboard)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(DT)
  library(readr)
  library(jsonlite)
  library(httr)
})

# ---------------------------
# 1) GitHub settings
# ---------------------------
GITHUB_USER   <- "fodiwuor"
GITHUB_REPO   <- "MSC_PROJECT"
GITHUB_BRANCH <- "main"
GITHUB_FOLDER <- "data"
DATA_FILE     <- "perfRshiny.csv"

# ---------------------------
# 2) GitHub helpers
# ---------------------------
gh_api_contents <- function() {
  api <- sprintf(
    "https://api.github.com/repos/%s/%s/contents/%s?ref=%s",
    GITHUB_USER, GITHUB_REPO, GITHUB_FOLDER, GITHUB_BRANCH
  )
  
  res <- httr::GET(
    api,
    httr::add_headers(
      "User-Agent" = "R-shiny",
      "Accept"     = "application/vnd.github+json"
    ),
    httr::timeout(30)
  )
  
  txt <- httr::content(res, as = "text", encoding = "UTF-8")
  if (httr::status_code(res) != 200) stop(txt)
  
  jsonlite::fromJSON(txt, flatten = TRUE)
}

gh_download_url <- function(items, filename) {
  hit <- items[items$name == filename, , drop = FALSE]
  if (nrow(hit) == 0) stop("CSV not found in GitHub repo.")
  hit$download_url[1]
}

download_to_temp <- function(url) {
  tf <- tempfile(fileext = ".csv")
  res <- httr::GET(
    url,
    httr::write_disk(tf, overwrite = TRUE),
    httr::timeout(60)
  )
  if (httr::status_code(res) != 200) stop("Download failed")
  tf
}

# ---------------------------
# 3) UI
# ---------------------------
ui <- dashboardPage(
  dashboardHeader(title = "Results Dashboard"),
  
  dashboardSidebar(
    sidebarMenu(
      menuItem("Results Table", tabName = "tbl", icon = icon("table"))
    ),
    
    tags$div(
      style = "padding:10px; font-size:12px; line-height:1.35;",
      tags$b("Authors:"), tags$br(),
      "Fredrick Orwa", tags$sup("*"), ", Charles Mutai, Ann Mwangi", tags$br(),
      tags$b("Affiliation:"), tags$br(),
      "Department of Mathematics, Physics, and Computing, Moi University, Kenya", tags$br(),
      tags$b("Emails:"), tags$br(),
      "orwafredrick95@gmail.com; charlimtai@gmail.com; annwsum@gmail.com", tags$br(),
      tags$b("*Corresponding author"), tags$br(),
      tags$b("GitHub:"), " ",
      tags$a(
        href = "https://github.com/fodiwuor/MSC_PROJECT/tree/main/data",
        "data folder",
        target = "_blank"
      )
    ),
    
    selectInput(
      "table_sel",
      "Select results table to view",
      choices = c(
        "Point estimate","Bias","Coverage",
        "Power","MSE","Standard error (SE)"
      ),
      selected = "Point estimate"
    ),
    
    # ✅ ADDED: "Filter by:" heading
    tags$div(
      style = "padding:6px 10px 0px 10px; font-weight:bold;",
      "Filter by:"
    ),
    
    uiOutput("rho_ui"),
    uiOutput("n_ui"),
    uiOutput("method_ui"),
    uiOutput("effect_ui")
  ),
  
  dashboardBody(
    tabItems(
      tabItem(
        tabName = "tbl",
        box(
          width = 12,
          title = textOutput("tbl_title"),
          DTOutput("tbl_out")
        )
      )
    )
  )
)

# ---------------------------
# 4) SERVER
# ---------------------------
server <- function(input, output, session) {
  
  gh_items <- reactiveVal(NULL)
  
  observeEvent(TRUE, {
    gh_items(gh_api_contents())
  }, once = TRUE)
  
  data_url <- reactive({
    req(gh_items())
    gh_download_url(gh_items(), DATA_FILE)
  })
  
  dat_all <- reactive({
    tf <- download_to_temp(data_url())
    d  <- read_csv(tf, show_col_types = FALSE, name_repair = "minimal")
    
    if ("true_val" %in% names(d)) {
      d$Estimand <- d$true_val
      d$true_val <- NULL
    }
    
    if ("rho" %in% names(d)) d$rho <- suppressWarnings(as.numeric(d$rho))
    if ("n"   %in% names(d)) d$n   <- suppressWarnings(as.numeric(d$n))
    
    d
  })
  
  output$rho_ui <- renderUI({
    d <- dat_all()
    req("rho" %in% names(d))
    choices <- sort(unique(d$rho))
    
    checkboxGroupInput(
      "rho_sel",
      "Autocorrelation (ρ)",
      choices = choices,
      selected = choices
    )
  })
  
  output$n_ui <- renderUI({
    d <- dat_all()
    req("n" %in% names(d))
    
    selectInput(
      "n_sel",
      "Length of series (n)",
      choices = c("All", sort(unique(d$n))),
      selected = "All"
    )
  })
  
  output$method_ui <- renderUI({
    d <- dat_all()
    req("Method" %in% names(d))
    
    selectInput(
      "m_sel",
      "Method",
      choices = c("All", sort(unique(d$Method))),
      selected = "All"
    )
  })
  
  output$effect_ui <- renderUI({
    d <- dat_all()
    req("EstimandScenario" %in% names(d))
    
    # ✅ CHANGED LABEL ONLY (still filters EstimandScenario)
    selectInput(
      "eff_sel",
      "Estimand scenario",
      choices = c("All", sort(unique(d$EstimandScenario))),
      selected = "All"
    )
  }) 
  
  dat_filt <- reactive({
    d <- dat_all()
    
    if (!is.null(input$rho_sel) && length(input$rho_sel) > 0) {
      d <- dplyr::filter(d, rho %in% as.numeric(input$rho_sel))
    }
    
    if (!is.null(input$n_sel) && input$n_sel != "All") {
      d <- dplyr::filter(d, n == as.numeric(input$n_sel))
    }
    
    if (!is.null(input$m_sel) && input$m_sel != "All") {
      d <- dplyr::filter(d, Method == input$m_sel)
    }
    
    if (!is.null(input$eff_sel) && input$eff_sel != "All") {
      d <- dplyr::filter(d, EstimandScenario == input$eff_sel)
    }
    
    d
  })
  
  pick_table_cols <- function(d, which) {
    base <- c("rho","n","EstimandScenario","Method")
    
    cols <- switch(
      which,
      "Point estimate"      = c(base,"Policy effect (Empirical SE)","Estimand"),
      "Bias"                = c(base,"Bias (MCSE)"),
      "Coverage"            = c(base,"Coverage 95% (MCSE)"),
      "Power"               = c(base,"Power (MCSE)"),
      "MSE"                 = c(base,"Mean square error (MCSE)"),
      "Standard error (SE)" = c(base,"Empirical SE (Model-based SE)"),
      base
    )
    
    dplyr::select(d, dplyr::any_of(cols))
  }
  
  output$tbl_title <- renderText(
    paste("Table:", input$table_sel)
  )
  
  output$tbl_out <- renderDT({
    d <- dat_filt()
    
    if (nrow(d) == 0) {
      return(
        DT::datatable(
          data.frame(
            Message = "No rows match your filters. Try selecting All for n/method/effect and all rho."
          ),
          rownames = FALSE
        )
      )
    }
    
    view_d <- pick_table_cols(d, input$table_sel)
    
    names(view_d)[names(view_d) == "rho"]              <- "Autocorrelation (ρ)"
    names(view_d)[names(view_d) == "n"]                <- "Length of series (n)"
    names(view_d)[names(view_d) == "EstimandScenario"] <- "Estimand scenario"
    names(view_d)[names(view_d) == "Policy effect (Empirical SE)"] <- "Policy effect estimate (Empirical SE)"
    
    DT::datatable(
      view_d,
      rownames = FALSE,
      options = list(pageLength = 25, scrollX = TRUE)
    )
  })
}

shinyApp(ui, server)
