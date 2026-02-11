library(shiny)
library(shinydashboard)
library(shinyjs)
library(shinyWidgets)
library(shinyjqui)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(plotly)

header <- dashboardHeader(
  title = "CORALIE",                      
  titleWidth = 600                        
)

sidebar <- dashboardSidebar(
  sidebarMenu(
    id = "main_tier",                     
    
    menuItem("Analysis Interface", 
             tabName = "tier1") 
    
  )
)

body <- dashboardBody(
  useShinyjs(),
  tags$style(HTML("
  /* Full-page loading overlay */
  #coralie-loading-overlay.loading-overlay {
    position: fixed;
    inset: 0;
    display: flex;
    align-items: center;
    justify-content: center;
    background: radial-gradient(circle at top,
      rgba(16, 19, 34, 0.45),
      rgba(5, 7, 18, 0.35) 60%,
      rgba(2, 3, 9, 0.25) 100%);
    z-index: 9999;
  }

  #coralie-loading-overlay.hidden {
    display: none;
  }
  
  #coralie-loading-bar {
  width: 220px;
  height: 4px;
  border-radius: 999px;
  background: rgba(129,140,248,0.25);
  overflow: hidden;
  margin-top: 6px;
}

#coralie-loading-bar-fill {
  height: 100%;
  width: 0%;
  border-radius: 999px;
  background: linear-gradient(90deg, #C2185B, #8E24AA, #1976D2);
  transition: width 0.18s ease-out;
}

  .coralie-loading-core {
    display: flex;
    flex-direction: column;
    align-items: center;
    justify-content: center;
    text-align: center;
    gap: 0.75rem;
    position: relative;
    z-index: 2;
  }

  /* Diamond logo block (rotated square) */
  .coralie-loading-logo {
    position: relative;
    width: 90px;
    height: 90px;
    border-radius: 20px;
    background: linear-gradient(135deg, #8E24AA, #C2185B, #1976D2);
    display: flex;
    align-items: center;
    justify-content: center;
    box-shadow: 0 10px 30px rgba(0,0,0,0.35);
    overflow: hidden;
    animation: coraliePulse 1.8s ease-in-out infinite;
  }

  .coralie-diamond {
    width: 70px;
    height: 70px;
    background: radial-gradient(circle at 20% 20%, #FFCDD2, #C2185B 35%, #8E24AA 80%);
    transform: rotate(45deg);
    position: relative;
    overflow: hidden;
  }

  /* Static cross to echo CORALIE pattern */
  .coralie-diamond::before,
  .coralie-diamond::after {
    content: '';
    position: absolute;
    background: rgba(255,255,255,0.13);
  }

  .coralie-diamond::before {
    width: 22%;
    height: 130%;
    left: 39%;
    top: -15%;
  }

  .coralie-diamond::after {
    height: 22%;
    width: 130%;
    top: 39%;
    left: -15%;
  }

  /* Moving scan bar inside diamond */
  .coralie-diamond-scan {
    position: absolute;
    width: 140%;
    height: 35%;
    top: -20%;
    left: -20%;
    background: linear-gradient(
      to bottom,
      rgba(255,255,255,0.0),
      rgba(255,255,255,0.35),
      rgba(255,255,255,0.0)
    );
    transform: rotate(-45deg);
    mix-blend-mode: screen;
    animation: coralieScan 2.2s ease-in-out infinite;
  }

  @keyframes coralieScan {
    0%   { transform: translateY(-60%) rotate(-45deg); opacity: 0.0; }
    15%  { opacity: 0.4; }
    50%  { transform: translateY(40%) rotate(-45deg); opacity: 0.7; }
    85%  { opacity: 0.4; }
    100% { transform: translateY(120%) rotate(-45deg); opacity: 0.0; }
  }

  /* Rotating inner highlight ring */
  .coralie-diamond-glow {
    position: absolute;
    inset: 10%;
    border-radius: 20px;
    border: 2px solid rgba(255,255,255,0.35);
    box-shadow:
      0 0 8px rgba(255,255,255,0.45),
      0 0 16px rgba(244,143,177,0.55);
    mix-blend-mode: screen;
    animation: coralieGlowRotate 3.5s linear infinite;
  }

  @keyframes coralieGlowRotate {
    0%   { transform: rotate(0deg) scale(0.96); opacity: 0.7; }
    40%  { transform: rotate(80deg) scale(1.02); opacity: 1.0; }
    80%  { transform: rotate(150deg) scale(0.98); opacity: 0.6; }
    100% { transform: rotate(180deg) scale(0.96); opacity: 0.7; }
  }

  /* Outer pulse stays as before */
  @keyframes coraliePulse {
    0%   { transform: scale(0.96); box-shadow: 0 8px 24px rgba(0,0,0,0.30); }
    50%  { transform: scale(1.03); box-shadow: 0 14px 36px rgba(0,0,0,0.45); }
    100% { transform: scale(0.96); box-shadow: 0 8px 24px rgba(0,0,0,0.30); }
  }

  #coralie-loading-text {
    color: #e0e7ff;
    font-size: 1.5rem;
    letter-spacing: 0.12em;
    text-transform: uppercase;
    text-shadow:
      0 0 4px rgba(129,140,248,0.9),
      0 0 12px rgba(167,139,250,0.9);
  }

  #coralie-loading-subtext {
    color: #c5cae9;
    font-size: 0.85rem;
    letter-spacing: 0.08em;
    text-transform: uppercase;
    opacity: 0.8;
  }

  /* Container for animated diagonal line segments */
  .coralie-stripes {
    position: fixed;
    inset: 0;
    pointer-events: none;
    overflow: hidden;
    z-index: 1;
  }

  .coralie-line {
    position: absolute;
    width: 180px;
    height: 2px;
    opacity: 0;
    transform-origin: 0 50%;  /* left end fixed */
    filter: blur(0.4px);
  }

  /* Brighter red/magenta lines: bottom-left -> top-right (45deg) */
  .coralie-line-blue {
    background: linear-gradient(90deg,
      rgba(233,30,99,0.0),
      rgba(255,255,255,0.9),
      rgba(233,30,99,0.95),
      rgba(244,143,177,0.0)
    );
    transform: rotate(45deg) scaleX(0);
    animation: coralieLineBlueDraw 2.4s ease-in-out infinite;
  }

  /* Brighter blue lines: top-left -> bottom-right (135deg) */
  .coralie-line-red {
    background: linear-gradient(90deg,
      rgba(33,150,243,0.0),
      rgba(255,255,255,0.9),
      rgba(33,150,243,0.95),
      rgba(129,212,250,0.0)
    );
    transform: rotate(135deg) scaleX(0);
    animation: coralieLineRedDraw 2.4s ease-in-out infinite;
  }

  /* Make some lines slightly thicker for variation */
  .coralie-line:nth-child(odd) {
    height: 3px;
  }

  /* Faster, more dynamic draw/fade */
  @keyframes coralieLineBlueDraw {
    0% {
      opacity: 0.0;
      transform: rotate(45deg)  scaleX(0.0);
    }
    20% {
      opacity: 0.8;
      transform: rotate(45deg)  scaleX(0.6);
    }
    50% {
      opacity: 0.6;
      transform: rotate(45deg)  scaleX(1.1);
    }
    80% {
      opacity: 0.2;
      transform: rotate(45deg)  scaleX(0.7);
    }
    100% {
      opacity: 0.0;
      transform: rotate(45deg)  scaleX(0.0);
    }
  }

  /* Red lines: move DOWN along 135deg while drawing */
  @keyframes coralieLineRedDraw {
    0% {
      opacity: 0.0;
      transform: rotate(135deg) scaleX(0.0);
    }
    20% {
      opacity: 0.8;
      transform: rotate(135deg) scaleX(0.6);
    }
    50% {
      opacity: 0.6;
      transform: rotate(135deg)  scaleX(1.1);
    }
    80% {
      opacity: 0.2;
      transform: rotate(135deg)  scaleX(0.7);
    }
    100% {
      opacity: 0.0;
      transform: rotate(135deg) scaleX(0.0);
    }
  }

  /* Positions + staggered delays for more apparent motion */
  .coralie-line:nth-child(1)  { left: 5%;   top: 82%; animation-delay: 0.0s; }
  .coralie-line:nth-child(2)  { left: 18%;  top: 64%; animation-delay: 0.3s; }
  .coralie-line:nth-child(3)  { left: 8%;   top: 32%; animation-delay: 0.6s; }
  .coralie-line:nth-child(4)  { left: 32%;  top: 88%; animation-delay: 0.9s; }
  .coralie-line:nth-child(5)  { left: 52%;  top: 72%; animation-delay: 1.2s; }
  .coralie-line:nth-child(6)  { left: 70%;  top: 22%; animation-delay: 0.4s; }
  .coralie-line:nth-child(7)  { left: 82%;  top: 46%; animation-delay: 0.8s; }
  .coralie-line:nth-child(8)  { left: 60%;  top: 12%; animation-delay: 1.1s; }
  .coralie-line:nth-child(9)  { left: 30%;  top: 18%; animation-delay: 1.5s; }
  .coralie-line:nth-child(10) { left: 46%;  top: 40%; animation-delay: 1.9s; }
")),
  
  tags$script(HTML("
  Shiny.addCustomMessageHandler('resize-tier1-bar', function() {
    var el = document.getElementById('tier1_subsystem_bar');
    if (!el) return;
    var rect = el.getBoundingClientRect();
    // Example: ensure minimum 400px but no more than 80% viewport
    var h = Math.max(400, Math.min(window.innerHeight * 0.8, rect.height));
    el.style.height = h + 'px';
  });
")),
  # JS handlers to show/hide + update text from server
  tags$script(HTML("
  Shiny.addCustomMessageHandler('coralie-toggle-loading', function(show) {
    var overlay = document.getElementById('coralie-loading-overlay');
    if (!overlay) return;
    if (show) {
      overlay.classList.remove('hidden');
    } else {
      overlay.classList.add('hidden');
    }
  });

  Shiny.addCustomMessageHandler('coralie-loading-text', function(msg) {
    var el = document.getElementById('coralie-loading-text');
    if (!el) return;
    el.textContent = msg;
  });

  Shiny.addCustomMessageHandler('coralie-loading-subtext', function(msg) {
    var el = document.getElementById('coralie-loading-subtext');
    if (!el) return;
    el.textContent = msg;
  });

  Shiny.addCustomMessageHandler('coralie-loading-progress', function(pct) {
    var bar = document.getElementById('coralie-loading-bar-fill');
    if (!bar) return;
    // clamp 0–100
    var v = Math.max(0, Math.min(100, pct));
    bar.style.width = v + '%';
  });
")),
  
  tags$style(HTML("
      body, .content-wrapper, .box-title, .sidebar-menu li a {
        font-family: Poppins, -apple-system, BlinkMacSystemFont, 'Segoe UI',
                     Roboto, 'Helvetica Neue', Arial, sans-serif;
      }

      /* CORALIE navbar + logo (already present) */
      .skin-purple .main-header .navbar {
        background: linear-gradient(135deg, #8E24AA, #C2185B, #1976D2);
        border: none;
      }
      .skin-purple .main-header .logo {
        background-color: #8E24AA;
        color: #FFFFFF;
        font-family: Poppins, -apple-system, BlinkMacSystemFont, 'Segoe UI',
                     Roboto, 'Helvetica Neue', Arial, sans-serif;
        font-weight: 600;
        font-size: 17px;
        letter-spacing: 0.03em;
        border: none;
      }
      .skin-purple .main-header .logo:hover {
        background-color: #7B1FA2;
      }

      /* xCheck-like background + box styling */
      .content-wrapper, .right-side {
        background-color: #f5f6fa;
      }

      .box {
        border-radius: 10px;
        border: 1px solid #e0e4f0;
        box-shadow: 0 4px 12px rgba(0,0,0,0.04);
        overflow: visible;
      }
      .box-header {
        border-bottom: 1px solid #e0e4f0;
        border-top-left-radius: 10px;
        border-top-right-radius: 10px;
        background-clip: padding-box;
      }
      .box-body {
        border-bottom-left-radius: 10px;
        border-bottom-right-radius: 10px;
        background-clip: padding-box;
      }

          :root {
        /* CORALIE-inspired primary colors */
        --coralie-primary:      #8E24AA;  /* left of header gradient */
        --coralie-primary-mid:  #C2185B;  /* magenta */
        --coralie-primary-dark: #6A1B9A;
      }

      /* Primary boxes pick up CORALIE colors */
      .skin-purple .box.box-primary {
        border-top-color: var(--coralie-primary-mid);
      }

      .skin-purple .box.box-primary .box-header {
        background: linear-gradient(135deg,
                   var(--coralie-primary),
                   var(--coralie-primary-mid));
        color: #ffffff;
      }

      .skin-purple .box.box-primary .box-header .box-title {
        color: #ffffff;
        letter-spacing: 0.03em;
      }

      .skin-purple .box.box-primary .box-body {
        background-color: #fbf5ff;  /* very light lavender/pink */
        color: #283046;
      }

      /* Buttons aligned with header gradient */
      .skin-purple .btn-primary {
        background: linear-gradient(135deg,
                     var(--coralie-primary),
                     var(--coralie-primary-mid));
        border-color: var(--coralie-primary-mid);
      }
      .skin-purple .btn-primary:hover,
      .skin-purple .btn-primary:focus {
        background: linear-gradient(135deg,
                     var(--coralie-primary-dark),
                     var(--coralie-primary-mid));
        border-color: var(--coralie-primary-dark);
      }

      .btn {
        border-radius: 20px;
        font-size: 12px;
        font-weight: 500;
        padding: 5px 14px;
      }

      .form-group {
        margin-bottom: 10px;
      }
      .control-label {
        font-weight: 500;
        font-size: 12px;
        color: #5f6473;
      }

      /* Fade-in-up animation (xCheck style) */
      .fade-in-up {
        opacity: 0;
        transform: translateY(10px);
        animation: fadeInUp 0.4s ease-out forwards;
      }
      @keyframes fadeInUp {
        from { opacity: 0; transform: translateY(10px); }
        to   { opacity: 1; transform: translateY(0); }
      }

      /* Fade-in for updated boxes */
      .fade-in-panel {
        animation: fadeInPanel 0.4s ease-in-out;
      }
      @keyframes fadeInPanel {
        from { opacity: 0; transform: translateY(4px); }
        to   { opacity: 1; transform: translateY(0); }
      }
    ")),
  tags$script(HTML("
      // Add fade-in-panel to any box containing an updated output
      document.onshinyvalue = function(e) {
        var el = document.getElementById(e.target.id);
        if (!el) return;
        var box = el.closest('.box');
        if (!box) return;
        box.classList.add('fade-in-panel');
        setTimeout(function() {
          box.classList.remove('fade-in-panel');
        }, 500);
      };
    ")),
  
  div(
    id = "coralie-loading-overlay",
    class = "loading-overlay hidden",
    div(
      class = "coralie-loading-core",
      div(
        class = "coralie-loading-logo",
        div(
          class = "coralie-diamond",
          div(class = "coralie-diamond-scan"),
          div(class = "coralie-diamond-glow")
        )
      ),
      div(id = "coralie-loading-text", "Preparing CORALIE analysis"),
      div(id = "coralie-loading-subtext", ""),
      div(
        id = "coralie-loading-bar",
        div(id = "coralie-loading-bar-fill")
      )
    ),
    div(
      class = "coralie-stripes",
      # 5 blue + 5 red lines, positions/delays via nth-child in CSS
      lapply(1:5, function(i) span(class = "coralie-line coralie-line-blue")),
      lapply(1:5, function(i) span(class = "coralie-line coralie-line-red"))
    )
  ),
  
  tabItems(
    tabItem(
      tabName = "tier1",
      class = "fade-in-up",
      fluidRow(
        box(
          width = 12, status = "primary", solidHeader = TRUE,
          title = "Fingerprint selection",
          shinyWidgets::pickerInput(
            inputId = "fingerprints_sel_1",
            label   = "Select two or more fingerprints",        # or "Select fingerprints"
            choices = NULL,        # populated server-side
            multiple = TRUE,
            options = list(
              `actions-box` = TRUE,          # adds Select All / Deselect All buttons
              `none-selected-text` = "No fingerprints selected"
            )
          ),
          br(),
          actionButton(
            inputId = "load_fingerprints",
            label   = "Select fingerprints"
          )
        )
      ),
      
      fluidRow(
        div(id = "tier1_reference_placeholder")
      ),
      fluidRow(
        div(id = "tier1_summary_placeholder")
      ),
      fluidRow(
        div(id = "tier1_analysis_placeholder")
      ),
      fluidRow(
        div(id = "tier1_corplot_grid_placeholder")
      ),
      fluidRow(
        div(id = "tier3_analysis_placeholder"))
    )
  )
)

ui <- dashboardPage(
  skin   = "purple",                      
  header = header,
  sidebar = sidebar,
  body   = body
)

server <- function(input, output, session) {
  fingerprints_dir <- "../../data/fingerprints"
  reference_dir   <- "../../data/reference_assays"
  
  fingerprints_list_1 <- reactiveVal(NULL)
  ref_box_inserted_1  <- reactiveVal(FALSE)
  reference_list_1    <- reactiveVal(NULL)
  
  fingerprint_display_names <- reactiveVal(NULL)
  reference_display_names   <- reactiveVal(NULL)
  
  fingerprint_display_order <- reactiveVal(NULL)
  reference_display_order   <- reactiveVal(NULL)
  
  summary_box_inserted_1    <- reactiveVal(FALSE)
  
  analysis_panel_inserted_1 <- reactiveVal(FALSE)
  analysis_done_1           <- reactiveVal(FALSE)
  tier1_cor_data            <- reactiveVal(NULL)
  
  tier2_selection <- reactiveVal(NULL)
  tier2_panel_inserted <- reactiveVal(FALSE)
  
  clicked_subsystem <- reactiveVal(NULL)
  tier3_box_inserted <- reactiveVal(FALSE)
  
  tier2_x_range <- reactiveVal(NULL)
  tier2_y_range <- reactiveVal(NULL)
  
  set_coralie_progress <- function(session, pct, subtext = NULL) {
    session$sendCustomMessage("coralie-loading-progress", round(pct))
    if (!is.null(subtext)) {
      session$sendCustomMessage("coralie-loading-subtext", subtext)
    }
  }
  
  # Populate choices from .rds files
  observe({
    files <- list.files(
      path = fingerprints_dir,
      pattern = "\\_.rds$",
      full.names = FALSE
    )
    stripped <- sub("\\_.rds$", "", files)
    
    shinyWidgets::updatePickerInput(
      session,
      inputId  = "fingerprints_sel_1",
      choices  = stripped,
      selected = NULL            # default: nothing selected
    )
  })
  
  # Keep the button disabled until >= 2 selected (optional but nice UX)
  observe({
    n_sel_1 <- length(input$fingerprints_sel_1)
    
    if (is.null(n_sel_1) || n_sel_1 < 2) {
      shinyjs::disable("load_fingerprints")
    } else {
      shinyjs::enable("load_fingerprints")
    }
  })
  
  human_react_meta<-read.csv("../../data/fingerprint_prep_objects/human_reaction_meta.csv", header = T, row.names = 1)
  mouse_react_meta<-read.csv("../../data/fingerprint_prep_objects/mouse_reaction_meta.csv", header = T, row.names = 1)
  react_meta_filt<-human_react_meta[intersect(rownames(human_react_meta), rownames(mouse_react_meta)),]
  
  # --- 3. Load fingerprints and, on first success, insert reference box ---
  observeEvent(input$load_fingerprints, {
    req(input$fingerprints_sel_1)
    req(length(input$fingerprints_sel_1) >= 2)
    
    session$sendCustomMessage("coralie-loading-text", "Loading fingerprints")
    session$sendCustomMessage("coralie-loading-subtext", "")
    session$sendCustomMessage("coralie-loading-progress", 0)
    session$sendCustomMessage("coralie-toggle-loading", TRUE)
    
    sel_names_1 <- input$fingerprints_sel_1
    paths <- file.path(
      fingerprints_dir,
      paste0(sel_names_1, "_.rds")
    )
    
    sizes <- file.info(paths)$size
    sizes[is.na(sizes) | sizes <= 0] <- 1
    total_size <- sum(sizes)
    cum_size <- 0
    
    fingerprints_1 <- list()
    for (i in seq_along(paths)) {
      p <- paths[i]
      nm <- sel_names_1[i]
      
      set_coralie_progress(
        session,
        100 * cum_size / total_size,
        sprintf("Loading fingerprint %d/%d: %s",
                i, length(paths), nm)
      )
      
      obj <- readRDS(p)
      fingerprints_1[[nm]] <- obj
      
      cum_size <- cum_size + sizes[i]
      set_coralie_progress(
        session,
        100 * cum_size / total_size
      )
    }
    
    names(fingerprints_1) <- sel_names_1
    fingerprints_list_1(fingerprints_1)
    
    raw_fp_names <- names(fingerprints_1)
    if (is.null(fingerprint_display_names()))
      fingerprint_display_names(raw_fp_names)
    if (is.null(fingerprint_display_order()))
      fingerprint_display_order(raw_fp_names)
    
    session$sendCustomMessage("coralie-toggle-loading", FALSE)
    session$sendCustomMessage("coralie-loading-subtext", "")
    
    if (!ref_box_inserted_1()) {
      insertUI(
        selector = "#tier1_reference_placeholder",
        where    = "afterBegin",
        ui = fluidRow(
          class = "fade-in-up",
          column(width = 12,
                 box(
                   width = 12, status = "primary", solidHeader = TRUE,
                   title = "Reference assay selection",
                   selectInput(
                     inputId  = "reference_sel_1",
                     label    = "Select one or more reference assays",
                     choices  = NULL,
                     multiple = TRUE
                   ),
                   actionButton(
                     inputId = "load_reference_1",
                     label   = "Select reference assay(s)"
                   )
                 )
          )
        ),
        immediate = TRUE
      )
      
      ref_box_inserted_1(TRUE)
      
      files_ref <- list.files(
        path       = reference_dir,
        pattern    = "\\.csv$",
        full.names = FALSE
      )
      stripped_ref <- sub("\\.csv$", "", files_ref)
      
      updateSelectInput(
        session,
        "reference_sel_1",
        choices  = stripped_ref,
        selected = NULL
      )
    }
  })
  
  observeEvent(input$load_reference_1, ignoreNULL = TRUE, {
    n_sel <- length(input$reference_sel_1)
    
    if (is.null(n_sel) || n_sel == 0) {
      showModal(
        modalDialog(
          title = "Selection required",
          "Please select at least one reference assay before proceeding.",
          easyClose = TRUE
        )
      )
      return(NULL)
    }
    
    session$sendCustomMessage("coralie-loading-text", "Loading reference assays")
    session$sendCustomMessage("coralie-loading-subtext", "")
    session$sendCustomMessage("coralie-loading-progress", 0)
    session$sendCustomMessage("coralie-toggle-loading", TRUE)
    
    sel_names <- input$reference_sel_1
    paths <- file.path(reference_dir, paste0(sel_names, ".csv"))
    
    sizes <- file.info(paths)$size
    sizes[is.na(sizes) | sizes <= 0] <- 1
    total_size <- sum(sizes)
    cum_size <- 0
    
    objs <- list()
    for (i in seq_along(paths)) {
      p  <- paths[i]
      nm <- sel_names[i]
      
      set_coralie_progress(
        session,
        100 * cum_size / total_size,
        sprintf("Loading reference assay %d/%d: %s",
                i, length(paths), nm)
      )
      
      df <- read.csv(p, stringsAsFactors = FALSE)
      objs[[nm]] <- df
      
      cum_size <- cum_size + sizes[i]
      set_coralie_progress(
        session,
        100 * cum_size / total_size
      )
    }
    
    names(objs) <- sel_names
    reference_list_1(objs)
    
    raw_ref_names <- names(objs)
    if (is.null(reference_display_names()))
      reference_display_names(raw_ref_names)
    if (is.null(reference_display_order()))
      reference_display_order(raw_ref_names)
    
    session$sendCustomMessage("coralie-toggle-loading", FALSE)
    session$sendCustomMessage("coralie-loading-subtext", "")
    
    ## Insert summary panel only once
    if (!summary_box_inserted_1()) {
      insertUI(
        selector = "#tier1_summary_placeholder",
        where    = "afterEnd",
        ui = fluidRow(
          class = "fade-in-up",
          column(
            width = 12,
            box(
              width = 12, status = "primary", solidHeader = TRUE,
              title = "Selection summary and labels",
              fluidRow(
                column(
                  width = 6,
                  h4("Fingerprints"),
                  orderInput(
                    inputId = "fingerprint_order",
                    label   = "Reorder names as desired",
                    items   = names(fingerprints_list_1()),
                    as_source = FALSE,
                    connect  = NULL
                  ),
                  uiOutput("fingerprint_rename_ui")
                ),
                column(
                  width = 6,
                  h4("Reference assays"),
                  orderInput(
                    inputId = "reference_order",
                    label   = "Reorder names as desired",
                    items   = names(reference_list_1()),
                    as_source = FALSE,
                    connect  = NULL
                  ),
                  uiOutput("reference_rename_ui")
                )
              ),
              br(),
              div(
                actionButton(
                  inputId = "begin_analysis_tier1",
                  label   = "Begin analysis"
                )
              )
            )
          )
        ),
        immediate = TRUE
      )
      summary_box_inserted_1(TRUE)
    } else {
      # Update orderInput lists if summary already exists
      shinyjqui::updateOrderInput(
        session, "fingerprint_order",
        items = names(fingerprints_list_1())
      )
      shinyjqui::updateOrderInput(
        session, "reference_order",
        items = names(reference_list_1())
      )
    }
  })
  
  observe({
    req(fingerprints_list_1())
    raw <- names(fingerprints_list_1())
    vals <- sapply(seq_along(raw), function(i) {
      id <- paste0("fp_name_", i)
      v  <- input[[id]]
      if (is.null(v) || v == "") raw[i] else v
    })
    fingerprint_display_names(vals)
  })
  
  observe({
    req(reference_list_1())
    raw <- names(reference_list_1())
    vals <- sapply(seq_along(raw), function(i) {
      id <- paste0("ref_name_", i)
      v  <- input[[id]]
      if (is.null(v) || v == "") raw[i] else v
    })
    reference_display_names(vals)
  })
  
  observeEvent(input$fingerprint_order, {
    req(fingerprints_list_1())
    # fingerprint_order should contain the raw names in desired order
    fingerprint_display_order(input$fingerprint_order)
  })
  
  observeEvent(input$reference_order, {
    req(reference_list_1())
    reference_display_order(input$reference_order)
  })
  
  fingerprints_for_display <- reactive({
    req(fingerprints_list_1())
    raw_list <- fingerprints_list_1()
    order    <- fingerprint_display_order()
    names_v  <- fingerprint_display_names()
    
    # keep object list intact, just create ordered view + labels
    ordered_list  <- raw_list[order]
    ordered_names <- names_v[match(names(ordered_list), names(raw_list))]
    list(objects = ordered_list, labels = ordered_names)
  })
  
  observeEvent(input$begin_analysis_tier1, {
    req(fingerprints_list_1(), reference_list_1())
    
    # reset state
    analysis_done_1(FALSE)
    tier1_cor_data(NULL)
    
    # loader
    session$sendCustomMessage("coralie-loading-text",
                              "Running CORALIE - Tier I Analysis")
    session$sendCustomMessage("coralie-loading-subtext", "")
    session$sendCustomMessage("coralie-loading-progress", 0)
    session$sendCustomMessage("coralie-toggle-loading", TRUE)
    
    # inputs
    fps_raw  <- fingerprints_list_1()
    refs_raw <- reference_list_1()
    
    fp_order <- fingerprint_display_order()
    fps <- if (!is.null(fp_order)) fps_raw[fp_order] else fps_raw
    
    ref_order <- reference_display_order()
    if (!is.null(ref_order)) {
      refs    <- refs_raw[ref_order]
      ref_ids <- ref_order
    } else {
      refs    <- refs_raw
      ref_ids <- names(refs_raw)
    }
    
    n_refs   <- length(ref_ids)
    fp_names <- names(fps)
    n_fp     <- length(fp_names)
    total_steps <- max(1, n_refs * n_fp)
    step <- 0L
    
    # containers
    # per reference: list(subsystem -> matrix(samples × fingerprints))
    ex_subsystems_all <- vector("list", length(ref_ids))
    names(ex_subsystems_all) <- ref_ids
    
    # also keep sample × fingerprint summary for the fingerprint corr heatmaps
    cor_tabs <- vector("list", length(ref_ids))
    names(cor_tabs) <- ref_ids
    
    # ---------- per reference assay ----------
    for (i in seq_along(ref_ids)) {
      ref_id <- ref_ids[i]
      assay  <- refs[[ref_id]]
      
      session$sendCustomMessage(
        "coralie-loading-subtext",
        sprintf("Reference %d/%d: %s", i, n_refs, ref_id)
      )
      
      # per fingerprint predictions by subsystem: fp -> data.frame(samples × subsystems)
      fp_subsystems <- list()
      
      for (j in seq_along(fp_names)) {
        fp_name    <- fp_names[j]
        fingerprint <- fps[[fp_name]]
        
        step <- step + 1L
        session$sendCustomMessage(
          "coralie-loading-subtext",
          sprintf("Reference %d/%d: %s – Fingerprint %d/%d: %s",
                  i, n_refs, ref_id, j, n_fp, fp_name)
        )
        session$sendCustomMessage(
          "coralie-loading-progress",
          round(100 * step / total_steps)
        )
        
        auc_list <- fingerprint$All$auc_primary
        rownames(auc_list) <- auc_list$index
        
        index        <- auc_list$index
        system_names <- auc_list[index, ]$subsystem
        final_cls    <- fingerprint$All$fingerprints_primary
        
        ex_subsystems_no_norm <- list()
        for (x in seq_along(index)) {
          idx        <- as.numeric(index[x])
          coef_names <- final_cls[[idx]]$coefnames
          
          newDat <- as.matrix(assay[, coef_names, drop = FALSE])
          colnames(newDat) <- coef_names
          
          pred_tot <- predict(
            final_cls[[idx]],
            newdata   = newDat,
            type      = "prob",
            na.action = na.pass
          )[, 1]
          
          ex_subsystems_no_norm[[system_names[x]]] <- pred_tot
        }
        
        ex_subsystems_no_norm <- as.data.frame(ex_subsystems_no_norm)
        fp_subsystems[[fp_name]] <- ex_subsystems_no_norm
      }
      
      # store ex_subsystems_no_norm per reference
      ex_subsystems_all[[ref_id]] <- fp_subsystems
      
      # sample × fingerprint “total score” for fingerprint correlation matrices
      total_scores <- lapply(fp_subsystems, function(df_sub) {
        score_raw <- rowMeans(df_sub, na.rm = TRUE)
        rng <- range(score_raw, na.rm = TRUE)
        if (diff(rng) == 0) {
          rep(0, length(score_raw))
        } else {
          (score_raw - rng[1]) / (rng[2] - rng[1])
        }
      })
      total_scores_df <- as.data.frame(total_scores)
      cor_tabs[[ref_id]] <- cor(total_scores_df, method = "pearson")
    }
    
    # ---------- subsystem-level average fingerprint correlations per reference ----------
    # For each reference and each subsystem:
    #   1) collect that subsystem’s predictions across fingerprints
    #   2) compute corr between fingerprints (samples × fingerprints)
    #   3) average off-diagonal correlations
    diag_stats_per_ref <- lapply(ref_ids, function(ref_id) {
      fp_subsystems <- ex_subsystems_all[[ref_id]]  # list(fp → df(samples × subsystems))
      if (is.null(fp_subsystems) || length(fp_subsystems) == 0) return(NULL)
      
      # set of all subsystems present in at least two fingerprints
      subs_by_fp <- lapply(fp_subsystems, colnames)
      all_subs   <- sort(unique(unlist(subs_by_fp)))
      
      # per subsystem: build sample × fingerprint matrix and compute mean pairwise corr
      sub_stats <- lapply(all_subs, function(sub) {
        # fingerprints that contain this subsystem
        fps_with_sub <- names(Filter(function(cols) sub %in% cols, subs_by_fp))
        if (length(fps_with_sub) < 2L) return(NULL)  # need at least 2 fingerprints
        
        # sample × fingerprint matrix for this subsystem
        mat <- do.call(
          cbind,
          lapply(fps_with_sub, function(fp_name) {
            fp_subsystems[[fp_name]][[sub]]
          })
        )
        colnames(mat) <- fps_with_sub
        
        # correlations between fingerprints for this subsystem
        corr <- suppressWarnings(cor(mat, use = "pairwise.complete.obs", method = "pearson"))
        if (!all(is.finite(corr))) {
          corr[!is.finite(corr)] <- NA_real_
        }
        
        # average off-diagonal
        off_diag <- corr[upper.tri(corr) | lower.tri(corr)]
        off_diag <- off_diag[is.finite(off_diag)]
        if (length(off_diag) == 0) return(NULL)
        
        data.frame(
          Subsystem        = sub,
          Average          = mean(off_diag),
          Count            = length(off_diag),
          Std_Dev          = if (length(off_diag) > 1) stats::sd(off_diag) else NA_real_,
          Reference_Assay  = ref_id,
          stringsAsFactors = FALSE
        )
      })
      
      sub_stats <- Filter(Negate(is.null), sub_stats)
      if (length(sub_stats) == 0) return(NULL)
      dplyr::bind_rows(sub_stats)
    })
    diag_stats_per_ref <- Filter(Negate(is.null), diag_stats_per_ref)
    
    if (length(diag_stats_per_ref) > 0) {
      diag_all <- dplyr::bind_rows(diag_stats_per_ref)
    } else {
      diag_all <- data.frame(
        Subsystem       = character(),
        Average         = numeric(),
        Count           = integer(),
        Std_Dev         = numeric(),
        Reference_Assay = character(),
        stringsAsFactors = FALSE
      )
    }
    
    # ---------- compare ref1 vs ref2 ----------
    if (length(ref_ids) >= 2 && nrow(diag_all) > 0) {
      ref1 <- ref_ids[1]
      ref2 <- ref_ids[2]
      
      diag_1 <- diag_all[diag_all$Reference_Assay == ref1, ]
      diag_2 <- diag_all[diag_all$Reference_Assay == ref2, ]
      
      merged_df <- merge(
        diag_1,
        diag_2,
        by = "Subsystem",
        suffixes = c("_1", "_2")
      )
      
      if (nrow(merged_df) > 0) {
        # normalize possible name variants to a standard set
        names(merged_df) <- sub("\\.x$", "_1", names(merged_df))
        names(merged_df) <- sub("\\.y$", "_2", names(merged_df))
        names(merged_df) <- sub("^Average$", "Average_1", names(merged_df)) # in case of no suffix
        names(merged_df) <- sub("^Average1$", "Average_1", names(merged_df))
        names(merged_df) <- sub("^Average2$", "Average_2", names(merged_df))
        names(merged_df) <- sub("^Std_Dev1$", "Std_Dev_1", names(merged_df))
        names(merged_df) <- sub("^Std_Dev2$", "Std_Dev_2", names(merged_df))
        names(merged_df) <- sub("^Count1$", "Count_1", names(merged_df))
        names(merged_df) <- sub("^Count2$", "Count_2", names(merged_df))
        
        A1 <- merged_df$Average_1
        A2 <- merged_df$Average_2
        S1 <- merged_df$Std_Dev_1
        S2 <- merged_df$Std_Dev_2
        N1 <- merged_df$Count_1
        N2 <- merged_df$Count_2
        
        merged_df$t_stat <- (A1 - A2) /
          sqrt((S1^2 / N1) + (S2^2 / N2))
        
        numerator <- (S1^2 / N1 + S2^2 / N2)^2
        denominator <- ((S1^2 / N1)^2 / pmax(N1 - 1, 1)) +
          ((S2^2 / N2)^2 / pmax(N2 - 1, 1))
        
        merged_df$df       <- numerator / denominator
        merged_df$p_value  <- 2 * stats::pt(-abs(merged_df$t_stat), df = merged_df$df)
        merged_df$abs_diff <- abs(A1 - A2)
        
        result_high_response <- subset(
          merged_df,
          p_value < 0.05 & t_stat > 0 &
            A1 > 0.2 & abs_diff > 0.05
        )
        result_low_response <- subset(
          merged_df,
          p_value < 0.05 & t_stat < 0 &
            A2 > 0.2 & abs_diff > 0.05
        )
        similar_high <- subset(
          merged_df,
          abs_diff < 0.05 & A1 > 0.2 & A2 > 0.2 & p_value > 0.05
        )
      } else {
        result_high_response <- merged_df[0, ]
        result_low_response  <- merged_df[0, ]
        similar_high         <- merged_df[0, ]
      }
    } else {
      result_high_response <- NULL
      result_low_response  <- NULL
      similar_high         <- NULL
    }
    
    # ---------- fingerprint-level correlation matrices (for existing plots) ----------
    cor_long <- do.call(
      rbind,
      lapply(ref_ids, function(ref_id) {
        mat <- cor_tabs[[ref_id]]
        as.data.frame(as.table(mat)) |>
          dplyr::rename(Fingerprint1 = Var1,
                        Fingerprint2 = Var2,
                        Correlation  = Freq) |>
          dplyr::mutate(Reference_Assay = ref_id)
      })
    )
    
    df_all <- do.call(
      rbind,
      lapply(ref_ids, function(ref_id) {
        cor_tab <- cor_tabs[[ref_id]]
        cor_df  <- as.data.frame(as.table(cor_tab))
        cor_df  <- cor_df[cor_df$Var1 != cor_df$Var2, ]
        cor_df  <- cor_df[!duplicated(t(apply(cor_df[, 1:2], 1, sort))), ]
        data.frame(
          Reference_Assay     = ref_id,
          Fingerprint_Pair    = paste(cor_df$Var1, cor_df$Var2, sep = "-"),
          Pearson_Correlation = cor_df$Freq,
          stringsAsFactors    = FALSE
        )
      })
    )
    
    # store results for downstream plots
    tier1_cor_data(list(
      ref_ids              = ref_ids,
      cor_tabs             = cor_tabs,
      df_all               = df_all,
      cor_long             = cor_long,
      diag_all             = diag_all,
      ex_subsystems_all    = ex_subsystems_all,
      result_high_response = result_high_response,
      result_low_response  = result_low_response,
      similar_high         = similar_high
    ))
    analysis_done_1(TRUE)
    
    # insert panels (unchanged)
    if (!analysis_panel_inserted_1()) {
      ref1_lab <- ref_ids[1]
      ref2_lab <- ref_ids[2]
      
      controls_col <- if (length(ref_ids) >= 2) {
        column(
          width = 3,
          radioButtons(
            "subsystem_bar_mode",
            "View:",
            choices = setNames(
              c("high", "low", "similar"),
              c(
                paste("Significantly higher in", ref1_lab),
                paste("Significantly higher in", ref2_lab),
                paste("Similar and high in", ref1_lab, "and", ref2_lab)
              )
            ),
            selected = "high"
          ),
          br(),
          selectInput(
            "subsystem_bar_sort",
            "Sort bars by:",
            choices = setNames(
              c("ref1", "ref2", "diff", "pval"),
              c(
                paste("Average in", ref1_lab),
                paste("Average in", ref2_lab),
                "Absolute difference",
                "Significance (p-value)"
              )
            ),
            selected = "ref1"
          )
        )
      } else {
        NULL
      }
      
      insertUI(
        selector = "#tier1_analysis_placeholder",
        where    = "afterEnd",
        ui = fluidRow(
          class = "fade-in-up",
          column(
            width = 12,
            box(
              width = 12, status = "primary", solidHeader = TRUE,
              title = "CORALIE - Tier I Correlation Plots",
              id    = "tier1_corrplot_grid_box", 
              plotly::plotlyOutput("tier1_corrplot_grid")
            ),
            box(
              width = 12, status = "primary", solidHeader = TRUE,
              title = "CORALIE - Tier I Changes in Correlation Across Reference Assays",
              plotly::plotlyOutput("tier1_corrplot_box")
            ),
            box(
              width = 12, status = "primary", solidHeader = TRUE,
              title = " CORALIE - Tier I Changes in Overall Subsystem Correlations Across Reference Assays",
              fluidRow(
                controls_col,
                column(
                  width = if (is.null(controls_col)) 12 else 9,
                  plotly::plotlyOutput("tier1_subsystem_bar", height = "70vh")
                )
              )
            )
          )
        ),
        immediate = TRUE
      )
      analysis_panel_inserted_1(TRUE)
    }
    
    session$sendCustomMessage("coralie-loading-progress", 100)
    session$sendCustomMessage("coralie-loading-subtext", "")
    session$sendCustomMessage("coralie-toggle-loading", FALSE)
  })
  
  
  
  tier1_corr_nrow <- reactiveVal(1)
  
  output$tier1_corrplot_grid <- plotly::renderPlotly({
    req(analysis_done_1())
    dat <- tier1_cor_data()
    req(dat)
    
    cor_long <- dat$cor_long
    
    p <- ggplot(cor_long, aes(
      x     = Fingerprint1,
      y     = Fingerprint2,
      fill  = Correlation,
      text  = paste0(
        "Ref: ", Reference_Assay,
        "<br>FP1: ", Fingerprint1,
        "<br>FP2: ", Fingerprint2,
        "<br>r = ", sprintf("%.3f", Correlation)
      )
    )) +
      geom_tile(color = "grey90") +
      scale_fill_gradient2(
        limits = c(-1, 1),
        low    = "#B2182B",
        mid    = "white",
        high   = "#2166AC"
      ) +
      facet_wrap(~ Reference_Assay) +
      coord_equal() +
      labs(x = NULL, y = NULL, fill = "r") +
      ggtitle("Correlation matrices of fingerprint predictions") +
      theme_minimal(base_size = 11) +
      theme(
        axis.text.x  = element_text(angle = 45, hjust = 1),
        panel.grid   = element_blank(),
        strip.text   = element_text(face = "bold"),
        plot.title   = element_text(hjust = 0.5, face = "bold")
      )
    
    plotly::ggplotly(
      p,
      tooltip = "text",
      source  = "tier1_corr_grid",
      height  = 350
    ) |>
      plotly::event_register("plotly_click") |>
      plotly::layout(
        legend = list(orientation = "h", x = 0.5, xanchor = "center"),
        margin = list(l = 60, r = 20, b = 80, t = 60)
      )
  })
  
  observeEvent(tier2_selection(), {
    sel <- tier2_selection()
    req(sel)
    
    if (!tier2_panel_inserted()) {
      insertUI(
        selector = "#tier1_analysis_placeholder",
        where    = "afterEnd",
        ui = fluidRow(
          class = "fade-in-up",
          id    = "tier2_panel_row",
          column(
            width = 12,
            box(
              width = 12, status = "primary", solidHeader = TRUE,
              title = "CORALIE - Tier II Correlation Analysis of Subsystem Models Between Two Fingerprints",
              # first row: full-width plot
              fluidRow(
                column(
                  width = 12,
                  plotlyOutput("tier2_subsystem_heatmap", height = "650px")
                )
              ),
              br(),
              # second row: controls below the plot
              fluidRow(
                column(
                  width = 12,
                  div(
                    style = "display: flex; gap: 20px; align-items: center; flex-wrap: wrap;",
                    div(
                      style = "min-width: 260px;",
                      div(
                        style = "min-width: 260px;",
                        sliderInput(
                          "tier2_corr_thresh",
                          "Correlation threshold:",
                          min   = -1,
                          max   =  1,
                          value =  0.5,
                          step  = 0.01
                        )
                      ),
                      div(
                        style = "min-width: 220px;",
                        radioButtons(
                          "tier2_corr_mode",
                          "Keep subsystems with:",
                          choices  = c(
                            "Any value ≥ threshold" = "gte",
                            "Any value ≤ threshold" = "lte"
                          ),
                          selected = "gte",
                          inline   = FALSE
                        )
                      )
                    ),
                    div(
                      style = "min-width: 220px;",
                      checkboxInput(
                        "tier2_run_ttest",
                        "Filter subsystems by t-test (p < 0.05)",
                        value = TRUE
                      )
                    ),
                    div(
                      style = "min-width: 260px;",
                      radioButtons(
                        "tier2_order_mode",
                        "Order features by:",
                        choices  = c(
                          "Name (alphabetical)"           = "name",
                          "Hierarchical clustering"       = "hclust",
                          "Optimize diagonal correlation" = "diag"
                        ),
                        selected = "hclust"
                      )
                    )
                  )
                )
              )
            )
          )
        ),
        immediate = TRUE
      )
      tier2_panel_inserted(TRUE)
    }
  })
  
  output$tier2_subsystem_heatmap <- plotly::renderPlotly({
    sel <- tier2_selection()
    req(sel)
    req(!is.null(sel$ref_id), !is.null(sel$fp1), !is.null(sel$fp2))
    
    ref_id <- sel$ref_id
    fp1    <- sel$fp1
    fp2    <- sel$fp2
    
    dat <- tier1_cor_data()
    req(dat, dat$ex_subsystems_all)
    ex_subsystems_all <- dat$ex_subsystems_all
    validate(need(ref_id %in% names(ex_subsystems_all),
                  "Selected reference assay not found in Tier I data."))
    subsystem_scores <- ex_subsystems_all[[ref_id]]
    validate(need(fp1 %in% names(subsystem_scores) && fp2 %in% names(subsystem_scores),
                  "Selected fingerprints not found in subsystem scores."))
    
    # containers from Tier I
    ex_subsystems_all <- dat$ex_subsystems_all   # list: ref -> list(fp -> df(samples × subsystems))
    ref_scores_list   <- dat$total_scores_all    # if you stored per-ref total_scores; if not, recompute below
    
    subsystem_scores <- ex_subsystems_all[[ref_id]]
    req(subsystem_scores[[fp1]], subsystem_scores[[fp2]])
    
    # get per-fingerprint scores (sample × 1) used in Tier I; recompute if not cached
    if (!is.null(ref_scores_list) && !is.null(ref_scores_list[[ref_id]])) {
      total_scores <- ref_scores_list[[ref_id]]
    } else {
      # fallback: mean of subsystem predictions per fingerprint
      total_scores <- lapply(subsystem_scores, function(df_sub) {
        score <- rowMeans(df_sub, na.rm = TRUE)
        rng <- range(score, na.rm = TRUE)
        if (diff(rng) == 0) rep(0, length(score)) else (score - rng[1]) / (rng[2] - rng[1])
      })
    }
    
    # Base correlation table between subsystems of fp1 and fp2
    cor_table_subsystem <- cor(subsystem_scores[[fp1]], subsystem_scores[[fp2]],
                               use = "pairwise.complete.obs", method = "pearson")
    cor_table_subsystem <- cor_table_subsystem[rowSums(!is.na(cor_table_subsystem)) > 0, , drop = FALSE]
    cor_table_subsystem <- cor_table_subsystem[, colSums(!is.na(cor_table_subsystem)) > 0, drop = FALSE]
    validate(need(nrow(cor_table_subsystem) > 0 && ncol(cor_table_subsystem) > 0,
                  "No overlapping subsystems between the selected fingerprints."))
    
    common_subs <- intersect(rownames(cor_table_subsystem), colnames(cor_table_subsystem))
    validate(need(length(common_subs) > 0,
                  "No common subsystems between the two fingerprints."))
    
    cor_table_subsystem <- cor_table_subsystem[common_subs, common_subs, drop = FALSE]
    
    # Optional t-test filtering
    if (isTRUE(input$tier2_run_ttest)) {
      # fingerprint 1
      out_1 <- ifelse(total_scores[[fp1]] > 0.5, 1, 0)
      sub_scores_1 <- data.frame(out = out_1, subsystem_scores[[fp1]])
      
      t_test_1 <- data.frame(
        sapply(
          names(subset(sub_scores_1, out == 1)[, -1, drop = FALSE]),
          function(col)
            t.test(
              subset(sub_scores_1, out == 1)[, -1][[col]],
              subset(sub_scores_1, out == 0)[, -1][[col]]
            )$p.value
        )
      )
      colnames(t_test_1)[1] <- "pval"
      fingerprint_features_1 <- rownames(subset(t_test_1, pval < 0.05))
      
      # fingerprint 2
      out_2 <- ifelse(total_scores[[fp2]] > 0.5, 1, 0)
      sub_scores_2 <- data.frame(out = out_2, subsystem_scores[[fp2]])
      
      t_test_2 <- data.frame(
        sapply(
          names(subset(sub_scores_2, out == 1)[, -1, drop = FALSE]),
          function(col)
            t.test(
              subset(sub_scores_2, out == 1)[, -1][[col]],
              subset(sub_scores_2, out == 0)[, -1][[col]]
            )$p.value
        )
      )
      colnames(t_test_2)[1] <- "pval"
      fingerprint_features_2 <- rownames(subset(t_test_2, pval < 0.05))
      
      # subset correlation table by significant subsystems in each fingerprint
      validate(need(length(fingerprint_features_1) > 0,
                    "No subsystems pass the t-test filter for the first fingerprint."))
      validate(need(length(fingerprint_features_2) > 0,
                    "No subsystems pass the t-test filter for the second fingerprint."))
      
      row_keep <- intersect(rownames(cor_table_subsystem), fingerprint_features_1)
      col_keep <- intersect(colnames(cor_table_subsystem), fingerprint_features_2)
      
      validate(need(length(row_keep) > 0 && length(col_keep) > 0,
                    "No overlapping significant subsystems between the two fingerprints."))
      
      cor_table_subsystem <- cor_table_subsystem[row_keep, col_keep, drop = FALSE]
    }
    
    # Apply correlation range filter from slider
    thresh <- input$tier2_corr_thresh %||% 0
    mode   <- input$tier2_corr_mode %||% "gte"
    
    row_min <- apply(cor_table_subsystem, 1, function(v) min(v, na.rm = TRUE))
    row_max <- apply(cor_table_subsystem, 1, function(v) max(v, na.rm = TRUE))
    col_min <- apply(cor_table_subsystem, 2, function(v) min(v, na.rm = TRUE))
    col_max <- apply(cor_table_subsystem, 2, function(v) max(v, na.rm = TRUE))
    
    if (mode == "gte") {
      row_keep <- row_max >= thresh
      col_keep <- col_max >= thresh
    } else { # "lte"
      row_keep <- row_min <= thresh
      col_keep <- col_min <= thresh
    }
    
    
    cor_display <- cor_table_subsystem[row_keep, col_keep, drop = FALSE]
    validate(need(nrow(cor_display) > 0 && ncol(cor_display) > 0,
                  "No subsystems meet the correlation threshold criteria."))
    
    # Drop rows/cols that are completely NA
    row_keep <- rowSums(!is.na(cor_display)) > 0
    col_keep <- colSums(!is.na(cor_display)) > 0
    cor_display <- cor_display[row_keep, col_keep, drop = FALSE]
    validate(need(nrow(cor_display) > 0 && ncol(cor_display) > 0,
                  "No subsystem pairs lie in the selected correlation range."))
    
    # Zoom‑dependent subsetting (uses ranges stored from relayout events)
    mat_full <- cor_display
    
    relayout <- plotly::event_data("plotly_relayout", source = "tier2_subsystem")
    
    if (!is.null(relayout) &&
        !is.null(relayout[["xaxis.range[0]"]]) &&
        !is.null(relayout[["xaxis.range[1]"]]) &&
        !is.null(relayout[["yaxis.range[0]"]]) &&
        !is.null(relayout[["yaxis.range[1]"]])) {
      
      x0 <- relayout[["xaxis.range[0]"]]
      x1 <- relayout[["xaxis.range[1]"]]
      y0 <- relayout[["yaxis.range[0]"]]
      y1 <- relayout[["yaxis.range[1]"]]
      
      # convert axis values to factor positions using the *current* labels
      x_levels <- colnames(mat_full)
      y_levels <- rownames(mat_full)
      
      # when you build plot_ly, x = x_levels, y = y_levels, so these are the same
      x_idx <- which(x_levels >= x0 & x_levels <= x1)
      y_idx <- which(y_levels >= y0 & y_levels <= y1)
      
      if (length(x_idx) > 0 && length(y_idx) > 0) {
        mat_zoom <- mat_full[y_idx, x_idx, drop = FALSE]
      } else {
        mat_zoom <- mat_full
      }
    } else {
      mat_zoom <- mat_full
    }
    
    # Ordering within the zoomed window
    order_mode <- input$tier2_order_mode %||% "hclust"
    
    mat_show <- mat_zoom
    
    x_vals <- colnames(mat_show)
    y_vals <- rownames(mat_show)
    
    x_range <- NULL
    y_range <- NULL
    if (!is.null(relayout)) {
      if (!is.null(relayout[["xaxis.range[0]"]]) && !is.null(relayout[["xaxis.range[1]"]])) {
        x_range <- c(relayout[["xaxis.range[0]"]], relayout[["xaxis.range[1]"]])
      }
      if (!is.null(relayout[["yaxis.range[0]"]]) && !is.null(relayout[["yaxis.range[1]"]])) {
        y_range <- c(relayout[["yaxis.range[0]"]], relayout[["yaxis.range[1]"]])
      }
    }
    
    if (order_mode == "name") {
      mat_show <- mat_show[
        order(rownames(mat_show)),
        order(colnames(mat_show)),
        drop = FALSE
      ]
    } else if (order_mode == "hclust") {
      if (nrow(mat_show) > 1) {
        row_ord <- hclust(dist(mat_show))$order
        mat_show <- mat_show[row_ord, , drop = FALSE]
      }
      if (ncol(mat_show) > 1) {
        col_ord <- hclust(dist(t(mat_show)))$order
        mat_show <- mat_show[, col_ord, drop = FALSE]
      }
    } else if (order_mode == "diag") {
      mat_show <- mat_zoom
      abs_cor <- abs(mat_show)
      
      nr <- nrow(abs_cor)
      nc <- ncol(abs_cor)
      
      # pad to square if needed
      if (nr > nc) {
        pad <- matrix(0, nrow = nr, ncol = nr - nc)
        abs_sq <- cbind(abs_cor, pad)
      } else if (nc > nr) {
        pad <- matrix(0, nrow = nc - nr, ncol = nc)
        abs_sq <- rbind(abs_cor, pad)
      } else {
        abs_sq <- abs_cor
      }
      
      assign_sq <- clue::solve_LSAP(abs_sq, maximum = TRUE)
      assign_idx <- as.integer(assign_sq)
      
      if (nr <= nc) {
        col_ord <- assign_idx[seq_len(nr)]
        row_ord <- seq_len(nr)
      } else {
        row_ord <- assign_idx[seq_len(nc)]
        col_ord <- seq_len(nc)
      }
      
      row_ord <- row_ord[row_ord <= nrow(mat_show)]
      col_ord <- col_ord[col_ord <= ncol(mat_show)]
      mat_show <- mat_show[row_ord, col_ord, drop = FALSE]
      
      # NEW: keep only diagonal where row and column names match
      common_names <- intersect(rownames(mat_show), colnames(mat_show))
      validate(need(length(common_names) > 0,
                    "No common subsystems to display on the diagonal."))
      
      mat_show <- mat_show[common_names, common_names, drop = FALSE]
    }
    
    ref_assay_name <- ref_id
    
    # Build Plotly heatmap with a source for relayout tracking
    p_plotly <- plotly::plot_ly(
      x      = colnames(mat_show),
      y      = rownames(mat_show),
      z      = mat_show,
      type   = "heatmap",
      colors = colorRamp(c("#B2182B", "white", "#2166AC")),
      zmin   = -1,
      zmax   =  1,
      source = "tier2_subsystem"
    ) |>
      plotly::colorbar(title = "r") |>
      plotly::layout(
        title = paste0(
          "Correlation of Effects Among Fingerprint Subsystem Features<br>",
          fp1, " vs. ", fp2,
          "<br>Reference Assay: ", ref_assay_name
        ),
        xaxis = list(
          title = fp2,
          range = if (!is.null(x_range)) x_range else NULL
        ),
        yaxis = list(
          title = fp1,
          range = if (!is.null(y_range)) y_range else NULL
        ),
        margin = list(l = 80, r = 20, b = 80, t = 80)
      )
    
    p_plotly
  })
  
  output$tier1_corrplot_box <- plotly::renderPlotly({
    req(analysis_done_1())
    dat    <- tier1_cor_data()
    req(dat)
    
    df_all  <- dat$df_all
    ref_ids <- dat$ref_ids
    
    df_all$row_id <- as.integer(as.factor(df_all$Fingerprint_Pair))
    df_all$Reference_Assay <- factor(df_all$Reference_Assay, levels = ref_ids)
    
    # ---- single-reference: simple boxplot, no pairwise stats ----
    if (length(ref_ids) < 2L) {
      p <- ggplot(df_all, aes(x = Reference_Assay, y = Pearson_Correlation)) +
        geom_boxplot(
          outlier.shape = NA,
          alpha         = 0.4,
          width         = 0.4,
          color         = "black",
          fill          = "#3182BD"
        ) +
        geom_jitter(
          width  = 0.1,
          size   = 1.6,
          alpha  = 0.6,
          color  = "#636363"
        ) +
        ggtitle("Correlation of fingerprint models within reference assay") +
        labs(x = NULL, y = "Pearson correlation") +
        theme_minimal(base_size = 12) +
        theme(
          plot.title       = element_text(face = "bold", hjust = 0.5),
          axis.text.x      = element_text(angle = 45, hjust = 1),
          panel.grid.minor = element_blank()
        )
      
      return(plotly::ggplotly(p, tooltip = c("x", "y")))
    }
    
    # ---- multi-reference: existing pairwise comparison logic ----
    base_ref <- ref_ids[1]
    
    pair_stats <- lapply(ref_ids[-1], function(ref2) {
      sub <- df_all[df_all$Reference_Assay %in% c(base_ref, ref2), ]
      pval <- tryCatch(
        stats::t.test(
          Pearson_Correlation ~ Reference_Assay,
          data = sub
        )$p.value,
        error = function(e) NA_real_
      )
      data.frame(
        group1 = base_ref,
        group2 = ref2,
        pval   = pval,
        stringsAsFactors = FALSE
      )
    })
    pair_stats <- do.call(rbind, pair_stats)
    
    pair_stats$p_signif <- dplyr::case_when(
      is.na(pair_stats$pval)  ~ "",
      pair_stats$pval < 0.001 ~ "***",
      pair_stats$pval < 0.01  ~ "**",
      pair_stats$pval < 0.05  ~ "*",
      TRUE                    ~ "ns"
    )
    
    y_base <- max(df_all$Pearson_Correlation, na.rm = TRUE)
    y_step <- 0.05 * diff(range(df_all$Pearson_Correlation, na.rm = TRUE))
    if (!is.finite(y_step) || y_step == 0) y_step <- 0.05
    
    pair_stats$y <- y_base + seq_len(nrow(pair_stats)) * y_step
    
    pair_stats$x1   <- as.numeric(factor(pair_stats$group1, levels = ref_ids))
    pair_stats$x2   <- as.numeric(factor(pair_stats$group2, levels = ref_ids))
    pair_stats$xmid <- (pair_stats$x1 + pair_stats$x2) / 2
    
    p <- ggplot(df_all, aes(x = Reference_Assay, y = Pearson_Correlation)) +
      geom_boxplot(
        outlier.shape = NA,
        alpha         = 0.4,
        width         = 0.35,
        color         = "black",
        fill          = "#3182BD"
      ) +
      geom_point(
        aes(group = row_id),
        position = position_jitter(width = 0.07),
        size     = 1.6,
        alpha    = 0.6,
        color    = "#636363"
      ) +
      geom_line(
        aes(group = row_id),
        color     = "#9E9E9E",
        alpha     = 0.4,
        linewidth = 0.25
      ) +
      geom_segment(
        data = pair_stats,
        aes(
          x    = x1,
          xend = x2,
          y    = y,
          yend = y
        ),
        inherit.aes = FALSE,
        linewidth   = 0.4,
        color       = "#444444"
      ) +
      geom_segment(
        data = pair_stats,
        aes(
          x    = x1,
          xend = x1,
          y    = y,
          yend = y - 0.015
        ),
        inherit.aes = FALSE,
        linewidth   = 0.4,
        color       = "#444444"
      ) +
      geom_segment(
        data = pair_stats,
        aes(
          x    = x2,
          xend = x2,
          y    = y,
          yend = y - 0.015
        ),
        inherit.aes = FALSE,
        linewidth   = 0.4,
        color       = "#444444"
      ) +
      geom_text(
        data = pair_stats,
        aes(
          x     = xmid,
          y     = y + 0.01,
          label = p_signif
        ),
        inherit.aes = FALSE,
        size        = 3,
        vjust       = 0
      ) +
      ggtitle("Changes in Correlations of Fingerprint Models Across Reference Assays") +
      labs(x = NULL, y = "Pearson correlation") +
      theme_minimal(base_size = 12) +
      theme(
        plot.title       = element_text(face = "bold", hjust = 0.5),
        axis.text.x      = element_text(angle = 45, hjust = 1),
        panel.grid.minor = element_blank()
      )
    
    plotly::ggplotly(p, tooltip = c("x", "y"))
  })
  
  observeEvent(plotly::event_data("plotly_click", source = "tier1_corr_grid"), {
    d <- plotly::event_data("plotly_click", source = "tier1_corr_grid")
    if (is.null(d) || nrow(d) == 0) return()

    dat <- tier1_cor_data()
    req(dat)
    cor_long <- dat$cor_long
    
    # ensure factors with known level order (same as in the plot)
    cor_long$Fingerprint1     <- factor(cor_long$Fingerprint1)
    cor_long$Fingerprint2     <- factor(cor_long$Fingerprint2)
    cor_long$Reference_Assay  <- factor(cor_long$Reference_Assay)
    
    # x/y are numeric positions along the factor axes
    x_idx <- round(d$x[1])
    y_idx <- round(d$y[1])
    ref_idx <- round(d$curveNumber[1])+1
    
    fp1_levels <- levels(cor_long$Fingerprint1)
    fp2_levels <- levels(cor_long$Fingerprint2)
    ref_levels <- levels(cor_long$Reference_Assay)
    
    # safety check
    if (x_idx < 1 || x_idx > length(fp1_levels)) return()
    if (y_idx < 1 || y_idx > length(fp2_levels)) return()
    
    fp1 <- fp1_levels[x_idx]
    fp2 <- fp2_levels[y_idx]
    ref <- ref_levels[ref_idx]
    
    # Fallback: if only one reference assay, use it
    if (is.null(ref) && length(levels(cor_long$Reference_Assay)) == 1L) {
      ref <- levels(cor_long$Reference_Assay)[1]
    }

    # Debug once
    # print(list(raw = d, fp1 = fp1, fp2 = fp2, ref = ref))
    
    if (is.null(ref) || is.null(fp1) || is.null(fp2) || fp1 == fp2) return()
    
    tier2_selection(list(
      ref_id = ref,
      fp1    = fp1,
      fp2    = fp2
    ))
  })
  
  
  observeEvent(plotly::event_data("plotly_relayout", source = "tier2_subsystem"), {
    ev <- plotly::event_data("plotly_relayout", source = "tier2_subsystem")
    if (is.null(ev)) return()
    
    # For heatmaps, ranges often appear as xaxis.range[0], xaxis.range[1], etc.
    xr <- ev[grep("^xaxis\\.range\\[", names(ev))]
    yr <- ev[grep("^yaxis\\.range\\[", names(ev))]
    
    if (length(xr) == 2) tier2_x_range(sort(unlist(xr)))
    if (length(yr) == 2) tier2_y_range(sort(unlist(yr)))
  })
  
  output$tier1_subsystem_bar <- plotly::renderPlotly({
    req(analysis_done_1())
    dat <- tier1_cor_data()
    req(dat)
    
    ref_ids  <- dat$ref_ids
    diag_all <- dat$diag_all
    
    # ---- single-reference branch ----
    if (length(ref_ids) == 1L || is.null(dat$result_high_response)) {
      validate(
        need(nrow(diag_all) > 0, "No subsystem correlation data to display.")
      )
      
      diag_single <- diag_all[diag_all$Reference_Assay == ref_ids[1], ]
      diag_single <- diag_single[order(-diag_single$Average), ]
      diag_single$Subsystem <- factor(diag_single$Subsystem, diag_single$Subsystem)
      
      p <- ggplot(diag_single, aes(x = Subsystem, y = Average, key = Subsystem)) +
        geom_col(
          fill = "#5C6BC0",
          width = 0.7
        ) +
        geom_hline(yintercept = 0, color = "grey85") +
        theme_minimal(base_size = 12) +
        theme(
          axis.text.x  = element_text(angle = 60, hjust = 1, vjust = 1, size = 10),
          axis.text.y  = element_text(size = 10, face = "bold"),
          axis.title.y = element_text(size = 11),
          panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
          panel.grid.minor = element_blank(),
          plot.title  = element_text(hjust = 0.5, face = "bold"),
          plot.margin = margin(t = 10, r = 10, b = 40, l = 80)
        ) +
        labs(
          x = NULL,
          y = "Average correlation",
          title = paste(
            "Average subsystem correlations for reference assay:",
            ref_ids[1]
          )
        )
      
      return(
        plotly::ggplotly(p, tooltip = c("y", "x","key"), source  = "tier1_subsystem_bar" ) |>
          plotly::layout(
            height = 550,
            margin = list(l = 80, r = 20, t = 60, b = 80)
          )
      )
    }
    
    # ---- multi-reference branch ----
    mode      <- input$subsystem_bar_mode
    sort_mode <- input$subsystem_bar_sort
    if (is.null(sort_mode) || length(sort_mode) != 1L) {
      sort_mode <- "ref1"
    }
    
    df <- switch(
      mode,
      high    = dat$result_high_response,
      low     = dat$result_low_response,
      similar = dat$similar_high,
      dat$result_high_response
    )
    validate(
      need(!is.null(df) && nrow(df) > 0,
           "No subsystems meeting the criteria for this view.")
    )
    
    # normalize names
    names(df) <- sub("\\.x$", "_1", names(df))
    names(df) <- sub("\\.y$", "_2", names(df))
    names(df) <- sub("^Average1$", "Average_1", names(df))
    names(df) <- sub("^Average2$", "Average_2", names(df))
    names(df) <- sub("^Std_Dev1$", "Std_Dev_1", names(df))
    names(df) <- sub("^Std_Dev2$", "Std_Dev_2", names(df))
    
    high_lab <- ref_ids[1]
    low_lab  <- ref_ids[2]
    
    # choose ordering metric based on user selection
    
    df <- df %>%
      dplyr::mutate(
        diff_value = Average_1 - Average_2,
        sig_score  = dplyr::case_when(
          is.na(p_value) ~ 0,
          TRUE           ~ -log10(p_value + 1e-16)
        )
      )
    
    if (sort_mode == "ref1") {
      df$order_value <- df$Average_1
    } else if (sort_mode == "ref2") {
      df$order_value <- df$Average_2
    } else if (sort_mode == "diff") {
      df$order_value <- abs(df$diff_value)
    } else if (sort_mode == "pval") {
      df$order_value <- df$sig_score
    } else {
      df$order_value <- df$Average_1
    }
    
    df <- df %>%
      dplyr::arrange(dplyr::desc(order_value))
    
    # highest at top
    df$Subsystem <- factor(df$Subsystem, levels = rev(unique(df$Subsystem)))
    
    metabolic_data <- df
    
    plot_data <- metabolic_data %>%
      tidyr::pivot_longer(
        cols = c(Average_1, Average_2),
        names_to  = "response_type",
        values_to = "average_response"
      ) %>%
      dplyr::mutate(
        std_dev = dplyr::if_else(
          response_type == "Average_1",
          Std_Dev_1,
          Std_Dev_2
        ),
        response_type = dplyr::recode(
          response_type,
          "Average_1" = high_lab,
          "Average_2" = low_lab
        )
      )
    
    # significance stars per subsystem (from df$p_value already computed upstream)
    sig_df <- metabolic_data %>%
      dplyr::mutate(
        p_signif = dplyr::case_when(
          is.na(p_value)  ~ "",
          p_value < 0.001 ~ "***",
          p_value < 0.01  ~ "**",
          p_value < 0.05  ~ "*",
          TRUE            ~ "ns"
        ),
        x_star = 1.05     # slightly beyond the max x (since your scale is roughly -1..1)
      )
    
    # standard bar aesthetic
    p <- ggplot(
      plot_data,
      aes(
        y    = Subsystem,
        x    = average_response,
        fill = response_type,
        key = Subsystem 
      )
    ) +
      geom_col(
        position = position_dodge(width = 0.7),
        width    = 0.6,
        color    = NA
      ) +
      geom_errorbar(
        aes(
          xmin = average_response - std_dev,
          xmax = average_response + std_dev
        ),
        position = position_dodge(width = 0.7),
        width    = 0.1,
        color    = "grey40",
        linewidth = 0.01
      ) +
      geom_vline(xintercept = 0, color = "grey90") +
      scale_fill_manual(
        values = setNames(
          c("#42A5F5", "#EF5350"),
          c(high_lab, low_lab)
        ),
        breaks = c(high_lab, low_lab),
        labels = c(high_lab, low_lab),
        name   = "Reference assay"
      ) +
      scale_x_continuous(expand = expansion(mult = c(0.01, 0.15))) +
      # significance stars, centered between the two bars
      geom_text(
        data = sig_df,
        aes(
          x     = x_star,
          y     = Subsystem,
          label = p_signif
        ),
        inherit.aes = FALSE,
        hjust       = 0,     # left-align text at x_star
        vjust       = 0.5,
        size        = 3.3
      ) +
      theme_minimal(base_size = 12) +
      theme(
        axis.text.y = element_text(hjust = 1, face = "bold", size = 10),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 11),
        panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
        panel.grid.major.x = element_line(color = "grey95", linewidth = 0.3),
        panel.grid.minor = element_blank(),
        legend.position = "top",
        legend.title = element_text(size = 11),
        legend.text  = element_text(size = 10),
        plot.title   = element_text(hjust = 0.5, face = "bold"),
        plot.margin  = margin(t = 10, r = 20, b = 40, l = 120)
      ) +
      labs(
        title = dplyr::case_when(
          mode == "high"    ~ paste("Subsystems higher in", high_lab),
          mode == "low"     ~ paste("Subsystems higher in", low_lab),
          mode == "similar" ~ paste("Subsystems similarly high in", high_lab, "and", low_lab),
          TRUE              ~ "Subsystem correlations"
        ),
        x = "Average Pearson correlation",
        y = "Subsystem"
      )
    
    plotly::ggplotly(
      p,
      tooltip = c("y", "x", "fill","key"),
      source = "tier1_subsystem_bar"
    ) |>
      plotly::layout(
        height = 650,
        margin = list(l = 120, r = 20, t = 70, b = 40)
      )
  })
  
  observeEvent(
    plotly::event_data("plotly_click", source = "tier1_subsystem_bar"),
    {
      ev <- plotly::event_data("plotly_click", source = "tier1_subsystem_bar")
      if (is.null(ev) || is.null(ev$key)) return()
      clicked_subsystem(as.character(ev$key))
    }
  )
  
  observeEvent(clicked_subsystem(), {
    req(clicked_subsystem())
    if (!tier3_box_inserted()) {
      insertUI(
        selector = "#tier3_analysis_placeholder",
        where    = "afterEnd",
        ui = fluidRow(
          class = "fade-in-up",
          column(
            width = 12,
            box(
              width = 12,
              status = "primary",
              solidHeader = TRUE,
              title = "CORALIE - Tier III Analysis of Reaction Activities Across Fingerprints and Reference Assays",
              plotly::plotlyOutput("tier3_reaction_plot", height = 600)
            )
          )
        )
      )
      tier3_box_inserted(TRUE)
    }
  })
  
  output$tier3_reaction_plot <- plotly::renderPlotly({
    req(analysis_done_1())
    req(clicked_subsystem())
    
    # session$sendCustomMessage("coralie-loading-text", "Running CORALIE - Tier III Analysis")
    # session$sendCustomMessage("coralie-toggle-loading", TRUE)
    
    dat <- tier1_cor_data()
    req(dat)
    
    ref_ids          <- dat$ref_ids
    ex_subsystems_all <- dat$ex_subsystems_all
    fps              <- fingerprints_list_1()
    req(fps)
    
    n_refs <- length(ref_ids)
    validate(
      need(n_refs >= 1, "Tier III: at least one reference assay is required.")
    )
    
    subsystem_name   <- clicked_subsystem()
    ref_assay_name_1 <- ref_ids[1]
    ref_assay_name_2 <- if (n_refs >= 2) ref_ids[2] else NULL
    
    # ------------------------------------------------------------------
    # 1) Collect reaction features used by models for this subsystem
    # ------------------------------------------------------------------
    react_list <- character(0)
    
    for (i in seq_along(fps)) {
      fingerprint <- fps[[i]]
      auc_list    <- fingerprint$All$auc_primary
      auc_list$subsystem <- make.names(auc_list$subsystem)
      
      if (subsystem_name %in% auc_list$subsystem) {
        model_num <- as.numeric(
          unlist(subset(auc_list, subsystem == subsystem_name)$index)
        )
        model     <- fingerprint$All$fingerprints_primary[[model_num]]
        features  <- model$coefnames
        react_list <- c(react_list, features)
      }
    }
    
    react_list <- unique(react_list)
    validate(
      need(length(react_list) > 0, "Tier III: no reaction features found for this subsystem.")
    )
    
    # ------------------------------------------------------------------
    # 2) Helper to build reaction prediction matrix for one reference
    # ------------------------------------------------------------------
    
    react_meta_filt$subsystem <- make.names(react_meta_filt$subsystem)
    
    build_react_matrix <- function(ref_id) {
      assay    <- reference_list_1()[[ref_id]]      # same data frame used in Tier I
      fp_names <- names(fps)
      
      react_pred_mat <- NULL
      sample_ids     <- assay[,1]
      
      for (fp_name in fp_names) {
        fingerprint <- fps[[fp_name]]
        auc_list    <- fingerprint$All$auc_primary
        auc_list$subsystem <- make.names(auc_list$subsystem)
        
        if (subsystem_name %in% auc_list$subsystem) {
          model_num <- as.numeric(
            unlist(subset(auc_list, subsystem == subsystem_name)$index)
          )
          model     <- fingerprint$All$fingerprints_primary[[model_num]]
          feats     <- intersect(model$coefnames, colnames(assay))
          if (length(feats) == 0) next
          
          new_dat <- as.matrix(assay[, feats, drop = FALSE])
          colnames(new_dat) <- feats
          
          pred <- predict(
            model,
            newdata   = new_dat,
            type      = "prob",
            na.action = na.pass
          )[ , 1]
          
          names(pred) <- sample_ids
          
          if (is.null(react_pred_mat)) {
            react_pred_mat <- cbind(pred)
            colnames(react_pred_mat) <- fp_name
          } else {
            react_pred_mat <- cbind(react_pred_mat, pred)
            colnames(react_pred_mat)[ncol(react_pred_mat)] <- fp_name
          }
        }
      }
      
      react_pred_mat
    }
    
    get_subsystem_react_matrix <- function(ref_id, subsystem_name) {
      ref_mat <- reference_list_1()[[ref_id]]  # samples × reactions
      validate(
        need(!is.null(ref_mat),
             paste("Tier III: no reaction matrix found for reference assay", ref_id))
      )
      
      sample_ids     <- ref_mat[,1]
      
      # Reactions in this subsystem (after make.names)
      subs_in_meta <- rownames(react_meta_filt)[
        react_meta_filt$subsystem == subsystem_name
      ]
      subs_in_meta <- intersect(subs_in_meta, colnames(ref_mat))
      
      validate(
        need(length(subs_in_meta) > 0,
             paste("Tier III: no reactions in react_meta_filt for subsystem", subsystem_name))
      )
      
      ref_sub <-ref_mat[, subs_in_meta, drop = FALSE]   # samples × reactions (raw values)
      
      keep_cols <- apply(ref_sub, 2, function(x) all(is.finite(x)))
      ref_sub   <- ref_sub[, keep_cols, drop = FALSE]
      
      validate(
        need(ncol(ref_sub) > 0,
             paste("Tier III: all reactions for subsystem", subsystem_name,
                   "contain missing or non-finite values."))
      )
      
      rownames(ref_sub) <- sample_ids
      ref_sub
    }
    
    # ------------------------------------------------------------------
    # 3) Correlations for reference assay 1 (always present)
    # ------------------------------------------------------------------
    react_pred_ref_1 <- build_react_matrix(ref_assay_name_1)
    validate(
      need(!is.null(react_pred_ref_1),
           "Tier III: no reaction predictions available for this subsystem in the first reference assay.")
    )
    
    subs_vec_ref_1 <- get_subsystem_react_matrix(ref_assay_name_1, subsystem_name)
    
    
    common_rows<-intersect(rownames(react_pred_ref_1),rownames(subs_vec_ref_1))
    
    react_pred_ref_1<-react_pred_ref_1[common_rows,]
    subs_vec_ref_1<-subs_vec_ref_1[common_rows,]
    
    corr_ref_1 <- cor(react_pred_ref_1, subs_vec_ref_1, use = "pairwise.complete.obs")
    react_cors_df_1 <- as.data.frame(corr_ref_1)
    
    keep_cols <- apply(react_cors_df_1, 2, function(x) all(is.finite(x)))
    react_cors_df_1   <- react_cors_df_1[, keep_cols, drop = FALSE]
    
    colnames(react_cors_df_1)<-react_meta_filt[colnames(react_cors_df_1),"RECON3D"]
    
    validate(
      need(nrow(react_cors_df_1) > 0,
           "Tier III: no valid correlations for the first reference assay.")
    )
    
    # ------------------------------------------------------------------
    # 4) If a second reference assay exists, compute its correlations
    # ------------------------------------------------------------------
    has_second_ref <- n_refs >= 2 && !is.null(ref_assay_name_2)
    
    if (has_second_ref) {
      react_pred_ref_2 <- build_react_matrix(ref_assay_name_2)
      validate(
        need(!is.null(react_pred_ref_2),
             "Tier III: no reaction predictions available for this subsystem in the second reference assay.")
      )
      
      sub_scores_ref_2 <- ex_subsystems_all[[ref_assay_name_2]]
      fp_names_sub_2   <- names(sub_scores_ref_2)
      validate(
        need(length(fp_names_sub_2) > 0,
             "Tier III: no subsystem scores found for the second reference assay.")
      )
      
      fp_target_2 <- fp_names_sub_2[1]
      subs_vec_ref_2 <- get_subsystem_react_matrix(ref_assay_name_2, subsystem_name)
      
      
      common_rows<-intersect(rownames(react_pred_ref_2),rownames(subs_vec_ref_2))
      
      react_pred_ref_2<-react_pred_ref_2[common_rows,]
      subs_vec_ref_2<-subs_vec_ref_2[common_rows,]
      
      corr_ref_2 <- cor(react_pred_ref_2, subs_vec_ref_2, use = "pairwise.complete.obs")
      react_cors_df_2 <- as.data.frame(corr_ref_2)
      
      keep_cols <- apply(react_cors_df_2, 2, function(x) all(is.finite(x)))
      react_cors_df_2   <- react_cors_df_2[, keep_cols, drop = FALSE]
      
      colnames(react_cors_df_2)<-react_meta_filt[colnames(react_cors_df_2),"RECON3D"]
      
      # Align reactions between ref 1 and ref 2
      common_cols <- intersect(colnames(react_cors_df_1), colnames(react_cors_df_2))
      
      react_cors_df_1 <- react_cors_df_1[, common_cols, drop = FALSE]
      react_cors_df_2 <- react_cors_df_2[, common_cols, drop = FALSE]
      
      
      validate(
        need(nrow(react_cors_df_1) > 0,
             "Tier III: no common reactions between the two reference assays.")
      )
    }
    # ------------------------------------------------------------------
    # 5) Build heatmaply plot(s)
    # ------------------------------------------------------------------
    
    # Decide whether to cluster rows
    do_cluster_cols <- TRUE
    if (ncol(react_cors_df_1) <= 2) {
      do_cluster_cols <- FALSE
    } else {
      # extra safety: check that dist() has no non-finite values
      d_test <- try(stats::dist(t(as.matrix(react_cors_df_1))), silent = TRUE)
      if (inherits(d_test, "try-error") || any(!is.finite(d_test))) {
        do_cluster_cols <- FALSE
      }
    }
    
    # always fine to cluster rows for this matrix; guard similarly if needed
    do_cluster_rows <- TRUE
    
    # Use ref 1 to derive row/column ordering
    hm_tmp <- heatmaply::heatmaply(
      react_cors_df_1,
      scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(
        low      = "red",
        high     = "blue",
        midpoint = 0,
        limits   = range(react_cors_df_1, na.rm = TRUE)
      ),
      main        = "Correlation ordering helper",
      Rowv        = do_cluster_rows,
      Colv        = do_cluster_cols,
      plot_method = "plotly"
    )
    
    
    rows    <- rev(hm_tmp$x$layout$yaxis2$ticktext)
    columns <- hm_tmp$x$layout$xaxis$ticktext
    
    if(is.null(rows)){
      p1 <- hm_tmp
      
      p1 <- plotly::layout(
        p1,
        yaxis = list(tickfont = list(size = 10, color = "black"))
      )
      
      p1$x$data[[1]]$colorscale <- list(
        list(0,  "red"),   # z = -1
        list(0.5, "white"),# z = 0
        list(1,  "blue")   # z = 1
      )
      p1$x$data[[1]]$zmin <- -1
      p1$x$data[[1]]$zmax <- 1
    }else{
      p1 <- heatmaply::heatmaply(
        react_cors_df_1[rows, columns, drop = FALSE],
        scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(
          low      = "red",
          high     = "blue",
          midpoint = 0,
          limits   = c(-1, 1)
        ),
        main = paste0(
          "<span style='font-size:12px;'>Correlation of Subsystem Reaction Activity and Fingerprint Prediction - ",
          subsystem_name, "<br>"
        ),
        Rowv        = FALSE,
        Colv        = FALSE,
        margins     = c(50, 50, 130, 50),
        plot_method = "plotly"
      )
      
      p1 <- plotly::layout(
        p1,
        yaxis = list(tickfont = list(size = 10, color = "black"))
      )
      
      p1$x$data[[1]]$colorscale <- list(
        list(0,  "red"),   # z = -1
        list(0.5, "white"),# z = 0
        list(1,  "blue")   # z = 1
      )
      p1$x$data[[1]]$zmin <- -1
      p1$x$data[[1]]$zmax <- 1
    }
    
    
    # If only one reference assay, return a single panel
    if (!has_second_ref) {
      p1 <- p1 |>
        plotly::layout(
          annotations = list(
            list(
              text   = ref_assay_name_1,
              x      = 0.5,
              y      = 1.02,
              xref   = "paper",
              yref   = "paper",
              xanchor = "center",
              yanchor = "bottom",
              showarrow = FALSE,
              font = list(size = 12, color = "black")
            )
          )
        )
      #session$sendCustomMessage("coralie-toggle-loading", FALSE)
      return(p1)
    }
    
    # Second panel
    p2 <- heatmaply::heatmaply(
      react_cors_df_2[rows, columns, drop = FALSE],
      scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(
        low      = "red",
        high     = "blue",
        midpoint = 0,
        limits   = c(-1, 1)
      ),
      main = paste0(
        "<span style='font-size:12px;'>Correlation of Subsystem Reaction Activity and Fingerprint Prediction - ",
        subsystem_name, "<br>"
      ),
      Rowv        = FALSE,
      Colv        = FALSE,
      margins     = c(50, 50, 130, 50),
      plot_method = "plotly"
    )
    
    p2 <- plotly::layout(
      p2,
      yaxis = list(tickfont = list(size = 10, color = "black"))
    )
    
    p2$x$data[[1]]$colorscale <- list(
      list(0,  "red"),
      list(0.5, "white"),
      list(1,  "blue")
    )
    p2$x$data[[1]]$zmin <- -1
    p2$x$data[[1]]$zmax <- 1
    
    # Two panels stacked
    fig <- plotly::subplot(
      p1, p2,
      nrows  = 2,
      shareX = TRUE,
      titleY = TRUE
    )
    
    fig <- fig |>
      plotly::layout(
        annotations = list(
          list(
            text   = ref_assay_name_1,
            x      = 0.5,     # center over top panel
            y      = 1.02,    # just above top panel
            xref   = "paper",
            yref   = "paper",
            xanchor = "center",
            yanchor = "bottom",
            showarrow = FALSE,
            font = list(size = 12, color = "black")
          ),
          list(
            text   = ref_assay_name_2,
            x      = 0.5,     # center over bottom panel
            y      = 0.48,    # just above bottom panel (adjust as needed)
            xref   = "paper",
            yref   = "paper",
            xanchor = "center",
            yanchor = "bottom",
            showarrow = FALSE,
            font = list(size = 12, color = "black")
          )
        )
      )
    #session$sendCustomMessage("coralie-toggle-loading", FALSE)
    fig
  })
  
}

shinyApp(ui, server)
