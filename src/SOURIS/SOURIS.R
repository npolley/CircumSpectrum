library(shiny)
library(shinydashboard)
library(shinyjs)
library(shinyWidgets)
library(shinyjqui)
library(plotly)
library(ggplot2)
library(dplyr)
library(tidyr)
library(uwot)

#----------------------------------------------------------
# UI
#----------------------------------------------------------

header <- dashboardHeader(
  title = "SOURIS",
  titleWidth = 600
)

sidebar <- dashboardSidebar(
  sidebarMenu(
    id = "main_tab",
    menuItem("Bulk Samples", tabName = "bulk", icon = icon("vial")),
    menuItem("Single-Cell", tabName = "sc",   icon = icon("table-cells-large"))
  )
)

body <- dashboardBody(
  useShinyjs(),
  tags$style(HTML("
  :root {
    --souris-navy:  #002B5C;
    --souris-mid:   #0D47A1;
    --souris-light: #90CAF9;
  }

  /* SOURIS header: replace CORALIE magenta wedge with blue-only gradient */
  .skin-purple .main-header .navbar {
    background: linear-gradient(
      135deg,
      var(--souris-navy),
      var(--souris-mid),
      var(--souris-light)
    );
    border: none;
  }

  .skin-purple .main-header .logo {
    background-color: var(--souris-navy);
    color: #FFFFFF;
    font-family: Poppins, -apple-system, BlinkMacSystemFont, 'Segoe UI',
                 Roboto, 'Helvetica Neue', Arial, sans-serif;
    font-weight: 600;
    font-size: 17px;
    letter-spacing: 0.03em;
    border: none;
  }

  .skin-purple .main-header .logo:hover {
    background-color: var(--souris-mid);
  }
")),
  
  tags$style(HTML("
    :root {
      --souris-navy:  #002B5C;
      --souris-mid:   #0D47A1;
      --souris-light: #90CAF9;
      --souris-bg:    #F4F6FB;
    }

    /* Global typography like CORALIE */
    body, .content-wrapper, .box-title, .sidebar-menu li a {
      font-family: Poppins, -apple-system, BlinkMacSystemFont, 'Segoe UI',
                   Roboto, 'Helvetica Neue', Arial, sans-serif;
    }

    /* Light, slightly blue page background */
    .content-wrapper, .right-side {
      background-color: var(--souris-bg);
    }

    /* Base box look (rounded, soft shadow) */
    .box {
      border-radius: 0px;
      border: 1px solid #e0e4f0;
      box-shadow: 0 4px 12px rgba(0, 0, 0, 0.04);
      overflow: visible;
    }

    .box-header {
      border-bottom: 1px solid #e0e4f0;
      border-top-left-radius: 0px;
      border-top-right-radius: 0px;
      background-clip: padding-box;
    }

    .box-body {
      border-bottom-left-radius: 0px;
      border-bottom-right-radius: 0px;
      background-clip: padding-box;
    }

    /* Primary boxes in SOURIS blue gradient (concordant with header) */
    .skin-purple .box.box-primary {
      border-top-color: var(--souris-mid);
    }

    .skin-purple .box.box-primary .box-header {
      background: linear-gradient(
        135deg,
        var(--souris-navy),
        var(--souris-mid)
      );
      color: #FFFFFF;
    }

    .skin-purple .box.box-primary .box-header .box-title {
      color: #FFFFFF;
      letter-spacing: 0.03em;
      font-weight: 600;
    }

    .skin-purple .box.box-primary .box-body {
      background-color: #E8F1FB;   /* light desaturated blue */
      color: #283046;
    }

    /* Buttons styled to match SOURIS gradient */
    .skin-purple .btn-primary {
      background: linear-gradient(
        135deg,
        var(--souris-navy),
        var(--souris-mid)
      );
      border-color: var(--souris-mid);
    }

    .skin-purple .btn-primary:hover,
    .skin-purple .btn-primary:focus {
      background: linear-gradient(
        135deg,
        var(--souris-mid),
        var(--souris-light)
      );
      border-color: var(--souris-mid);
    }

    .btn {
      border-radius: 20px;
      font-size: 12px;
      font-weight: 500;
      padding: 5px 14px;
    }

    /* Form controls as in CORALIE */
    .form-group {
      margin-bottom: 10px;
    }

    .control-label {
      font-weight: 500;
      font-size: 12px;
      color: #5f6473;
    }

    /* Fade-in animation for new panels (same as CORALIE) */
    .fade-in-up {
      opacity: 0;
      transform: translateY(10px);
      animation: fadeInUp 0.4s ease-out forwards;
    }
    @keyframes fadeInUp {
      from { opacity: 0; transform: translateY(10px); }
      to   { opacity: 1; transform: translateY(0); }
    }

    .fade-in-panel {
      animation: fadeInPanel 0.4s ease-in-out;
    }
    @keyframes fadeInPanel {
      from { opacity: 0; transform: translateY(4px); }
      to   { opacity: 1; transform: translateY(0); }
    }
  ")),
  tags$style(HTML("
  /* SOURIS full-page loading overlay */
  #souris-loading-overlay.loading-overlay {
    position: fixed;
    inset: 0;
    display: flex;
    align-items: center;
    justify-content: center;
    background: radial-gradient(
      circle at top,
      rgba(0, 43, 92, 0.70),     /* navy core */
      rgba(13, 71, 161, 0.50) 60%,  /* mid blue band */
      rgba(0, 20, 40, 0.40)   100%  /* deep outer */
    );
    z-index: 9999;
  }

  #souris-loading-overlay.hidden {
    display: none;
  }

  /* Central loading core, tuned to SOURIS blues */
  .souris-loading-core {
    display: flex;
    flex-direction: column;
    align-items: center;
    justify-content: center;
    text-align: center;
    gap: 0.75rem;
    position: relative;
    z-index: 2;
  }

  #souris-loading-text {
    color: #E3F2FD;
    font-size: 1.5rem;
    letter-spacing: 0.12em;
    text-transform: uppercase;
    text-shadow:
      0 0 4px rgba(144, 202, 249, 0.9),
      0 0 12px rgba(13, 71, 161, 0.9);
  }

  #souris-loading-subtext {
    color: #BBDEFB;
    font-size: 0.85rem;
    letter-spacing: 0.08em;
    text-transform: uppercase;
    opacity: 0.85;
  }

  /* SOURIS loading bar */
  #souris-loading-bar {
    width: 220px;
    height: 4px;
    border-radius: 999px;
    background: rgba(144, 202, 249, 0.25);
    overflow: hidden;
    margin-top: 6px;
  }

  #souris-loading-bar-fill {
    height: 100%;
    width: 0;
    border-radius: 999px;
    background: linear-gradient(
      90deg,
      var(--souris-navy),
      var(--souris-mid),
      var(--souris-light)
    );
    transition: width 0.18s ease-out;
  }
")),tags$style(HTML("
  :root {
    --souris-navy:  #002B5C;
    --souris-mid:   #0D47A1;
    --souris-light: #90CAF9;
  }

  /* Base hex: now with subtle vibration */
  .souris-logo-hex {
    position: relative;
    width: 110px;   /* was 80px */
    height: 110px;  /* was 80px */
    margin-bottom: 18px;
    background: radial-gradient(circle at 30% 25%,
                                #0E1C33,
                                var(--souris-mid) 45%,
                                var(--souris-navy) 85%);
    clip-path: polygon(
      25% 5%, 75% 5%,
      100% 50%, 75% 95%,
      25% 95%, 0% 50%
    );
    display: flex;
    align-items: center;
    justify-content: center;
    box-shadow: 0 14px 30px rgba(0,0,0,0.38);
    overflow: visible;
    animation: sourisHexVibrate 1400ms ease-in-out infinite;
  }

  @keyframes sourisHexVibrate {
    0%   { transform: translate(0px, 0px) scale(1.00); }
    20%  { transform: translate(0.8px, -0.6px) scale(1.01); }
    40%  { transform: translate(-0.7px, 0.5px) scale(1.015); }
    50%  { transform: translate(0.4px, -0.4px) scale(1.02); }
    60%  { transform: translate(-0.6px, 0.7px) scale(1.015); }
    80%  { transform: translate(0.3px, -0.3px) scale(1.01); }
    100% { transform: translate(0px, 0px) scale(1.00); }
  }

  /* Outer glowing hexagon shell */
  .souris-logo-hex-glow {
    position: absolute;
    inset: -8%;
    clip-path: polygon(
      25% 5%, 75% 5%,
      100% 50%, 75% 95%,
      25% 95%, 0% 50%
    );
    border: 2px solid rgba(144,202,249,0.0);
    box-shadow:
      0 0 0 rgba(144,202,249,0.0),
      0 0 0 rgba(13,71,161,0.0);
    mix-blend-mode: screen;
    pointer-events: none;
    animation: sourisHexGlowRing 2.2s ease-in-out infinite;
  }

  @keyframes sourisHexGlowRing {
    0% {
      border-color: rgba(144,202,249,0.0);
      box-shadow: 0 0 0 rgba(144,202,249,0.0);
      opacity: 0.0;
    }
    40% {
      border-color: rgba(144,202,249,0.45);
      box-shadow: 0 0 14px rgba(144,202,249,0.7);
      opacity: 0.8;
    }
    50% {
      border-color: rgba(144,202,249,0.9);
      box-shadow:
        0 0 22px rgba(144,202,249,0.95),
        0 0 40px rgba(13,71,161,0.85),
        0 0 60px rgba(255,255,255,0.75);
      opacity: 1.0;
    }
    65% {
      border-color: rgba(144,202,249,0.5);
      box-shadow: 0 0 16px rgba(144,202,249,0.7);
      opacity: 0.85;
    }
    100% {
      border-color: rgba(144,202,249,0.0);
      box-shadow: 0 0 0 rgba(144,202,249,0.0);
      opacity: 0.0;
    }
  }

  /* Rotating inner hex outline */
    .souris-logo-hex::before {
    content: \"\";
    position: absolute;
    inset: 10%;
    clip-path: polygon(
      25% 5%, 75% 5%,
      100% 50%, 75% 95%,
      25% 95%, 0% 50%
    );
    border: 2px solid rgba(255,255,255,0.85);
    box-shadow:
      0 0 12px rgba(255,255,255,0.7),
      0 0 26px rgba(144,202,249,0.7);
    mix-blend-mode: screen;
    transform-origin: 50% 50%;
    animation: sourisHexOutlineSpin 4s linear infinite;
  }

  @keyframes sourisHexOutlineSpin {
    0%   { transform: rotate(0deg);   opacity: 0.9; }
    50%  { transform: rotate(180deg); opacity: 1.0; }
    100% { transform: rotate(360deg); opacity: 0.9; }
  }

  /* Multiple rotating spokes: three layers with phase offsets */

  .souris-logo-hex::after {
    content: \"\";
    position: absolute;
    inset: 22%;
    border-radius: 999px;
    background: linear-gradient(
      to right,
      transparent 0,
      transparent 49%,
      rgba(255,255,255,0.95) 50%,
      transparent 51%,
      transparent 100%
    );
    mix-blend-mode: screen;
    transform-origin: 50% 50%;
    animation: sourisSpokeSpinMain 1.9s ease-in-out infinite;
  }

  /* extra spokes using child pseudo elements via shadows */
  .souris-logo-hex::after {
    box-shadow:
      0 0 0 0 rgba(255,255,255,0),         /* base handled by animation */
      0 0 0 0 rgba(255,255,255,0);         /* placeholder */
  }

  /* Additional spoke layers using inner elements */
  .souris-logo-hex span.souris-spoke {
    position: absolute;
    inset: 26%;
    border-radius: 999px;
    background: linear-gradient(
      to right,
      transparent 0,
      transparent 49%,
      rgba(255,255,255,0.85) 50%,
      transparent 51%,
      transparent 100%
    );
    mix-blend-mode: screen;
    transform-origin: 50% 50%;
    pointer-events: none;
  }

  .souris-logo-hex span.souris-spoke.spoke-1 {
    animation: sourisSpokeSpinMain 1.9s ease-in-out infinite;
  }
  .souris-logo-hex span.souris-spoke.spoke-2 {
    animation: sourisSpokeSpinOffset 2.1s ease-in-out infinite;
  }
  .souris-logo-hex span.souris-spoke.spoke-3 {
    animation: sourisSpokeSpinOffset2 2.3s ease-in-out infinite;
  }

  @keyframes sourisSpokeSpinMain {
    0% {
      transform: rotate(0deg) scale(1);
      box-shadow: 0 0 0 0 rgba(255,255,255,0.0);
      opacity: 0.7;
    }
    35% {
      transform: rotate(40deg) scale(1.02);
      box-shadow: 0 0 8px 1px rgba(255,255,255,0.4);
      opacity: 0.85;
    }
    45% {
      transform: rotate(43deg) scale(1.05);
      box-shadow:
        0 0 16px 3px rgba(255,255,255,0.8),
        0 0 26px 6px rgba(144,202,249,0.7);
      opacity: 0.95;
    }
    50% {
      transform: rotate(45deg) scale(1.08);
      box-shadow:
        0 0 24px 6px rgba(255,255,255,1.0),
        0 0 40px 12px rgba(144,202,249,0.95),
        0 0 60px 18px rgba(13,71,161,0.85);
      opacity: 1.0;
    }
    55% {
      transform: rotate(47deg) scale(1.04);
      box-shadow:
        0 0 14px 3px rgba(255,255,255,0.75),
        0 0 26px 8px rgba(144,202,249,0.55);
      opacity: 0.9;
    }
    75% {
      transform: rotate(70deg) scale(1.01);
      box-shadow: 0 0 6px 1px rgba(255,255,255,0.3);
      opacity: 0.8;
    }
    100% {
      transform: rotate(90deg) scale(1);
      box-shadow: 0 0 0 0 rgba(255,255,255,0.0);
      opacity: 0.7;
    }
  }

  @keyframes sourisSpokeSpinOffset {
    0%   { transform: rotate(30deg)  scale(1);   opacity: 0.6; }
    35%  { transform: rotate(70deg)  scale(1.02); opacity: 0.8; }
    50%  { transform: rotate(75deg)  scale(1.06); opacity: 0.95; }
    100% { transform: rotate(120deg) scale(1);   opacity: 0.6; }
  }

  @keyframes sourisSpokeSpinOffset2 {
    0%   { transform: rotate(-30deg)  scale(1);   opacity: 0.6; }
    35%  { transform: rotate(10deg)   scale(1.02); opacity: 0.8; }
    50%  { transform: rotate(15deg)   scale(1.06); opacity: 0.95; }
    100% { transform: rotate(60deg)   scale(1);   opacity: 0.6; }
  }
")),
  
  # Loading overlay container
  div(
    id = "souris-loading-overlay",
    class = "loading-overlay hidden",
    div(
      class = "souris-loading-core",
      div(
        class = "souris-logo-hex",
        div(class = "souris-logo-hex-glow"),
        tags$span(class = "souris-spoke spoke-1"),
        tags$span(class = "souris-spoke spoke-2"),
        tags$span(class = "souris-spoke spoke-3")
      ),
      div(id = "souris-loading-text",   "Preparing SOURIS analysis"),
      div(id = "souris-loading-subtext", ""),
      div(
        id = "souris-loading-bar",
        div(id = "souris-loading-bar-fill")
      )
    )
  ),
  
  # JS handlers for SOURIS loading overlay
  tags$script(HTML("
    Shiny.addCustomMessageHandler('souris-toggle-loading', function(show) {
      var overlay = document.getElementById('souris-loading-overlay');
      if (!overlay) return;
      if (show) overlay.classList.remove('hidden'); else overlay.classList.add('hidden');
    });
    Shiny.addCustomMessageHandler('souris-loading-text', function(msg) {
      var el = document.getElementById('souris-loading-text');
      if (!el) return;
      el.textContent = msg;
    });
    Shiny.addCustomMessageHandler('souris-loading-subtext', function(msg) {
      var el = document.getElementById('souris-loading-subtext');
      if (!el) return;
      el.textContent = msg;
    });
  ")),
  
  tabItems(
    
    #------------------------------------------------------
    # Bulk Samples tab (CORALIE-like panels, 1 fingerprint)
    #------------------------------------------------------
    tabItem(
      tabName = "bulk",
      class   = "fade-in-up",
      
      # Fingerprint selection
      fluidRow(
        box(
          width = 12, status = "primary", solidHeader = TRUE,
          title = "Fingerprint selection",
          shinyWidgets::pickerInput(
            inputId = "souris_fingerprints_sel",
            label   = "Select one fingerprint",
            choices = NULL,          # populated server-side
            multiple = FALSE,        # <<< only one fingerprint allowed
            options = list(
              `actions-box`        = FALSE,
              `none-selected-text` = "No fingerprint selected"
            )
          ),
          br(),
          actionButton(
            inputId = "souris_load_fingerprints",
            label   = "Select fingerprint"
          )
        )
      ),
      
      # Placeholders for dynamic UI (same pattern as CORALIE)
      fluidRow(div(id = "bulk_reference_placeholder")),
      fluidRow(div(id = "bulk_summary_placeholder")),
      fluidRow(div(id = "bulk_analysis_placeholder")),
      fluidRow(div(id = "bulk_corplotgrid_placeholder")),
      fluidRow(div(id = "bulk_tier3_placeholder"))
    ),
    
    #------------------------------------------------------
    # Single-Cell tab (stub)
    #------------------------------------------------------
    tabItem(
      tabName = "sc",
      class   = "fade-in-up",
      fluidRow(
        box(
          width = 12, status = "primary", solidHeader = TRUE,
          title = "Single-Cell (coming soon)",
          p("Single-cell panels will be implemented here.")
        )
      )
    )
  )
)

ui <- dashboardPage(
  skin   = "purple",
  header = header,
  sidebar = sidebar,
  body   = body
)

#----------------------------------------------------------
# SERVER
#----------------------------------------------------------

server <- function(input, output, session) {
  
  # paths analogous to CORALIE
  fingerprints_dir <- "../../data/fingerprints"
  reference_dir    <- "../../data/reference_assays"
  
  # reactive state
  souris_fingerprints_list <- reactiveVal(NULL)
  souris_reference_list    <- reactiveVal(NULL)
  
  ref_box_inserted_bulk    <- reactiveVal(FALSE)
  summary_box_inserted     <- reactiveVal(FALSE)
  analysis_panel_inserted  <- reactiveVal(FALSE)
  analysis_done_bulk       <- reactiveVal(FALSE)
  
  df_all_predictions <- reactiveVal(NULL)
  df_all_scores <- reactiveVal(NULL)
  
  detail_panel_inserted <- reactiveVal(FALSE)
  clicked_subsystem     <- reactiveVal(NULL)
  
  reaction_order <- reactiveVal(NULL)
  
  # helper for loading overlay
  set_souris_progress <- function(text = NULL, subtext = NULL, pct = NULL, show = NULL) {
    if (!is.null(text))
      session$sendCustomMessage("souris-loading-text", text)
    if (!is.null(subtext))
      session$sendCustomMessage("souris-loading-subtext", subtext)
    if (!is.null(pct))
      session$sendCustomMessage("souris-loading-progress", round(pct))
    if (!is.null(show))
      session$sendCustomMessage("souris-toggle-loading", show)
  }
  
  auc_two_groups <- function(pred, grp) {
    # grp: factor with 2 levels, length(pred)
    ord <- order(pred)
    pred <- pred[ord]
    grp  <- grp[ord]
    y <- as.integer(grp == levels(grp)[2])  # 0/1
    n1 <- sum(y == 0)
    n2 <- sum(y == 1)
    if (n1 == 0 || n2 == 0) return(NA_real_)
    # Mann–Whitney U / AUC
    ranks <- rank(pred)
    R2 <- sum(ranks[y == 1])
    U  <- R2 - n2 * (n2 + 1) / 2
    U / (n1 * n2)
  }
  
  human_react_meta<-read.csv("../../data/fingerprint_prep_objects/human_reaction_meta.csv", header = T, row.names = 1)
  mouse_react_meta<-read.csv("../../data/fingerprint_prep_objects/mouse_reaction_meta.csv", header = T, row.names = 1)
  react_meta_filt<-human_react_meta[intersect(rownames(human_react_meta), rownames(mouse_react_meta)),]
  
  get_reaction_preds <- function(ref_id, subsel) {
    # 1) Get the reference assay matrix
    ref_mat <- souris_reference_list()[[ref_id]]
    validate(need(!is.null(ref_mat),
                  paste("No reaction matrix found for reference assay", ref_id, ".")))
    
    # 2) Select reactions whose subsystem matches subsel
    #    (make.names if you normalized names earlier)
    subs_react <- rownames(react_meta_filt)[react_meta_filt$subsystem == subsel]
    subs_react <- intersect(subs_react, colnames(ref_mat))
    validate(need(length(subs_react) > 0,
                  paste("No reactions in react_meta_filt for subsystem", subsel, "in", ref_id, ".")))
    
    # 3) Return the reference matrix restricted to those reaction columns
    mat <- as.matrix(ref_mat[, subs_react, drop = FALSE])
    rownames(mat) <- rownames(ref_mat)  # samples
    mat
  }
  
  observe({
    files <- list.files(path = fingerprints_dir, pattern = "\\.rds$", full.names = FALSE)
    stripped <- sub("\\_.rds$", "", files)
    shinyWidgets::updatePickerInput(
      session,
      inputId = "souris_fingerprints_sel",
      choices = stripped,
      selected = NULL
    )
  })
  
  # enable button only if one fingerprint selected
  observe({
    n_sel <- length(input$souris_fingerprints_sel)
    if (is.null(n_sel) || n_sel != 1) {
      shinyjs::disable("souris_load_fingerprints")
    } else {
      shinyjs::enable("souris_load_fingerprints")
    }
  })
  
  #--------------------------------------------------------
  # Load fingerprint (single) and insert reference box
  #--------------------------------------------------------
  observeEvent(input$souris_load_fingerprints, {
    req(input$souris_fingerprints_sel)
    validate(need(length(input$souris_fingerprints_sel) == 1,
                  "Please select exactly one fingerprint."))
    
    set_souris_progress("Loading fingerprint", "", TRUE)
    
    sel_name <- input$souris_fingerprints_sel
    path     <- file.path(fingerprints_dir, paste0(sel_name, "_.rds"))
    
    size      <- file.info(path)$size
    totalSize <- ifelse(is.na(size), 1, size)
    cumSize   <- 0
    
    
    obj <- readRDS(path)
    cumSize <- totalSize
    set_souris_progress(
      subtext = sprintf("Loading %s", sel_name),
      pct     = 100 * cumSize / totalSize
    )
    
    souris_fingerprints_list(setNames(list(obj), sel_name))
    
    set_souris_progress(subtext = "", pct = 100, show = FALSE)
    
    set_souris_progress(show = FALSE)
    
    # Insert reference assay selection box once
    if (!ref_box_inserted_bulk()) {
      insertUI(
        selector = "#bulk_reference_placeholder",
        where    = "afterBegin",
        ui = fluidRow(
          class = "fade-in-up",
          column(
            width = 12,
            box(
              width = 12, status = "primary", solidHeader = TRUE,
              title = "Reference assay selection",
              selectInput(
                inputId = "souris_reference_sel",
                label   = "Select one or more reference assays",
                choices = NULL,
                multiple = TRUE
              ),
              actionButton(
                inputId = "souris_load_reference",
                label   = "Select reference assay(s)"
              )
            )
          )
        ),
        immediate = TRUE
      )
      ref_box_inserted_bulk(TRUE)
    }
  })
  
  #--------------------------------------------------------
  # Populate reference choices and load reference assays
  #--------------------------------------------------------
  observe({
    req(ref_box_inserted_bulk())
    files_ref  <- list.files(path = reference_dir, pattern = "\\.csv$", full.names = FALSE)
    stripped   <- sub("\\.csv$", "", files_ref)

    updateSelectInput(
      session,
      "souris_reference_sel",
      choices  = stripped,
      selected = NULL
    )
  })
  
  observeEvent(input$souris_load_reference, {
    n_sel <- length(input$souris_reference_sel)
    if (is.null(n_sel) || n_sel == 0) {
      showModal(modalDialog(
        title = "Selection required",
        "Please select at least one reference assay before proceeding.",
        easyClose = TRUE
      ))
      return(NULL)
    }
    
    set_souris_progress("Loading reference assays", "", pct = 0, show = TRUE)
    
    sel_names <- input$souris_reference_sel
    paths     <- file.path(reference_dir, paste0(sel_names, ".csv"))
    
    sizes <- file.info(paths)$size
    sizes[is.na(sizes) | sizes == 0] <- 1
    totalSize <- sum(sizes)
    cumSize   <- 0
    
    objs <- vector("list", length(paths))
    for (i in seq_along(paths)) {
      p  <- paths[i]
      nm <- sel_names[i]
      
      set_souris_progress(
        subtext = sprintf("Loading reference assay %d/%d: %s", i, length(paths), nm),
        pct     = 100 * cumSize / totalSize
      )
      
      df <- read.csv(p, stringsAsFactors = FALSE)
      objs[[i]] <- df
      cumSize <- cumSize + sizes[i]
    }
    
    names(objs) <- sel_names
    souris_reference_list(objs)
    
    set_souris_progress(show = FALSE)
    
    # Here you can insert summary / Tier I / Tier II / Tier III panels
    # for SOURIS bulk, using the same structure as CORALIE but restricted
    # to one fingerprint (souris_fingerprints_list()).
    
    # Example: insert a placeholder analysis panel once
    if (!analysis_panel_inserted()) {
      insertUI(
        selector = "#bulk_analysis_placeholder",
        where    = "afterEnd",
        ui = tagList(
          fluidRow(
            class = "fade-in-up",
            column(
              width = 6,
              box(
                width = 12, status = "primary", solidHeader = TRUE,
                title = "Total Fingerprint Scores by Reference Assay",
                plotly::plotlyOutput("souris_violin_fp", height = "500px")
              )
            ),
            column(
              width = 6,
              box(
                width = 12, status = "primary", solidHeader = TRUE,
                title = "Subsystem AUC Performance Distinguishing Reference Assays",
                div(
                  style = "height: 500px; overflow-y: auto;",
                  plotly::plotlyOutput("souris_auc_bar", height = "800px")
                )
              )
            )
          ),
          fluidRow(
            class = "fade-in-up",
            column(
              width = 12,
              box(
                width = 12, status = "primary", solidHeader = TRUE,
                id    = "souris_umap_box",
                title = "UMAP of samples coloured by reference assay",
                plotly::plotlyOutput("souris_umap", height = "500px")
              )
            )
          )
        ),
        immediate = TRUE
      )
      analysis_panel_inserted(TRUE)
    }
  })
  
  output$souris_violin_fp <- plotly::renderPlotly({
    req(souris_fingerprints_list())
    req(souris_reference_list())
    req(input$souris_reference_sel)
    
    set_souris_progress(
      text    = "Running subsystem models",
      subtext = "Computing fingerprint predictions across reference assays",
      show    = TRUE
    )
    
    on.exit(set_souris_progress(show = FALSE), add = TRUE)
    
    validate(
      need(length(input$souris_reference_sel) > 0,
           "Select at least one reference assay to view predictions."),
      need(length(souris_fingerprints_list()) >= 1,
           "Fingerprint object not loaded.")
    )
    
    fp_obj <- souris_fingerprints_list()[[1]]
    validate(
      need(!is.null(fp_obj$All), "Fingerprint object has no 'All' component."),
      need(!is.null(fp_obj$All$fingerprints_primary),
           "No primary fingerprint classifiers found.")
    )
    
    fingerprint <- fp_obj$All
    final_classifiers <- fingerprint$fingerprints_primary
    auc_list          <- fingerprint$auc_primary
    
    rownames(auc_list) <- auc_list$index
    index        <- auc_list$index
    system_names <- auc_list[index, ]$subsystem
    
    ref_ids    <- input$souris_reference_sel
    all_violin <- list()
    all_scores <- list()
    
    normalize01 <- function(v) {
      rng <- range(v, na.rm = TRUE)
      if (diff(rng) == 0 || !is.finite(diff(rng))) return(rep(0.5, length(v)))
      (v - rng[1]) / (rng[2] - rng[1])
    }
    
    for (ref_id in ref_ids) {
      ref_assay <- souris_reference_list()[[ref_id]]
      validate(need(!is.null(ref_assay),
                    paste("Reference assay", ref_id, "not found.")))
      
      ex_subsystems_no_norm <- vector("list", length(index))
      names(ex_subsystems_no_norm) <- system_names
      
      for (x in seq_along(index)) {
        idx  <- as.numeric(index[x])
        if (is.null(final_classifiers[[idx]])) next
        
        coef  <- final_classifiers[[idx]]$coefnames
        feats <- intersect(coef, colnames(ref_assay))
        if (!length(feats)) next
        
        newDat <- as.matrix(ref_assay[, feats, drop = FALSE])
        colnames(newDat) <- feats
        
        pred_tot <- predict(
          final_classifiers[[idx]],
          newdata   = newDat,
          type      = "prob",
          na.action = na.pass
        )[ , 1]
        
        ex_subsystems_no_norm[[x]] <- pred_tot
      }
      
      ex_subsystems_no_norm <- ex_subsystems_no_norm[!vapply(ex_subsystems_no_norm, is.null, logical(1))]
      
      if (!length(ex_subsystems_no_norm)) next  # skip empty ref
      
      df_ref <- lapply(names(ex_subsystems_no_norm), function(subsys) {
        data.frame(
          Reference  = ref_id,
          Subsystem  = subsys,
          Prediction = ex_subsystems_no_norm[[subsys]],
          stringsAsFactors = FALSE
        )
      })
      all_violin[[ref_id]] <- dplyr::bind_rows(df_ref)
      
      mat <- data.frame(ex_subsystems_no_norm, check.names = FALSE)
      mat <- as.data.frame(lapply(mat, normalize01))
      
      score <- rowMeans(mat, na.rm = TRUE)
      score <- (score - min(score, na.rm = TRUE)) /
        (max(score, na.rm = TRUE) - min(score, na.rm = TRUE))
      
      df_score <- data.frame(
        Sample    = seq_along(score),
        Reference = ref_id,
        Score     = score,
        stringsAsFactors = FALSE
      )
      all_scores[[ref_id]] <- df_score
    }
    
    if (!length(all_violin)) {
      validate(need(FALSE, "No valid subsystem predictions could be computed for the selected reference assays."))
    }
    
    df_all    <- dplyr::bind_rows(all_violin)
    df_scores <- dplyr::bind_rows(all_scores)
    
    df_all_predictions(df_all)
    df_all_scores(df_scores)
    
    validate(need(nrow(df_scores) > 0,
                  "No aggregated fingerprint scores could be computed for the selected reference assays."))
    
    # Set factor order so the first level is the first selected reference assay
    df_scores$Reference <- factor(df_scores$Reference,
                                  levels = unique(df_scores$Reference))
    
    ref_levels <- levels(df_scores$Reference)
    baseline   <- ref_levels[1]
    
    p_to_stars <- function(p) {
      if (is.na(p)) return("ns")
      if (p < 0.0001) return("****")
      if (p < 0.001)  return("***")
      if (p < 0.01)   return("**")
      if (p < 0.05)   return("*")
      "ns"
    }
    
    tests <- lapply(ref_levels[-1], function(ref) {
      x <- df_scores$Score[df_scores$Reference == baseline]
      y <- df_scores$Score[df_scores$Reference == ref]
      p <- stats::wilcox.test(x, y, exact = FALSE)$p.value
      data.frame(
        group1 = baseline,
        group2 = ref,
        p      = p,
        p_lab  = p_to_stars(p)
      )
    })
    stat_df <- dplyr::bind_rows(tests)
    
    # if only one reference selected, no stars
    if (nrow(stat_df) == 0) {
      y_max <- max(df_scores$Score, na.rm = TRUE)
      
      p <- ggplot2::ggplot(
        df_scores,
        ggplot2::aes(x = Reference, y = Score, fill = Reference)
      ) +
        ggplot2::geom_violin(alpha = 0.6, scale = "width", trim = FALSE) +
        ggplot2::geom_boxplot(width = 0.12, outlier.size = 0.5,
                              alpha = 0.9, color = "grey20") +
        ggplot2::geom_jitter(
          width = 0.08, size = 1.2, alpha = 0.7,
          color = "black"
        ) +
        ggplot2::labs(
          x = "Reference assay",
          y = "Aggregated fingerprint score",
          title = "Aggregated fingerprint scores per reference assay"
        ) +
        ggplot2::theme_minimal(base_size = 11) +
        ggplot2::theme(
          axis.text.x  = ggplot2::element_text(angle = 60, hjust = 1, vjust = 1),
          plot.title   = ggplot2::element_text(hjust = 0.5, face = "bold"),
          legend.position = "none"
        )
      
      set_souris_progress(
        subtext = "Rendering subsystem score plots",
        pct     = 100
      )
      
      return(plotly::ggplotly(p, tooltip = c("x", "y")))
    }
    
    # if we have at least one comparison, compute y positions and add stars
    y_max <- max(df_scores$Score, na.rm = TRUE)
    stat_df$y_pos <- seq(y_max * 1.05, length.out = nrow(stat_df), by = y_max * 0.05)
    
    p <- ggplot2::ggplot(
      df_scores,
      ggplot2::aes(x = Reference, y = Score, fill = Reference)
    ) +
      ggplot2::geom_violin(alpha = 0.6, scale = "width", trim = FALSE) +
      ggplot2::geom_boxplot(width = 0.12, outlier.size = 0.5,
                            alpha = 0.9, color = "grey20") +
      ggplot2::geom_jitter(
        width = 0.08, size = 1.2, alpha = 0.7,
        color = "black"
      ) +
      ggplot2::geom_segment(
        data = stat_df,
        ggplot2::aes(
          x    = as.numeric(group1),
          xend = as.numeric(group2),
          y    = y_pos,
          yend = y_pos
        ),
        inherit.aes = FALSE
      ) +
      ggplot2::geom_text(
        data = stat_df,
        ggplot2::aes(
          x = (as.numeric(group1) + as.numeric(group2)) / 2,
          y = y_pos,
          label = p_lab
        ),
        vjust = -0.3,
        inherit.aes = FALSE
      ) +
      ggplot2::labs(
        x = "Reference assay",
        y = "Aggregated fingerprint score",
        title = "Aggregated fingerprint scores per reference assay"
      ) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(
        axis.text.x  = ggplot2::element_text(angle = 60, hjust = 1, vjust = 1),
        plot.title   = ggplot2::element_text(hjust = 0.5, face = "bold"),
        legend.position = "none"
      )
    
    plotly::ggplotly(p, tooltip = c("x", "y"))
  })
  
  output$souris_auc_bar <- plotly::renderPlotly({
    df_all <- df_all_predictions()
    req(df_all)
    
    # need 2+ references
    refs <- unique(df_all$Reference)
    validate(need(length(refs) >= 2,
                  "Select at least two reference assays to compute subsystem AUC."))
    
    ref1 <- refs[1]
    ref2 <- refs[2]
    
    # keep only first two references
    df_pair <- df_all[df_all$Reference %in% c(ref1, ref2), ]
    df_pair$Reference <- factor(df_pair$Reference, levels = c(ref1, ref2))
    
    # compute AUC per subsystem
    auc_df <- df_pair |>
      dplyr::group_by(Subsystem) |>
      dplyr::summarise(
        AUC = auc_two_groups(Prediction, Reference),
        .groups = "drop"
      )
    
    auc_df <- dplyr::filter(auc_df, !is.na(AUC))
    validate(need(nrow(auc_df) > 0,
                  "No valid AUC values could be computed for subsystem models."))
    
    p_auc <- ggplot2::ggplot(
      auc_df,
      ggplot2::aes(x = reorder(Subsystem, AUC), y = AUC, customdata = Subsystem)
    ) +
      ggplot2::geom_col(fill = "#0D47A1") +
      ggplot2::coord_flip() +
      ggplot2::geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
      ggplot2::labs(
        x = "Subsystem",
        y = "AUC (first vs second reference)",
        title = paste("Subsystem AUC:", ref1, "vs", ref2)
      ) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
      )
    
    plotly::ggplotly(p_auc, tooltip = c("x", "y"), source = "souris_auc")
  })
  
  output$souris_umap <- plotly::renderPlotly({
    df_all <- df_all_predictions()
    req(df_all)
    
    validate(need(nrow(df_all) >= 3,
                  "Not enough subsystem predictions to compute a UMAP embedding."))
    
    # Each point: subsystem in a given reference.
    # Features: prediction vector across samples.
    # We need one row per (Subsystem, Reference) and columns = samples.
    df_all$SampleID <- ave(df_all$Prediction, df_all$Reference, df_all$Subsystem,
                           FUN = seq_along)
    
    # wide: rows = Subsystem-Reference, cols = SampleID
    wide_mat <- df_all |>
      dplyr::select(Reference, Subsystem, SampleID, Prediction) |>
      tidyr::pivot_wider(
        id_cols  = c(Reference, Subsystem),
        names_from  = SampleID,
        values_from = Prediction
      )
    
    feat_mat <- as.matrix(wide_mat[, !(names(wide_mat) %in% c("Reference", "Subsystem")), drop = FALSE])
    
    # remove rows with all NA
    keep <- rowSums(!is.na(feat_mat)) > 0
    feat_mat <- feat_mat[keep, , drop = FALSE]
    wide_mat <- wide_mat[keep, , drop = FALSE]
    
    validate(need(nrow(feat_mat) >= 3,
                  "Not enough subsystem points with valid predictions for UMAP."))
    
    # simple NA imputation by column means
    if (anyNA(feat_mat)) {
      cm <- apply(feat_mat, 2, function(x) if (all(is.na(x))) 0 else mean(x, na.rm = TRUE))
      for (j in seq_len(ncol(feat_mat))) {
        nas <- is.na(feat_mat[, j])
        if (any(nas)) feat_mat[nas, j] <- cm[j]
      }
    }
    
    # n_neighbors must be < nrow
    nn <- min(15, max(2, nrow(feat_mat) - 1))
    
    emb <- uwot::umap(
      feat_mat,
      n_neighbors = nn,
      min_dist    = 0.3,
      metric      = "euclidean"
    )
    
    umap_df <- data.frame(
      UMAP1     = emb[, 1],
      UMAP2     = emb[, 2],
      Reference = wide_mat$Reference,
      Subsystem = wide_mat$Subsystem,
      stringsAsFactors = FALSE
    )
    
    p_umap <- ggplot2::ggplot(
      umap_df,
      ggplot2::aes(
        x     = UMAP1,
        y     = UMAP2,
        color = Reference,
        text  = paste0("Subsystem: ", Subsystem, "<br>Reference: ", Reference)
      )
    ) +
      ggplot2::geom_point(alpha = 0.8, size = 2) +
      ggplot2::labs(
        x = "UMAP1",
        y = "UMAP2",
        title = "UMAP of subsystem models coloured by reference assay"
      ) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
      )
    
    plotly::ggplotly(p_umap, tooltip = c("text"))
  })
  
  observeEvent(plotly::event_data("plotly_click", source = "souris_auc"), {
    d <- plotly::event_data("plotly_click", source = "souris_auc")

    if (is.null(d)) return()
    
    subsel <- as.character(d$customdata)  # subsystem name from bar label
    clicked_subsystem(subsel)
    
    set_souris_progress(
      text    = "Analyzing subsystem reactions",
      subtext = sprintf("Computing reaction-level metrics for '%s'", subsel),
      show    = TRUE
    )
    
    if (!detail_panel_inserted()) {
      insertUI(
        selector = "#bulk_tier3_placeholder",
        where    = "afterEnd",
        ui = tagList(
          fluidRow(
            class = "fade-in-up",
            column(
              width = 6,
              box(
                width = 12, status = "primary", solidHeader = TRUE,
                title = "Selected subsystem prediction by reference assay",
                plotly::plotlyOutput("souris_violin_subsys", height = "400px")
              )
            ),
            column(
              width = 6,
              box(
                width = 12, status = "primary", solidHeader = TRUE,
                title = "Reaction-level AUC for selected subsystem",
                plotly::plotlyOutput("souris_auc_reactions", height = "400px")
              )
            )
          ),fluidRow(
            class = "fade-in-up",
            box(
              width = 12, status = "primary", solidHeader = TRUE,
              title = "UMAP of reactions within selected subsystem",
              plotly::plotlyOutput("souris_umap_reactions", height = "500px")
            )
          ),fluidRow(
            class = "fade-in-up",
            box(
              width = 12, status = "primary", solidHeader = TRUE,
              title = "Reaction predictions by reference assay",
              plotly::plotlyOutput("souris_violin_reactions", height = "500px")
            )
          )
        ),
        immediate = TRUE
      )
      detail_panel_inserted(TRUE)
    }
    set_souris_progress(subtext = "Rendering subsystem detail plots")
  })
  
  output$souris_violin_subsys <- plotly::renderPlotly({
    df_all <- df_all_predictions()
    subsel <- clicked_subsystem()
    req(df_all, subsel)

    df_sub <- df_all[df_all$Subsystem == subsel, ]
    View(df_sub)
    validate(need(nrow(df_sub) > 0,
                  "No predictions found for the selected subsystem."))
    
    p_sub <- ggplot2::ggplot(
      df_sub,
      ggplot2::aes(x = Reference, y = Prediction, fill = Reference)
    ) +
      ggplot2::geom_violin(alpha = 0.6, scale = "width", trim = FALSE) +
      ggplot2::geom_boxplot(width = 0.12, outlier.size = 0.5,
                            alpha = 0.9, color = "grey20") +
      ggplot2::geom_jitter(width = 0.08, size = 1.2, alpha = 0.7,
                           color = "black") +
      ggplot2::labs(
        x = "Reference assay",
        y = "Subsystem prediction probability",
        title = paste("Subsystem:", subsel)
      ) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(
        axis.text.x  = ggplot2::element_text(angle = 60, hjust = 1, vjust = 1),
        plot.title   = ggplot2::element_text(hjust = 0.5, face = "bold"),
        legend.position = "none"
      )
    
    plotly::ggplotly(p_sub, tooltip = c("x", "y"))
  })
  
  output$souris_auc_reactions <- plotly::renderPlotly({
    subsel <- clicked_subsystem()
    req(subsel)
    
    refs <- input$souris_reference_sel
    validate(need(length(refs) >= 2,
                  "Select at least two reference assays to compute reaction-level AUC."))
    ref1 <- refs[1]
    ref2 <- refs[2]
    
    # reaction matrices: rows = samples, cols = reactions
    mat1 <- get_reaction_preds(ref1, subsel)
    mat2 <- get_reaction_preds(ref2, subsel)
    validate(need(!is.null(mat1) && !is.null(mat2),
                  "No reaction predictions available for the selected subsystem."))
    
    # intersect reactions and align samples as in your Tier III code
    common_rxn <- intersect(colnames(mat1), colnames(mat2))
    validate(need(length(common_rxn) > 0,
                  "No overlapping reactions for this subsystem in the two reference assays."))
    
    mat1 <- mat1[, common_rxn, drop = FALSE]
    mat2 <- mat2[, common_rxn, drop = FALSE]
    
    # build long df: each row = sample-reaction, with ref label
    df_rxn <- rbind(
      data.frame(Reference = ref1,
                 Sample    = rownames(mat1),
                 as.data.frame(mat1, check.names = FALSE)),
      data.frame(Reference = ref2,
                 Sample    = rownames(mat2),
                 as.data.frame(mat2, check.names = FALSE))
    )
    
    # compute AUC per reaction across the two references
    rxn_names <- common_rxn
    auc_rxn <- sapply(rxn_names, function(r) {
      preds <- df_rxn[[r]]
      grp   <- factor(df_rxn$Reference, levels = c(ref1, ref2))
      auc_two_groups(preds, grp)
    })
    
    auc_df <- data.frame(
      Reaction = rxn_names,
      AUC      = auc_rxn,
      stringsAsFactors = FALSE
    )
    auc_df <- dplyr::filter(auc_df, !is.na(AUC))
    validate(need(nrow(auc_df) > 0,
                  "No valid reaction-level AUC values could be computed."))
    
    rxn_order <- auc_df$Reaction
    reaction_order(rxn_order)
    
    p_rxn <- ggplot2::ggplot(
      auc_df,
      ggplot2::aes(x = reorder(Reaction, AUC), y = AUC)
    ) +
      ggplot2::geom_col(fill = "#0D47A1") +
      ggplot2::coord_flip() +
      ggplot2::geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
      ggplot2::labs(
        x = "Reaction",
        y = "AUC (first vs second reference)",
        title = paste("Reaction AUC for subsystem", subsel, ":", ref1, "vs", ref2)
      ) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
      )
    
    plotly::ggplotly(p_rxn, tooltip = c("x", "y"))
  })
  
  output$souris_umap_reactions <- plotly::renderPlotly({
    subsel <- clicked_subsystem()
    req(subsel)
    
    set_souris_progress(
      text    = "Embedding reactions",
      subtext = sprintf("UMAP on reactions in subsystem '%s'", subsel),
      show    = TRUE
    )
    on.exit(set_souris_progress(show = FALSE), add = TRUE)
    
    refs <- input$souris_reference_sel
    validate(need(length(refs) >= 2,
                  "Select at least two reference assays to compute reaction-level UMAP."))
    
    ref1 <- refs[1]
    ref2 <- refs[2]
    
    mat1 <- get_reaction_preds(ref1, subsel)  # samples × reactions
    mat2 <- get_reaction_preds(ref2, subsel)
    validate(need(!is.null(mat1) && !is.null(mat2),
                  "No reaction predictions available for the selected subsystem."))
    
    # intersect reactions and align samples by row position
    common_rxn <- intersect(colnames(mat1), colnames(mat2))
    validate(need(length(common_rxn) >= 3,
                  "Not enough common reactions for UMAP."))
    mat1 <- mat1[, common_rxn, drop = FALSE]
    mat2 <- mat2[, common_rxn, drop = FALSE]
    
    # Build long: each row = reaction-reference, features = samples
    df_rxn <- rbind(
      data.frame(
        Reaction  = common_rxn,
        Reference = ref1,
        t(mat1),                    # rows: reactions, cols: samples
        check.names = FALSE
      ),
      data.frame(
        Reaction  = common_rxn,
        Reference = ref2,
        t(mat2),
        check.names = FALSE
      )
    )
    
    feat_mat <- as.matrix(df_rxn[, !(names(df_rxn) %in% c("Reaction", "Reference")), drop = FALSE])
    
    # remove rows with all NA
    keep <- rowSums(!is.na(feat_mat)) > 0
    feat_mat <- feat_mat[keep, , drop = FALSE]
    df_rxn   <- df_rxn[keep, , drop = FALSE]
    validate(need(nrow(feat_mat) >= 3,
                  "Not enough reaction points with valid values for UMAP."))
    
    # impute remaining NA by column mean
    if (anyNA(feat_mat)) {
      cm <- apply(feat_mat, 2, function(x) if (all(is.na(x))) 0 else mean(x, na.rm = TRUE))
      for (j in seq_len(ncol(feat_mat))) {
        nas <- is.na(feat_mat[, j])
        if (any(nas)) feat_mat[nas, j] <- cm[j]
      }
    }
    
    nn <- min(15, max(2, nrow(feat_mat) - 1))
    
    emb <- uwot::umap(
      feat_mat,
      n_neighbors = nn,
      min_dist    = 0.3,
      metric      = "euclidean"
    )
    
    umap_df <- data.frame(
      UMAP1    = emb[, 1],
      UMAP2    = emb[, 2],
      Reaction = df_rxn$Reaction,
      Reference = df_rxn$Reference,
      stringsAsFactors = FALSE
    )
    
    p <- ggplot2::ggplot(
      umap_df,
      ggplot2::aes(
        x     = UMAP1,
        y     = UMAP2,
        color = Reference,
        text  = paste0("Reaction: ", Reaction, "<br>Reference: ", Reference)
      )
    ) +
      ggplot2::geom_point(alpha = 0.8, size = 2) +
      ggplot2::labs(
        x = "UMAP1",
        y = "UMAP2",
        title = paste("UMAP of reactions in subsystem", subsel)
      ) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
      )
    
    plotly::ggplotly(p, tooltip = c("text"))
  })
  
  output$souris_violin_reactions <- plotly::renderPlotly({
    subsel <- clicked_subsystem()
    req(subsel)
    
    refs <- input$souris_reference_sel
    validate(need(length(refs) >= 1,
                  "Select at least one reference assay to view reaction predictions."))
    
    # long df of raw reaction values from reference assays
    df_rxn_long <- do.call(
      rbind,
      lapply(refs, function(ref_id) {
        mat <- get_reaction_preds(ref_id, subsel)  # samples × reactions (raw assay values)
        validate(need(!is.null(mat),
                      paste("No reactions for subsystem", subsel, "in", ref_id)))
        if (!nrow(mat) || !ncol(mat)) return(NULL)
        
        df <- as.data.frame(mat, check.names = FALSE)
        df$Sample <- rownames(mat)
        df_long <- tidyr::pivot_longer(
          df,
          cols      = -Sample,
          names_to  = "Reaction",
          values_to = "Value"        # raw reaction value from reference assay
        )
        df_long$Reference <- ref_id
        df_long
      })
    )
    
    rxn_order <- reaction_order()
    if (!is.null(rxn_order)) {
      rxn_levels <- intersect(rxn_order, unique(df_rxn_long$Reaction))
    } else {
      rxn_levels <- sort(unique(df_rxn_long$Reaction))
    }
    df_rxn_long$Reaction  <- factor(df_rxn_long$Reaction, levels = rxn_levels)
    df_rxn_long$Reference <- factor(df_rxn_long$Reference,
                                    levels = unique(df_rxn_long$Reference))
    
    p <- ggplot2::ggplot(
      df_rxn_long,
      ggplot2::aes(x = Reaction, y = Value, fill = Reference)
    ) +
      # drop geom_violin() entirely
      ggplot2::geom_boxplot(
        outlier.size = 0.4,
        alpha   = 0.9,
        color   = "grey20",
        position = ggplot2::position_dodge2(width = 0.9, padding = 0.2)  # more spacing
      ) +
      ggplot2::geom_jitter(
        ggplot2::aes(text = paste0(
          "Reaction: ", Reaction,
          "<br>Reference: ", Reference,
          "<br>Sample: ", Sample,
          "<br>Value: ", sprintf('%.3f', Value)
        )),
        size     = 0.6,
        alpha    = 0.6,
        color    = "black",
        position = ggplot2::position_jitterdodge(
          jitter.width = 0.15,
          dodge.width  = 0.9
        )
      ) +
      ggplot2::labs(
        x = "Reaction",
        y = "Reaction value (reference assay)",
        title = paste("Per-reaction values in subsystem", subsel)
      ) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(
        axis.text.x  = ggplot2::element_text(angle = 60, hjust = 1, vjust = 1),
        plot.title   = ggplot2::element_text(hjust = 0.5, face = "bold"),
        legend.position = "top"
      )
    
    plotly::ggplotly(p, tooltip = c("text"))
  })
}

shinyApp(ui, server)