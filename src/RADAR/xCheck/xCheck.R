library(caTools)
library(matrixTests)
library(ggplot2)
library(plotly)
library(shinydashboard)
library(reactable)
library(cowplot)
library(shinyWidgets)
library(data.table)
library(MatrixGenerics)

reactMeta<-read.csv("../../../data/RADAR_xCheck_cohort/human_reaction_meta.csv")
colnames(reactMeta)[1]<-"metabolite"

cohort_registry <- read.csv("../../../data/RADAR_xCheck_cohort/xCheck_cohort_registry.csv", stringsAsFactors = FALSE)

experiment_registry <- read.csv("../../../data/RADAR_xCheck_experimental_assay/xCheck_experiment_registry.csv", stringsAsFactors = FALSE)

header <- dashboardHeader(
  title = "RADAR | xCheck",
  titleWidth = 600
)

body <- dashboardBody(
  #Cosmetic CSS and javascript elements
  tags$head(
    tags$style(HTML("
  .fade-in-up {
    opacity: 0;
    transform: translateY(10px);
    animation: fadeInUp 0.4s ease-out forwards;
  }

  @keyframes fadeInUp {
    from {
      opacity: 0;
      transform: translateY(10px);
    }
    to {
      opacity: 1;
      transform: translateY(0);
    }
  }
")),
    tags$style(HTML("
    #global-loading-text {
      color: #e0e7ff; /* bright, cool lavender-blue */
      font-size: 2rem;
      letter-spacing: 0.03em;
      text-shadow:
        0 0 4px rgba(129, 140, 248, 0.8),
        0 0 10px rgba(167, 139, 250, 0.7);
    }
  ")),
    tags$style(HTML("
    #global-loading-overlay.loading-overlay {
      position: fixed;
      inset: 0;
      display: flex;
      align-items: center;
      justify-content: center;
      background: radial-gradient(
          circle at top,
          rgba(16, 19, 34, 0.45),
          rgba(5, 7, 18, 0.35) 60%,
          rgba(2, 3, 9, 0.25) 100%
        );
      z-index: 9999;
    }

    .loading-core {
      display: flex;
      flex-direction: column;   /* stack logo above text */
      align-items: center;      /* center horizontally */
      justify-content: center;  /* center vertically within overlay */
      text-align: center;
      gap: 0.75rem;
      position: relative;
      z-index: 2;
    }

    ..sparkles {
      position: fixed;
      inset: 0;
      pointer-events: none;
      overflow: hidden;
      z-index: 1;
    }

    .sparkle {
      position: absolute;
      border-radius: 999px;
      background: radial-gradient(circle, #ffffff, #a78bfa);
      box-shadow:
        0 0 8px rgba(129, 140, 248, 0.9),
        0 0 18px rgba(167, 139, 250, 0.8);
      opacity: 0;
      animation: sparkleTwinkle 2.8s ease-in-out infinite;
    }

    @keyframes sparkleTwinkle {
      0%, 100% {
        opacity: 0;
        transform: scale(0.3);
      }
      35% {
        opacity: 1;
        transform: scale(1);
      }
      70% {
        opacity: 0.2;
        transform: scale(0.4);
      }
    }

    .sparkle:nth-child(1)  { top: 8%;   left: 12%; width: 6px;  height: 6px;  animation-delay: 0.1s; }
    .sparkle:nth-child(2)  { top: 22%;  left: 78%; width: 7px;  height: 7px;  animation-delay: 0.5s; }
    .sparkle:nth-child(3)  { top: 40%;  left: 18%; width: 5px;  height: 5px;  animation-delay: 0.9s; }
    .sparkle:nth-child(4)  { top: 68%;  left: 30%; width: 8px;  height: 8px;  animation-delay: 1.3s; }
    .sparkle:nth-child(5)  { top: 82%;  left: 70%; width: 6px;  height: 6px;  animation-delay: 1.7s; }
    .sparkle:nth-child(6)  { top: 15%;  left: 50%; width: 9px;  height: 9px;  animation-delay: 0.3s; }
    .sparkle:nth-child(7)  { top: 30%;  left: 35%; width: 4px;  height: 4px;  animation-delay: 0.8s; }
    .sparkle:nth-child(8)  { top: 55%;  left: 85%; width: 7px;  height: 7px;  animation-delay: 1.1s; }
    .sparkle:nth-child(9)  { top: 73%;  left: 12%; width: 10px; height: 10px; animation-delay: 1.6s; }
    .sparkle:nth-child(10) { top: 10%;  left: 88%; width: 5px;  height: 5px;  animation-delay: 0.7s; }
    .sparkle:nth-child(11) { top: 48%;  left: 60%; width: 9px;  height: 9px;  animation-delay: 1.0s; }
    .sparkle:nth-child(12) { top: 33%;  left: 5%;  width: 6px;  height: 6px;  animation-delay: 1.4s; }
    .sparkle:nth-child(13) { top: 63%;  left: 48%; width: 8px;  height: 8px;  animation-delay: 1.8s; }
    .sparkle:nth-child(14) { top: 25%;  left: 92%; width: 4px;  height: 4px;  animation-delay: 0.2s; }
    .sparkle:nth-child(15) { top: 88%;  left: 40%; width: 9px;  height: 9px;  animation-delay: 1.9s; }
    .sparkle:nth-child(16) { top: 5%;   left: 30%; width: 7px;  height: 7px;  animation-delay: 0.4s; }
    .sparkle:nth-child(17) { top: 52%;  left: 8%;  width: 6px;  height: 6px;  animation-delay: 1.2s; }
    .sparkle:nth-child(18) { top: 60%;  left: 92%; width: 10px; height: 10px; animation-delay: 1.5s; }
    .sparkle:nth-child(19) { top: 18%;  left: 65%; width: 5px;  height: 5px;  animation-delay: 0.6s; }
    .sparkle:nth-child(20) { top: 78%;  left: 55%; width: 8px;  height: 8px;  animation-delay: 1.1s; }
  ")),
    tags$style(HTML("
  .btn {
    border-radius: 20px;
    font-size: 12px;
    font-weight: 500;
    padding: 5px 14px;
  }
  .btn-primary {
    background-color: var(--radar-primary);
    border-color: var(--radar-primary);
  }
  .btn-primary:hover, .btn-primary:focus {
    background-color: var(--radar-primary-dark);
    border-color: var(--radar-primary-dark);
  }
")),
    tags$style(HTML("
  .form-group {
    margin-bottom: 10px;
  }
  .control-label {
    font-weight: 500;
    font-size: 12px;
    color: #5f6473;
  }
  .selectize-control.single .selectize-input,
  .vscomp-toggle-button {
    border-radius: 6px;
    border-color: #d3d7e5;
    min-height: 32px;
  }
")),
    tags$style(HTML("
  body, .content-wrapper, .box-title, .sidebar-menu > li > a {
    font-family: 'Poppins', -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif;
  }
  .box-title {
    font-weight: 500;
    font-size: 15px;
  }
  .main-header .logo {
    font-size: 18px;
    letter-spacing: 0.03em;
  }
")),
    tags$style(HTML("
  .content-wrapper, .right-side {
    background-color: #f5f6fa;
  }
")),
    tags$style(HTML("
    .box {
    border-radius: 10px;
    border: 1px solid #e0e4f0;
    box-shadow: 0 4px 12px rgba(0,0,0,0.04);
    overflow: visible;  /* keep dropdowns visible */
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
")),
    tags$style(HTML("
    :root {
      --radar-primary: #5c6bc0;      /* solid purple */
      --radar-primary-dark: #3f51b5; /* hover */
    }

    /* Primary boxes (status = 'primary') */
    .skin-purple .box.box-primary {
      border-top-color: var(--radar-primary);
    }
    .skin-purple .box.box-primary > .box-header {
      background-color: var(--radar-primary);
      color: #ffffff;
    }
    .skin-purple .box.box-primary > .box-header .box-title {
      color: #ffffff;
    }
    .skin-purple .box.box-primary > .box-body {
      background-color: #ffffff;
      color: #283046;
    }

    /* Primary buttons */
    .skin-purple .btn-primary {
      background-color: var(--radar-primary);
      border-color: var(--radar-primary);
    }
    .skin-purple .btn-primary:hover,
    .skin-purple .btn-primary:focus {
      background-color: var(--radar-primary-dark);
      border-color: var(--radar-primary-dark);
    }
  ")),
      tags$link(
        rel = "stylesheet",
        href = "https://fonts.googleapis.com/css2?family=Poppins:wght@500;600&display=swap"
      ),
      tags$style(HTML("
    .skin-purple .main-header .logo {
      font-family: 'Poppins', -apple-system, BlinkMacSystemFont, 'Segoe UI',
                   Roboto, 'Helvetica Neue', Arial, sans-serif;
      font-weight: 600;
      font-size: 100px;
      letter-spacing: 0.03em;
    }
  ")),
    # Global styles (loader + fade-in)
    tags$style(HTML("
  /* Apply gradient only to the navbar bar across the page */
  .skin-purple .main-header .navbar {
    background: linear-gradient(135deg, #8e24aa, #3f51b5, #29b6f6);
    border: none;
  }

  /* Make the logo tile a flat color that matches the gradient at the left */
  .skin-purple .main-header .logo {
    background-color: #8e24aa;
    color: #ffffff;
    font-weight: 600;
    font-size: 17px;
    border: none;
  }

  .skin-purple .main-header .logo:hover {
    background-color: #7b1fa2;
  }

  .box.box-primary {
    border-top-color = #8e24aa;
  }
  .box.box-primary > .box-header {
    background: #ffffff;
    color: #333;
  }
")),
    tags$style(HTML("
      .loading-overlay {
        position: fixed;
        top: 0; left: 0; right: 0; bottom: 0;
        background: rgba(255,255,255,0.85);
        z-index: 3000;
        display: none;
        align-items: center;
        justify-content: center;
        flex-direction: column;
        font-size: 18px;
        color: #555;
      }
 .loading-logo {
    position: relative;
    width: 80px;
    height: 80px;
    border-radius: 24px;
    background: linear-gradient(135deg, #8e24aa, #3f51b5, #29b6f6);
    display: flex;
    align-items: center;
    justify-content: center;
    box-shadow: 0 10px 30px rgba(0,0,0,0.3);
  }

  .loading-logo svg {
      display: block;
      width: 72px;
      height: 72px;
      margin: 0 auto;           /* ensure center in its div */
    }

  .hex-spiral {
    fill: none;
    stroke: #ffffff;
    stroke-width: 4.5;
    stroke-linecap: round;
    stroke-linejoin: round;
    stroke-dasharray: 260;
    stroke-dashoffset: 260;
    animation: drawHexSpiral 1.5s ease-in-out infinite;
  }

  @keyframes drawHexSpiral {
    0%   { stroke-dashoffset: 260; opacity: 0.0; }
    10%  { opacity: 1.0; }
    65%  { stroke-dashoffset: 0;   opacity: 1.0; }
    100% { stroke-dashoffset: 0;   opacity: 0.0; }
  }

      /* Fade-in animation for panels */
      .fade-in-panel {
        animation: fadeInPanel 0.4s ease-in-out;
      }
      @keyframes fadeInPanel {
        from { opacity: 0; transform: translateY(4px); }
        to   { opacity: 1; transform: translateY(0); }
      }

      /* Dropdowns above normal content */
      .vscomp-wrapper,
      .vscomp-wrapper .vscomp-dropdown {
        z-index: 2000 !important;
      }
      .selectize-control .selectize-dropdown {
        z-index: 2000 !important;
      }

      /* (optional) disable pointer events class if you use it elsewhere */
      .radar-loading-disabled {
        pointer-events: none !important;
      }
    ")),

    # Loader + fade-in JS handlers
    tags$script(HTML("
      Shiny.addCustomMessageHandler('toggle-loading', function(show) {
        var overlay = document.getElementById('global-loading-overlay');
        if (overlay) {
          overlay.style.display = show ? 'flex' : 'none';
        }
      });

      Shiny.addCustomMessageHandler('loading-text', function(msg) {
        var el = document.getElementById('global-loading-text');
        if (!el) return;
        el.textContent = msg;
      });

      // Fade-in any box containing an updated output
      $(document).on('shiny:value', function(e) {
        var $box = $('#' + e.target.id).closest('.box');
        if ($box.length) {
          $box.addClass('fade-in-panel');
          setTimeout(function() {
            $box.removeClass('fade-in-panel');
          }, 500);
        }
      });
    ")),

    # Highlight subsystem handler (plot1 tick bold)
    tags$script(HTML("
      Shiny.addCustomMessageHandler('highlight-subsystem', function(idx) {
        var ticks = $('.shiny-plot-output#plot1').find('text');
        ticks.css('font-weight', 'normal');
        if (idx > 0 && idx <= ticks.length) {
          $(ticks[idx-1]).css('font-weight', 'bold');
        }
      });
    "))
  ),

  # Full-page overlay element
  tags$div(
    id = "global-loading-overlay",
    class = "loading-overlay",

    tags$div(
      class = "loading-core",
      # logo on top
      tags$div(
        class = "loading-logo",
        tags$svg(
          viewBox = "0 0 100 100",
          tags$path(
            class = "hex-spiral",
            d = paste(
              "M 50 14",
              "L 78 30",
              "L 78 56",
              "L 50 74",
              "L 22 56",
              "L 22 30",
              "L 46 22",
              "L 70 34",
              "L 70 52",
              "L 50 64",
              "L 30 52",
              "L 30 34",
              "L 50 26",
              "L 62 38",
              "L 62 48",
              "L 50 54",
              sep = " "
            )
          )
        )
      ),
      # text under logo
      tags$div(
        id = "global-loading-text",
        "Preparing RADAR | xCheck analysis…"
      )
    ),

    tags$div(
      class = "sparkles",
      lapply(1:20, function(i) tags$span(class = "sparkle"))
    )
  ),
  tabItems(tabItem(
    tabName = "obs",
  fluidRow(
    style = "margin-top: 10px;",
    column(
      width = 12,
      box(
        width = 12,
        title = "Select Cohort",
        solidHeader = TRUE, status = "primary",
        box(
          width = 4,
          selectInput(
            label   = "Select Cohort",
            inputId = "cohort_select",
            # value = folder_name, label = "cohort_name - disease"
            choices = stats::setNames(
              cohort_registry$folder_name,
              paste0(cohort_registry$cohort_name, " - ", cohort_registry$disease)
            ),
            selected = cohort_registry$folder_name[1]
          )
        )
      )
    )
  ),

  fluidRow(
    class = "fade-in-up",
    style = "margin-top: 10px;",
    column(
      width = 12,
      box(
        width = 12,
        title = "Filter by Stratification Attributes",
        solidHeader = TRUE,
        status = "primary",

        # Stratification type
        box(
          width = 4,
          selectInput(
            label   = "Stratification type (optional)",
            inputId = "stratType",
            choices = c(
              "Three-Class Continuous (25th / 75th percentile)" = 1,
              "Two-Class Continuous (median split)"             = 2
            ),
            selected = 1
          )
        ),

        # Clinical stratification
        box(
          width = 4,
          virtualSelectInput(
            label   = "Clinical stratification (leave empty for no stratification)",
            inputId = "clinStrat",
            choices = NULL,
            showValueAsTags = TRUE,
            search = TRUE,
            multiple = TRUE,
            placeholder = "No clinical stratification (optional)"
          )
        ),

        # Gene stratification
        box(
          width = 4,
          virtualSelectInput(
            label   = "Gene expression stratification (leave empty for no stratification)",
            inputId = "geneStrat",
            choices = NULL,
            showValueAsTags = TRUE,
            search = TRUE,
            multiple = TRUE,
            placeholder = "No gene stratification (optional)"
          )
        ),

        # Info text + submit button
        fluidRow(
          style = "margin-top: 10px;",
          column(
            width = 8,
            tags$div(
              style = "padding: 10px 15px; color: #555;",
              tags$strong("Note: "),
              "If you leave both stratification fields empty, the full cohort will be used (no stratification)."
            )
          ),
          column(
            width = 4,
            div(
              style = "text-align: right; padding-right: 15px;",
              actionButton("submitStrat", "Apply stratification")
            )
          )
        )
      ),
      tags$div(id = "stratification")
    )
  ),

  conditionalPanel(
    condition = "input.submitReport > 0 || input.submitReportDefault > 0",

    # Row 1: plot + guide + downloads
    fluidRow(
      class = "fade-in-up",
      style = "margin-top: 10px;",
      column(
        width = 12,
        box(
          width = 12, status = "primary", solidHeader = TRUE,
          title = textOutput("plot1Title"),
          plotlyOutput("plot1_obs", height = 650),

          tags$div(
            style = "margin-top: 10px; padding: 10px 12px; border-radius: 6px;
                   background-color: #e8ebff;
                   border: 1px solid #d1d7ff;
                   font-size: 13px; line-height: 1.4; color: #273043;",
            tags$div(
              style = "display: flex; align-items: flex-start; gap: 8px;",
              tags$span(
                "\u2139",
                style = "font-size: 18px; line-height: 1;
                       color: var(--radar-primary);"
              ),
              tags$div(
                tags$div(
                  "How to use this plot:",
                  style = "font-weight: 600; margin-bottom: 3px;"
                ),
                tags$ul(
                  style = "padding-left: 18px; margin: 0;",
                  tags$li("Enlarge regions of interest using the zoom and pan tools."),
                  tags$li("Use the box or lasso selection tools in the Plotly toolbar to select reactions."),
                  tags$li("Click on selection to analyse differential flux across the selected comparisons.")
                )
              )
            )
          ),

          downloadLink("dlPlot1", "Download Plot as PDF"),
          div(
            style = "margin-top: 10px;",
            downloadButton(
              "reacts_obs",
              "Download Fingerprint Preparation Data",
              style = "width: 100%;"
            )
          )
        )
      )
    ),

    # Row 2: outer/inner selection
    fluidRow(
      class = "fade-in-up",
      column(
        width = 4,
        style = "display: flex;",  # make the column stretch
        box(
          title = textOutput("outerTitle"),
          solidHeader = TRUE, status = "primary",
          width = 12,
          style = "flex: 1;",      # box fills the column height
          checkboxGroupInput(
            "outerSelect",
            label = NULL,
            choices = c(),
            selected = NULL
          )
        )
      ),

      column(
        width = 4,
        style = "display: flex;",
        box(
          title = textOutput("innerTitle"),
          solidHeader = TRUE, status = "primary",
          width = 12,
          style = "flex: 1;",
          checkboxGroupInput(
            "innerSelect",
            label = NULL,
            choices = c("high", "baseline", "low")
          )
        )
      )
    )
  ),

  conditionalPanel(
    condition = "output.selected_points_n_obs > 0",

    fluidRow(
      class = "fade-in-up",
      column(
        style = "display:none;",
        width = 12,
        box(
          width = 12,
          selectizeInput(
            "boxplot",
            "Select Metabolic Flux to Analyze",
            choices = c("")
          )
        )
      )
    ),

    fluidRow(
      class = "fade-in-up",
      column(
        width = 12,
        box(
          width = 12,
          reactableOutput("info_obs")
        )
      )
    ),

    fluidRow(
      class = "fade-in-up",
      column(
        width = 12,
        box(
          width = 12,
          plotOutput("plot2_obs", height = 500),
          actionLink("invert", "Invert Boxplot | "),
          downloadLink("dlplot2_obs", "Download Plot as PDF")
        )
      )
    )
  )
  ),
  tabItem(
    tabName = "exp",

    h2("Metabolic Fingerprint Designer for Experiments"),

    # 1) Experiment + stratification
    fluidRow(
      column(
        width = 12,
        box(
          width = 12,
          title = "Select Experiment and Stratification Attributes",
          solidHeader = TRUE,
          status = "primary",
          column(
            width = 6,
            selectInput(
              label   = "Select Experiment",
              inputId = "experiment_select",              ### CHANGED ID
              choices = stats::setNames(
                experiment_registry$folder_name,
                paste0(experiment_registry$experiment_name, " | ", experiment_registry$treatment)
              ),
              selected = experiment_registry$folder_name[1]
            )
          ),
          column(
            width = 6,
            virtualSelectInput(
              label   = "Select Stratification",
              inputId = "strat_exp",                  ### CHANGED ID
              choices = c(),
              showValueAsTags = TRUE,
              search  = TRUE,
              multiple = TRUE
            )
          ),
          fluidRow(
            column(
              width = 12,
              actionButton("submitStrat_exp", "Submit"),   ### CHANGED ID
                            ### CHANGED ID
            )
          )
        )
      )
    ),div(id = "stratification_exp"),
  conditionalPanel(
    condition = "input.submitReport_exp > 0 || input.submitReportDefault_exp > 0",
    # 2) Plots 1 & 2, downloads
    fluidRow(
      class = "fade-in-up",
      style = "margin-top: 10px;",
      column(
        width = 12,
        box(
          width = 12, status = "primary", solidHeader = TRUE,
          title = textOutput("plot1Title_exp"),
          plotlyOutput("plot1_exp", height = 650),
          
          tags$div(
            style = "margin-top: 10px; padding: 10px 12px; border-radius: 6px;
                   background-color: #e8ebff;
                   border: 1px solid #d1d7ff;
                   font-size: 13px; line-height: 1.4; color: #273043;",
            tags$div(
              style = "display: flex; align-items: flex-start; gap: 8px;",
              tags$span(
                "\u2139",
                style = "font-size: 18px; line-height: 1;
                       color: var(--radar-primary);"
              ),
              tags$div(
                tags$div(
                  "How to use this plot:",
                  style = "font-weight: 600; margin-bottom: 3px;"
                ),
                tags$ul(
                  style = "padding-left: 18px; margin: 0;",
                  tags$li("Enlarge regions of interest using the zoom and pan tools."),
                  tags$li("Use the box or lasso selection tools in the Plotly toolbar to select reactions."),
                  tags$li("Click on selection to analyse differential flux across the selected comparisons.")
                )
              )
            )
          ),
          
          downloadLink("dlPlot1", "Download Plot as PDF"),
          div(
            style = "margin-top: 10px;",
            downloadButton(
              "reacts_exp",
              "Download Fingerprint Preparation Data",
              style = "width: 100%;"
            )
          )
        )
      )
    )
    ),

  conditionalPanel(
    condition = "output.selected_points_n_exp > 0",
    
    fluidRow(
      class = "fade-in-up",
      column(
        style = "display:none;",
        width = 12,
        box(
          width = 12,
          selectizeInput(
            "boxplot",
            "Select Metabolic Flux to Analyze",
            choices = c("")
          )
        )
      )
    ),
    
    fluidRow(
      class = "fade-in-up",
      column(
        width = 12,
        box(
          width = 12,
          reactableOutput("info_exp")
        )
      )
    ),
    
    fluidRow(
      class = "fade-in-up",
      column(
        width = 12,
        box(
          width = 12,
          plotOutput("plot2_exp", height = 500),
          actionLink("invert", "Invert Boxplot | "),
          downloadLink("dlplot2_exp", "Download Plot as PDF")
        )
      )
    )
  )
  ))


)


sidebar <- dashboardSidebar(  ### ADDED (replaces disable=TRUE)
  sidebarMenu(
    id = "main_tabs",
    menuItem(
      "Observational Studies",
      tabName = "obs",
      icon = icon("project-diagram")
    ),
    menuItem(
      "Experiments",
      tabName = "exp",
      icon = icon("flask")
    )
  )
)

ui<-dashboardPage(
  skin = "purple",
  header = header,
  sidebar = sidebar,
  body = body
)

server <-function(input,output,session){
  cohort_name<-reactiveVal(NULL)
  cohort_flux<-reactiveVal(NULL)
  experiment_flux<-reactiveVal(NULL)
  gene<-reactiveVal(NULL)
  clinical<-reactiveVal(NULL)
  meta_clinical<-reactiveVal(NULL)

  choices_gene_2_class<-reactiveVal(NULL)
  choices_gene_3_class<-reactiveVal(NULL)
  choices_clinical_2_class<-reactiveVal(NULL)
  choices_clinical_3_class<-reactiveVal(NULL)

  type_strat<-reactiveVal(NULL)
  clinical_strat<-reactiveVal(NULL)
  gene_strat<-reactiveVal(NULL)

  type_outer<-reactiveVal(NULL)
  outerVarFor<-reactiveVal(NULL)
  outerVarRev<-reactiveVal(NULL)

  type_inner<-reactiveVal(NULL)
  inner<-reactiveVal(NULL)

  upper_label_rv <- reactiveVal("upper")
  lower_label_rv <- reactiveVal("lower")

  metaFinal_out <- reactiveVal(NULL)

  errors_list<-reactiveVal(NULL)

  strat_ui_inserted  <- reactiveVal(FALSE)
  outer_ui_inserted  <- reactiveVal(FALSE)
  inner_ui_inserted  <- reactiveVal(FALSE)
  summary_ui_inserted <- reactiveVal(FALSE)

  metabs_tab_cohort <- reactiveVal(NULL)
  metabs_tab_exp <- reactiveVal(NULL)
  metabs_tab <- reactiveVal(NULL)

  analysis_ready <- reactiveVal(FALSE)
  analysis_ready_exp <- reactiveVal(FALSE)

  experiment_name   <- reactiveVal(NULL)
  experiment_flux   <- reactiveVal(NULL)
  experiment_meta     <- reactiveVal(NULL)
  
  exp_strat_ui_inserted  <- reactiveVal(FALSE)
  exp_outer_ui_inserted  <- reactiveVal(FALSE)
  exp_inner_ui_inserted  <- reactiveVal(FALSE)
  exp_summary_ui_inserted <- reactiveVal(FALSE)
  
  strat_exp <- reactiveVal(NULL)
  metaFinal_out <- reactiveVal(NULL)
  outer_exp            <- reactiveVal(NULL)
  inner_exp            <- reactiveVal(NULL)
  top_exp              <- reactiveVal(NULL)
  bottom_exp           <- reactiveVal(NULL)
  innerVars_exp        <- reactiveVal(NULL)
  error_exp            <- reactiveVal(NULL)
  customAUC_exp        <- reactiveVal(NULL)
  customFC_exp         <- reactiveVal(NULL)
  customFDR_exp        <- reactiveVal(NULL)

  radar_theme <- theme_minimal(base_family = "Poppins") +
    theme(
      panel.grid.major = element_line(color = "#e2e6f3", size = 0.3),
      panel.grid.minor = element_blank(),
      axis.title       = element_text(size = 11, color = "#444"),
      axis.text        = element_text(size = 9, color = "#555"),
      strip.background = element_rect(fill = "#f0f2fa", color = NA),
      strip.text       = element_text(size = 9, face = "bold"),
      plot.title       = element_text(size = 12, face = "bold", hjust = 0.5),
      legend.position  = "none"
    )

  run_analysis <- function(AUC_cutoff, FC_cutoff, FDR_cutoff) {

  # ---- 0) Global checks ----
  if (length(errors_list()) > 0) {
    showModal(modalDialog(
      title = "Cannot run analysis, paste-errors present",
      paste(errors_list(), collapse = "\n"),
      easyClose = TRUE
    ))
    return(invisible(NULL))
  }

  metaFinal <- metaFinal_out()

  if (is.null(metaFinal) || nrow(metaFinal) == 0) {
    showModal(modalDialog(
      title = "No filtered cohort",
      "Please configure stratification, outer and inner comparisons, then run the pre-analysis summary before starting the analysis.",
      easyClose = TRUE
    ))
    return(invisible(NULL))
  }

  # ---- 1) Basic objects ----
  flux_mat <- if (input$main_tabs == "exp") experiment_flux() else cohort_flux()
  fluxfull <- cbind(metaFinal, flux_mat[rownames(metaFinal), ])
  fluxcomp <- fluxfull

  outer <- colnames(metaFinal)[1]
  inner <- colnames(metaFinal)[2]

  outer_for_vec <- if (!is.null(outerVarFor())) outerVarFor()[[1]] else character(0)
  outer_rev_vec <- if (!is.null(outerVarRev())) outerVarRev()[[1]] else character(0)
  outer_all     <- c(outer_for_vec, outer_rev_vec)

  parse_outer <- function(x) {
    parts <- strsplit(x, " - ", fixed = TRUE)[[1]]
    if (length(parts) >= 2) {
      list(var = parts[1], level = parts[2])
    } else {
      list(var = x, level = x)
    }
  }

  outer_comp_name <- if (length(outer_all) > 0) parse_outer(outer_all[1])$var else outer_exp()
  upper_label     <- if (length(outer_for_vec) > 0) parse_outer(outer_for_vec[1])$level else top_exp()
  lower_label     <- if (length(outer_rev_vec) > 0) parse_outer(outer_rev_vec[1])$level else bottom_exp()

  output$outerTitle <- renderText(outer)
  output$innerTitle <- renderText(inner)
  output$plot1Title <- renderText(
    paste0(
      "Significant Reactions Across All Subsystems - ",
      outer_comp_name, " (", upper_label, " v. ", lower_label, ")"
    )
  )
  
  output$plot1Title_exp <- renderText(
    paste0(
      "Significant Reactions Across All Subsystems - ",
      outer_exp(), " (", top_exp(), " v. ", bottom_exp(), ")"
    )
  )

  # ---- 2) ROC / t‑test / FDR loop per inner group ----
  outer_select <- c("upper", "lower")
  inner_opts   <- unique(metaFinal[, inner])

  metab_names <- c()
  data_filt   <- fluxcomp

  for (i in seq_along(inner_opts)) {

    data_filt_var <- subset(data_filt, inner == inner_opts[i])

    AUC_col <- caTools::colAUC(
      X       = data_filt_var[, -(1:ncol(metaFinal)), drop = FALSE],
      y       = factor(data_filt_var[, outer]),
      plotROC = FALSE
    )

    AUC_col <- as.numeric(AUC_col[1, ])

    flux_idx <- seq_len(ncol(data_filt_var))[(ncol(metaFinal) + 1):ncol(data_filt_var)]

    upper_mat <- subset(data_filt_var, outer == "upper")[, flux_idx, drop = FALSE]
    lower_mat <- subset(data_filt_var, outer == "lower")[, flux_idx, drop = FALSE]

    # guard: if no rows or no columns, skip this inner level
    if (nrow(upper_mat) == 0 || nrow(lower_mat) == 0 || ncol(upper_mat) == 0) {
      next
    }

    logFC <- log2(
      colMeans(upper_mat, na.rm = TRUE) /
        colMeans(lower_mat, na.rm = TRUE)
    )

    short <- as.data.frame(logFC[which(abs(logFC) > FC_cutoff & AUC_col >= AUC_cutoff)])

    t_col <- matrixTests::col_t_welch(
      x = as.matrix(subset(data_filt_var, outer == "upper")[, rownames(short), drop = FALSE]),
      y = as.matrix(subset(data_filt_var, outer == "lower")[, rownames(short), drop = FALSE])
    )

    final <- cbind(t_col$mean.y, t_col$mean.x, short, t_col$statistic, t_col$pvalue)
    colnames(final) <- c("var_1_mean", "var_2_mean", "log2FC", "t_welch_statistic", "pval")
    final <- subset(final, subset = pval < 0.05)

    final$fdr <- p.adjust(final$pval, method = "fdr")
    final     <- subset(final, subset = fdr < FDR_cutoff)

    metab_names <- c(metab_names, rownames(final))
  }

  metab_names <- unique(metab_names)

  # ---- 3) Build full test list for all inner_opts ----
  test_list <- data.frame()

  for (i in seq_along(inner_opts)) {
    nm <- c("outer", metab_names)
    data_filt_var <- subset(data_filt, subset = inner == inner_opts[i])
    data_filt_var <- data_filt_var[, nm]

    t_col <- matrixTests::col_t_welch(
      as.matrix(subset(data_filt_var, subset = outer == "upper")[, metab_names]),
      as.matrix(subset(data_filt_var, subset = outer == "lower")[, metab_names])
    )

    final <- as.data.frame(cbind(
      metabolite       = metab_names,
      inner_comparison = inner_opts[i],
      var_1_mean       = as.numeric(t_col$mean.y),
      var_2_mean       = as.numeric(t_col$mean.x),
      t_welch_statistic = as.numeric(t_col$statistic),
      pval             = as.numeric(t_col$pvalue)
    ))

    colnames(final) <- c("metabolite", "inner_comparison",
                         "var_1_mean", "var_2_mean",
                         "t_welch_statistic", "pval")

    final$sign_t_log_pval <- ifelse(
      final$t_welch_statistic < 0,
      -(-log10(as.numeric(final$pval))),
      -log10(as.numeric(final$pval))
    )

    test_list <- rbind(test_list, final)
  }

  # ---- 4) Merge with metadata, set sig, inner_comparison factor, etc. ----
  metabstab <- merge(test_list, reactMeta, by = "metabolite", all.x = TRUE, all.y = FALSE)

  upper_label_rv(upper_label)
  lower_label_rv(lower_label)

  metabstab$significant_category <- ifelse(
    metabstab$sign_t_log_pval > 1.301, upper_label,
    ifelse(metabstab$sign_t_log_pval < -1.301, lower_label, "ns")
  )

  # Ensure inner_comparison exists and has correct levels
  if(!is.null(input$innerType)){
  metabstab$inner_comparison <- switch(
    as.numeric(input$innerType),
    { # 1: three‑class continuous
      factor(metabstab$inner_comparison, levels = c("high", "baseline", "low"))
    },
    { # 2: two‑class continuous
      factor(metabstab$inner_comparison, levels = c("high", "low"))
    },
    { # 3: binary
      factor(metabstab$inner_comparison, levels = unique(metabstab$inner_comparison))
    },
    { # 4: categorical (or fallback)
      factor(metabstab$inner_comparison, levels = unique(metabstab$inner_comparison))
    }
  )}

  if (is.null(inner()) || inner() == "" || as.numeric(input$innerType) == 0) {
    metabstab$inner_comparison <- "All"
  }

  updateCheckboxGroupInput(
    session, "innerSelect",
    choices  = unique(data_filt[, inner]),
    selected = unique(data_filt[, inner])
  )
  # ---- 5) Store for downstream reactives ----

  metabs_tab(metabstab)  # or use reactiveVal instead of <<-
}

  observeEvent(input$cohort_select, {
    req(input$cohort_select)

    session$sendCustomMessage("toggle-loading", TRUE)
    session$sendCustomMessage("loading-text", "Reading cohort registry…")
    on.exit({
      session$sendCustomMessage("toggle-loading", FALSE)
      session$sendCustomMessage("loading-text", "Cohort loaded successfully")
    }, add = TRUE)

    row <- subset(cohort_registry, folder_name == input$cohort_select)

    cohort_name(paste0(row$cohort_name, " - ", row$disease))

    prefix <- row$folder_name

    flux_file          <- paste0("../../../data/RADAR_xCheck_cohort/",file.path(prefix, paste0(prefix, "_flux.csv")))
    gene_file          <- paste0("../../../data/RADAR_xCheck_cohort/",file.path(prefix, paste0(prefix, "_norm_filtered.csv")))
    clinical_file      <- paste0("../../../data/RADAR_xCheck_cohort/",file.path(prefix, paste0(prefix, "_clinical.csv")))
    meta_clinical_file <- paste0("../../../data/RADAR_xCheck_cohort/",file.path(prefix, paste0(prefix, "_meta_clinical.csv")))

    session$sendCustomMessage("loading-text", "Loading flux matrices…")
    flux <- as.data.frame(as.matrix(data.table::fread(flux_file), rownames = 1))

    session$sendCustomMessage("loading-text", "Loading gene expression data…")
    gene_mat <- as.data.frame(as.matrix(data.table::fread(gene_file, header = TRUE), rownames = 1))

    session$sendCustomMessage("loading-text", "Loading clinical annotations…")
    clinical_df <- read.csv(clinical_file, header = TRUE, row.names = 1)
    meta_clin_df <- read.csv(meta_clinical_file, header = TRUE)

    session$sendCustomMessage("loading-text", "Finalizing cohort setup…")
    cohort_flux(flux)
    gene(gene_mat)
    clinical(clinical_df)
    meta_clinical(meta_clin_df)


  })

  observeEvent(list(input$cohort_select,input$stratType),{

      req(gene(), clinical(), meta_clinical())

    session$sendCustomMessage("toggle-loading", TRUE)
    session$sendCustomMessage("loading-text", "Preparing stratification options…")
    on.exit({
      session$sendCustomMessage("toggle-loading", FALSE)
      session$sendCustomMessage("loading-text", "Stratifications loaded successfully")
    }, add = TRUE)

      genes <- colnames(gene())
      n_genes <- length(genes)

      session$sendCustomMessage("loading-text", "Building 2-class gene stratifications…")
      results_genes_2 <- vector("list", n_genes)
      names(results_genes_2) <- genes

      for (i in seq_along(genes)) {
        results_genes_2[[i]] <- c(
          paste0(genes[i], " - high"),
          paste0(genes[i], " - low")
        )
      }
      choices_gene_2_class(results_genes_2)

      session$sendCustomMessage("loading-text", "Building 3-class gene stratifications…")
      results_genes_3 <- vector("list", n_genes)
      names(results_genes_3) <- genes
      for (i in seq_along(genes)) {
        results_genes_3[[i]] <- c(
          paste0(genes[i], " - high"),
          paste0(genes[i], " - baseline"),
          paste0(genes[i], " - low")
        )
      }
      choices_gene_3_class(results_genes_3)

      # Clinical stratifications
      session$sendCustomMessage("loading-text", "Building clinical stratifications…")
      clinicalVars <- colnames(clinical())
      n_clin <- length(clinicalVars)

      results_clinical_2 <- list()
      results_clinical_3 <- list()

      for (i in seq_along(clinicalVars)) {
        var <- clinicalVars[i]
        if (var %in% subset(meta_clinical(), type == "continuous")$variable) {
          results_clinical_2[[var]] <- c(paste0(var, " - high"),
                                         paste0(var, " - low"))
          results_clinical_3[[var]] <- c(paste0(var, " - high"),
                                         paste0(var, " - baseline"),
                                         paste0(var, " - low"))
        } else {
          vals <- na.omit(unique(clinical()[[var]]))
          labels <- paste0(var, " - ", vals)
          results_clinical_2[[var]] <- labels
          results_clinical_3[[var]] <- labels
        }
      }

      choices_clinical_2_class(results_clinical_2)
      choices_clinical_3_class(results_clinical_3)

      session$sendCustomMessage("loading-text", "Updating stratification controls…")
      # Final UI updates (coarse progress to 1.0)
      if (as.numeric(input$stratType) == 1) {
        updateVirtualSelect("clinStrat", choices = choices_clinical_3_class())
        updateVirtualSelect("geneStrat", choices = choices_gene_3_class())
      } else {
        updateVirtualSelect("clinStrat", choices = choices_clinical_2_class())
        updateVirtualSelect("geneStrat", choices = choices_gene_2_class())
      }
  })

  observeEvent(input$submitStrat, {
    req(input$submitStrat > 0)
    if (!isTRUE(outer_ui_inserted())){
      insertUI(
        selector = "#stratification",
        where = "afterEnd",
        immediate = TRUE,
        ui = fluidRow(class = "fade-in-up",column(width=6,box(width=12,title="Select Outer Comparison",solidHeader = TRUE, status="primary",
                                         box(width = 12, selectInput(label="Select Outer Comparison Type", inputId = "outerType", choices = c("Three-Class Continuous (from 25th and 75th percentiles)"=1,"Two-Class Continuous (from median)"=2, "Two-Class Binary (yes/no, positive/negative)"=3, "Categorical"=4),selected=1)),
                                         box(width = 12, selectizeInput(label="Select Variable Type", inputId = "outerVarType", choices = c("Clinical Data"=1,"Gene Expression"=2),selected=1)),
                                         box(width = 6, virtualSelectInput(label="Select Forward Variable(s)", inputId = "outerVarFor", choices = colnames(clinical),showValueAsTags = TRUE, search = TRUE, multiple = TRUE)),
                                         box(width = 6, virtualSelectInput(label="Select Inverted Variable(s)", inputId = "outerVarRev", choices = colnames(clinical),showValueAsTags = TRUE, search = TRUE, multiple = TRUE)),
                                         fluidRow(column(width=12,actionButton("submitOuter", "Submit outer settings"),tags$div(id = 'outerButton'))))),tags$div(id = 'outer')

        )
      )
      outer_ui_inserted(TRUE)
    }
    type_strat(as.numeric(input$stratType))
    clinical_strat(list(input$clinStrat))
    gene_strat(list(input$geneStrat))
  })

  observeEvent(input$outerType,{
    req(input$outerType)
    if(input$outerType == 1 | input$outerType == 2){
      updateSelectizeInput(session, "outerVarType", choices = c("Clinical Data"=1,"Gene Expression"=2),selected=1)
    }else if(input$outerType == 3){
      updateSelectizeInput(session, "outerVarType", choices = c("Clinical Data"=1),selected=1)
    }else{
      updateSelectizeInput(session, "outerVarType", choices = c("Clinical Data"=1),selected=1)
    }
  })

  observeEvent(list(input$cohort_select,input$outerType,input$outerVarType),{
    req(input$outerType)
    session$sendCustomMessage("toggle-loading", TRUE)
    session$sendCustomMessage("loading-text", "Reading cohort registry…")
    on.exit({
      session$sendCustomMessage("toggle-loading", FALSE)
      session$sendCustomMessage("loading-text", "Cohort loaded successfully")
    }, add = TRUE)
    session$sendCustomMessage("loading-text", "Updating outer parameter choices…")
    if(input$outerType == 1){
      if(input$outerVarType == 1){
        updateVirtualSelect("outerVarFor", label = "Select Upper Variable(s)", choices = choices_clinical_3_class()[subset(meta_clinical(), type == "continuous")$variable])
        updateVirtualSelect("outerVarRev", label = "Select Lower Variable(s)", choices = choices_clinical_3_class()[subset(meta_clinical(), type == "continuous")$variable])
      }else{
        updateVirtualSelect("outerVarFor", label = "Select Upper Variable(s)", choices = choices_gene_3_class())
        updateVirtualSelect("outerVarRev", label = "Select Lower Variable(s)", choices = choices_gene_3_class())
      }
    }else if(input$outerType == 2){
      if(input$outerVarType == 1){
        updateVirtualSelect("outerVarFor", label = "Select Upper Variable(s)", choices = choices_clinical_2_class()[subset(meta_clinical(), type == "continuous")$variable])
        updateVirtualSelect("outerVarRev", label = "Select Lower Variable(s)", choices = choices_clinical_2_class()[subset(meta_clinical(), type == "continuous")$variable])
      }else{
        updateVirtualSelect("outerVarFor", label = "Select Upper Variable(s)", choices = choices_gene_2_class())
        updateVirtualSelect("outerVarRev", label = "Select Lower Variable(s)", choices = choices_gene_2_class())
      }
    }else if(input$outerType == 3){
      updateVirtualSelect("outerVarFor", label = "Select Upper Variable(s)", choices = choices_clinical_3_class()[subset(meta_clinical(), type == "binary")$variable])
      updateVirtualSelect("outerVarRev", label = "Select Lower Variable(s)", choices = choices_clinical_3_class()[subset(meta_clinical(), type == "binary")$variable])

    }else{
      updateVirtualSelect("outerVarFor", label = "Select Upper Variable","outerVarType", choices = choices_clinical_3_class()[subset(meta_clinical(), type == "categorical")$variable])
      updateVirtualSelect("outerVarRev", label = "Select Lower Variable","outerVarType", choices = choices_clinical_3_class()[subset(meta_clinical(), type == "categorical")$variable])
    }

  })

  observeEvent(input$submitOuter, {
    req(input$submitOuter > 0)

    if (is.null(input$outerVarFor) && is.null(input$outerVarRev)) {
      showModal(
        modalDialog(
          title = "Missing outer comparison",
          "At least one outer comparison needs to be selected for upper/lower.",
          easyClose = TRUE,
          footer = modalButton("OK")
        )
      )
      return()  # do nothing else
    }

    if (!isTRUE(inner_ui_inserted())){
      insertUI(
        selector  = "#outer",
        where     = "beforeEnd",
        immediate = TRUE,
        ui = column(
          class = "fade-in-up",
          width = 6,
          box(
            width = 12,
            title = "Select Inner Comparison (optional)",
            solidHeader = TRUE,
            status = "primary",

            # Inner comparison type, with an explicit "None" option
            box(
              width = 12,
              selectInput(
                label   = "Inner comparison type",
                inputId = "innerType",
                choices = c(
                  "None (no inner comparison)"                            = 0,
                  "Three-Class Continuous (from 25th and 75th percentiles)" = 1,
                  "Two-Class Continuous (from median)"                      = 2,
                  "Two-Class Binary (yes/no, positive/negative)"           = 3,
                  "Categorical"                                            = 4
                ),
                selected = 0
              )
            ),

            # Variable type – disabled / ignored when innerType = 0
            box(
              width = 12,
              selectInput(
                label   = "Inner variable type",
                inputId = "innerVarType",
                choices = c("Clinical Data" = 1, "Gene Expression" = 2),
                selected = 1
              )
            ),

            # Inner variable – optional, with a clear 'none' placeholder
            box(
              width = 12,
              virtualSelectInput(
                label       = "Inner variable (leave empty for none)",
                inputId     = "innerVar",
                choices     = choices_clinical_3_class(),
                showValueAsTags = TRUE,
                search      = TRUE,
                multiple    = FALSE,
                placeholder = "No inner comparison (default)"
              )
            ),

            fluidRow(
              column(
                width = 12,
                actionButton("submitInner", "Submit inner settings")
              )
            )
          )
        ),
        tags$div(id = "inner")
      )
      inner_ui_inserted(TRUE)
    }
    type_outer(input$outerType)
    outerVarFor(list(input$outerVarFor))
    outerVarRev(list(input$outerVarRev))
  })

  observeEvent(input$innerType,{
    req(input$innerType)
    if(input$innerType == 1 | input$innerType == 2){
      updateSelectizeInput(session, "innerVarType", choices = c("Clinical Data"=1,"Gene Expression"=2),selected=1)
    }else if(input$innerType == 3){
      updateSelectizeInput(session, "innerVarType", choices = c("Clinical Data"=1),selected=1)
    }else{
      updateSelectizeInput(session, "innerVarType", choices = c("Clinical Data"=1),selected=1)
    }
  })

  observeEvent(list(input$cohort_select,input$innerType,input$innerVarType),{
    req(input$innerType)
    if(input$innerType == 1 | input$innerType == 2){
      if(input$innerVarType == 1){
        updateVirtualSelect("innerVar", label = "Select Variable", choices = subset(meta_clinical(), type == "continuous")$variable)
      }else{
        updateVirtualSelect("innerVar", label = "Select Variable", choices = colnames(gene()))
      }
    }else if(input$innerType == 3){
      updateVirtualSelect("innerVar", label = "Select Variable", choices = subset(meta_clinical(), type == "binary")$variable)
    }else{
      updateVirtualSelect("innerVar", label = "Select Variable", choices = subset(meta_clinical(), type == "categorical")$variable)
    }
  })

  observeEvent(input$submitInner, {
    req(input$submitInner > 0)
    if (!isTRUE(summary_ui_inserted())){
      insertUI(
        selector = "#outer",
        where = "afterEnd",
        immediate = TRUE,
        ui = fluidRow(column(width=12,box(width=12,title="RADAR | xCheck Pre-Analysis Summary",solidHeader = TRUE, status="primary",
                                          box(width = 6, verbatimTextOutput('reportText')),
                                          box(width = 6, noUiSliderInput(inputId = "customAUC", label = "AUC:", min = 0.6, max = 1.0, value = 0.65, tooltip = TRUE, step = 0.01),
                                              noUiSliderInput(inputId = "customFC", label = "Log2FC:", min = 0, max = 5, value = 1, tooltip = TRUE, step = 0.01),
                                              noUiSliderInput(inputId = "customFDR", label = "FDR:", min = 0, max = 1, value = 0.05, tooltip = TRUE, step = 0.01)),
                                          fluidRow(column(width=12,actionButton("submitReport", "Begin Analysis (set parameters)")),column(width=12,actionButton("submitReportDefault", "Begin Analysis (p < 0.05 only)")))))

        )
      )
      summary_ui_inserted(TRUE)
    }
      type_inner(input$innerType)
      inner(input$innerVar)

  })

  observeEvent(
    list(input$submitStrat, input$submitOuter, input$submitInner),
    {
      # Require that stratification, outer and inner settings exist
      req(type_strat(), type_outer(), type_inner())

      warnings <- character(0)
      errors   <- character(0)

      ########################
      # 1. STRATIFIED COHORT #
      ########################

      if (is.null(clinical_strat()[[1]]) && is.null(gene_strat()[[1]])) {

        # No stratification: keep all samples
        strat_names <- rownames(gene())

      } else if (!is.null(clinical_strat()[[1]]) && is.null(gene_strat()[[1]])) {

        # Clinical-only stratification
        clinical_keys <- lapply(clinical_strat(), function(x) {
          split_elements <- unlist(strsplit(x, " - "))
          split_elements[seq_along(split_elements) %% 2 == 1]
        })
        clinical_values <- lapply(clinical_strat(), function(x) {
          split_elements <- unlist(strsplit(x, " - "))
          split_elements[seq_along(split_elements) %% 2 == 0]
        })

        dictionary <- vector("list", length(unique(clinical_keys[[1]])))
        names(dictionary) <- unique(clinical_keys[[1]])
        for (i in seq_along(clinical_keys[[1]])) {
          key   <- clinical_keys[[1]][i]
          value <- clinical_values[[1]][i]
          dictionary[[key]] <- c(dictionary[[key]], value)
        }

        strat_df <- as.data.frame(clinical()[, unique(clinical_keys[[1]])])
        colnames(strat_df) <- unique(clinical_keys[[1]])
        rownames(strat_df) <- rownames(clinical())

        for (i in seq_len(ncol(strat_df))) {
          if (colnames(strat_df)[i] %in% subset(meta_clinical(), type == "continuous")$variable) {
            if (type_strat() == 1) {
              q25 <- quantile(strat_df[, i], 0.25, na.rm = TRUE)
              q75 <- quantile(strat_df[, i], 0.75, na.rm = TRUE)
              strat_df[, i] <- ifelse(
                strat_df[, i] > q75, "high",
                ifelse(strat_df[, i] >= q25, "baseline", "low")
              )
            } else {
              q50 <- quantile(strat_df[, i], 0.5, na.rm = TRUE)
              strat_df[, i] <- ifelse(strat_df[, i] > q50, "high", "low")
            }
          }
        }

        for (i in seq_along(clinical_keys[[1]])) {
          strat_df <- subset(
            strat_df,
            subset = strat_df[, clinical_keys[[1]][i]] %in% dictionary[[clinical_keys[[1]][i]]]
          )
        }
        strat_names <- rownames(strat_df)

      } else if (is.null(clinical_strat()[[1]]) && !is.null(gene_strat()[[1]])) {

        # Gene-only stratification
        gene_keys <- lapply(gene_strat(), function(x) {
          split_elements <- unlist(strsplit(x, " - "))
          split_elements[seq_along(split_elements) %% 2 == 1]
        })
        gene_values <- lapply(gene_strat(), function(x) {
          split_elements <- unlist(strsplit(x, " - "))
          split_elements[seq_along(split_elements) %% 2 == 0]
        })

        dictionary <- vector("list", length(unique(gene_keys[[1]])))
        names(dictionary) <- unique(gene_keys[[1]])
        for (i in seq_along(gene_keys[[1]])) {
          key   <- gene_keys[[1]][i]
          value <- gene_values[[1]][i]
          dictionary[[key]] <- c(dictionary[[key]], value)
        }

        strat_df <- as.data.frame(gene()[, unique(gene_keys[[1]])])
        rownames(strat_df) <- rownames(gene())
        colnames(strat_df) <- unique(gene_keys[[1]])

        for (i in seq_len(ncol(strat_df))) {
          if (type_strat() == 1) {
            q25 <- quantile(strat_df[, i], 0.25, na.rm = TRUE)
            q75 <- quantile(strat_df[, i], 0.75, na.rm = TRUE)
            strat_df[, i] <- ifelse(
              strat_df[, i] > q75, "high",
              ifelse(strat_df[, i] >= q25, "baseline", "low")
            )
          } else {
            q50 <- quantile(strat_df[, i], 0.5, na.rm = TRUE)
            strat_df[, i] <- ifelse(strat_df[, i] > q50, "high", "low")
          }
        }

        for (i in seq_along(gene_keys[[1]])) {
          strat_df <- subset(
            strat_df,
            subset = strat_df[, gene_keys[[1]][i]] %in% dictionary[[gene_keys[[1]][i]]]
          )
        }
        strat_names <- rownames(strat_df)

      } else {

        # Clinical + gene stratification
        clinical_keys <- lapply(clinical_strat(), function(x) {
          split_elements <- unlist(strsplit(x, " - "))
          split_elements[seq_along(split_elements) %% 2 == 1]
        })
        clinical_values <- lapply(clinical_strat(), function(x) {
          split_elements <- unlist(strsplit(x, " - "))
          split_elements[seq_along(split_elements) %% 2 == 0]
        })
        gene_keys <- lapply(gene_strat(), function(x) {
          split_elements <- unlist(strsplit(x, " - "))
          split_elements[seq_along(split_elements) %% 2 == 1]
        })
        gene_values <- lapply(gene_strat(), function(x) {
          split_elements <- unlist(strsplit(x, " - "))
          split_elements[seq_along(split_elements) %% 2 == 0]
        })

        dictionary_clinical <- vector("list", length(unique(clinical_keys[[1]])))
        names(dictionary_clinical) <- unique(clinical_keys[[1]])
        for (i in seq_along(clinical_keys[[1]])) {
          key   <- clinical_keys[[1]][i]
          value <- clinical_values[[1]][i]
          dictionary_clinical[[key]] <- c(dictionary_clinical[[key]], value)
        }

        dictionary_gene <- vector("list", length(unique(gene_keys[[1]])))
        names(dictionary_gene) <- unique(gene_keys[[1]])
        for (i in seq_along(gene_keys[[1]])) {
          key   <- gene_keys[[1]][i]
          value <- gene_values[[1]][i]
          dictionary_gene[[key]] <- c(dictionary_gene[[key]], value)
        }

        strat_df_clinical <- as.data.frame(clinical()[, unique(clinical_keys[[1]])])
        rownames(strat_df_clinical) <- rownames(clinical())
        colnames(strat_df_clinical) <- unique(clinical_keys[[1]])

        for (i in seq_len(ncol(strat_df_clinical))) {
          if (colnames(strat_df_clinical)[i] %in% subset(meta_clinical(), type == "continuous")$variable) {
            if (type_strat() == 1) {
              q25 <- quantile(strat_df_clinical[, i], 0.25, na.rm = TRUE)
              q75 <- quantile(strat_df_clinical[, i], 0.75, na.rm = TRUE)
              strat_df_clinical[, i] <- ifelse(
                strat_df_clinical[, i] > q75, "high",
                ifelse(strat_df_clinical[, i] >= q25, "baseline", "low")
              )
            } else {
              q50 <- quantile(strat_df_clinical[, i], 0.5, na.rm = TRUE)
              strat_df_clinical[, i] <- ifelse(strat_df_clinical[, i] > q50, "high", "low")
            }
          }
        }

        strat_df_gene <- as.data.frame(gene()[, unique(gene_keys[[1]])])
        rownames(strat_df_gene) <- rownames(gene())
        colnames(strat_df_gene) <- unique(gene_keys[[1]])

        for (i in seq_len(ncol(strat_df_gene))) {
          if (type_strat() == 1) {
            q25 <- quantile(strat_df_gene[, i], 0.25, na.rm = TRUE)
            q75 <- quantile(strat_df_gene[, i], 0.75, na.rm = TRUE)
            strat_df_gene[, i] <- ifelse(
              strat_df_gene[, i] > q75, "high",
              ifelse(strat_df_gene[, i] >= q25, "baseline", "low")
            )
          } else {
            q50 <- quantile(strat_df_gene[, i], 0.5, na.rm = TRUE)
            strat_df_gene[, i] <- ifelse(strat_df_gene[, i] > q50, "high", "low")
          }
        }

        strat_df <- cbind(strat_df_clinical, strat_df_gene)
        dictionary <- merge(dictionary_clinical, dictionary_gene, all = TRUE)

        for (i in seq_along(clinical_keys[[1]])) {
          strat_df <- subset(
            strat_df,
            subset = strat_df[, clinical_keys[[1]][i]] %in% dictionary[[clinical_keys[[1]][i]]]
          )
        }
        for (i in seq_along(gene_keys[[1]])) {
          strat_df <- subset(
            strat_df,
            subset = strat_df[, gene_keys[[1]][i]] %in% dictionary[[gene_keys[[1]][i]]]
          )
        }
        strat_names <- rownames(strat_df)
      }

      if (length(strat_names) < 2) {
        errors <- c(errors, "Stratified Cohort Must Contain at Least Two Samples")
        metaFinal_out(NULL)
        errors_list(errors)
        output$reportText <- renderText(paste(
          "Cohort:", cohort_name(), "\n",
          "Errors:", paste(errors, collapse = "; "), "\n"
        ))
        return()
      }

      ##############################
      # 2. OUTER DF AND OUTER FLAG #
      ##############################

      outer_df <- cbind(clinical(), gene())[strat_names, , drop = FALSE]

      upper_keys <- lapply(outerVarFor(), function(x) {
        split_elements <- unlist(strsplit(x, " - "))
        split_elements[seq_along(split_elements) %% 2 == 1]
      })
      upper_values <- lapply(outerVarFor(), function(x) {
        split_elements <- unlist(strsplit(x, " - "))
        split_elements[seq_along(split_elements) %% 2 == 0]
      })
      lower_keys <- lapply(outerVarRev(), function(x) {
        split_elements <- unlist(strsplit(x, " - "))
        split_elements[seq_along(split_elements) %% 2 == 1]
      })
      lower_values <- lapply(outerVarRev(), function(x) {
        split_elements <- unlist(strsplit(x, " - "))
        split_elements[seq_along(split_elements) %% 2 == 0]
      })

      total_keys <- unique(c(upper_keys[[1]], lower_keys[[1]]))
      outer_df   <- as.data.frame(outer_df[, total_keys, drop = FALSE])
      rownames(outer_df) <- strat_names
      colnames(outer_df) <- total_keys

      dictionary_upper <- vector("list", length(unique(upper_keys[[1]])))
      names(dictionary_upper) <- unique(upper_keys[[1]])
      for (i in seq_along(upper_keys[[1]])) {
        key   <- upper_keys[[1]][i]
        value <- upper_values[[1]][i]
        dictionary_upper[[key]] <- c(dictionary_upper[[key]], value)
      }

      dictionary_lower <- vector("list", length(unique(lower_keys[[1]])))
      names(dictionary_lower) <- unique(lower_keys[[1]])
      for (i in seq_along(lower_keys[[1]])) {
        key   <- lower_keys[[1]][i]
        value <- lower_values[[1]][i]
        dictionary_lower[[key]] <- c(dictionary_lower[[key]], value)
      }

      # Define continuous clinical variables once
      continuous_clinical <- if ("continuous" %in% meta_clinical()[, "type"]) {
        subset(meta_clinical(), type == "continuous")$variable
      } else {
        character(0)
      }

      # Discretize outer_df variables and map to high/baseline/low or high/low
      for (i in seq_len(ncol(outer_df))) {
        var_name <- colnames(outer_df)[i]

        if (var_name %in% continuous_clinical || var_name %in% colnames(gene())) {

          x_num <- suppressWarnings(as.numeric(outer_df[, i]))
          if (all(is.na(x_num))) {
            next
          }

          if (type_outer() == 1) {
            q25 <- quantile(x_num, 0.25, na.rm = TRUE)
            q75 <- quantile(x_num, 0.75, na.rm = TRUE)
            outer_df[, i] <- ifelse(
              x_num > q75, "high",
              ifelse(x_num >= q25, "baseline", "low")
            )
          } else {
            q50 <- quantile(x_num, 0.5, na.rm = TRUE)
            outer_df[, i] <- ifelse(x_num > q50, "high", "low")
          }
        }
      }

      # Map to upper/lower using dictionaries
      for (i in seq_along(total_keys)) {
        key <- total_keys[i]
        if (!is.null(dictionary_upper[[key]])) {
          outer_df[, i] <- ifelse(outer_df[, i] %in% dictionary_upper[[key]], "upper", outer_df[, i])
        }
        if (!is.null(dictionary_lower[[key]])) {
          outer_df[, i] <- ifelse(outer_df[, i] %in% dictionary_lower[[key]], "lower", outer_df[, i])
        }
      }

      outer_df <- na.omit(outer_df)
      if (nrow(outer_df) == 0) {
        errors <- c(errors, "No samples remain after outer discretization")
        metaFinal_out(NULL)
        errors_list(errors)
        output$reportText <- renderText(paste(
          "Cohort:", cohort_name(), "\n",
          "Errors:", paste(errors, collapse = "; "), "\n"
        ))
        return()
      }

      # Determine rows that are consistently upper or lower
      v1 <- do.call(pmin, c(outer_df, na.rm = TRUE))
      v2 <- do.call(pmax, c(outer_df, na.rm = TRUE))
      outer_label <- ifelse(v1 == v2, v1, NA_character_)

      valid_rows <- which(outer_label %in% c("upper", "lower"))
      if (length(valid_rows) == 0) {
        errors <- c(errors, "No samples remain after applying outer comparison filters")
        metaFinal_out(NULL)
        errors_list(errors)
        output$reportText <- renderText(paste(
          "Cohort:", cohort_name(), "\n",
          "Errors:", paste(errors, collapse = "; "), "\n"
        ))
        return()
      }

      outer_df    <- outer_df[valid_rows, , drop = FALSE]
      outer_label <- outer_label[valid_rows]

      #############################
      # 3. INNER DF FROM SCRATCH  #
      #############################

      inner_df <- data.frame(
        outer = outer_label,
        row.names = rownames(outer_df),
        stringsAsFactors = FALSE
      )

      inner_choices <- cbind(clinical(), gene())[rownames(inner_df), , drop = FALSE]

      if (inner() == "" || type_inner() == 0) {
        inner_df$inner <- "All"
      } else {
        if (!inner() %in% colnames(inner_choices)) {
          errors <- c(errors, paste0("Inner variable '", inner(), "' not found in current cohort"))
          metaFinal_out(NULL)
          errors_list(errors)
          output$reportText <- renderText(paste(
            "Cohort:", cohort_name(), "\n",
            "Errors:", paste(errors, collapse = "; "), "\n"
          ))
          return()
        }
        inner_df$inner <- inner_choices[, inner(), drop = TRUE]

        if (inner() %in% continuous_clinical || inner() %in% colnames(gene())) {
          if (type_inner() == 1) {
            q25 <- quantile(as.numeric(inner_df$inner), 0.25, na.rm = TRUE)
            q75 <- quantile(as.numeric(inner_df$inner), 0.75, na.rm = TRUE)
            inner_df$inner <- ifelse(
              as.numeric(inner_df$inner) > q75, "high",
              ifelse(as.numeric(inner_df$inner) >= q25, "baseline", "low")
            )
          } else if (type_inner() == 2) {
            q50 <- quantile(as.numeric(inner_df$inner), 0.5, na.rm = TRUE)
            inner_df$inner <- ifelse(as.numeric(inner_df$inner) > q50, "high", "low")
          }
          # type_inner() 3/4: leave categorical/binary as-is
        }
      }

      inner_df <- na.omit(inner_df)
      if (nrow(inner_df) == 0) {
        errors <- c(errors, "No samples remain after applying inner comparison filters")
        metaFinal_out(NULL)
        errors_list(errors)
        output$reportText <- renderText(paste(
          "Cohort:", cohort_name(), "\n",
          "Errors:", paste(errors, collapse = "; "), "\n"
        ))
        return()
      }

      ##################################
      # 4. STORE META AND UPDATE TEXT  #
      ##################################

      metaFinal_out(inner_df)
      errors_list(errors)

      clinical_strat_text <- paste(clinical_strat(), sep = ", ")
      gene_strat_text     <- paste(gene_strat(),     sep = ", ")
      outer_upper_text    <- paste(outerVarFor(),    sep = ", ")
      outer_lower_text    <- paste(outerVarRev(),    sep = ", ")

      output$reportText <- renderText({
        paste(
          paste0("Cohort: ", cohort_name(), "\n"),
          paste0("Patient Stratifications (Clinical): ", clinical_strat_text, "\n"),
          paste0("Patient Stratifications (Gene): ", gene_strat_text, "\n"),
          paste0("Outer Comparisons (Upper): ", outer_upper_text, "\n"),
          paste0("Outer Comparisons (Lower): ", outer_lower_text, "\n"),
          paste0("Inner Comparisons: ", inner(), "\n"),
          paste0("Minimum AUC: ", input$customAUC, "\n"),
          paste0("Minimum Log2FC: ", input$customFC, "\n"),
          paste0("Maximum FDR: ", input$customFDR, "\n"),
          paste0("Errors: ", paste(errors, collapse = "; "), "\n"),
          paste0("Report Generated at: ", Sys.time())
        )
      })
    }
  )



    observeEvent(input$submitReport, {
      session$sendCustomMessage("toggle-loading", TRUE)
      session$sendCustomMessage("loading-text", "Reading cohort registry…")
      on.exit({
        session$sendCustomMessage("toggle-loading", FALSE)
        session$sendCustomMessage("loading-text", "Cohort loaded successfully")
      }, add = TRUE)
      session$sendCustomMessage("loading-text", "Running xCheck analysis…")
      run_analysis(
        AUC_cutoff = input$customAUC,
        FC_cutoff  = input$customFC,
        FDR_cutoff = input$customFDR
      )
      analysis_ready(TRUE)
    })

    observeEvent(input$submitReportDefault, {
      session$sendCustomMessage("toggle-loading", TRUE)
      session$sendCustomMessage("loading-text", "Reading cohort registry…")
      on.exit({
        session$sendCustomMessage("toggle-loading", FALSE)
        session$sendCustomMessage("loading-text", "Cohort loaded successfully")
      }, add = TRUE)
      session$sendCustomMessage("loading-text", "Running xCheck analysis…")
      run_analysis(
        AUC_cutoff = 0,
        FC_cutoff  = 0,
        FDR_cutoff = 0.05
      )
      analysis_ready(TRUE)
    })

    observeEvent(input$submitReport_exp, {
      req(input$main_tabs == "exp")
      if (is.null(input$outerType_exp)) {
        showModal(
          modalDialog(
            title = "Missing outer comparison",
            "An outer comparison parameter must be selected.",
            easyClose = TRUE,
            footer = modalButton("OK")
          )
        )
        return()  # do nothing else
      }
      if (is.null(input$outerVarTop_exp) && is.null(input$outerVarBottom_exp)) {
        showModal(
          modalDialog(
            title = "Missing outer comparison",
            "At least one outer comparison needs to be selected for upper/lower.",
            easyClose = TRUE,
            footer = modalButton("OK")
          )
        )
        return()  # do nothing else
      }
      session$sendCustomMessage("toggle-loading", TRUE)
      session$sendCustomMessage("loading-text", "Reading experiment registry…")
      on.exit({
        session$sendCustomMessage("toggle-loading", FALSE)
        session$sendCustomMessage("loading-text", "Experiment loaded successfully")
      }, add = TRUE)
      session$sendCustomMessage("loading-text", "Running xCheck analysis…")
      run_analysis(
        AUC_cutoff = input$customAUC_exp,
        FC_cutoff  = input$customFC_exp,
        FDR_cutoff = input$customFDR_exp
      )
      analysis_ready_exp(TRUE)
    })
    
    observeEvent(input$submitReportDefault_exp, {
      req(input$main_tabs == "exp")
      if (is.null(input$outerType_exp)) {
        showModal(
          modalDialog(
            title = "Missing outer comparison",
            "An outer comparison parameter must be selected.",
            easyClose = TRUE,
            footer = modalButton("OK")
          )
        )
        return()  # do nothing else
      }
      if (is.null(input$outerVarTop_exp) && is.null(input$outerVarBottom_exp)) {
        showModal(
          modalDialog(
            title = "Missing outer comparison",
            "At least one outer comparison needs to be selected for upper/lower.",
            easyClose = TRUE,
            footer = modalButton("OK")
          )
        )
        return()  # do nothing else
      }
      session$sendCustomMessage("toggle-loading", TRUE)
      session$sendCustomMessage("loading-text", "Reading cohort registry…")
      on.exit({
        session$sendCustomMessage("toggle-loading", FALSE)
        session$sendCustomMessage("loading-text", "Cohort loaded successfully")
      }, add = TRUE)
      session$sendCustomMessage("loading-text", "Running xCheck analysis…")
      run_analysis(
        AUC_cutoff = 0,
        FC_cutoff  = 0,
        FDR_cutoff = 0.05
      )
      analysis_ready_exp(TRUE)
    })
    
    observeEvent(list(analysis_ready(),analysis_ready_exp()), {
      req(analysis_ready() || analysis_ready_exp())
      up_lab   <- upper_label_rv()
      down_lab <- lower_label_rv()

      choices_named <- c(
        setNames(up_lab,   paste0("Significant - ", up_lab)),
        setNames("ns",     "Not significant"),
        setNames(down_lab, paste0("Significant - ", down_lab))
      )

      updateCheckboxGroupInput(
        session, "outerSelect",
        choices  = as.list(choices_named),
        selected = unname(choices_named)
      )
    })


    metabs_fin <- reactive({
      req(analysis_ready() || analysis_ready_exp())
      req(exists("metabs_tab"), !is.null(metabs_tab()))

      if (!"inner_comparison" %in% colnames(metabs_tab())) {
        return(metabs_tab())
      }

      if (nrow(metabs_tab()) == 0) {
        return(metabs_tab())
      }

      inner_levels <- unique(metabs_tab()[,"inner_comparison", drop = TRUE])

      if (length(inner_levels) == 1) {
        subset(metabs_tab(), significant_category %in% input$outerSelect)
      } else {
        subset(
          metabs_tab(),
          inner_comparison %in% input$innerSelect & significant_category %in% input$outerSelect
        )
      }
    })

    plot1_combined <- reactive({
      req(isTRUE(analysis_ready()) || isTRUE(analysis_ready_exp()))
      
      mf         <- metabs_fin()       # experimental metabs
      up_lab     <- upper_label_rv()
      down_lab   <- lower_label_rv()
      
      
      inner_levels <- unique(mf$inner_comparison)
      
      if (length(inner_levels) <= 1) {
        ggplot(
          mf,
          aes(
            x   = reorder(subsystem, -sign_t_log_pval),
            y   = sign_t_log_pval,
            key = metabolite
          )
        ) +
          geom_point(
            aes(color = significant_category, fill = significant_category),
            size = 1.4
          ) +
          geom_hline(yintercept =  1.3, size = .05) +
          geom_hline(yintercept = -1.3, size = .05) +
          xlab("subsystem") +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
          scale_color_manual(
            values = c(
              setNames("#1f77b4", up_lab),
              setNames("#d62728", down_lab),
              "ns" = "grey60"
            ),
            drop = FALSE
          ) +
          scale_fill_manual(
            values = c(
              setNames("#1f77b4", up_lab),
              setNames("#d62728", down_lab),
              "ns" = "grey80"
            ),
            drop = FALSE
          )
        
      } else {
        mf$subsystem_base <- factor(mf$subsystem)
        n_inner           <- length(inner_levels)
        
        base_x  <- as.numeric(mf$subsystem_base)
        offsets <- seq(-0.3, 0.3, length.out = n_inner)
        names(offsets) <- inner_levels
        
        mf$x_pos <- base_x + offsets[as.character(mf$inner_comparison)]
        mf <- mf[order(mf$subsystem, mf$metabolite, mf$x_pos), ]
        
        ggplot(
          mf,
          aes(
            x   = x_pos,
            y   = sign_t_log_pval,
            key = metabolite
          )
        ) +
          geom_hline(yintercept =  1.3, size = .05) +
          geom_hline(yintercept = -1.3, size = .05) +
          geom_line(
            aes(group = interaction(subsystem, metabolite)),
            linewidth = 0.3,
            color = "grey70",
            alpha = 0.6
          ) +
          geom_point(
            aes(color = inner_comparison),
            size = 1.8
          ) +
          scale_x_continuous(
            breaks = sort(unique(base_x)),
            labels = levels(mf$subsystem_base)
          ) +
          xlab("subsystem") +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
          scale_color_brewer(palette = "Set1", drop = FALSE)
      }
    })

    selected_points_plot1_obs <- reactive({
      req(analysis_ready())
      
      # if plot not drawn yet, or not on the right tab, abort
      #if (is.null(output$plot1)) return(NULL)
      
      sel <- event_data("plotly_selected", source = "p1_obs")
      if (is.null(sel) || nrow(sel) == 0) return(NULL)
      
      mf <- metabs_fin()
      as.data.frame(mf[mf$metabolite %in% sel$key, , drop = FALSE])
    })

    selected_points_plot1_exp <- reactive({
      req(analysis_ready_exp())
      
      # if plot not drawn yet, or not on the right tab, abort
      #if (is.null(output$plot1)) return(NULL)
      
      sel <- event_data("plotly_selected", source = "p1_exp")
      if (is.null(sel) || nrow(sel) == 0) return(NULL)
      
      mf <- metabs_fin()
      as.data.frame(mf[mf$metabolite %in% sel$key, , drop = FALSE])
    })
    
    output$selected_points_n_obs <- reactive({
      
      sp <- selected_points_plot1_obs()

      # handle NULL or unexpected types safely
      if (is.null(sp)) return(0)
      if (is.data.frame(sp) || is.matrix(sp)) {
        return(nrow(sp))
      }
      if (is.list(sp) && !is.data.frame(sp)) {
        # some plotly selections come as a list; adjust as needed
        return(length(sp$x))  # or another appropriate element
      }

      0
    })
    
    output$selected_points_n_exp <- reactive({
      sp <- selected_points_plot1_exp()
      
      # handle NULL or unexpected types safely
      if (is.null(sp)) return(0)
      if (is.data.frame(sp) || is.matrix(sp)) {
        return(nrow(sp))
      }
      if (is.list(sp) && !is.data.frame(sp)) {
        # some plotly selections come as a list; adjust as needed
        return(length(sp$x))  # or another appropriate element
      }
      
      0
    })
    
    observe({
      outputOptions(output, "selected_points_n_obs", suspendWhenHidden = FALSE)
      outputOptions(output, "selected_points_n_exp", suspendWhenHidden = FALSE)
    })

    output$plot1_obs <- renderPlotly({
      req(input$main_tabs == "obs")
      req(analysis_ready())
      p <- plot1_combined()
      ggplotly(p, source = "p1_obs", tooltip = c("x", "y", "metabolite")) |>
        event_register("plotly_selected")
    })
    
    output$plot1_exp <- renderPlotly({
      req(input$main_tabs == "exp")
      req(analysis_ready_exp())
      p <- plot1_combined()
      ggplotly(p, source = "p1_exp", tooltip = c("x", "y", "metabolite")) |>
        event_register("plotly_selected")
    })

    plot1pdf <- reactive({
      ggplot(metabs_fin(), aes(x=as.numeric(t_welch_statistic), y=reorder(subsystem, -as.numeric(t_welch_statistic)), color=significant_category, fill = significant_category)) + geom_point(size = 1, position = position_dodge(width=0.5))+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+geom_vline(xintercept=1.3, size=.05)+geom_vline(xintercept=-1.3, size=.05)+ylab("subsystem")+xlab("log(pval) with sign of t-statistic")+ggtitle(paste0("Significant Reactions Across All Subsystems - (",upper_label_rv()," v. ",lower_label_rv(),")"))
    })

    output$reacts_obs <- downloadHandler(
      filename = "(enter name for fingerprint preparation object).rds",
      content = function(file) {
        saveRDS(list(metabs_fin(),metaFinal_out(),data.frame(AUC=input$customAUC,Log2FC=input$customFC,FDR=input$customFDR),"cohort",list(upper_label_rv(),lower_label_rv()),list(unique(metabs_fin()$inner_comparison))),file)
      }
    )

    output$reacts_exp <- downloadHandler(
      filename = "(enter name for fingerprint preparation object).rds",
      content = function(file) {
        saveRDS(list(metabs_fin(),metaFinal_out(),data.frame(AUC=input$customAUC,Log2FC=input$customFC,FDR=input$customFDR),"experiment",list(top_exp(),bottom_exp()),list(unique(metabs_fin()$inner_comparison))),file)
      }
    )
    
    systems<-reactive({
      if(length(unique(metabs_fin()[,"inner_comparison"])) == 1){
        ggbld <- ggplot_build(ggplot(metabs_fin(), aes(x=reorder(subsystem, -sign_t_log_pval), y=sign_t_log_pval, color=sig, fill = sig)) + geom_point(size = 1, position = position_dodge(width=0.5))+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+geom_hline(yintercept=1.3, size=.05)+geom_hline(yintercept=-1.3, size=.05)+xlab("subsystem"))
        ggbld$layout$panel_params[[1]]$x$limits
      }else{
        ggbld <- ggplot_build(ggplot(metabs_fin(), aes(x=reorder(subsystem, -sign_t_log_pval), y=sign_t_log_pval, color=inner_comparison, fill = inner_comparison)) + geom_point(size = 1, position = position_dodge(width=0.5))+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+geom_hline(yintercept=1.3, size=.05)+geom_hline(yintercept=-1.3, size=.05)+xlab("subsystem"))
        ggbld$layout$panel_params[[1]]$x$limits
      }

    })


    output$dlPlot1 <- downloadHandler(
      filename = function() {
        paste('all_subsystems.pdf', sep='')
      },
      content = function(file) {
        withProgress(message = 'Generating PDF', value = 0, {
          incProgress(0.3, detail = "Preparing data")

          pdf(file, width = 11, height = 8.5)  # Standard letter size

          # Plot or write your content to the PDF here
          # For example:
          plot(plot1pdf())
          # Add more plots or text as needed

          incProgress(0.3, detail = "Creating PDF")

          # Close the PDF device
          dev.off()

          incProgress(0.4, detail = "Finalizing")
        })
      },
      contentType = "application/pdf"
    )

    ## 1) Subsystem index from clicking plot1 (Plotly)
    selected_subsystem_index <- reactive({
      req(input$main_tabs == "obs")
      click <- event_data("plotly_click", source = "p1_obs")
      req(click)
      # if subsystem is on x as discrete factor, pointNumber+1 indexes systems()
      idx <- click$pointNumber[1] + 1
      idx
    })
    
    selected_subsystem_index <- reactive({
      req(input$main_tabs == "exp")
      click <- event_data("plotly_click", source = "p1_exp")
      req(click)
      # if subsystem is on x as discrete factor, pointNumber+1 indexes systems()
      idx <- click$pointNumber[1] + 1
      idx
    })

      # Table of brushed reactions
    selected_flux <- reactiveVal(NULL)

    output$info_obs <- renderReactable({
      df <- selected_points_plot1_obs()
      df$var_1_mean      <- signif(as.numeric(df$var_1_mean), 4)
      df$var_2_mean      <- signif(as.numeric(df$var_2_mean), 4)
      df$sign_t_log_pval <- signif(as.numeric(df$sign_t_log_pval), 4)
      df$pval <- signif(as.numeric(df$pval), 4)
      df$t_welch_statistic <- signif(as.numeric(df$t_welch_statistic), 4)

      reactable::reactable(
        df,
        striped    = TRUE,
        highlight  = TRUE,
        compact    = TRUE,
        selection  = "single",
        onClick    = "select"   # clicking a row selects it
      )
    })
    
    output$info_exp <- renderReactable({
      df <- selected_points_plot1_exp()
      df$var_1_mean      <- signif(as.numeric(df$var_1_mean), 4)
      df$var_2_mean      <- signif(as.numeric(df$var_2_mean), 4)
      df$sign_t_log_pval <- signif(as.numeric(df$sign_t_log_pval), 4)
      df$pval <- signif(as.numeric(df$pval), 4)
      df$t_welch_statistic <- signif(as.numeric(df$t_welch_statistic), 4)
      
      reactable::reactable(
        df,
        striped    = TRUE,
        highlight  = TRUE,
        compact    = TRUE,
        selection  = "single",
        onClick    = "select"   # clicking a row selects it
      )
    })

    observeEvent(reactable::getReactableState("info_obs", "selected"), {
      sel <- reactable::getReactableState("info_obs", "selected")
      df  <- selected_points_plot1_obs()
      if (!is.null(sel) && length(sel) == 1) {
        flux_name <- df$metabolite[sel]
        selected_flux(flux_name)
        updateSelectizeInput(session, "boxplot", selected = flux_name)
      }
    })
    
    observeEvent(reactable::getReactableState("info_exp", "selected"), {
      sel <- reactable::getReactableState("info_exp", "selected")
      df  <- selected_points_plot1_exp()
      if (!is.null(sel) && length(sel) == 1) {
        flux_name <- df$metabolite[sel]
        selected_flux(flux_name)
        updateSelectizeInput(session, "boxplot", selected = flux_name)
      }
    })

      # Download brushed data as CSV
      data_obs <- reactive({
        as.data.frame(selected_points_plot1_obs())
      })
      
      data_exp <- reactive({
        as.data.frame(selected_points_plot1_exp())
      })


      output$data_obs <- downloadHandler(
        filename = function() {
          paste("_selected_reactions_", outer, ".csv", sep = "")
        },
        content = function(file) {
          write.csv(data_obs(), file, row.names = FALSE)
        }
      )
      
      output$data_exp <- downloadHandler(
        filename = function() {
          paste("_selected_reactions_", outer, ".csv", sep = "")
        },
        content = function(file) {
          write.csv(data_exp(), file, row.names = FALSE)
        }
      )

      # Update boxplot choices from brushed points
      observe({
        df <- selected_points_plot1_obs()
        updateSelectizeInput(
          session, "boxplot",
          choices = unique(df$metabolite)
        )
      })
      
      observe({
        df <- selected_points_plot1_exp()
        updateSelectizeInput(
          session, "boxplot",
          choices = unique(df$metabolite)
        )
      })


    plot2_obs <- reactive({
      if (input$boxplot == "") {
        return(ggplot())
      }

      metaFinal <- metaFinal_out()
      req(metaFinal)

      session$sendCustomMessage("toggle-loading", TRUE)
      session$sendCustomMessage("loading-text", "Reading cohort registry…")
      on.exit({
        session$sendCustomMessage("toggle-loading", FALSE)
        session$sendCustomMessage("loading-text", "Cohort loaded successfully")
      }, add = TRUE)
      session$sendCustomMessage("loading-text", "Loading plots…")

      outer_var <- colnames(metaFinal)[1]
      inner_var <- colnames(metaFinal)[2]

      outer_type <- as.numeric(input$outerType)
      has_inner  <- !(is.null(inner()) || inner() == "" || type_inner() == 0)
      inner_name <- if (has_inner) inner() else "None"

      # Base data
      flux_mat <- flux_mat <- if (input$main_tabs == "exp") experiment_flux() else cohort_flux()                # your flux matrix reactiveVal
      req(!is.null(flux_mat))

      data_filt <- cbind(
        metaFinal,
        flux_mat[rownames(metaFinal), , drop = FALSE]
      )

      df <- data_filt[, c(colnames(metaFinal), input$boxplot), drop = FALSE]
      colnames(df)[1] <- "outer_group"
      colnames(df)[2] <- "inner"
      # ---------- Parse outer selections for labels ----------

      outer_for_vec <- if (!is.null(outerVarFor())) outerVarFor()[[1]] else character(0)
      outer_rev_vec <- if (!is.null(outerVarRev())) outerVarRev()[[1]] else character(0)
      outer_all     <- c(outer_rev_vec,outer_for_vec)

      parse_outer <- function(x) {
        parts <- strsplit(x, " - ", fixed = TRUE)[[1]]
        if (length(parts) >= 2) list(var = parts[1], level = parts[2]) else list(var = x, level = x)
      }

      # Default labels
      cont_name       <- NULL
      outer_level_map <- NULL
      outer_levels    <- sort(unique(df$outer_group))

      if (length(outer_all) > 0) {
        parsed <- lapply(outer_all, parse_outer)
        cont_name <- parsed[[1]]$var                    # e.g. "FLT3.ITD"
        for_levels <- vapply(outer_for_vec, function(z) parse_outer(z)$level, character(1))
        rev_levels <- vapply(outer_rev_vec, function(z) parse_outer(z)$level, character(1))
        outer_levels <- unique(c(for_levels, rev_levels))

        # Map original outer_group values to pretty labels:
        # assume outer_group encodes upper/lower in your metaFinal,
        # but here we simply relabel based on order of levels we derived
        if (length(outer_levels) >= 2) {
          outer_level_map <- setNames(outer_levels, c("upper", "lower")[seq_along(outer_levels)])
        }
      }

      # Apply outer label mapping to df$outer_group (if available)
      df$outer_group <- factor(as.character(df$outer_group))
      if (!is.null(outer_level_map)) {
        # Map factor levels "upper"/"lower" (or whatever you used) to user-friendly labels
        df$outer_group <- factor(
          outer_level_map[as.character(df$outer_group)],
          levels = rev(outer_levels)
        )
      }

      # ---------- BOX PLOTS (original behavior, now with pretty x labels) ----------
      flux_col <- input$boxplot

      pval_to_stars <- function(p) {
        if (is.na(p)) return("ns")
        if (p < 0.001) return("***")
        if (p < 0.01)  return("**")
        if (p < 0.05)  return("*")
        "ns"
      }

      if (!has_inner) {
        # global test between the two outer groups
        sub <- df[!is.na(df$outer_group) & !is.na(df[[flux_col]]), ]
        if (nrow(sub) >= 2 && nlevels(sub$outer_group) == 2) {
          p_val <- tryCatch(
            {
              t.test(df[[flux_col]] ~ df$outer_group)$p.value
            },
            error = function(e) NA_real_
          )
        } else {
          p_val <- NA_real_
        }
        star_label <- pval_to_stars(p_val)

        # y position for the star: slightly above max
        y_max <- max(sub[[flux_col]], na.rm = TRUE)
        y_star <- y_max + 0.05 * diff(range(sub[[flux_col]], na.rm = TRUE))

        p_box <- ggplot(sub, aes(x = outer_group, y = .data[[flux_col]])) +
          geom_boxplot(
            width  = 0.5,
            outlier.shape = NA,
            fill   = "#e8edff",
            color  = "#5c6bc0",
            alpha  = 0.7
          ) +
          geom_jitter(
            width  = 0.12,
            alpha  = 0.5,
            size   = 1.4,
            color  = "#3949ab"
          ) +
          # add one star spanning both boxes
          geom_segment(aes(x = 1, xend = 2, y = y_star, yend = y_star),
                       size = 0.3, color = "#757575") +
          geom_text(aes(x = 1.5, y = y_star, label = star_label),
                    vjust = -0.4, size = 3.2, color = "#424242") +
          radar_theme +
          xlab(if (is.null(cont_name)) "Outer group" else cont_name) +
          ylab(flux_col)

      } else {
        # per-inner tests
        inner_levels <- sort(unique(df$inner))
        # compute a small data frame with positions and stars
        star_df <- do.call(rbind, lapply(inner_levels, function(lv) {
          sub <- df[df$inner == lv & !is.na(df$outer_group) & !is.na(df[[flux_col]]), ]
          if (nrow(sub) < 2 || nlevels(sub$outer_group) < 2) {
            return(NULL)
          }
          p_val <- tryCatch(
            {
              t.test(sub[[flux_col]] ~ sub$outer_group)$p.value
            },
            error = function(e) NA_real_
          )
          star_label <- pval_to_stars(p_val)
          y_max <- max(sub[[flux_col]], na.rm = TRUE)
          y_star <- y_max + 0.05 * diff(range(sub[[flux_col]], na.rm = TRUE))
          data.frame(
            inner = lv,
            x_start = 1,
            x_end   = 2,
            x_mid   = 1.5,
            y_star  = y_star,
            label   = star_label,
            stringsAsFactors = FALSE
          )
        }))

        p_box <- ggplot(df, aes(x = outer_group, y = .data[[flux_col]])) +
          geom_boxplot(
            width  = 0.5,
            outlier.shape = NA,
            fill   = "#e8edff",
            color  = "#5c6bc0",
            alpha  = 0.7
          ) +
          geom_jitter(
            width  = 0.12,
            alpha  = 0.5,
            size   = 1.4,
            color  = "#3949ab"
          ) +
          facet_grid(. ~ inner) + radar_theme +
          xlab(if (is.null(cont_name)) "Outer group" else cont_name) +
          ylab(flux_col)

        if (!is.null(star_df) && nrow(star_df) > 0) {
          p_box <- p_box +
            geom_segment(
              data = star_df,
              aes(x = x_start, xend = x_end, y = y_star, yend = y_star),
              inherit.aes = FALSE
            ) +
            geom_text(
              data = star_df,
              aes(x = x_mid, y = y_star, label = label),
              vjust = -0.3,
              inherit.aes = FALSE
            )
        }
      }

      # ---------- IF OUTER NOT CONTINUOUS-DERIVED: ONLY BOX PLOTS ----------

      if (!(outer_type %in% c(1, 2))) {
        return(p_box)
      }

      # ---------- RECONSTRUCT RAW CONTINUOUS OUTER VARIABLE ----------

      if (is.null(cont_name)) {
        return(p_box)
      }

      all_covars <- cbind(clinical(), gene())
      if (!cont_name %in% colnames(all_covars)) {
        return(p_box)
      }

      vals  <- all_covars[rownames(metaFinal), cont_name]
      x_num <- suppressWarnings(as.numeric(vals))
      if (all(is.na(x_num))) {
        return(p_box)
      }

      df$outer_cont <- x_num[match(rownames(df), rownames(metaFinal))]

      # ---------- CORRELATION PLOTS (outer_cont vs flux, r annotated) ----------

      flux_col <- input$boxplot

      if (!has_inner) {
        sub <- df[!is.na(df$outer_cont) & !is.na(df[[flux_col]]), ]
        if (nrow(sub) == 0) {
          return(p_box)
        }
        r_val <- suppressWarnings(cor(sub$outer_cont, sub[[flux_col]], use = "complete.obs"))
        p_corr <- ggplot(sub, aes(x = outer_cont, y = .data[[flux_col]])) +
          geom_point(alpha = 0.55, size = 1.5, color = "#3949ab") +
          geom_smooth(method = "lm", se = FALSE, color = "#ef6c00", size = 0.7) +
          xlab(cont_name) +
          ylab(flux_col) + radar_theme +
          ggtitle(
            paste0("Correlation across full cohort (r = ",
                   sprintf("%.2f", r_val), ")")
          )
      } else {
        inner_levels <- sort(unique(df$inner))
        corr_plots <- lapply(inner_levels, function(lv) {
          sub <- df[df$inner == lv & !is.na(df$outer_cont) & !is.na(df[[flux_col]]), ]
          if (nrow(sub) == 0) {
            return(ggplot() + ggtitle(paste("No data for", lv)))
          }
          r_val <- suppressWarnings(cor(sub$outer_cont, sub[[flux_col]], use = "complete.obs"))
          ggplot(sub, aes(x = outer_cont, y = .data[[flux_col]])) +
            geom_point(alpha = 0.55, size = 1.5, color = "#3949ab") +
            geom_smooth(method = "lm", se = FALSE, color = "#ef6c00", size = 0.7) +
            xlab(cont_name) +
            ylab(flux_col) +
            ggtitle(
              paste0(inner_name, " = ", lv, " (r = ", sprintf("%.2f", r_val), ")")
            )
        })
        p_corr <- cowplot::plot_grid(plotlist = corr_plots, ncol = 1)
      }

      cowplot::plot_grid(p_box, p_corr, ncol = 2, rel_widths = c(2, 1))
    })
    
    plot2_exp <- reactive({
      if (input$boxplot == "") {
        return(ggplot())
      }
      
      metaFinal <- metaFinal_out()
      req(metaFinal)
      
      session$sendCustomMessage("toggle-loading", TRUE)
      session$sendCustomMessage("loading-text", "Reading cohort registry…")
      on.exit({
        session$sendCustomMessage("toggle-loading", FALSE)
        session$sendCustomMessage("loading-text", "Cohort loaded successfully")
      }, add = TRUE)
      session$sendCustomMessage("loading-text", "Loading plots…")
      
      outer_var <- colnames(metaFinal)[1]
      inner_var <- colnames(metaFinal)[2]
      
      outer_type <- as.numeric(input$outerType)
      has_inner  <- !(is.null(inner()) || inner() == "" || type_inner() == 0)
      inner_name <- if (has_inner) inner() else "None"
      
      # Base data
      flux_mat <- flux_mat <- if (input$main_tabs == "exp") experiment_flux() else cohort_flux()                # your flux matrix reactiveVal
      req(!is.null(flux_mat))
      
      data_filt <- cbind(
        metaFinal,
        flux_mat[rownames(metaFinal), , drop = FALSE]
      )
      
      df <- data_filt[, c(colnames(metaFinal), input$boxplot), drop = FALSE]
      colnames(df)[1] <- "outer_group"
      colnames(df)[2] <- "inner"
      # ---------- Parse outer selections for labels ----------
      
      outer_for_vec <- if (!is.null(top_exp())) top_exp()[[1]] else character(0)
      outer_rev_vec <- if (!is.null(bottom_exp())) bottom_exp()[[1]] else character(0)
      outer_all     <- c(outer_rev_vec,outer_for_vec)
      
      parse_outer <- function(x) {
        parts <- strsplit(x, " - ", fixed = TRUE)[[1]]
        if (length(parts) >= 2) list(var = parts[1], level = parts[2]) else list(var = x, level = x)
      }
      
      # Default labels
      cont_name       <- NULL
      outer_level_map <- NULL
      outer_levels    <- sort(unique(df$outer_group))
      
      if (length(outer_all) > 0) {
        parsed <- lapply(outer_all, parse_outer)
        cont_name <- outer_exp()                  # e.g. "FLT3.ITD"
        for_levels <- vapply(outer_for_vec, function(z) parse_outer(z)$level, character(1))
        rev_levels <- vapply(outer_rev_vec, function(z) parse_outer(z)$level, character(1))
        outer_levels <- unique(c(for_levels, rev_levels))
        
        # Map original outer_group values to pretty labels:
        # assume outer_group encodes upper/lower in your metaFinal,
        # but here we simply relabel based on order of levels we derived
        if (length(outer_levels) >= 2) {
          outer_level_map <- setNames(outer_levels, c("upper", "lower")[seq_along(outer_levels)])
        }
      }
      
      # Apply outer label mapping to df$outer_group (if available)
      df$outer_group <- factor(as.character(df$outer_group))
      if (!is.null(outer_level_map)) {
        # Map factor levels "upper"/"lower" (or whatever you used) to user-friendly labels
        df$outer_group <- factor(
          outer_level_map[as.character(df$outer_group)],
          levels = rev(outer_levels)
        )
      }
      
      # ---------- BOX PLOTS (original behavior, now with pretty x labels) ----------
      flux_col <- input$boxplot
      
      pval_to_stars <- function(p) {
        if (is.na(p)) return("ns")
        if (p < 0.001) return("***")
        if (p < 0.01)  return("**")
        if (p < 0.05)  return("*")
        "ns"
      }
      
      if (!has_inner) {
        # global test between the two outer groups
        sub <- df[!is.na(df$outer_group) & !is.na(df[[flux_col]]), ]
        if (nrow(sub) >= 2 && nlevels(sub$outer_group) == 2) {
          p_val <- tryCatch(
            {
              t.test(df[[flux_col]] ~ df$outer_group)$p.value
            },
            error = function(e) NA_real_
          )
        } else {
          p_val <- NA_real_
        }
        star_label <- pval_to_stars(p_val)
        
        # y position for the star: slightly above max
        y_max <- max(sub[[flux_col]], na.rm = TRUE)
        y_star <- y_max + 0.05 * diff(range(sub[[flux_col]], na.rm = TRUE))
        
        p_box <- ggplot(sub, aes(x = outer_group, y = .data[[flux_col]])) +
          geom_boxplot(
            width  = 0.5,
            outlier.shape = NA,
            fill   = "#e8edff",
            color  = "#5c6bc0",
            alpha  = 0.7
          ) +
          geom_jitter(
            width  = 0.12,
            alpha  = 0.5,
            size   = 1.4,
            color  = "#3949ab"
          ) +
          # add one star spanning both boxes
          geom_segment(aes(x = 1, xend = 2, y = y_star, yend = y_star),
                       size = 0.3, color = "#757575") +
          geom_text(aes(x = 1.5, y = y_star, label = star_label),
                    vjust = -0.4, size = 3.2, color = "#424242") +
          radar_theme +
          xlab(if (is.null(cont_name)) "Outer group" else cont_name) +
          ylab(flux_col)
        
      } else {
        # per-inner tests
        inner_levels <- sort(unique(df$inner))
        # compute a small data frame with positions and stars
        star_df <- do.call(rbind, lapply(inner_levels, function(lv) {
          sub <- df[df$inner == lv & !is.na(df$outer_group) & !is.na(df[[flux_col]]), ]
          if (nrow(sub) < 2 || nlevels(sub$outer_group) < 2) {
            return(NULL)
          }
          p_val <- tryCatch(
            {
              t.test(sub[[flux_col]] ~ sub$outer_group)$p.value
            },
            error = function(e) NA_real_
          )
          star_label <- pval_to_stars(p_val)
          y_max <- max(sub[[flux_col]], na.rm = TRUE)
          y_star <- y_max + 0.05 * diff(range(sub[[flux_col]], na.rm = TRUE))
          data.frame(
            inner = lv,
            x_start = 1,
            x_end   = 2,
            x_mid   = 1.5,
            y_star  = y_star,
            label   = star_label,
            stringsAsFactors = FALSE
          )
        }))
        
        p_box <- ggplot(df, aes(x = outer_group, y = .data[[flux_col]])) +
          geom_boxplot(
            width  = 0.5,
            outlier.shape = NA,
            fill   = "#e8edff",
            color  = "#5c6bc0",
            alpha  = 0.7
          ) +
          geom_jitter(
            width  = 0.12,
            alpha  = 0.5,
            size   = 1.4,
            color  = "#3949ab"
          ) +
          facet_grid(. ~ inner) + radar_theme +
          xlab(if (is.null(cont_name)) "Outer group" else cont_name) +
          ylab(flux_col)
        
        if (!is.null(star_df) && nrow(star_df) > 0) {
          p_box <- p_box +
            geom_segment(
              data = star_df,
              aes(x = x_start, xend = x_end, y = y_star, yend = y_star),
              inherit.aes = FALSE
            ) +
            geom_text(
              data = star_df,
              aes(x = x_mid, y = y_star, label = label),
              vjust = -0.3,
              inherit.aes = FALSE
            )
        }
      }
      
      p_box
    })
    
    output$plot2_obs<- renderPlot({
      plot2_obs()
    })
    
    output$plot2_exp<- renderPlot({
      plot2_exp()
    })
    

    output$dlplot2_obs <- downloadHandler(
      filename = function() {
        paste(input$boxplot,'_',label, '.pdf', sep='')
      },
      content = function(file) {
        ggsave(file, plot2_obs(), width = 12, height = 8, dpi = 400, units = "in")
      }
    )
    
    output$dlplot2_exp <- downloadHandler(
      filename = function() {
        paste(input$boxplot,'_',label, '.pdf', sep='')
      },
      content = function(file) {
        ggsave(file, plot2_exp(), width = 12, height = 8, dpi = 400, units = "in")
      }
    )

    observeEvent(input$experiment_select, {
      req(input$experiment_select)
      analysis_ready_exp(FALSE)
      session$sendCustomMessage("toggle-loading", TRUE)
      session$sendCustomMessage("loading-text", "Reading experiment registry…")
      on.exit({
        session$sendCustomMessage("toggle-loading", FALSE)
        session$sendCustomMessage("loading-text", "Experiment loaded successfully")
      }, add = TRUE)

      row <- subset(experiment_registry, folder_name == input$experiment_select)
      experiment_name(paste0(row$experiment_name, " | ", row$treatment))
      
      prefix <- row$folder_name

      # repopulate stratification choices
      flux_file          <- paste0("../../../data/RADAR_xCheck_experimental_assay/",prefix, "/", prefix, "_flux.csv")
      meta_experiment_file <- paste0("../../../data/RADAR_xCheck_experimental_assay/",prefix, "/", prefix, "_meta.csv")
      
      session$sendCustomMessage("loading-text", "Loading flux matrices…")
      flux <- as.data.frame(as.matrix(data.table::fread(flux_file), rownames = 1))
      
      session$sendCustomMessage("loading-text", "Loading experiment metadata…")
      meta <- read.csv(meta_experiment_file, header = TRUE, row.names = 1)
      
      session$sendCustomMessage("loading-text", "Finalizing experiment setup…")
      experiment_flux(flux)
      experiment_meta(meta)
      
      result <- lapply(experiment_meta(), function(col) {
        list(unique(col))
      })
      
      updateVirtualSelect("strat_exp", choices = result,selected = unique(unlist(experiment_meta())))

    })

    # Apply stratification
    observeEvent(list(input$submitStrat_exp), {
      req(input$submitStrat_exp > 0)
      analysis_ready_exp(FALSE)
      session$sendCustomMessage("toggle-loading", TRUE)
      session$sendCustomMessage("loading-text", "Reading experiment registry…")
      on.exit({
        session$sendCustomMessage("toggle-loading", FALSE)
        session$sendCustomMessage("loading-text", "Experiment loaded successfully")
      }, add = TRUE)
      session$sendCustomMessage("loading-text", "Loading experimental parameters…")
      if (!isTRUE(exp_outer_ui_inserted())) {
        insertUI(
          selector = "#stratification_exp",
          where = "afterEnd",
          immediate = TRUE,
          ui = tags$div(
            id = "exp_outer_inner_row",
            fluidRow(
            class = "fade-in-up",
            column(
              width = 6,
              box(
                width = 12,
                title = "Select Outer Comparison",
                solidHeader = TRUE,
                status = "primary",
                selectizeInput(
                  "outerType_exp", "Select Outer Comparison",
                  choices = character(0)
                ),
                column(
                  width = 6,
                  virtualSelectInput(
                    "outerVarTop_exp", "Select Top Variable",
                    choices = list(), multiple = TRUE,
                    showValueAsTags = TRUE, search = TRUE
                  )
                ),
                column(
                  width = 6,
                  virtualSelectInput(
                    "outerVarBottom_exp", "Select Bottom Variable",
                    choices = list(), multiple = TRUE,
                    showValueAsTags = TRUE, search = TRUE
                  )
                ),
                fluidRow(
                  column(
                    width = 12,
                    actionButton("submitOuter_exp", "Submit"),
                    div(id = "outer_exp")
                  )
                )
              ),
              tags$div(id = "outer_exp")   # anchor inside left column
            ),
            column(
              width = 6,
              tags$div(id = "inner_exp_col")  # empty right column anchor
            )
          )
        )
        )
      }
      exp_outer_ui_inserted(TRUE)
      matching_rows <- apply(experiment_meta(), 1, function(row) all(row %in% input$strat_exp))
      filtered_df <- experiment_meta()[matching_rows, , drop = FALSE]
      metaFinal_out(filtered_df)

      outer_exp(NULL)
      inner_exp(NULL)
      top_exp(NULL)
      bottom_exp(NULL)
      innerVars_exp(NULL)
      error_exp(NULL)

      # Insert outer box only once
      

      updateSelectizeInput(
        session, "outerType_exp",
        choices = colnames(filtered_df)[colnames(filtered_df) != "none"],
        selected = NULL
      )
    })

    # Outer selection updates
    observeEvent(input$outerType_exp, {
      req(metaFinal_out())
      meta <- metaFinal_out()
      updateVirtualSelect("outerVarTop_exp",
                          choices = unique(meta[[input$outerType_exp]]))
      updateVirtualSelect("outerVarBottom_exp",
                          choices = unique(meta[[input$outerType_exp]]))
    })

    observeEvent(input$submitOuter_exp, {
      req(input$submitOuter_exp)
      analysis_ready_exp(FALSE)
      inner_exp(NULL)
      innerVars_exp(NULL)
      error_exp(NULL)
      
      if (is.null(input$outerVarTop_exp) && is.null(input$outerVarBottom_exp)) {
        showModal(
          modalDialog(
            title = "Missing outer comparison",
            "At least one outer comparison needs to be selected for upper/lower.",
            easyClose = TRUE,
            footer = modalButton("OK")
          )
        )
        return()  # do nothing else
      }

      outer_exp(input$outerType_exp)
      top_exp(input$outerVarTop_exp)
      bottom_exp(input$outerVarBottom_exp)

      # Insert inner UI only once
      if (!isTRUE(exp_inner_ui_inserted())) {
        insertUI(
          selector = "#inner_exp_col",
          where = "beforeEnd",
          immediate = TRUE,
          ui = 
            box(
              width = 12,
              title = "Select Inner Comparison",
              solidHeader = TRUE,
              status = "primary",
              selectizeInput(
                "innerType_exp", "Select Inner Comparison",
                choices = character(0)
              ),
              virtualSelectInput(
                "innerVar_exp", "Select Variables",
                choices = list(), multiple = TRUE,
                showValueAsTags = TRUE, search = TRUE
              ),
              fluidRow(
                column(12, actionButton("submitInner_exp", "Submit")),
              )
            )

        )
      }
      exp_inner_ui_inserted(TRUE)
      meta <- metaFinal_out()
      updateSelectizeInput(
        session, "innerType_exp",
        choices = colnames(meta)[colnames(meta) != outer_exp()],
        selected = NULL
      )
    })

    observeEvent(list(input$submitOuter_exp, input$innerType_exp), {
      req(input$submitOuter_exp, input$innerType_exp)
      meta <- metaFinal_out()
      updateVirtualSelect(
        "innerVar_exp",
        choices = unique(meta[[input$innerType_exp]]),
        selected = unique(meta[[input$innerType_exp]])
      )
    })

    observeEvent(input$submitInner_exp, {
      req(input$submitInner_exp)
      analysis_ready_exp(FALSE)
      inner_exp(input$innerType_exp)
      innerVars_exp(input$innerVar_exp)
      error_exp(NULL)

      # Insert summary box
      if(!isTRUE(exp_summary_ui_inserted())){
      insertUI(
        selector = "#exp_outer_inner_row",
        where = "afterEnd",
        immediate = TRUE,
        ui = fluidRow(
          column(
            width = 12,
            box(
              width = 12,
              title = "RADAR xCheck Pre-Analysis Summary",
              solidHeader = TRUE,
              status = "primary",
              column(
                width = 6,
                verbatimTextOutput("reportText_exp")     ### CHANGED
              ),
              column(
                width = 6,
                noUiSliderInput("customAUC_exp", "AUC",   ### CHANGED
                                min = 0.6, max = 1.0, value = 0.65, step = 0.01),
                noUiSliderInput("customFC_exp", "Log2FC", ### CHANGED
                                min = 0, max = 5, value = 1, step = 0.01),
                noUiSliderInput("customFDR_exp", "FDR",   ### CHANGED
                                min = 0, max = 1, value = 0.05, step = 0.01)
              ),
              fluidRow(
                column(
                  width = 12,
                  actionButton("submitReport_exp", "Begin Analysis (set parameters)"),
                  actionButton("submitReportDefault_exp", "Begin Analysis (p < 0.05)")
                )
              )
            )
          )
        )
      )
      }
      exp_summary_ui_inserted(TRUE)
    })

    # Summary text
    observeEvent(
      list(input$submitOuter_exp, input$submitInner_exp),
      {
        req(input$submitInner_exp)
        warnings <- character(0)
        errors   <- character(0)
        
        ###################################
        # 1. STRATIFIED EXPERIMENT SAMPLES
        ###################################
        
        # experiment_meta() should be a data.frame with per‑sample annotations
        # and rownames matching rows of experiment_flux()
        meta_exp <- experiment_meta()
        
        # Here you implement whatever logic you use to select samples
        # based on input$strat_exp (analogous to clinical_strat/gene_strat).
        # The result should be a character vector of kept sample IDs:
        strat_names_exp <- rownames(meta_exp)   # <- placeholder: keep all for now
        
        if (length(strat_names_exp) < 2) {
          errors <- c(errors, "Stratified experiment must contain at least two samples")
          metaFinal_out(NULL)
          errors_list(errors)
          output$reportText_exp <- renderText(paste(
            "Experiment:", experiment_name(), "\n",
            "Errors:", paste(errors, collapse = "; "), "\n"
          ))
          return()
        }
        
        ################################
        # 2. OUTER LABELS FOR EXPERIMENT
        ################################
        
        # Build a 1‑column data.frame with outer labels "upper"/"lower"
        # using your experiment‑specific UI state (outer_exp(), top_exp(), bottom_exp()).
        # For simplicity, assume outer_exp() is a single variable name in meta_exp,
        # and top_exp()/bottom_exp() are sets of levels considered "upper"/"lower".
        
        outer_var  <- outer_exp()                  # e.g. "treatment"
        top_levels <- top_exp()                    # e.g. c("DrugA")
        bot_levels <- bottom_exp()                 # e.g. c("Control")
        
        outer_df_exp <- meta_exp[strat_names_exp, outer_var, drop = FALSE]
        colnames(outer_df_exp) <- outer_var
        
        # Map to "upper"/"lower"
        outer_label <- ifelse(
          outer_df_exp[[outer_var]] %in% top_levels,  "upper",
          ifelse(outer_df_exp[[outer_var]] %in% bot_levels, "lower", NA_character_)
        )
        
        valid_rows <- which(!is.na(outer_label))
        if (length(valid_rows) == 0) {
          errors <- c(errors, "No experiment samples remain after outer comparison filters")
          metaFinal_out(NULL)
          errors_list(errors)
          output$reportText_exp <- renderText(paste(
            "Experiment:", experiment_name(), "\n",
            "Errors:", paste(errors, collapse = "; "), "\n"
          ))
          return()
        }
        
        outer_label <- outer_label[valid_rows]
        sample_ids  <- strat_names_exp[valid_rows]
        
        ################################
        # 3. INNER LABELS FOR EXPERIMENT
        ################################
        
        inner_df_exp <- data.frame(
          outer = outer_label,
          row.names = sample_ids,
          stringsAsFactors = FALSE
        )
        
        if (is.null(inner_exp()) || inner_exp() == "" || length(innerVars_exp()) == 0) {
          # no inner comparison
          inner_df_exp$inner <- "All"
        } else {
          inner_var <- inner_exp()   # e.g. "timepoint"
          vals <- meta_exp[sample_ids, inner_var, drop = TRUE]
          inner_df_exp$inner <- vals
          
          # Optionally discretize continuous inner variables,
          # analogous to the cohort logic (q25/q75 or median)
          # depending on your inner comparison type for experiments.
        }
        
        # Drop any rows with missing outer/inner
        inner_df_exp <- na.omit(inner_df_exp)
        if (nrow(inner_df_exp) == 0) {
          errors <- c(errors, "No experiment samples remain after inner comparison filters")
          metaFinal_out(NULL)
          errors_list(errors)
          output$reportText_exp <- renderText(paste(
            "Experiment:", experiment_name(), "\n",
            "Errors:", paste(errors, collapse = "; "), "\n"
          ))
          return()
        }
        
        #########################
        # 4. STORE FOR ANALYSIS #
        #########################
        
        metaFinal_out(inner_df_exp)
        errors_list(errors)
        
        output$reportText_exp <- renderText({
          paste(
            paste0("Experiment: ", experiment_name()),
            paste0("Outer Category: ", outer_exp()),
            paste0("Outer Variables: ", paste(top_exp(), collapse = ", "),
                   " (top) vs ", paste(bottom_exp(), collapse = ", "), " (bottom)"),
            paste0("Inner Category: ", inner_exp()),
            paste0("Inner Variables: ", paste(innerVars_exp(), collapse = ", ")),
            paste0("Minimum AUC: ", input$customAUC_exp),
            paste0("Minimum Log2FC: ", input$customFC_exp),
            paste0("Maximum FDR: ", input$customFDR_exp),
            paste0("Errors: ", if (length(errors) == 0) "NONE" else paste(errors, collapse = "; ")),
            paste0("Report Generated at: ", Sys.time()),
            sep = "\n"
          )
        })
        
      }
    )
 
    observeEvent(input$main_tabs, {
      if (input$main_tabs == "obs") {
        ## Leaving experimental tab → reset experimental state
        
        analysis_ready_exp(FALSE)
        
      } else if (input$main_tabs == "exp") {
        ## Leaving cohort tab → reset cohort state
        
        analysis_ready(FALSE)
        
      }
    })
}

shinyApp(ui,server)
