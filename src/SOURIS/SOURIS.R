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
library(ggpubr)
library(Seurat)
library(caret)
library(cowplot)
library(Seurat)

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
    menuItem("Single-Cell", tabName = "sc", icon = icon("table-cells-large"))
  )
)

body <- dashboardBody(
  useShinyjs(),
  tags$style(HTML("
  :root {
    --souris-navy:  #002B5C;
    --souris-mid:   #0D47A1;
    --souris-light: #90CAF9;
    --souris-bg:    #F4F6FB;
  }
  .souris-plot-fade {
    opacity: 0;
    transform: translateY(6px);
    transition: opacity 1.2s ease-out, transform 1.2s ease-out;
  }
  .souris-plot-fade.souris-plot-visible {
    opacity: 1;
    transform: translateY(0);
    animation: souris-flash 0.8s ease-out 1;
  }
  @keyframes souris-flash {
    0%, 100% {
      box-shadow: 0 0 0 rgba(144,202,249,0);
      transform: translateY(0) scale(1);
    }
    40%, 60% {
      box-shadow: 0 0 16px rgba(144,202,249,0.9);
      transform: translateY(0) scale(1.01);
    }
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
    border: none !important;
    box-shadow: 0 10px 28px rgba(0, 0, 0, 0.1);
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
    background-color: #E8F1FB;
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

  /* SOURIS full-page loading overlay */
  #souris-loading-overlay.loading-overlay {
    position: fixed;
    inset: 0;
    display: flex;
    align-items: center;
    justify-content: center;
    background: radial-gradient(
      circle at top,
      rgba(0, 43, 92, 0.70),
      rgba(13, 71, 161, 0.50) 60%,
      rgba(0, 20, 40, 0.40)   100%
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
  #souris-auc-hint {
    position: absolute;
    top: 8px;
    right: 12px;
    z-index: 10;
    pointer-events: none;
    opacity: 0;
    transform: translateY(-4px);
    transition: opacity 0.5s ease, transform 0.5s ease;
  }

  #souris-auc-hint .souris-auc-pill {
    display: inline-flex;
    align-items: center;
    gap: 6px;
    padding: 4px 10px;
    border-radius: 999px;
    background: rgba(13, 71, 161, 0.92); /* SOURIS mid blue */
    box-shadow: 0 4px 12px rgba(0, 0, 0, 0.18);
    color: #E3F2FD;
    font-size: 11px;
    letter-spacing: 0.04em;
    text-transform: uppercase;
  }

  #souris-auc-hint .souris-auc-dot {
    width: 6px;
    height: 6px;
    border-radius: 999px;
    background: #BBDEFB;
    box-shadow: 0 0 8px rgba(187, 222, 251, 0.9);
    animation: sourisPulseDot 1.4s ease-in-out infinite;
  }

  #souris-auc-hint .souris-auc-text {
    white-space: nowrap;
  }

  @keyframes sourisPulseDot {
    0%   { transform: scale(1);   opacity: 0.7; }
    50%  { transform: scale(1.6); opacity: 1.0; }
    100% { transform: scale(1);   opacity: 0.7; }
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

  /* Base hex: subtle vibration */
  .souris-logo-hex {
    position: relative;
    width: 110px;
    height: 110px;
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
    box-shadow:
      0 0 0 0 rgba(255,255,255,0),
      0 0 0 0 rgba(255,255,255,0);
  }

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

  /* Atmospheric statistical glyphs behind the hex */
  .souris-loading-atmosphere {
    position: absolute;
    inset: 0;
    z-index: 0;
    pointer-events: none;
    overflow: visible;
  }

  .souris-loading-atmosphere .souris-stat-glyph {
    position: absolute;
    width: 40px;
    height: 90px;
    opacity: 0;
    transform-origin: 50% 100%;
    animation: sourisStatAssemble 5s ease-in-out infinite;
  }

 .souris-stat-box {
    position: absolute;
    left: 14px;
    width: 12px;
    height: 40px;
    top: 25px;
    border-radius: 4px;
    border: 1px solid rgba(144, 202, 249, 0.35);
    background: linear-gradient(
      to bottom,
      rgba(144, 202, 249, 0.12),
      rgba(13, 71, 161, 0.26)
    );
  }

   .souris-stat-line {
    position: absolute;
    left: 20px;
    top: 8px;
    width: 2px;
    height: 75px;
    background: linear-gradient(
      to bottom,
      rgba(227, 242, 253, 0.3),
      rgba(13, 71, 161, 0.85)
    );
    box-shadow: 0 0 6px rgba(144, 202, 249, 0.4);
    border-radius: 999px;
  }

  .souris-stat-circle {
    position: absolute;
    width: 6px;
    height: 6px;
    border-radius: 999px;
    background: radial-gradient(circle, #E3F2FD 0%, #90CAF9 60%, rgba(0,43,92,0.0) 100%);
    box-shadow: 0 0 8px rgba(144, 202, 249, 0.7);
  }

/* Positions: corners, edges, inner ring */
   .souris-stat-glyph.g1  { top:  3%;  left:  6%;  animation-delay: 0.0s; }
  .souris-stat-glyph.g2  { top: 12%;  left: 38%; animation-delay: 0.3s; }
  .souris-stat-glyph.g3  { top:  7%;  left: 72%; animation-delay: 0.6s; }
  .souris-stat-glyph.g4  { top: 18%;  right: 4%;  animation-delay: 0.9s; }

  .souris-stat-glyph.g5  { top: 26%;  left: 14%; animation-delay: 0.5s; }
  .souris-stat-glyph.g6  { top: 20%;  left: 55%; animation-delay: 0.8s; }
  .souris-stat-glyph.g7  { top: 32%;  right: 22%; animation-delay: 1.1s; }

  .souris-stat-glyph.g8  { top: 40%;  left: 8%;  animation-delay: 0.7s; }
  .souris-stat-glyph.g9  { top: 44%;  left: 48%; animation-delay: 1.0s; }
  .souris-stat-glyph.g10 { top: 38%;  right: 10%; animation-delay: 1.3s; }

  .souris-stat-glyph.g11 { top: 56%;  left: 22%; animation-delay: 0.9s; }
  .souris-stat-glyph.g12 { top: 60%;  left: 68%; animation-delay: 1.2s; }
  .souris-stat-glyph.g13 { top: 52%;  right: 16%; animation-delay: 1.5s; }

  .souris-stat-glyph.g14 { bottom: 26%; left: 10%; animation-delay: 1.1s; }
  .souris-stat-glyph.g15 { bottom: 18%; left: 46%; animation-delay: 1.4s; }
  .souris-stat-glyph.g16 { bottom: 22%; right: 20%; animation-delay: 1.7s; }

  .souris-stat-glyph.g17 { bottom: 8%;  left: 18%; animation-delay: 1.3s; }
  .souris-stat-glyph.g18 { bottom: 5%;  right: 6%;  animation-delay: 1.6s; }

  /* Subtle scale/rotation variation */
  .souris-stat-glyph.g2  { transform: scale(0.85) rotate(-6deg); }
  .souris-stat-glyph.g3  { transform: scale(0.9)  rotate(4deg);  }
  .souris-stat-glyph.g4  { transform: scale(0.8)  rotate(-3deg); }
  .souris-stat-glyph.g5  { transform: scale(0.95) rotate(5deg);  }
  .souris-stat-glyph.g6  { transform: scale(0.9)  rotate(-5deg); }
  .souris-stat-glyph.g7  { transform: scale(0.85) rotate(3deg);  }
  .souris-stat-glyph.g8  { transform: scale(0.8)  rotate(-4deg); }
  .souris-stat-glyph.g9  { transform: scale(0.9)  rotate(2deg);  }
  .souris-stat-glyph.g10 { transform: scale(0.85) rotate(-2deg); }
  .souris-stat-glyph.g11 { transform: scale(0.9)  rotate(4deg);  }
  .souris-stat-glyph.g12 { transform: scale(0.8)  rotate(-5deg); }
  .souris-stat-glyph.g13 { transform: scale(0.9)  rotate(3deg);  }
  .souris-stat-glyph.g14 { transform: scale(0.85) rotate(-4deg); }
  .souris-stat-glyph.g15 { transform: scale(0.8)  rotate(2deg);  }
  .souris-stat-glyph.g16 { transform: scale(0.9)  rotate(-3deg); }
  .souris-stat-glyph.g17 { transform: scale(0.85) rotate(5deg);  }
  .souris-stat-glyph.g18 { transform: scale(0.8)  rotate(-2deg); }

  @keyframes sourisStatAssemble {
    0% {
      opacity: 0;
      transform: translateY(22px) scale(0.85);
      filter: blur(4px);
    }
    20% {
      opacity: 0.25;
      transform: translateY(10px) scale(0.95);
      filter: blur(2px);
    }
    40% {
      opacity: 0.6;
      transform: translateY(0) scale(1.0);
      filter: blur(0.6px);
    }
    70% {
      opacity: 0.6;
      transform: translateY(-5px) scale(1.03);
      filter: blur(1.2px);
    }
    100% {
      opacity: 0;
      transform: translateY(-16px) scale(0.9);
      filter: blur(4px);
    }
  }
")),

  # Loading overlay container
  div(
    id = "souris-loading-overlay", class = "loading-overlay hidden",

    # Full-screen atmospheric glyph layer
    div(
      class = "souris-loading-atmosphere",
      div(
        class = "souris-stat-glyph g1",
        div(class = "souris-stat-line"),
        div(class = "souris-stat-box"),
        div(class = "souris-stat-circle", style = "top: 4px;  left: 17px;"),
        div(class = "souris-stat-circle", style = "top: 70px; left: 26px;")
      ),
      div(
        class = "souris-stat-glyph g2",
        div(class = "souris-stat-line"),
        div(class = "souris-stat-box", style = "top: 30px; height: 32px;"),
        div(class = "souris-stat-circle", style = "top: 10px; left: 22px;")
      ),
      div(
        class = "souris-stat-glyph g3",
        div(class = "souris-stat-line"),
        div(class = "souris-stat-box", style = "top: 22px; height: 45px;"),
        div(class = "souris-stat-circle", style = "top: 6px;  left: 18px;"),
        div(class = "souris-stat-circle", style = "top: 78px; left: 15px;")
      ),
      div(
        class = "souris-stat-glyph g4",
        div(class = "souris-stat-line"),
        div(class = "souris-stat-box", style = "top: 26px; height: 36px;"),
        div(class = "souris-stat-circle", style = "top: 8px;  left: 21px;")
      ),
      div(
        class = "souris-stat-glyph g5",
        div(class = "souris-stat-line"),
        div(class = "souris-stat-box", style = "top: 20px; height: 42px;"),
        div(class = "souris-stat-circle", style = "top: 4px;  left: 18px;"),
        div(class = "souris-stat-circle", style = "top: 72px; left: 23px;")
      ),
      div(
        class = "souris-stat-glyph g6",
        div(class = "souris-stat-line"),
        div(class = "souris-stat-box", style = "top: 28px; height: 30px;"),
        div(class = "souris-stat-circle", style = "top: 12px; left: 20px;")
      ),
      div(
        class = "souris-stat-glyph g7",
        div(class = "souris-stat-line"),
        div(class = "souris-stat-box", style = "top: 24px; height: 40px;"),
        div(class = "souris-stat-circle", style = "top: 6px;  left: 19px;"),
        div(class = "souris-stat-circle", style = "top: 76px; left: 24px;")
      ),
      div(
        class = "souris-stat-glyph g8",
        div(class = "souris-stat-line"),
        div(class = "souris-stat-box", style = "top: 32px; height: 34px;"),
        div(class = "souris-stat-circle", style = "top: 14px; left: 21px;")
      ),
      div(
        class = "souris-stat-glyph g9",
        div(class = "souris-stat-line"),
        div(class = "souris-stat-box", style = "top: 26px; height: 44px;"),
        div(class = "souris-stat-circle", style = "top: 8px;  left: 19px;"),
        div(class = "souris-stat-circle", style = "top: 74px; left: 16px;")
      ),
      div(
        class = "souris-stat-glyph g10",
        div(class = "souris-stat-line"),
        div(class = "souris-stat-box", style = "top: 30px; height: 32px;"),
        div(class = "souris-stat-circle", style = "top: 12px; left: 21px;")
      ),
      div(
        class = "souris-stat-glyph g11",
        div(class = "souris-stat-line"),
        div(class = "souris-stat-box", style = "top: 22px; height: 46px;"),
        div(class = "souris-stat-circle", style = "top: 6px;  left: 18px;"),
        div(class = "souris-stat-circle", style = "top: 78px; left: 25px;")
      ),
      div(
        class = "souris-stat-glyph g12",
        div(class = "souris-stat-line"),
        div(class = "souris-stat-box", style = "top: 28px; height: 38px;"),
        div(class = "souris-stat-circle", style = "top: 10px; left: 20px;")
      ),
      div(
        class = "souris-stat-glyph g13",
        div(class = "souris-stat-line"),
        div(class = "souris-stat-box", style = "top: 24px; height: 42px;"),
        div(class = "souris-stat-circle", style = "top: 8px;  left: 19px;"),
        div(class = "souris-stat-circle", style = "top: 76px; left: 23px;")
      ),
      div(
        class = "souris-stat-glyph g14",
        div(class = "souris-stat-line"),
        div(class = "souris-stat-box", style = "top: 30px; height: 34px;"),
        div(class = "souris-stat-circle", style = "top: 12px; left: 21px;")
      ),
      div(
        class = "souris-stat-glyph g15",
        div(class = "souris-stat-line"),
        div(class = "souris-stat-box", style = "top: 22px; height: 40px;"),
        div(class = "souris-stat-circle", style = "top: 6px;  left: 18px;"),
        div(class = "souris-stat-circle", style = "top: 74px; left: 24px;")
      ),
      div(
        class = "souris-stat-glyph g16",
        div(class = "souris-stat-line"),
        div(class = "souris-stat-box", style = "top: 28px; height: 36px;"),
        div(class = "souris-stat-circle", style = "top: 10px; left: 20px;")
      ),
      div(
        class = "souris-stat-glyph g17",
        div(class = "souris-stat-line"),
        div(class = "souris-stat-box", style = "top: 26px; height: 42px;"),
        div(class = "souris-stat-circle", style = "top: 8px;  left: 19px;"),
        div(class = "souris-stat-circle", style = "top: 76px; left: 22px;")
      ),
      div(
        class = "souris-stat-glyph g18",
        div(class = "souris-stat-line"),
        div(class = "souris-stat-box", style = "top: 30px; height: 34px;"),
        div(class = "souris-stat-circle", style = "top: 12px; left: 21px;")
      )
    ),

    # Centered loading core on top
    div(
      class = "souris-loading-core",
      div(
        class = "souris-logo-hex",
        div(class = "souris-logo-hex-glow"),
        tags$span(class = "souris-spoke spoke-1"),
        tags$span(class = "souris-spoke spoke-2"),
        tags$span(class = "souris-spoke spoke-3")
      ),
      div(id = "souris-loading-text", "Preparing SOURIS analysis"),
      div(id = "souris-loading-subtext", ""),
      div(
        id = "souris-loading-bar",
        div(id = "souris-loading-bar-fill")
      )
    )
  ),

  # JS handlers for SOURIS loading overlay
  tags$script(HTML("
  Shiny.addCustomMessageHandler('souris-auc-header-subtitle', function(msg) {
    var el = document.getElementById('souris-auc-subtitle');
    if (!el) return;
    el.textContent = msg;
  });
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
    Shiny.addCustomMessageHandler('souris-loading-progress', function(pct) {
      var fill = document.getElementById('souris-loading-bar-fill');
      if (!fill) return;
      var v = Math.max(0, Math.min(100, pct));
      fill.style.width = v + '%';
    });
    Shiny.addCustomMessageHandler('souris-show-plot', function(id) {
    var el = document.getElementById(id);
    if (!el) return;
    el.classList.add('souris-plot-visible');
  });

    function sourisShowAucHint() {
    var wrapper = document.getElementById('souris-auc-wrapper');
    var hint    = document.getElementById('souris-auc-hint');
    if (!wrapper || !hint) return;

    // Only create once
    if (!hint.hasChildNodes()) {
      var pill  = document.createElement('div');
      pill.className = 'souris-auc-pill';

      var dot = document.createElement('div');
      dot.className = 'souris-auc-dot';

      var txt = document.createElement('div');
      txt.className = 'souris-auc-text';
      txt.textContent = 'Click a bar for reaction details';

      pill.appendChild(dot);
      pill.appendChild(txt);
      hint.appendChild(pill);
    }

    // Fade in
    hint.style.opacity   = '1';
    hint.style.transform = 'translateY(0px)';

    // Auto-fade out after a few seconds
    setTimeout(function() {
      hint.style.opacity   = '0';
      hint.style.transform = 'translateY(-4px)';
    }, 4500);
  }

  // Expose hook to Shiny
  Shiny.addCustomMessageHandler('souris-show-auc-hint', function(msg) {
    sourisShowAucHint();
  });
  ")),
  tabItems(

    #------------------------------------------------------
    # Bulk Samples tab (CORALIE-like panels, 1 fingerprint)
    #------------------------------------------------------
    tabItem(
      tabName = "bulk",
      class = "fade-in-up",

      # Fingerprint selection
      fluidRow(
        box(
          width = 12, status = "primary", solidHeader = TRUE,
          title = "Fingerprint selection",
          shinyWidgets::pickerInput(
            inputId = "souris_fingerprints_sel",
            label = "Select one fingerprint",
            choices = NULL, # populated server-side
            multiple = FALSE, # <<< only one fingerprint allowed
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
      class = "fade-in-up",

      # Fingerprint selection (same as bulk)
      fluidRow(
        box(
          width = 12, status = "primary", solidHeader = TRUE,
          title = "Fingerprint selection",
          shinyWidgets::pickerInput(
            inputId = "souris_fingerprints_sel_sc",
            label = "Select one fingerprint",
            choices = NULL,
            multiple = FALSE,
            options = list(
              `actions-box` = FALSE,
              `none-selected-text` = "No fingerprint selected"
            )
          ),
          br(),
          actionButton(
            inputId = "souris_load_fingerprints_sc",
            label = "Select fingerprint"
          )
        )
      ),

      # Placeholder where the reference selection box will be inserted
      fluidRow(div(id = "sc_reference_placeholder")),

      # Placeholder where the single-cell object selector will be inserted
      fluidRow(div(id = "sc_scobj_placeholder")),

      # (Later: analysis panels for single-cell)
      fluidRow(div(id = "sc_analysis_placeholder")),
      fluidRow(div(id = "sc_tier3_placeholder"))
    )
  )
)

ui <- dashboardPage(
  skin = "purple",
  header = header,
  sidebar = sidebar,
  body = body
)

#----------------------------------------------------------
# SERVER
#----------------------------------------------------------

server <- function(input, output, session) {
  # paths analogous to CORALIE
  fingerprints_dir <- "../../data/fingerprints"
  reference_dir <- "../../data/reference_assays"
  sc_dir <- "../../data/sc_objects"

  # reactive state
  souris_fingerprints_list <- reactiveVal(NULL)
  souris_reference_list <- reactiveVal(NULL)

  ref_box_inserted_bulk <- reactiveVal(FALSE)
  summary_box_inserted <- reactiveVal(FALSE)
  analysis_panel_inserted <- reactiveVal(FALSE)
  analysis_done_bulk <- reactiveVal(FALSE)

  df_all_predictions <- reactiveVal(NULL)
  df_all_scores <- reactiveVal(NULL)

  detail_panel_inserted <- reactiveVal(FALSE)
  clicked_subsystem <- reactiveVal(NULL)

  reaction_order <- reactiveVal(NULL)

  subsystem_done <- reactiveVal(FALSE)
  reactions_done <- reactiveVal(FALSE)

  sc_ref_box_inserted <- reactiveVal(FALSE)
  sc_scobj_box_inserted <- reactiveVal(FALSE)
  souris_sc_object <- reactiveVal(NULL)

  sc_analysis_panel_inserted <- reactiveVal(FALSE)

  sc_fingerprint_scores <- reactiveVal(NULL)

  clicked_subsystem_sc <- reactiveVal(NULL)
  sc_detail_panel_inserted <- reactiveVal(FALSE)

  clicked_reaction_sc <- reactiveVal(NULL)
  sc_reaction_feature_panel_inserted <- reactiveVal(FALSE)

  # helper for loading overlay
  set_souris_progress <- function(text = NULL, subtext = NULL, pct = NULL, show = NULL) {
    if (!is.null(text)) {
      session$sendCustomMessage("souris-loading-text", text)
    }
    if (!is.null(subtext)) {
      session$sendCustomMessage("souris-loading-subtext", subtext)
    }
    if (!is.null(pct)) {
      session$sendCustomMessage("souris-loading-progress", round(pct))
    }
    if (!is.null(show)) {
      session$sendCustomMessage("souris-toggle-loading", show)
    }
  }

  auc_two_groups <- function(pred, grp) {
    # grp: factor with 2 levels, length(pred)
    ord <- order(pred)
    pred <- pred[ord]
    grp <- grp[ord]
    y <- as.integer(grp == levels(grp)[2]) # 0/1
    n1 <- sum(y == 0)
    n2 <- sum(y == 1)
    if (n1 == 0 || n2 == 0) {
      return(NA_real_)
    }
    # Mann–Whitney U / AUC
    ranks <- rank(pred)
    R2 <- sum(ranks[y == 1])
    U <- R2 - n2 * (n2 + 1) / 2
    U / (n1 * n2)
  }

  human_react_meta <- read.csv("../../data/fingerprint_prep_objects/human_reaction_meta.csv", header = T, row.names = 1)
  mouse_react_meta <- read.csv("../../data/fingerprint_prep_objects/mouse_reaction_meta.csv", header = T, row.names = 1)
  react_meta_filt <- human_react_meta[intersect(rownames(human_react_meta), rownames(mouse_react_meta)), ]

  get_recon_names <- function(reaction_ids) {
    # reaction_ids are current column names (internal IDs)
    common_ids <- intersect(reaction_ids, rownames(react_meta_filt))
    recon <- as.character(react_meta_filt[common_ids, "RECON3D"])
    # fall back to internal ID if RECON3D is NA/empty
    recon[is.na(recon) | recon == ""] <- common_ids[is.na(recon) | recon == ""]
    out <- setNames(recon, common_ids)
    out
  }

  get_reaction_preds <- function(ref_id, subsel) {
    # 1) Get the reference assay matrix
    ref_mat <- souris_reference_list()[[ref_id]]
    validate(need(
      !is.null(ref_mat),
      paste("No reaction matrix found for reference assay", ref_id, ".")
    ))

    # 2) Select reactions whose subsystem matches subsel
    #    (make.names if you normalized names earlier)
    subs_react <- rownames(react_meta_filt)[react_meta_filt$subsystem == subsel]
    subs_react <- intersect(subs_react, colnames(ref_mat))
    validate(need(
      length(subs_react) > 0,
      paste("No reactions in react_meta_filt for subsystem", subsel, "in", ref_id, ".")
    ))

    # 3) Return the reference matrix restricted to those reaction columns
    mat <- as.matrix(ref_mat[, subs_react, drop = FALSE])
    rownames(mat) <- rownames(ref_mat) # samples
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

    shinyWidgets::updatePickerInput(
      session,
      inputId = "souris_fingerprints_sel_sc", # single-cell
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

  observe({
    nsel <- length(input$souris_fingerprints_sel_sc)
    if (is.null(nsel) || nsel != 1) {
      shinyjs::disable("souris_load_fingerprints_sc")
    } else {
      shinyjs::enable("souris_load_fingerprints_sc")
    }
  })

  #--------------------------------------------------------
  # Load fingerprint (single) and insert reference box
  #--------------------------------------------------------
  observeEvent(input$souris_load_fingerprints, {
    req(input$souris_fingerprints_sel)
    validate(need(
      length(input$souris_fingerprints_sel) == 1,
      "Please select exactly one fingerprint."
    ))

    set_souris_progress("Loading fingerprint", "", TRUE)

    sel_name <- input$souris_fingerprints_sel
    path <- file.path(fingerprints_dir, paste0(sel_name, "_.rds"))

    size <- file.info(path)$size
    total_Size <- ifelse(is.na(size), 1, size)
    cum_Size <- 0


    obj <- readRDS(path)
    cum_Size <- total_Size
    set_souris_progress(
      subtext = sprintf("Loading %s", sel_name),
      pct     = 100 * cum_Size / total_Size
    )

    souris_fingerprints_list(setNames(list(obj), sel_name))

    set_souris_progress(subtext = "", pct = 100, show = FALSE)

    set_souris_progress(show = FALSE)

    # Insert reference assay selection box once
    if (!ref_box_inserted_bulk()) {
      insertUI(
        selector = "#bulk_reference_placeholder",
        where = "afterBegin",
        ui = fluidRow(
          class = "fade-in-up",
          column(
            width = 12,
            box(
              width = 12, status = "primary", solidHeader = TRUE,
              title = "Reference assay selection",
              selectInput(
                inputId = "souris_reference_sel",
                label = "Select exactly two reference assays",
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

  observeEvent(input$souris_load_fingerprints_sc, {
    req(input$souris_fingerprints_sel_sc)
    validate(
      need(
        length(input$souris_fingerprints_sel_sc) == 1,
        "Please select exactly one fingerprint."
      )
    )

    set_souris_progress("Loading fingerprint", NULL, NULL, TRUE)

    sel_name <- input$souris_fingerprints_sel_sc
    path <- file.path(fingerprints_dir, paste0(sel_name, "_.rds"))

    size <- file.info(path)$size
    total_Size <- ifelse(is.na(size) || size <= 0, 1, size)
    cum_Size <- 0

    obj <- readRDS(path)
    cum_Size <- total_Size

    set_souris_progress(
      subtext = sprintf("Loading %s", sel_name),
      pct     = 100 * cum_Size / total_Size
    )

    # reuse the same fingerprint list (overwrites with SC choice, which is fine)
    souris_fingerprints_list(setNames(list(obj), sel_name))

    souris_sc_object(NULL)
    sc_fingerprint_scores(NULL)
    df_all_predictions(NULL)
    df_all_scores(NULL)

    set_souris_progress(subtext = "", pct = 100, show = FALSE)

    # Insert SC reference assay selection box once
    if (!isTRUE(sc_ref_box_inserted())) {
      insertUI(
        selector = "#sc_reference_placeholder",
        where = "afterBegin",
        ui = fluidRow(
          class = "fade-in-up",
          column(
            width = 12,
            box(
              width = 12, status = "primary", solidHeader = TRUE,
              title = "Reference assay selection",
              selectInput(
                inputId = "souris_reference_sel_sc",
                label = "Select exactly two reference assays",
                choices = NULL,
                multiple = TRUE
              ),
              actionButton(
                inputId = "souris_load_reference_sc",
                label = "Select reference assays"
              )
            )
          )
        ),
        immediate = TRUE
      )
      sc_ref_box_inserted(TRUE)
    }
  })

  #--------------------------------------------------------
  # Populate reference choices and load reference assays
  #--------------------------------------------------------
  observe({
    req(ref_box_inserted_bulk())
    files_ref <- list.files(path = reference_dir, pattern = "\\.csv$", full.names = FALSE)
    stripped <- sub("\\.csv$", "", files_ref)

    updateSelectInput(
      session,
      "souris_reference_sel",
      choices  = stripped,
      selected = NULL
    )
  })

  observe({
    req(sc_ref_box_inserted())
    files_ref <- list.files(path = reference_dir, pattern = "\\.csv$", full.names = FALSE)
    stripped <- sub("\\.csv$", "", files_ref)

    updateSelectInput(
      session,
      "souris_reference_sel_sc",
      choices  = stripped,
      selected = NULL
    )
  })

  observeEvent(input$souris_load_reference, {
    n_sel <- length(input$souris_reference_sel)
    if (is.null(n_sel) || n_sel != 2) {
      showModal(modalDialog(
        title = "Selection required",
        "Please select exactly two reference assays before proceeding.",
        easyClose = TRUE
      ))
      return(NULL)
    }

    set_souris_progress("Loading reference assays", "", pct = 0, show = TRUE)

    sel_names <- input$souris_reference_sel
    paths <- file.path(reference_dir, paste0(sel_names, ".csv"))

    sizes <- file.info(paths)$size
    sizes[is.na(sizes) | sizes == 0] <- 1
    total_Size <- sum(sizes)
    cum_Size <- 0

    objs <- vector("list", length(paths))
    for (i in seq_along(paths)) {
      p <- paths[i]
      nm <- sel_names[i]

      set_souris_progress(
        subtext = sprintf("Loading reference assay %d/%d: %s", i, length(paths), nm),
        pct = 100 * cum_Size / total_Size, show = TRUE
      )

      df <- read.csv(p, stringsAsFactors = FALSE)
      objs[[i]] <- df
      cum_Size <- cum_Size + sizes[i]
    }

    names(objs) <- sel_names
    souris_reference_list(objs)

    set_souris_progress(
      text    = "Running subsystem models",
      subtext = "Computing fingerprint predictions across reference assays",
      show    = TRUE
    )

    on.exit(set_souris_progress(show = FALSE), add = TRUE)

    validate(
      need(
        length(input$souris_reference_sel) > 0,
        "Select at least one reference assay to view predictions."
      ),
      need(
        length(souris_fingerprints_list()) >= 1,
        "Fingerprint object not loaded."
      )
    )

    fp_obj <- souris_fingerprints_list()[[1]]

    validate(
      need(!is.null(fp_obj$All), "Fingerprint object has no 'All' component."),
      need(
        !is.null(fp_obj$All$fingerprints_primary),
        "No primary fingerprint classifiers found."
      )
    )

    fingerprint <- fp_obj$All
    final_classifiers <- fingerprint$fingerprints_primary
    auc_list <- fingerprint$auc_primary

    rownames(auc_list) <- auc_list$index
    index <- auc_list$index
    system_names <- auc_list[index, ]$subsystem

    ref_ids <- input$souris_reference_sel
    all_violin <- list()
    all_scores <- list()

    normalize01 <- function(v) {
      rng <- range(v, na.rm = TRUE)
      if (diff(rng) == 0 || !is.finite(diff(rng))) {
        return(rep(0.5, length(v)))
      }
      (v - rng[1]) / (rng[2] - rng[1])
    }

    total_Size <- sum(sizes)
    cum_Size <- 0

    for (ref_id in ref_ids) {
      ref_assay <- souris_reference_list()[[ref_id]]
      validate(need(
        !is.null(ref_assay),
        paste("Reference assay", ref_id, "not found.")
      ))

      ex_subsystems_no_norm <- vector("list", length(index))
      names(ex_subsystems_no_norm) <- system_names

      for (x in seq_along(index)) {
        set_souris_progress(
          subtext = sprintf(
            "Subsystem %d/%d: %s",
            x, length(system_names), system_names[x]
          ),
          pct = 100 * cum_Size / (length(index) * 2), show = TRUE
        )
        cum_Size <- cum_Size + 1
        idx <- as.numeric(index[x])
        if (is.null(final_classifiers[[idx]])) next

        coef <- final_classifiers[[idx]]$coefnames
        feats <- intersect(coef, colnames(ref_assay))
        if (!length(feats)) next

        newDat <- as.matrix(ref_assay[, feats, drop = FALSE])
        colnames(newDat) <- feats

        keep_cols <- colSums(!is.na(newDat)) > 0
        newDat <- newDat[, keep_cols, drop = FALSE]
        if (!ncol(newDat)) next # nothing usable

        keep_rows <- rowSums(is.na(newDat)) == 0
        newDat <- newDat[keep_rows, , drop = FALSE]
        if (!nrow(newDat)) next

        pred_tot <- predict(
          final_classifiers[[idx]],
          newdata   = newDat,
          type      = "prob",
          na.action = na.pass
        )[, 1]

        ex_subsystems_no_norm[[x]] <- pred_tot
      }

      ex_subsystems_no_norm <- ex_subsystems_no_norm[!vapply(ex_subsystems_no_norm, is.null, logical(1))]

      if (!length(ex_subsystems_no_norm)) next # skip empty ref

      df_ref <- lapply(names(ex_subsystems_no_norm), function(subsys) {
        data.frame(
          Reference = ref_id,
          Subsystem = subsys,
          Prediction = ex_subsystems_no_norm[[subsys]],
          stringsAsFactors = FALSE
        )
      })
      all_violin[[ref_id]] <- dplyr::bind_rows(df_ref)

      mat <- data.frame(ex_subsystems_no_norm, check.names = FALSE)
      mat <- as.data.frame(lapply(mat, normalize01))

      score <- rowMeans(mat, na.rm = TRUE)
      score <- (score - min(score, na.rm = TRUE)) / (max(score, na.rm = TRUE) - min(score, na.rm = TRUE))

      df_score <- data.frame(
        Sample = seq_along(score),
        Reference = ref_id,
        Score = score,
        stringsAsFactors = FALSE
      )
      all_scores[[ref_id]] <- df_score
    }

    if (!length(all_violin)) {
      validate(need(FALSE, "No valid subsystem predictions could be computed for the selected reference assays."))
    }

    df_all <- dplyr::bind_rows(all_violin)
    df_scores <- dplyr::bind_rows(all_scores)

    df_all_predictions(df_all)
    df_all_scores(df_scores)

    # Example: insert a placeholder analysis panel once
    if (!analysis_panel_inserted()) {
      insertUI(
        selector = "#bulk_analysis_placeholder",
        where = "afterEnd",
        ui = tagList(
          fluidRow(
            class = "fade-in-up",
            column(
              width = 6,
              box(
                width = 12, status = "primary", solidHeader = TRUE,
                title = NULL,
                div(
                  style = "display:flex; flex-direction:column; height:800px;",
                  # Sticky header
                  div(
                    id = "souris-fp-header",
                    style = "
        position: sticky;
        top: 0;
        z-index: 5;
        padding: 6px 10px;
        background: linear-gradient(135deg, #002B5C, #0D47A1);
        color: #E3F2FD;
        display: flex;
        align-items: center;
        justify-content: space-between;
        gap: 8px;
        box-shadow: 0 2px 6px rgba(0,0,0,0.15);
      ",
                    tags$div(
                      style = "display:flex; flex-direction:column; gap:2px;",
                      tags$div(
                        style = "font-size:12px; font-weight:600; letter-spacing:0.06em; text-transform:uppercase;",
                        "Total fingerprint scores"
                      ),
                      tags$div(
                        id = "souris-fp-subtitle",
                        style = "font-size:11px; color:#BBDEFB;",
                        "Per reference assay · Higher scores indicate stronger fingerprint signal"
                      )
                    ),
                    tags$div(
                      style = "display:flex; align-items:center; gap:4px; font-size:11px;",
                      tags$div(
                        style = "width:6px; height:6px; border-radius:999px; background:#BBDEFB; box-shadow:0 0 8px rgba(187,222,251,0.9);"
                      ),
                      tags$div("Interactive")
                    )
                  ),
                  # Scrollable body
                  div(
                    class = "souris-plot-fade",
                    id = "souris-violin-fp-wrapper",
                    style = "flex:1 1 auto; overflow-y:auto; padding:4px 0 0 0;",
                    plotly::plotlyOutput("souris_violin_fp", height = "100%")
                  )
                )
              )
            ),
            column(
              width = 6,
              box(
                width = 12, status = "primary", solidHeader = TRUE,
                title = NULL,
                div(
                  style = "display:flex; flex-direction:column; height:800px;",
                  # Sticky header
                  div(
                    id = "souris-auc-header",
                    style = "
        position: sticky;
        top: 0;
        z-index: 5;
        padding: 6px 10px;
        background: linear-gradient(135deg, #002B5C, #0D47A1);
        color: #E3F2FD;
        display: flex;
        align-items: center;
        justify-content: space-between;
        gap: 8px;
        box-shadow: 0 2px 6px rgba(0,0,0,0.15);
      ",
                    tags$div(
                      style = "display:flex; flex-direction:column; gap:2px;",
                      tags$div(
                        style = "font-size:12px; font-weight:600; letter-spacing:0.06em; text-transform:uppercase;",
                        "Subsystem AUC"
                      ),
                      tags$div(
                        id = "souris-auc-subtitle",
                        style = "font-size:11px; color:#BBDEFB;",
                        # filled/update from server via JS
                        "Select two assays · Click a bar for reaction-level details"
                      )
                    ),
                    tags$div(
                      style = "display:flex; align-items:center; gap:4px; font-size:11px;",
                      tags$div(
                        style = "width:6px; height:6px; border-radius:999px; background:#BBDEFB; box-shadow:0 0 8px rgba(187,222,251,0.9);"
                      ),
                      tags$div("Interactive")
                    )
                  ),
                  # Scrollable body
                  div(
                    class = "souris-plot-fade",
                    id = "souris-auc-bar-wrapper",
                    style = "flex:1 1 auto; overflow-y:auto; padding:4px 0 0 0;",
                    plotly::plotlyOutput("souris_auc_bar", height = "100%")
                  )
                )
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
    req(df_all_scores())
    df_scores <- df_all_scores()

    validate(need(
      nrow(df_scores) > 0,
      "No aggregated fingerprint scores could be computed for the selected reference assays."
    ))

    # Set factor order so the first level is the first selected reference assay
    df_scores$Reference <- factor(df_scores$Reference,
      levels = unique(df_scores$Reference)
    )

    ref_levels <- levels(df_scores$Reference)
    baseline <- ref_levels[1]


    p_to_stars <- function(p) {
      if (is.na(p)) {
        return("ns")
      }
      if (p < 0.0001) {
        return("****")
      }
      if (p < 0.001) {
        return("***")
      }
      if (p < 0.01) {
        return("**")
      }
      if (p < 0.05) {
        return("*")
      }
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

    # if we have at least one comparison, compute y positions and add stars
    y_max <- max(df_scores$Score, na.rm = TRUE)
    stat_df$y_pos <- seq(y_max * 1.05, length.out = nrow(stat_df), by = y_max * 0.05)

    ref1 <- ref_levels[1]
    ref2 <- ref_levels[2]

    auc_total <- auc_two_groups(
      pred = df_scores$Score,
      grp  = df_scores$Reference
    )
    auc_label <- sprintf("AUC = %.3f", auc_total)

    p <- ggplot2::ggplot(
      df_scores,
      ggplot2::aes(x = Reference, y = Score, fill = Reference)
    ) +
      ggplot2::geom_violin(
        alpha  = 0.5,
        scale  = "width",
        trim   = FALSE,
        color  = NA
      ) +
      ggplot2::geom_boxplot(
        width = 0.14,
        outlier.size = 0.6,
        outlier.alpha = 0.5,
        fill = "white",
        color = "#455A64",
        alpha = 0.9,
        size = 0.3
      ) +
      ggplot2::geom_jitter(
        width  = 0.06,
        size   = 0.8,
        alpha  = 0.35,
        color  = "#263238"
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
        y = "Total fingerprint score",
        title = auc_label
      ) +
      ggplot2::scale_fill_manual(
        values = c(
          "#1565C0", # ref 1
          "#42A5F5", # ref 2
          "#90CAF9" # extra levels, if any
        )
      ) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(
        # transparent / very light background to blend with box
        panel.background = ggplot2::element_rect(fill = NA, colour = NA),
        plot.background = ggplot2::element_rect(fill = NA, colour = NA),
        panel.grid.major.y = ggplot2::element_line(color = "#E0E4F0", linewidth = 0.25),
        panel.grid.major.x = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(
          angle = 45, hjust = 1, vjust = 1,
          color = "#283046",
          size = 9
        ),
        axis.text.y = ggplot2::element_text(color = "#283046", size = 9),
        axis.title.x = ggplot2::element_text(color = "#5f6473", size = 10, margin = ggplot2::margin(t = 6)),
        axis.title.y = ggplot2::element_text(color = "#5f6473", size = 10, margin = ggplot2::margin(r = 6)),
        plot.title = ggplot2::element_text(
          hjust = 0.5,
          face = "bold",
          color = "#283046",
          size = 11,
          margin = ggplot2::margin(b = 6)
        ),
        legend.position = "none"
      )
    session$sendCustomMessage("souris-show-plot", "souris-violin-fp-wrapper")
    plt <- plotly::ggplotly(p, tooltip = c("x", "y"))
    plt <- plotly::config(
      plt,
      toImageButtonOptions = list(
        format   = "svg",   # this makes the camera button download SVG
        filename = "souris_violin_fp"
      )
    )
    
    plt
  })

  output$souris_auc_bar <- plotly::renderPlotly({
    req(df_all_predictions())
    df_all <- df_all_predictions()

    refs <- unique(df_all$Reference)
    validate(need(
      length(refs) >= 2,
      "Select at least two reference assays to compute subsystem AUC."
    ))

    ref1 <- refs[1]
    ref2 <- refs[2]

    df_pair <- df_all[df_all$Reference %in% c(ref1, ref2), ]
    df_pair$Reference <- factor(df_pair$Reference, levels = c(ref1, ref2))

    auc_df <- df_pair |>
      dplyr::group_by(Subsystem) |>
      dplyr::summarise(
        AUC = auc_two_groups(Prediction, Reference),
        .groups = "drop"
      ) |>
      dplyr::filter(!is.na(AUC)) |>
      dplyr::arrange(AUC)

    validate(need(
      nrow(auc_df) > 0,
      "No valid AUC values could be computed for subsystem models."
    ))


    session$sendCustomMessage("souris-show-plot", "souris-auc-bar-wrapper")
    plt<-plotly::plot_ly(
      data = auc_df,
      x = ~AUC,
      y = ~ factor(Subsystem, levels = Subsystem),
      type = "bar",
      orientation = "h",
      source = "souris_auc",
      customdata = ~Subsystem,
      hoverinfo = "text",
      text = ~ paste0(
        "Subsystem: ", Subsystem,
        "<br>AUC: ", round(AUC, 3)
      ),
      marker = list(
        color = "rgba(21,101,192,0.85)",
        line  = list(color = "rgba(21,101,192,1)", width = 1)
      )
    ) |>
      plotly::layout(
        xaxis = list(
          title = "AUC",
          zeroline = FALSE,
          gridcolor = "#E0E4F0",
          tickfont = list(size = 10, color = "#283046"),
          titlefont = list(size = 11, color = "#5f6473")
        ),
        yaxis = list(
          title = "Subsystem",
          automargin = TRUE,
          tickfont = list(size = 10, color = "#283046"),
          titlefont = list(size = 11, color = "#5f6473"),
          categoryorder = "array",
          categoryarray = auc_df$Subsystem
        ),
        title = list(text = NULL),
        shapes = list(
          list(
            type = "line",
            x0 = 0.5, x1 = 0.5,
            y0 = -0.5, y1 = nrow(auc_df) + 0.5,
            line = list(color = "#EF5350", dash = "dash", width = 1)
          )
        ),
        bargap = 0.35,
        plot_bgcolor = "rgba(0,0,0,0)",
        paper_bgcolor = "rgba(0,0,0,0)",
        hovermode = "closest"
      )
    plt <- plotly::config(
      plt,
      toImageButtonOptions = list(
        format   = "svg",
        filename = "souris_auc_bar"
      )
    )
    
    plt
  })

  observeEvent(input$souris_load_reference_sc, {
    nsel <- length(input$souris_reference_sel_sc)
    if (is.null(nsel) || nsel != 2) {
      showModal(
        modalDialog(
          title = "Selection required",
          "Please select exactly two reference assays before proceeding.",
          easyClose = TRUE
        )
      )
      return(NULL)
    }

    set_souris_progress("Loading reference assays", NULL, pct = 0, show = TRUE)

    sel_names <- input$souris_reference_sel_sc
    paths <- file.path(reference_dir, paste0(sel_names, ".csv"))
    sizes <- file.info(paths)$size
    sizes[is.na(sizes) | sizes <= 0] <- 1
    total_Size <- sum(sizes)
    cum_Size <- 0

    objs <- vector("list", length(paths))

    for (i in seq_along(paths)) {
      p <- paths[i]
      nm <- sel_names[i]

      set_souris_progress(
        subtext = sprintf(
          "Loading reference assay %d/%d: %s",
          i, length(paths), nm
        ),
        pct = 100 * cum_Size / total_Size, show = TRUE
      )

      df <- read.csv(p, stringsAsFactors = FALSE)
      objs[[i]] <- df
      cum_Size <- cum_Size + sizes[i]
    }

    names(objs) <- sel_names
    souris_reference_list(objs)

    fpobj <- souris_fingerprints_list()[[1]]
    validate(
      need(!is.null(fpobj$All), "Fingerprint object has no All component."),
      need(
        !is.null(fpobj$All$fingerprints_primary),
        "No primary fingerprint classifiers found."
      )
    )

    set_souris_progress("Calculating fingerprint scores", NULL, pct = 0, show = TRUE)

    fingerprint <- fpobj$All
    finalclassifiers <- fingerprint$fingerprints_primary
    auclist <- fingerprint$auc_primary
    rownames(auclist) <- auclist$index
    index <- auclist$index
    systemnames <- auclist$subsystem
    refids <- sel_names

    allviolin <- list()
    allscores <- list()
    
    cum_Size <- 0
    
    for (refid in refids) {
      refassay <- souris_reference_list()[[refid]]
      validate(
        need(
          !is.null(refassay),
          paste("Reference assay", refid, "not found.")
        )
      )

      exsubsystems_nonorm <- vector("list", length(index))
      names(exsubsystems_nonorm) <- systemnames

      i <- 1
      for (x in seq_along(index)) {
        idx <- as.numeric(index[x])

        set_souris_progress(
          subtext = sprintf(
            "Subsystem %d/%d: %s",
            i, length(systemnames), systemnames[i]
          ),
          pct = 100 * cum_Size / (length(index)*2), show = TRUE
        )
        cum_Size <- cum_Size + 1
        
        if (is.null(finalclassifiers[[idx]])) next

        coef <- finalclassifiers[[idx]]$coefnames
        feats <- intersect(coef, colnames(refassay))
        if (!length(feats)) next

        newDat <- as.matrix(refassay[, feats, drop = FALSE])
        colnames(newDat) <- feats

        keepcols <- colSums(!is.na(newDat)) > 0
        newDat <- newDat[, keepcols, drop = FALSE]
        if (!ncol(newDat)) next

        keeprows <- rowSums(is.na(newDat)) == 0
        newDat <- newDat[keeprows, , drop = FALSE]
        if (!nrow(newDat)) next

        predtot <- predict(
          finalclassifiers[[idx]],
          newdata   = newDat,
          type      = "prob",
          na.action = na.pass
        )[, 1]
        exsubsystems_nonorm[[x]] <- predtot

        i <- i + 1
      }

      exsubsystems_nonorm <- exsubsystems_nonorm[
        !vapply(exsubsystems_nonorm, is.null, logical(1))
      ]
      if (!length(exsubsystems_nonorm)) next

      dfref <- lapply(
        names(exsubsystems_nonorm),
        function(subsys) {
          data.frame(
            Sample = souris_reference_list()[[refid]][["X"]],
            Reference = refid,
            Subsystem = subsys,
            Prediction = exsubsystems_nonorm[[subsys]],
            stringsAsFactors = FALSE
          )
        }
      )
      allviolin[[refid]] <- dplyr::bind_rows(dfref)

      mat <- data.frame(exsubsystems_nonorm, check.names = FALSE)

      score <- rowMeans(mat, na.rm = TRUE)

      dfscore <- data.frame(
        Sample = souris_reference_list()[[refid]][["X"]],
        Reference = refid,
        Score = score,
        stringsAsFactors = FALSE
      )
      allscores[[refid]] <- dfscore
    }

    dfall <- dplyr::bind_rows(allviolin)
    dfscores <- dplyr::bind_rows(allscores)

    dfscores$Score <- (dfscores$Score - min(dfscores$Score, na.rm = TRUE)) /
      (max(dfscores$Score, na.rm = TRUE) - min(dfscores$Score, na.rm = TRUE))

    df_all_predictions(dfall)
    df_all_scores(dfscores)

    ## Attach cell IDs now so they can match Seurat later
    dfscores$Cell <- unlist(
      lapply(souris_reference_list(), function(df) df[["X"]]),
      use.names = FALSE
    )

    sc_fingerprint_scores(dfscores)
    set_souris_progress(show = FALSE)

    ## After reference loading is complete, insert SC object selection panel
    if (!isTRUE(sc_scobj_box_inserted())) {
      insertUI(
        selector = "#sc_scobj_placeholder",
        where = "afterBegin",
        ui = fluidRow(
          class = "fade-in-up",
          column(
            width = 12,
            box(
              width = 12, status = "primary", solidHeader = TRUE,
              title = "Single-cell object selection",
              shinyWidgets::pickerInput(
                inputId = "souris_scobject_sel",
                label = "Select single-cell object",
                choices = NULL,
                multiple = FALSE,
                options = list(
                  `actions-box` = FALSE,
                  `none-selected-text` = "No single-cell object selected"
                )
              ),
              actionButton(
                inputId = "souris_load_scobject",
                label = "Load single-cell object"
              )
            )
          )
        ),
        immediate = TRUE
      )
      sc_scobj_box_inserted(TRUE)

      ## Populate SC object choices from .rds files in sc_dir
      files_sc <- list.files(path = sc_dir, pattern = "\\.rds$", full.names = FALSE)
      stripped_sc <- sub("\\.rds$", "", files_sc)

      shinyWidgets::updatePickerInput(
        session,
        inputId = "souris_scobject_sel",
        choices = stripped_sc,
        selected = NULL
      )
    }
  })

  observeEvent(input$souris_load_scobject, {
    req(input$souris_scobject_sel)
    selname <- input$souris_scobject_sel
    path <- file.path(sc_dir, paste0(selname, ".rds"))

    set_souris_progress("Loading single-cell object", NULL, pct = 0, show = TRUE)
    on.exit(set_souris_progress(show = FALSE), add = TRUE)

    obj_sc <- readRDS(path)

    refs <- souris_reference_list()
    validate(
      need(
        !is.null(refs) && length(refs) > 0,
        "Reference assays must be loaded before loading the single-cell object."
      )
    )
    refids <- unique(unlist(lapply(refs, function(df) df[[1]])))

    sc_meta_row_names <- rownames(obj_sc@meta.data)
    common_ids <- intersect(refids, sc_meta_row_names)

    if (length(common_ids) == 0) {
      showModal(
        modalDialog(
          title = "No matching samples",
          paste(
            "No IDs from the first column of the reference assays were found",
            "in the single-cell object metadata rownames. Please select a compatible object."
          ),
          easyClose = TRUE
        )
      )
      souris_sc_object(NULL)
      return(NULL)
    }

    obj_sc <- obj_sc[, sc_meta_row_names %in% common_ids]

    df_sc_scores <- sc_fingerprint_scores()
    if (is.null(df_sc_scores)) {
      showModal(
        modalDialog(
          title = "Scores not available",
          "Please run the single-cell fingerprint analysis (load SC reference assays) before loading the single-cell object.",
          easyClose = TRUE
        )
      )
      return(NULL)
    }

    df_sc_scores <- df_sc_scores[match(colnames(obj_sc), df_sc_scores$Cell), ]

    obj_sc$SOURIS_fingerprint_score <- df_sc_scores$Score
    obj_sc$SOURIS_reference_group <- df_sc_scores$Reference

    ## Just store the SC object in memory; no metadata changes yet
    souris_sc_object(setNames(list(obj_sc), selname))

    ## Insert SC analysis panel (violin + AUC + FeaturePlot shell) only now
    if (!isTRUE(sc_analysis_panel_inserted())) {
      insertUI(
        selector = "#sc_analysis_placeholder",
        where = "afterEnd",
        ui = tagList(
          fluidRow(
            class = "fade-in-up",
            column(
              width = 6,
              box(
                width = 12, status = "primary", solidHeader = TRUE,
                title = NULL,
                div(
                  style = "display:flex; flex-direction:column; height:800px;",
                  # Sticky header (can reuse same text or tweak)
                  div(
                    id = "souris-fp-sc-header",
                    style = "
          position: sticky;
          top: 0;
          z-index: 5;
          padding: 6px 10px;
          background: linear-gradient(135deg, #002B5C, #0D47A1);
          color: #E3F2FD;
          display: flex;
          align-items: center;
          justify-content: space-between;
          gap: 8px;
          box-shadow: 0 2px 6px rgba(0,0,0,0.15);
        ",
                    tags$div(
                      style = 'display:flex; flex-direction:column; gap:2px;',
                      tags$div(
                        style = 'font-size:12px; font-weight:600; letter-spacing:0.06em; text-transform:uppercase;',
                        "Total fingerprint scores (SC)"
                      ),
                      tags$div(
                        id = "souris-fp-sc-subtitle",
                        style = 'font-size:11px; color:#BBDEFB;',
                        "Per reference assay · Higher scores indicate stronger fingerprint signal"
                      )
                    ),
                    tags$div(
                      style = "display:flex; align-items:center; gap:4px; font-size:11px;",
                      tags$div(
                        style = "width:6px; height:6px; border-radius:999px; background:#BBDEFB; box-shadow:0 0 8px rgba(187,222,251,0.9);"
                      ),
                      tags$div("Interactive")
                    )
                  ),
                  # Scrollable body with fade-in plot
                  div(
                    class = "souris-plot-fade",
                    id    = "souris-violin-fp-sc-wrapper",
                    style = "flex:1 1 auto; overflow-y:auto; padding:4px 0 0 0;",
                    plotly::plotlyOutput("souris_violin_fp_sc", height = "100%")
                  )
                )
              )
            ),
            column(
              width = 6,
              box(
                width = 12, status = "primary", solidHeader = TRUE,
                title = NULL,
                div(
                  style = "display:flex; flex-direction:column; height:800px;",
                  # Sticky header
                  div(
                    id = "souris-auc-sc-header",
                    style = "
          position: sticky;
          top: 0;
          z-index: 5;
          padding: 6px 10px;
          background: linear-gradient(135deg, #002B5C, #0D47A1);
          color: #E3F2FD;
          display: flex;
          align-items: center;
          justify-content: space-between;
          gap: 8px;
          box-shadow: 0 2px 6px rgba(0,0,0,0.15);
        ",
                    tags$div(
                      style = 'display:flex; flex-direction:column; gap:2px;',
                      tags$div(
                        style = 'font-size:12px; font-weight:600; letter-spacing:0.06em; text-transform:uppercase;',
                        "Subsystem AUC (SC)"
                      ),
                      tags$div(
                        id = "souris-auc-sc-subtitle",
                        style = 'font-size:11px; color:#BBDEFB;',
                        "Select two assays · Click a bar to open reaction-level analyses"
                      )
                    ),
                    tags$div(
                      style = "display:flex; align-items:center; gap:4px; font-size:11px;",
                      tags$div(
                        style = "width:6px; height:6px; border-radius:999px; background:#BBDEFB; box-shadow:0 0 8px rgba(187,222,251,0.9);"
                      ),
                      tags$div("Interactive")
                    )
                  ),
                  # Scrollable body with fade-in plot
                  div(
                    class = "souris-plot-fade",
                    id    = "souris-auc-bar-sc-wrapper",
                    style = "flex:1 1 auto; overflow-y:auto; padding:4px 0 0 0;",
                    plotly::plotlyOutput("souris_auc_bar_sc", height = "100%")
                  )
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
                title = NULL,
                div(
                  style = "display:flex; flex-direction:column; height:550px;",
                  div(
                    id = "souris-sc-feature-fp-header",
                    style = "
            position: sticky;
            top: 0;
            z-index: 5;
            padding: 6px 10px;
            background: linear-gradient(135deg, #002B5C, #0D47A1);
            color: #E3F2FD;
            display: flex;
            align-items: center;
            justify-content: space-between;
            gap: 8px;
            box-shadow: 0 2px 6px rgba(0,0,0,0.15);
          ",
                    tags$div(
                      style = 'display:flex; flex-direction:column; gap:2px;',
                      tags$div(
                        style = 'font-size:12px; font-weight:600; letter-spacing:0.06em; text-transform:uppercase;',
                        "Single-cell fingerprint scores FeaturePlot"
                      ),
                      tags$div(
                        style = 'font-size:11px; color:#BBDEFB;',
                        "Total fingerprint scores per reference group · single-cell"
                      )
                    ),
                    div(
                      style = "display:flex; align-items:center; justify-content:flex-end;",
                      downloadButton(
                        "download_souris_sc_featureplot",
                        "Download PDF",
                        style = "
                background: rgba(227,242,253,0.12);
                border: 1px solid #BBDEFB;
                color: #E3F2FD;
                padding: 2px 10px;
                font-size: 10px;
                border-radius: 999px;
              "
                      )
                    )
                  ),
                  div(
                    style = "flex:1 1 auto; overflow-y:auto; padding:4px 0 0 0;",
                    plotOutput(
                      "souris_sc_featureplot",
                      height = "100%", width = "100%"
                    )
                  )
                )
              )
            )
          )
        ),
        immediate = TRUE
      )
      sc_analysis_panel_inserted(TRUE)
    }
  })


  output$souris_violin_fp_sc <- plotly::renderPlotly({
    req(df_all_scores())
    dfscores <- df_all_scores()
    
    validate(
      need(
        nrow(dfscores) > 0,
        "No aggregated fingerprint scores could be computed for the selected reference assays."
      )
    )
    
    # factor order like bulk
    dfscores$Reference <- factor(dfscores$Reference,
                                 levels = unique(dfscores$Reference))
    reflevels <- levels(dfscores$Reference)
    baseline  <- reflevels[1]
    
    # Wilcoxon vs baseline for stars
    p_to_stars <- function(p) {
      if (is.na(p)) return("ns")
      if (p < 0.0001) return("****")
      if (p < 0.001)  return("***")
      if (p < 0.01)   return("**")
      if (p < 0.05)   return("*")
      "ns"
    }
    
    tests <- lapply(reflevels[-1], function(ref) {
      x <- dfscores$Score[dfscores$Reference == baseline]
      y <- dfscores$Score[dfscores$Reference == ref]
      p <- stats::wilcox.test(x, y, exact = FALSE)$p.value
      data.frame(group1 = baseline, group2 = ref,
                 p = p, plab = p_to_stars(p))
    })
    statdf <- dplyr::bind_rows(tests)
    
    ymax <- max(dfscores$Score, na.rm = TRUE)
    if (nrow(statdf) > 0) {
      statdf$ypos <- seq(ymax * 1.05,
                         length.out = nrow(statdf),
                         by = ymax * 0.05)
    }
    
    validate(
      need(
        length(reflevels) == 2,
        "Select exactly two reference assays to compute AUC."
      )
    )
    ref1 <- reflevels[1]
    ref2 <- reflevels[2]
    
    auctotal <- auc_two_groups(
      pred = dfscores$Score,
      grp  = dfscores$Reference
    )
    auclabel <- sprintf("AUC = %.3f", auctotal)
    
    p <- ggplot2::ggplot(
      dfscores,
      ggplot2::aes(x = Reference, y = Score, fill = Reference)
    ) +
      ggplot2::geom_violin(
        alpha = 0.5,
        scale = "width",
        trim  = FALSE,
        color = NA
      ) +
      ggplot2::geom_boxplot(
        width         = 0.14,
        outlier.size  = 0.6,
        outlier.alpha = 0.5,
        fill          = "white",
        color         = "#455A64",
        alpha         = 0.9,
        size          = 0.3
      ) +
      ggplot2::geom_jitter(
        width = 0.06,
        size  = 0.8,
        alpha = 0.35,
        color = "#263238"
      ) +
      # significance segments + labels (if any)
      ggplot2::geom_segment(
        data        = statdf,
        ggplot2::aes(
          x    = as.numeric(group1),
          xend = as.numeric(group2),
          y    = ypos,
          yend = ypos
        ),
        inherit.aes = FALSE
      ) +
      ggplot2::geom_text(
        data        = statdf,
        ggplot2::aes(
          x     = (as.numeric(group1) + as.numeric(group2)) / 2,
          y     = ypos,
          label = plab
        ),
        vjust       = -0.3,
        inherit.aes = FALSE,
        size        = 3
      ) +
      ggplot2::scale_fill_manual(
        values = c(
          "#1565C0", # ref 1
          "#42A5F5", # ref 2
          "#90CAF9" # extra levels, if any
        )
      ) +
      ggplot2::labs(
        x     = "Reference assay",
        y     = "Total fingerprint score",
        title = auclabel
      ) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(
        panel.background   = ggplot2::element_rect(fill = NA, colour = NA),
        plot.background    = ggplot2::element_rect(fill = NA, colour = NA),
        panel.grid.major.y = ggplot2::element_line(color = "#E0E4F0", linewidth = 0.25),
        panel.grid.major.x = ggplot2::element_blank(),
        panel.grid.minor   = ggplot2::element_blank(),
        axis.text.x  = ggplot2::element_text(
          angle = 45, hjust = 1, vjust = 1,
          color = "#283046",
          size  = 9
        ),
        axis.text.y  = ggplot2::element_text(color = "#283046", size = 9),
        axis.title.x = ggplot2::element_text(color = "#5f6473", size = 10, margin = ggplot2::margin(t = 6)),
        axis.title.y = ggplot2::element_text(color = "#5f6473", size = 10, margin = ggplot2::margin(r = 6)),
        plot.title   = ggplot2::element_text(
          hjust  = 0.5,
          face   = "bold",
          color  = "#283046",
          size   = 11,
          margin = ggplot2::margin(b = 6)
        ),
        legend.position = "none"
      )
    
    session$sendCustomMessage("souris-show-plot", "souris-violin-fp-sc-wrapper")
    
    plt <- plotly::ggplotly(p, tooltip = c("x", "y"))
    plt <- plotly::config(
      plt,
      toImageButtonOptions = list(
        format   = "svg",   # this makes the camera button download SVG
        filename = "souris_violin_fp_sc"
      )
    )
    
    plt
  })




  output$souris_auc_bar_sc <- plotly::renderPlotly({
    df_all <- df_all_predictions()
    req(df_all)
    
    refs <- unique(df_all$Reference)
    validate(need(
      length(refs) >= 2,
      "Select at least two reference assays to compute subsystem AUC."
    ))
    
    ref1 <- refs[1]
    ref2 <- refs[2]
    
    df_pair <- df_all[df_all$Reference %in% c(ref1, ref2), ]
    df_pair$Reference <- factor(df_pair$Reference, levels = c(ref1, ref2))
    
    auc_df <- df_pair |>
      dplyr::group_by(Subsystem) |>
      dplyr::summarise(
        AUC = auc_two_groups(Prediction, Reference),
        .groups = "drop"
      ) |>
      dplyr::filter(!is.na(AUC)) |>
      dplyr::arrange(AUC)
    
    validate(need(
      nrow(auc_df) > 0,
      "No valid AUC values could be computed for subsystem models."
    ))
    
    pauc <- ggplot2::ggplot(
      auc_df,
      ggplot2::aes(
        x         = reorder(Subsystem, AUC),
        y         = AUC,
        customdata = Subsystem
      )
    ) +
      ggplot2::geom_col(
        fill = "rgba(21,101,192,0.85)"
      ) +
      ggplot2::coord_flip() +
      ggplot2::geom_hline(
        yintercept = 0.5,
        linetype   = "dashed",
        color      = "#EF5350",
        linewidth  = 0.6
      ) +
      ggplot2::labs(
        x     = "Subsystem",
        y     = "AUC (first vs second reference)"
      ) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(
        panel.background   = ggplot2::element_rect(fill = NA, colour = NA),
        plot.background    = ggplot2::element_rect(fill = NA, colour = NA),
        panel.grid.major.y = ggplot2::element_line(color = "#E0E4F0", linewidth = 0.25),
        panel.grid.major.x = ggplot2::element_blank(),
        panel.grid.minor   = ggplot2::element_blank(),
        axis.text.x  = ggplot2::element_text(color = "#283046", size = 9),
        axis.text.y  = ggplot2::element_text(color = "#283046", size = 9),
        axis.title.x = ggplot2::element_text(color = "#5f6473", size = 10, margin = ggplot2::margin(t = 6)),
        axis.title.y = ggplot2::element_text(color = "#5f6473", size = 10, margin = ggplot2::margin(r = 6)),
        plot.title   = ggplot2::element_text(
          hjust  = 0.5,
          face   = "bold",
          color  = "#283046",
          size   = 11,
          margin = ggplot2::margin(b = 6)
        ),
        legend.position = "none"
      )
    
    session$sendCustomMessage("souris-show-plot", "souris-auc-bar-sc-wrapper")
    
    plt <- plotly::ggplotly(pauc, tooltip = c("x", "y"), source = "souris_auc_sc")
    plt <- plotly::config(
      plt,
      toImageButtonOptions = list(
        format   = "svg",   # this makes the camera button download SVG
        filename = "souris_auc_bar_sc"
      )
    )
    
    plt
  })

  observeEvent(
    plotly::event_data("plotly_click", source = "souris_auc_sc"),
    {
      d <- plotly::event_data("plotly_click", source = "souris_auc_sc")
      if (is.null(d)) {
        return()
      }

      subsel <- as.character(d$customdata) # subsystem name from SC AUC bar
      clicked_subsystem_sc(subsel)

      subsystem_done(FALSE)
      reactions_done(FALSE)

      set_souris_progress(
        text    = "Analyzing subsystem reactions (single-cell)",
        subtext = sprintf("Computing reaction-level metrics for '%s'", subsel),
        show    = TRUE
      )

      if (!sc_detail_panel_inserted()) {
        insertUI(
          selector = "#sc_tier3_placeholder", # add this div in SC UI
          where = "afterEnd",
          ui = tagList(
            fluidRow(
              class = "fade-in-up",
              column(
                width = 6,
                box(
                  width = 12, status = "primary", solidHeader = TRUE,
                  title = NULL,
                  div(
                    style = "display:flex; flex-direction:column; height:800px;",
                    # Sticky header
                    div(
                      id = "souris-subsys-sc-header",
                      style = "
        position: sticky;
        top: 0;
        z-index: 5;
        padding: 6px 10px;
        background: linear-gradient(135deg, #002B5C, #0D47A1);
        color: #E3F2FD;
        display: flex;
        align-items: center;
        justify-content: space-between;
        gap: 8px;
        box-shadow: 0 2px 6px rgba(0,0,0,0.15);
      ",
                      tags$div(
                        style = 'display:flex; flex-direction:column; gap:2px;',
                        tags$div(
                          style = 'font-size:12px; font-weight:600; letter-spacing:0.06em; text-transform:uppercase;',
                          "Subsystem prediction (SC)"
                        ),
                        tags$div(
                          id = "souris-subsys-sc-subtitle",
                          style = 'font-size:11px; color:#BBDEFB;',
                          "Per reference assay · Higher scores indicate stronger subsystem model signal"
                        )
                      ),
                      tags$div(
                        style = "display:flex; align-items:center; gap:4px; font-size:11px;",
                        tags$div(
                          style = "width:6px; height:6px; border-radius:999px; background:#BBDEFB; box-shadow:0 0 8px rgba(187,222,251,0.9);"
                        ),
                        tags$div("Interactive")
                      )
                    ),
                    # Scrollable body with fade-in plot
                    div(
                      class = "souris-plot-fade",
                      id    = "souris-violin-subsys-sc-wrapper",
                      style = "flex:1 1 auto; overflow-y:auto; padding:4px 0 0 0;",
                      plotly::plotlyOutput("souris_violin_subsys_sc", height = "100%")
                    )
                  )
                )
              ),
              column(
                width = 6,
                box(
                  width = 12, status = "primary", solidHeader = TRUE,
                  title = NULL,
                  div(
                    style = "display:flex; flex-direction:column; height:800px;",
                    # Sticky header
                    div(
                      id = "souris-auc-reactions-sc-header",
                      style = "
        position: sticky;
        top: 0;
        z-index: 5;
        padding: 6px 10px;
        background: linear-gradient(135deg, #002B5C, #0D47A1);
        color: #E3F2FD;
        display: flex;
        align-items: center;
        justify-content: space-between;
        gap: 8px;
        box-shadow: 0 2px 6px rgba(0,0,0,0.15);
      ",
                      tags$div(
                        style = 'display:flex; flex-direction:column; gap:2px;',
                        tags$div(
                          style = 'font-size:12px; font-weight:600; letter-spacing:0.06em; text-transform:uppercase;',
                          "Reaction AUC (SC)"
                        ),
                        tags$div(
                          id = "souris-auc-reactions-sc-subtitle",
                          style = 'font-size:11px; color:#BBDEFB;',
                          "Per reaction · Click reaction bar to view Seurat FeaturePlot"
                        )
                      ),
                      tags$div(
                        style = "display:flex; align-items:center; gap:4px; font-size:11px;",
                        tags$div(
                          style = "width:6px; height:6px; border-radius:999px; background:#BBDEFB; box-shadow:0 0 8px rgba(187,222,251,0.9);"
                        ),
                        tags$div("Interactive")
                      )
                    ),
                    # Scrollable body with fade-in plot
                    div(
                      class = "souris-plot-fade",
                      id    = "souris-auc-reactions-sc-wrapper",
                      style = "flex:1 1 auto; overflow-y:auto; padding:4px 0 0 0;",
                      plotly::plotlyOutput("souris_auc_reactions_sc", height = "100%")
                    )
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
                  title = NULL,
                  div(
                    style = "display:flex; flex-direction:column; height:500px;",
                    div(
                      id = "souris-violin-reactions-sc-header",
                      style = "
        position: sticky;
        top: 0;
        z-index: 5;
        padding: 6px 10px;
        background: linear-gradient(135deg, #002B5C, #0D47A1);
        color: #E3F2FD;
        display: flex;
        align-items: center;
        justify-content: space-between;
        gap: 8px;
        box-shadow: 0 2px 6px rgba(0,0,0,0.15);
      ",
                      tags$div(
                        style = 'display:flex; flex-direction:column; gap:2px;',
                        tags$div(
                          style = 'font-size:12px; font-weight:600; letter-spacing:0.06em; text-transform:uppercase;',
                          "Reaction distributions (SC)"
                        ),
                        tags$div(
                          id = "souris-violin-reactions-sc-subtitle",
                          style = 'font-size:11px; color:#BBDEFB;',
                          "Per reaction · violins and boxplots per reference assay"
                        )
                      ),
                      tags$div(
                        style = "display:flex; align-items:center; gap:4px; font-size:11px;",
                        tags$div(
                          style = "width:6px; height:6px; border-radius:999px; background:#BBDEFB; box-shadow:0 0 8px rgba(187,222,251,0.9);"
                        ),
                        tags$div("Interactive")
                      )
                    ),
                    div(
                      class = "souris-plot-fade",
                      id    = "souris-violin-reactions-sc-wrapper",
                      style = "flex:1 1 auto; overflow-y:auto; padding:4px 0 0 0;",
                      plotly::plotlyOutput("souris_violin_reactions_sc", height = "100%")
                    )
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
                  title = NULL,
                  div(
                    style = "display:flex; flex-direction:column; height:500px;",
                    div(
                      id = "souris-sc-feature-subsystem-header",
                      style = "
            position: sticky;
            top: 0;
            z-index: 5;
            padding: 6px 10px;
            background: linear-gradient(135deg, #002B5C, #0D47A1);
            color: #E3F2FD;
            display: flex;
            align-items: center;
            justify-content: space-between;
            gap: 8px;
            box-shadow: 0 2px 6px rgba(0,0,0,0.15);
          ",
                      tags$div(
                        style = 'display:flex; flex-direction:column; gap:2px;',
                        tags$div(
                          style = 'font-size:12px; font-weight:600; letter-spacing:0.06em; text-transform:uppercase;',
                          uiOutput("souris_sc_featureplot_subsystem_title")
                        ),
                        tags$div(
                          style = 'font-size:11px; color:#BBDEFB;',
                          "Subsystem scores per reference group · single-cell"
                        )
                      )
                    ),
                    div(
                      style = "flex:1 1 auto; overflow-y:auto; padding:4px 0 0 0;",
                      plotOutput(
                        "souris_sc_featureplot_subsystem",
                        height = "100%", width = "100%"
                      )
                    )
                  )
                )
              ),
              div(id = "sc_reaction_featureplot_placeholder")
            )
          ),
          immediate = TRUE
        )
        sc_detail_panel_inserted(TRUE)
      }

      set_souris_progress(subtext = "Rendering subsystem detail plots (single-cell)")
    }
  )


  observeEvent(plotly::event_data("plotly_click", source = "souris_auc"), {
    d <- plotly::event_data("plotly_click", source = "souris_auc")

    if (is.null(d)) {
      return()
    }

    subsel <- as.character(d$customdata) # subsystem name from bar label
    clicked_subsystem(subsel)

    subsystem_done(FALSE)
    reactions_done(FALSE)

    set_souris_progress(
      text    = "Analyzing subsystem reactions",
      subtext = sprintf("Computing reaction-level metrics for '%s'", subsel),
      show    = TRUE
    )

    if (!detail_panel_inserted()) {
      insertUI(
        selector = "#bulk_tier3_placeholder",
        where = "afterEnd",
        ui = tagList(
          fluidRow(
            class = "fade-in-up",
            column(
              width = 6,
              box(
                width = 12, status = "primary", solidHeader = TRUE,
                title = NULL,
                div(
                  style = "display:flex; flex-direction:column; height:800px;",
                  # Sticky header
                  div(
                    id = "souris-subsys-header",
                    style = "
        position: sticky;
        top: 0;
        z-index: 5;
        padding: 6px 10px;
        background: linear-gradient(135deg, #002B5C, #0D47A1);
        color: #E3F2FD;
        display: flex;
        align-items: center;
        justify-content: space-between;
        gap: 8px;
        box-shadow: 0 2px 6px rgba(0,0,0,0.15);
      ",
                    tags$div(
                      style = "display:flex; flex-direction:column; gap:2px;",
                      tags$div(
                        style = "font-size:12px; font-weight:600; letter-spacing:0.06em; text-transform:uppercase;",
                        "Subsystem prediction"
                      ),
                      tags$div(
                        id = "souris-subsys-subtitle",
                        style = "font-size:11px; color:#BBDEFB;",
                        "Per reference assay · Higher scores indicate stronger subsystem model signal"
                      )
                    ),
                    tags$div(
                      style = "display:flex; align-items:center; gap:4px; font-size:11px;",
                      tags$div(
                        style = "width:6px; height:6px; border-radius:999px; background:#BBDEFB; box-shadow:0 0 8px rgba(187,222,251,0.9);"
                      ),
                      tags$div("Interactive")
                    )
                  ),
                  # Scrollable body with fade-in plot
                  div(
                    class = "souris-plot-fade",
                    id = "souris-violin-subsys-wrapper",
                    style = "flex:1 1 auto; overflow-y:auto; padding:4px 0 0 0;",
                    plotly::plotlyOutput("souris_violin_subsys", height = "100%")
                  )
                )
              )
            ),
            column(
              width = 6,
              box(
                width = 12, status = "primary", solidHeader = TRUE,
                title = NULL,
                div(
                  style = "display:flex; flex-direction:column; height:800px;",
                  # Sticky header
                  div(
                    id = "souris-auc-reactions-header",
                    style = "
        position: sticky;
        top: 0;
        z-index: 5;
        padding: 6px 10px;
        background: linear-gradient(135deg, #002B5C, #0D47A1);
        color: #E3F2FD;
        display: flex;
        align-items: center;
        justify-content: space-between;
        gap: 8px;
        box-shadow: 0 2px 6px rgba(0,0,0,0.15);
      ",
                    tags$div(
                      style = "display:flex; flex-direction:column; gap:2px;",
                      tags$div(
                        style = "font-size:12px; font-weight:600; letter-spacing:0.06em; text-transform:uppercase;",
                        "Reaction AUC"
                      ),
                      tags$div(
                        id = "souris-auc-reactions-subtitle",
                        style = "font-size:11px; color:#BBDEFB;",
                        "Per reaction · AUC between the two references"
                      )
                    ),
                    tags$div(
                      style = "display:flex; align-items:center; gap:4px; font-size:11px;",
                      tags$div(
                        style = "width:6px; height:6px; border-radius:999px; background:#BBDEFB; box-shadow:0 0 8px rgba(187,222,251,0.9);"
                      ),
                      tags$div("Interactive")
                    )
                  ),
                  # Scrollable body with fade-in plot
                  div(
                    class = "souris-plot-fade",
                    id = "souris-auc-reactions-wrapper",
                    style = "flex:1 1 auto; overflow-y:auto; padding:4px 0 0 0;",
                    plotly::plotlyOutput("souris_auc_reactions", height = "100%")
                  )
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
                title = NULL,
                div(
                  style = "display:flex; flex-direction:column; height:500px;",
                  div(
                    id = "souris-violin-reactions-header",
                    style = "
        position: sticky;
        top: 0;
        z-index: 5;
        padding: 6px 10px;
        background: linear-gradient(135deg, #002B5C, #0D47A1);
        color: #E3F2FD;
        display: flex;
        align-items: center;
        justify-content: space-between;
        gap: 8px;
        box-shadow: 0 2px 6px rgba(0,0,0,0.15);
      ",
                    tags$div(
                      style = 'display:flex; flex-direction:column; gap:2px;',
                      tags$div(
                        style = 'font-size:12px; font-weight:600; letter-spacing:0.06em; text-transform:uppercase;',
                        "Reaction distributions"
                      ),
                      tags$div(
                        id = "souris-violin-reactions-subtitle",
                        style = 'font-size:11px; color:#BBDEFB;',
                        "Per reaction · boxplots per reference assay"
                      )
                    ),
                    tags$div(
                      style = "display:flex; align-items:center; gap:4px; font-size:11px;",
                      tags$div(
                        style = "width:6px; height:6px; border-radius:999px; background:#BBDEFB; box-shadow:0 0 8px rgba(187,222,251,0.9);"
                      ),
                      tags$div("Interactive")
                    )
                  ),
                  div(
                    class = "souris-plot-fade",
                    id    = "souris-violin-reactions-wrapper",
                    style = "flex:1 1 auto; overflow-y:auto; padding:4px 0 0 0;",
                    plotly::plotlyOutput("souris_violin_reactions", height = "100%")
                  )
                )
              )
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

    validate(need(
      nrow(df_sub) > 0,
      "No predictions found for the selected subsystem."
    ))

    # scale prediction 0–1 like the main violin
    rng <- range(df_sub$Prediction, na.rm = TRUE)
    if (diff(rng) == 0 || !is.finite(diff(rng))) {
      df_sub$Prediction_scaled <- 0.5
    } else {
      df_sub$Prediction_scaled <- (df_sub$Prediction - rng[1]) / (rng[2] - rng[1])
    }

    # Make sure Reference is a factor with the current order
    df_sub$Reference <- factor(df_sub$Reference, levels = unique(df_sub$Reference))
    ref_levels <- levels(df_sub$Reference)
    validate(need(
      length(ref_levels) == 2,
      "Subsystem AUC is defined here only for exactly two reference assays."
    ))

    ref1 <- ref_levels[1]
    ref2 <- ref_levels[2]

    # AUC of subsystem prediction between the two references
    auc_sub <- auc_two_groups(
      pred = df_sub$Prediction,
      grp  = df_sub$Reference
    )
    auc_label <- sprintf("AUC = %.3f", auc_sub)

    p_sub <- ggplot2::ggplot(
      df_sub,
      ggplot2::aes(x = Reference, y = Prediction_scaled, fill = Reference)
    ) +
      ggplot2::geom_violin(
        alpha = 0.5,
        scale = "width",
        trim  = FALSE,
        color = NA
      ) +
      ggplot2::geom_boxplot(
        width         = 0.14,
        outlier.size  = 0.6,
        outlier.alpha = 0.5,
        fill          = "white",
        color         = "#455A64",
        alpha         = 0.9,
        size          = 0.3
      ) +
      ggplot2::geom_jitter(
        width = 0.06,
        size  = 0.8,
        alpha = 0.35,
        color = "#263238"
      ) +
      ggplot2::labs(
        x     = "Reference assay",
        y     = "Subsystem prediction probability",
        title = paste("Subsystem:", subsel, "|", auc_label)
      ) +
      ggplot2::scale_fill_manual(
        values = c(
          "#1565C0", # ref 1
          "#42A5F5" # ref 2 (extra levels, if any will recycle)
        )
      ) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(
        panel.background = ggplot2::element_rect(fill = NA, colour = NA),
        plot.background = ggplot2::element_rect(fill = NA, colour = NA),
        panel.grid.major.y = ggplot2::element_line(color = "#E0E4F0", linewidth = 0.25),
        panel.grid.major.x = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(
          angle = 45, hjust = 1, vjust = 1,
          color = "#283046",
          size = 9
        ),
        axis.text.y = ggplot2::element_text(color = "#283046", size = 9),
        axis.title.x = ggplot2::element_text(color = "#5f6473", size = 10, margin = ggplot2::margin(t = 6)),
        axis.title.y = ggplot2::element_text(color = "#5f6473", size = 10, margin = ggplot2::margin(r = 6)),
        plot.title = ggplot2::element_text(
          hjust  = 0.5,
          face   = "bold",
          color  = "#283046",
          size   = 11,
          margin = ggplot2::margin(b = 6)
        ),
        legend.position = "none"
      )

    subsystem_done(TRUE)
    if (subsystem_done() && reactions_done()) {
      set_souris_progress(show = FALSE)
    }
    session$sendCustomMessage("souris-show-plot", "souris-violin-subsys-wrapper")

    plt <- plotly::ggplotly(p_sub, tooltip = c("x", "y"))
    plt <- plotly::config(
      plt,
      toImageButtonOptions = list(
        format   = "svg",   # this makes the camera button download SVG
        filename = "souris_violin_subsys"
      )
    )
    
    plt
  })

  output$souris_violin_subsys_sc <- plotly::renderPlotly({
    df_all <- df_all_predictions()
    subsel <- clicked_subsystem_sc()
    req(df_all, subsel)
    
    df_sub <- df_all[df_all$Subsystem == subsel, ]
    
    validate(need(
      nrow(df_sub) > 0,
      "No predictions found for the selected subsystem."
    ))
    
    # scale prediction 0–1 like the main violin
    rng <- range(df_sub$Prediction, na.rm = TRUE)
    if (diff(rng) == 0 || !is.finite(diff(rng))) {
      df_sub$Prediction_scaled <- 0.5
    } else {
      df_sub$Prediction_scaled <- (df_sub$Prediction - rng[1]) / (rng[2] - rng[1])
    }
    
    # Make sure Reference is a factor with the current order
    df_sub$Reference <- factor(df_sub$Reference, levels = unique(df_sub$Reference))
    ref_levels <- levels(df_sub$Reference)
    validate(need(
      length(ref_levels) == 2,
      "Subsystem AUC is defined here only for exactly two reference assays."
    ))
    
    ref1 <- ref_levels[1]
    ref2 <- ref_levels[2]
    
    # AUC of subsystem prediction between the two references
    auc_sub <- auc_two_groups(
      pred = df_sub$Prediction,
      grp  = df_sub$Reference
    )
    auc_label <- sprintf("AUC = %.3f", auc_sub)
    
    p_sub <- ggplot2::ggplot(
      df_sub,
      ggplot2::aes(x = Reference, y = Prediction_scaled, fill = Reference)
    ) +
      ggplot2::geom_violin(
        alpha = 0.5,
        scale = "width",
        trim  = FALSE,
        color = NA
      ) +
      ggplot2::geom_boxplot(
        width         = 0.14,
        outlier.size  = 0.6,
        outlier.alpha = 0.5,
        fill          = "white",
        color         = "#455A64",
        alpha         = 0.9,
        size          = 0.3
      ) +
      ggplot2::geom_jitter(
        width = 0.06,
        size  = 0.8,
        alpha = 0.35,
        color = "#263238"
      ) +
      ggplot2::labs(
        x     = "Reference assay",
        y     = "Subsystem prediction probability",
        title = paste("Subsystem:", subsel, "|", auc_label)
      ) +
      ggplot2::scale_fill_manual(
        values = c(
          "#1565C0",  # ref 1
          "#42A5F5"   # ref 2 (recycles if extra)
        )
      ) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(
        panel.background   = ggplot2::element_rect(fill = NA, colour = NA),
        plot.background    = ggplot2::element_rect(fill = NA, colour = NA),
        panel.grid.major.y = ggplot2::element_line(color = "#E0E4F0", linewidth = 0.25),
        panel.grid.major.x = ggplot2::element_blank(),
        panel.grid.minor   = ggplot2::element_blank(),
        axis.text.x  = ggplot2::element_text(
          angle = 45, hjust = 1, vjust = 1,
          color = "#283046",
          size  = 9
        ),
        axis.text.y  = ggplot2::element_text(color = "#283046", size = 9),
        axis.title.x = ggplot2::element_text(color = "#5f6473", size = 10, margin = ggplot2::margin(t = 6)),
        axis.title.y = ggplot2::element_text(color = "#5f6473", size = 10, margin = ggplot2::margin(r = 6)),
        plot.title   = ggplot2::element_text(
          hjust  = 0.5,
          face   = "bold",
          color  = "#283046",
          size   = 11,
          margin = ggplot2::margin(b = 6)
        ),
        legend.position = "none"
      )
    
    subsystem_done(TRUE)
    if (subsystem_done() && reactions_done()) {
      set_souris_progress(show = FALSE)
    }
    
    session$sendCustomMessage("souris-show-plot", "souris-violin-subsys-sc-wrapper")
    
    plt <- plotly::ggplotly(p_sub, tooltip = c("x", "y"))
    plt <- plotly::config(
      plt,
      toImageButtonOptions = list(
        format   = "svg",   # this makes the camera button download SVG
        filename = "souris_violin_subsys_sc"
      )
    )
    
    plt
  })

  output$souris_auc_reactions <- plotly::renderPlotly({
    subsel <- clicked_subsystem()
    req(subsel)

    refs <- input$souris_reference_sel
    validate(need(
      length(refs) >= 2,
      "Select at least two reference assays to compute reaction-level AUC."
    ))
    ref1 <- refs[1]
    ref2 <- refs[2]

    # reaction matrices: rows = samples, cols = reactions
    mat1 <- get_reaction_preds(ref1, subsel)
    mat2 <- get_reaction_preds(ref2, subsel)
    validate(need(
      !is.null(mat1) && !is.null(mat2),
      "No reaction predictions available for the selected subsystem."
    ))

    # intersect reactions and align samples
    common_rxn <- intersect(colnames(mat1), colnames(mat2))
    validate(need(
      length(common_rxn) > 0,
      "No overlapping reactions for this subsystem in the two reference assays."
    ))

    mat1 <- mat1[, common_rxn, drop = FALSE]
    mat2 <- mat2[, common_rxn, drop = FALSE]

    # build long df: each row = sample-reaction, with ref label
    df_rxn <- rbind(
      data.frame(
        Reference = ref1,
        Sample = rownames(mat1),
        as.data.frame(mat1, check.names = FALSE)
      ),
      data.frame(
        Reference = ref2,
        Sample = rownames(mat2),
        as.data.frame(mat2, check.names = FALSE)
      )
    )

    # compute AUC per reaction across the two references
    rxn_names <- common_rxn
    auc_rxn <- sapply(rxn_names, function(r) {
      preds <- df_rxn[[r]]
      grp <- factor(df_rxn$Reference, levels = c(ref1, ref2))
      auc_two_groups(preds, grp)
    })

    auc_df <- data.frame(
      Reaction = rxn_names,
      AUC = auc_rxn,
      stringsAsFactors = FALSE
    )
    auc_df <- dplyr::filter(auc_df, !is.na(AUC))
    validate(need(
      nrow(auc_df) > 0,
      "No valid reaction-level AUC values could be computed."
    ))

    recon_map <- get_recon_names(auc_df$Reaction)
    auc_df$ReactionLabel <- recon_map[auc_df$Reaction]

    # same visual style as subsystem AUC bar
    p_rxn <- ggplot2::ggplot(
      auc_df,
      ggplot2::aes(
        x         = reorder(ReactionLabel, AUC),
        y         = AUC,
        customdata = ReactionLabel
      )
    ) +
      ggplot2::geom_col(
        fill = "rgba(21,101,192,0.85)"  # same blue as subsystem AUC (approx)
      ) +
      ggplot2::coord_flip() +
      ggplot2::geom_hline(
        yintercept = 0.5,
        linetype   = "dashed",
        color      = "#EF5350",
        linewidth  = 0.6
      ) +
      ggplot2::labs(
        x     = "Reaction (RECON3D)",
        y     = "AUC (first vs second reference)",
        title = paste("Reaction AUC for subsystem", subsel)
      ) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(
        panel.background   = ggplot2::element_rect(fill = NA, colour = NA),
        plot.background    = ggplot2::element_rect(fill = NA, colour = NA),
        panel.grid.major.y = ggplot2::element_line(color = "#E0E4F0", linewidth = 0.25),
        panel.grid.major.x = ggplot2::element_blank(),
        panel.grid.minor   = ggplot2::element_blank(),
        axis.text.x  = ggplot2::element_text(color = "#283046", size = 9),
        axis.text.y  = ggplot2::element_text(color = "#283046", size = 9),
        axis.title.x = ggplot2::element_text(color = "#5f6473", size = 10, margin = ggplot2::margin(t = 6)),
        axis.title.y = ggplot2::element_text(color = "#5f6473", size = 10, margin = ggplot2::margin(r = 6)),
        plot.title   = ggplot2::element_text(
          hjust  = 0.5,
          face   = "bold",
          color  = "#283046",
          size   = 11,
          margin = ggplot2::margin(b = 6)
        ),
        legend.position = "none"
      )

    # fade-in trigger for the reactions AUC panel
    session$sendCustomMessage("souris-show-plot", "souris-auc-reactions-wrapper")

    plt <- plotly::ggplotly(p_rxn, tooltip = c("x", "y"))
    plt <- plotly::config(
      plt,
      toImageButtonOptions = list(
        format   = "svg",   # this makes the camera button download SVG
        filename = "souris_auc_reactions"
      )
    )
    
    plt
  })

  output$souris_auc_reactions_sc <- plotly::renderPlotly({
    subsel <- clicked_subsystem_sc()
    req(subsel)
    
    refs <- input$souris_reference_sel_sc
    validate(need(
      length(refs) >= 2,
      "Select at least two reference assays to compute reaction-level AUC."
    ))
    ref1 <- refs[1]
    ref2 <- refs[2]
    
    # reaction matrices: rows = samples, cols = reactions
    mat1 <- get_reaction_preds(ref1, subsel)
    mat2 <- get_reaction_preds(ref2, subsel)
    validate(need(
      !is.null(mat1) && !is.null(mat2),
      "No reaction predictions available for the selected subsystem."
    ))
    
    # intersect reactions and align samples
    common_rxn <- intersect(colnames(mat1), colnames(mat2))
    validate(need(
      length(common_rxn) > 0,
      "No overlapping reactions for this subsystem in the two reference assays."
    ))
    
    mat1 <- mat1[, common_rxn, drop = FALSE]
    mat2 <- mat2[, common_rxn, drop = FALSE]
    
    # build long df: each row = sample-reaction, with ref label
    df_rxn <- rbind(
      data.frame(
        Reference = ref1,
        Sample    = rownames(mat1),
        as.data.frame(mat1, check.names = FALSE)
      ),
      data.frame(
        Reference = ref2,
        Sample    = rownames(mat2),
        as.data.frame(mat2, check.names = FALSE)
      )
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
    validate(need(
      nrow(auc_df) > 0,
      "No valid reaction-level AUC values could be computed."
    ))
    
    recon_map <- get_recon_names(auc_df$Reaction)
    auc_df$ReactionLabel <- recon_map[auc_df$Reaction]
    
    p_rxn <- ggplot2::ggplot(
      auc_df,
      ggplot2::aes(
        x         = reorder(ReactionLabel, AUC),
        y         = AUC,
        customdata = ReactionLabel
      )
    ) +
      ggplot2::geom_col(
        fill = "rgba(21,101,192,0.85)"
      ) +
      ggplot2::coord_flip() +
      ggplot2::geom_hline(
        yintercept = 0.5,
        linetype   = "dashed",
        color      = "#EF5350",
        linewidth  = 0.6
      ) +
      ggplot2::labs(
        x     = "Reaction (RECON3D)",
        y     = "AUC (first vs second reference)",
        title = paste("Reaction AUC for subsystem", subsel, "|", ref1, "vs", ref2)
      ) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(
        panel.background   = ggplot2::element_rect(fill = NA, colour = NA),
        plot.background    = ggplot2::element_rect(fill = NA, colour = NA),
        panel.grid.major.y = ggplot2::element_line(color = "#E0E4F0", linewidth = 0.25),
        panel.grid.major.x = ggplot2::element_blank(),
        panel.grid.minor   = ggplot2::element_blank(),
        axis.text.x  = ggplot2::element_text(color = "#283046", size = 9),
        axis.text.y  = ggplot2::element_text(color = "#283046", size = 9),
        axis.title.x = ggplot2::element_text(color = "#5f6473", size = 10, margin = ggplot2::margin(t = 6)),
        axis.title.y = ggplot2::element_text(color = "#5f6473", size = 10, margin = ggplot2::margin(r = 6)),
        plot.title   = ggplot2::element_text(
          hjust  = 0.5,
          face   = "bold",
          color  = "#283046",
          size   = 11,
          margin = ggplot2::margin(b = 6)
        ),
        legend.position = "none"
      )
    
    session$sendCustomMessage("souris-show-plot", "souris-auc-reactions-sc-wrapper")
    
    plt <- plotly::ggplotly(p_rxn, tooltip = c("x", "y"), source = "souris_auc_rxn_sc")
    plt <- plotly::config(
      plt,
      toImageButtonOptions = list(
        format   = "svg",   # this makes the camera button download SVG
        filename = "souris_auc_reactions_sc"
      )
    )
    
    plt
  })

  output$souris_violin_reactions <- plotly::renderPlotly({
    subsel <- clicked_subsystem()
    req(subsel)

    refs <- input$souris_reference_sel
    validate(need(
      length(refs) >= 2,
      "Select at least two reference assays to compare reaction predictions."
    ))

    ref1 <- refs[1]
    ref2 <- refs[2]

    # 1) Build long df of raw reaction values from reference assays
    df_rxn_long <- do.call(
      rbind,
      lapply(refs, function(ref_id) {
        mat <- get_reaction_preds(ref_id, subsel) # samples × reactions
        validate(need(
          !is.null(mat),
          paste("No reactions for subsystem", subsel, "in", ref_id)
        ))
        if (!nrow(mat) || !ncol(mat)) {
          return(NULL)
        }

        df <- as.data.frame(mat, check.names = FALSE)
        df$Sample <- rownames(mat)
        df_long <- tidyr::pivot_longer(
          df,
          cols      = -Sample,
          names_to  = "Reaction",
          values_to = "Value"
        )
        df_long$Reference <- ref_id
        df_long
      })
    )

    validate(need(
      nrow(df_rxn_long) > 0,
      "No reaction values available to plot."
    ))

    df_rxn_long$Reference <- factor(df_rxn_long$Reference,
      levels = unique(df_rxn_long$Reference)
    )

    recon_map <- get_recon_names(unique(df_rxn_long$Reaction))
    df_rxn_long$ReactionLabel <- recon_map[df_rxn_long$Reaction]

    # 2) Per-reaction Wilcoxon between ref1 and ref2 (for ordering + stars)
    df_pair <- df_rxn_long[df_rxn_long$Reference %in% c(ref1, ref2), ]

    stats_df <- df_pair |>
      dplyr::group_by(ReactionLabel) |>
      dplyr::summarise(
        p_value = tryCatch(
          stats::wilcox.test(Value ~ Reference, exact = FALSE)$p.value,
          error = function(e) NA_real_
        ),
        .groups = "drop"
      ) |>
      dplyr::mutate(
        p_value = ifelse(is.na(p_value), 1, p_value),
        sig_score = -log10(p_value + 1e-12)
      )

    rxn_levels <- stats_df |>
      dplyr::arrange(dplyr::desc(sig_score)) |>
      dplyr::pull(ReactionLabel)

    df_rxn_long$ReactionLabel <- factor(df_rxn_long$ReactionLabel, levels = rxn_levels)

    # 3) Star labels per reaction
    p_to_stars <- function(p) {
      if (is.na(p)) {
        return("ns")
      }
      if (p < 0.0001) {
        return("****")
      }
      if (p < 0.001) {
        return("***")
      }
      if (p < 0.01) {
        return("**")
      }
      if (p < 0.05) {
        return("*")
      }
      "ns"
    }

    max_y <- df_rxn_long |>
      dplyr::group_by(ReactionLabel) |>
      dplyr::summarise(y_max = max(Value, na.rm = TRUE), .groups = "drop")

    star_df <- stats_df |>
      dplyr::inner_join(max_y, by = "ReactionLabel") |>
      dplyr::mutate(
        label = vapply(p_value, p_to_stars, character(1)),
        x_num = as.numeric(factor(ReactionLabel, levels = levels(df_rxn_long$ReactionLabel))),
        y_pos = y_max * 1.08
      )

    
    df_rxn_long$ReactionLabel <- factor(df_rxn_long$ReactionLabel, levels = rxn_levels)
    df_rxn_long$Reference     <- factor(df_rxn_long$Reference, levels = c(ref1, ref2))
    
    rxn_index <- as.numeric(df_rxn_long$ReactionLabel)
    offset <- ifelse(
      df_rxn_long$Reference == ref1, -0.2,
      ifelse(df_rxn_long$Reference == ref2, 0.2, 0)
    )
    
    p <- ggplot2::ggplot(
      df_rxn_long,
      ggplot2::aes(
        y    = Value,
        fill = Reference
      )
    ) +
      # Violin exactly at the boxplot x positions
      ggplot2::geom_violin(
        ggplot2::aes(
          x     = rxn_index + offset,
          group = interaction(ReactionLabel, Reference)
        ),
        alpha  = 0.5,
        scale  = "width",
        trim   = FALSE,
        color  = NA,
        width  = 0.05  
      ) +
      # Boxplot on top, same x mapping
      ggplot2::geom_boxplot(
        ggplot2::aes(
          x     = rxn_index + offset,
          group = interaction(ReactionLabel, Reference)
        ),
        outlier.size = 0.4,
        alpha        = 0.9,
        color        = "#455A64",
        width        = 0.35
      ) +
      ggplot2::geom_segment(
        data = star_df,
        ggplot2::aes(
          x    = x_num - 0.25,
          xend = x_num + 0.25,
          y    = y_pos,
          yend = y_pos
        ),
        inherit.aes = FALSE,
        linewidth   = 0.4
      ) +
      ggplot2::geom_text(
        data = star_df,
        ggplot2::aes(
          x     = x_num,
          y     = y_pos,
          label = label
        ),
        vjust        = -0.2,
        inherit.aes  = FALSE,
        size         = 3
      ) +
      ggplot2::scale_x_continuous(
        breaks = seq_along(levels(df_rxn_long$ReactionLabel)),
        labels = levels(df_rxn_long$ReactionLabel)
      ) +
      ggplot2::scale_fill_manual(
        values = setNames(
          c("#1565C0", "#42A5F5"),
          c(ref1,      ref2)
        )
      ) +
      ggplot2::labs(
        x = "Reaction (RECON3D)",
        y = "Reaction value (reference assay)",
        title = paste0(
          "Per-reaction values in subsystem ", subsel,
          " | ", ref1, " vs ", ref2
        )
      ) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(
        panel.background   = ggplot2::element_rect(fill = NA, colour = NA),
        plot.background    = ggplot2::element_rect(fill = NA, colour = NA),
        panel.grid.major.y = ggplot2::element_line(color = "#E0E4F0", linewidth = 0.25),
        panel.grid.major.x = ggplot2::element_blank(),
        panel.grid.minor   = ggplot2::element_blank(),
        axis.text.x  = ggplot2::element_text(
          angle = 45, hjust = 1, vjust = 1,
          color = "#283046",
          size  = 9
        ),
        axis.text.y  = ggplot2::element_text(color = "#283046", size = 9),
        axis.title.x = ggplot2::element_text(color = "#5f6473", size = 10, margin = ggplot2::margin(t = 6)),
        axis.title.y = ggplot2::element_text(color = "#5f6473", size = 10, margin = ggplot2::margin(r = 6)),
        plot.title   = ggplot2::element_text(
          hjust  = 0.5,
          face   = "bold",
          color  = "#283046",
          size   = 11,
          margin = ggplot2::margin(b = 6)
        ),
        legend.position = "none"
      )

    reactions_done(TRUE)
    if (subsystem_done() && reactions_done()) {
      set_souris_progress(show = FALSE)
    }

    session$sendCustomMessage("souris-show-plot", "souris-violin-reactions-wrapper")
    
    plt <- plotly::ggplotly(p, tooltip = c("text"))
    plt <- plotly::config(
      plt,
      toImageButtonOptions = list(
        format   = "svg",   # this makes the camera button download SVG
        filename = "souris_violin_reactions"
      )
    )
    
    plt
  })

  output$souris_violin_reactions_sc <- plotly::renderPlotly({
    subsel <- clicked_subsystem_sc()
    req(subsel)
    
    refs <- input$souris_reference_sel_sc
    validate(need(
      length(refs) >= 2,
      "Select at least two reference assays to compare reaction predictions."
    ))
    ref1 <- refs[1]
    ref2 <- refs[2]
    
    # 1) Long df of raw reaction values (unchanged)
    df_rxn_long <- do.call(
      rbind,
      lapply(refs, function(ref_id) {
        mat <- get_reaction_preds(ref_id, subsel)
        validate(need(
          !is.null(mat),
          paste("No reactions for subsystem", subsel, "in", ref_id)
        ))
        if (!nrow(mat) || !ncol(mat)) return(NULL)
        
        df <- as.data.frame(mat, check.names = FALSE)
        df$Sample <- rownames(mat)
        df_long <- tidyr::pivot_longer(
          df,
          cols      = -Sample,
          names_to  = "Reaction",
          values_to = "Value"
        )
        df_long$Reference <- ref_id
        df_long
      })
    )
    
    validate(need(
      nrow(df_rxn_long) > 0,
      "No reaction values available to plot."
    ))
    
    # Join react_meta_filt metadata using safe matching
    idx <- match(df_rxn_long$Reaction, rownames(react_meta_filt))
    meta_df <- react_meta_filt[
      idx,
      c("RECON3D", "equation", "name", "compartments", "geneAssociation"),
      drop = FALSE
    ]
    
    missing_meta <- is.na(idx)
    if (any(missing_meta)) {
      meta_df$RECON3D[missing_meta] <- df_rxn_long$Reaction[missing_meta]
    }
    meta_df$name[is.na(meta_df$name)]               <- ""
    meta_df$equation[is.na(meta_df$equation)]       <- ""
    meta_df$compartments[is.na(meta_df$compartments)] <- ""
    meta_df$geneAssociation[is.na(meta_df$geneAssociation)] <- ""
    
    df_rxn_long <- cbind(df_rxn_long, meta_df)
    
    df_rxn_long$Reference <- factor(
      df_rxn_long$Reference,
      levels = unique(df_rxn_long$Reference)
    )
    
    recon_map <- get_recon_names(unique(df_rxn_long$Reaction))
    df_rxn_long$ReactionLabel <- recon_map[df_rxn_long$Reaction]
    
    # 2) Per-reaction Wilcoxon (ordering + stars) (unchanged)
    df_pair <- df_rxn_long[df_rxn_long$Reference %in% c(ref1, ref2), ]
    
    stats_df <- df_pair |>
      dplyr::group_by(ReactionLabel) |>
      dplyr::summarise(
        p_value = tryCatch(
          stats::wilcox.test(Value ~ Reference, exact = FALSE)$p.value,
          error = function(e) NA_real_
        ),
        .groups = "drop"
      ) |>
      dplyr::mutate(
        p_value  = ifelse(is.na(p_value), 1, p_value),
        sig_score = -log10(p_value + 1e-12)
      )
    
    rxn_levels <- stats_df |>
      dplyr::arrange(dplyr::desc(sig_score)) |>
      dplyr::pull(ReactionLabel)
    
    df_rxn_long$ReactionLabel <- factor(df_rxn_long$ReactionLabel,
                                        levels = rxn_levels)
    df_rxn_long$Reference <- factor(df_rxn_long$Reference,
                                    levels = c(ref1, ref2))
    
    # 3) Star labels per reaction (unchanged)
    p_to_stars <- function(p) {
      if (is.na(p))    return("ns")
      if (p < 0.0001) return("****")
      if (p < 0.001)  return("***")
      if (p < 0.01)   return("**")
      if (p < 0.05)   return("*")
      "ns"
    }
    
    max_y <- df_rxn_long |>
      dplyr::group_by(ReactionLabel) |>
      dplyr::summarise(y_max = max(Value, na.rm = TRUE), .groups = "drop")
    
    star_df <- stats_df |>
      dplyr::inner_join(max_y, by = "ReactionLabel") |>
      dplyr::mutate(
        label = vapply(p_value, p_to_stars, character(1)),
        x_num = as.numeric(factor(ReactionLabel,
                                  levels = levels(df_rxn_long$ReactionLabel))),
        y_pos = y_max * 1.08
      )
    
    # Build x positions and hover text per row
    rxn_index <- as.numeric(df_rxn_long$ReactionLabel)
    offset <- ifelse(
      df_rxn_long$Reference == ref1, -0.2,
      ifelse(df_rxn_long$Reference == ref2, 0.2, 0)
    )
    x_pos <- rxn_index + offset
    
    hover_text <- paste0(
      "RECON3D: ", df_rxn_long$RECON3D, "\n",
      "Name: ", df_rxn_long$name, "\n",
      "Equation: ", df_rxn_long$equation, "\n",
      "Compartments: ", df_rxn_long$compartments, "\n",
      "GPR: ", df_rxn_long$geneAssociation, "\n",
      "Reference: ", df_rxn_long$Reference, "\n",
      "Value: ", signif(df_rxn_long$Value, 4)
    )
    
    # Build plotly directly: one violin trace per reference
    plt <- plotly::plot_ly()
    
    for (ref in c(ref1, ref2)) {
      df_sub <- df_rxn_long[df_rxn_long$Reference == ref, , drop = FALSE]
      x_sub  <- x_pos[df_rxn_long$Reference == ref]
      txt    <- hover_text[df_rxn_long$Reference == ref]
      
      plt <- plt |>
        plotly::add_trace(
          type = "violin",
          x    = x_sub,
          y    = df_sub$Value,
          text = txt,
          hoverinfo = "text",
          name = ref,
          legendgroup = ref,
          scalemode = "width",
          spanmode  = "hard",
          opacity   = 0.5,
          line = list(color = "rgba(0,0,0,0)"),
          fillcolor = if (ref == ref1) "#1565C0" else "#42A5F5",
          showlegend = FALSE
        )
    }
    
    # Add boxplot layer as separate traces if you still want them
    for (ref in c(ref1, ref2)) {
      df_sub <- df_rxn_long[df_rxn_long$Reference == ref, , drop = FALSE]
      x_sub  <- x_pos[df_rxn_long$Reference == ref]
      txt    <- hover_text[df_rxn_long$Reference == ref]
      
      plt <- plt |>
        plotly::add_boxplot(
          x    = x_sub,
          y    = df_sub$Value,
          text = txt,
          hoverinfo = "text",
          name = ref,
          marker = list(
            color = if (ref == ref1) "#1565C0" else "#42A5F5",
            size  = 2
          ),
          line   = list(color = "#455A64"),
          fillcolor = "rgba(255,255,255,0.9)",
          boxpoints = FALSE,
          showlegend = FALSE
        )
    }
    
    # Add significance segments and labels as shapes/annotations
    for (i in seq_len(nrow(star_df))) {
      plt <- plt |>
        plotly::add_segments(
          x    = star_df$x_num[i] - 0.25,
          xend = star_df$x_num[i] + 0.25,
          y    = star_df$y_pos[i],
          yend = star_df$y_pos[i],
          inherit = FALSE,
          showlegend = FALSE,
          line = list(color = "black", width = 0.4)
        ) |>
        plotly::add_text(
          x = star_df$x_num[i],
          y = star_df$y_pos[i],
          text = star_df$label[i],
          showlegend = FALSE,
          textposition = "top center"
        )
    }
    
    plt <- plt |>
      plotly::layout(
        xaxis = list(
          title = "Reaction (RECON3D)",
          tickmode = "array",
          tickvals = seq_along(levels(df_rxn_long$ReactionLabel)),
          ticktext = levels(df_rxn_long$ReactionLabel)
        ),
        yaxis = list(
          title = "Reaction value (reference assay)"
        ),
        title = list(
          text = paste0(
            "Per-reaction values in subsystem ", subsel,
            " | ", ref1, " vs ", ref2
          )
        ),
        hovermode = "closest",
        paper_bgcolor = "rgba(0,0,0,0)",
        plot_bgcolor  = "rgba(0,0,0,0)"
      )
    
    reactions_done(TRUE)
    if (subsystem_done() && reactions_done()) {
      set_souris_progress(show = FALSE)
    }
    
    session$sendCustomMessage(
      "souris-show-plot",
      "souris-violin-reactions-sc-wrapper"
    )
    
    plt
  })

  output$souris_sc_featureplot <- renderPlot({
    req(souris_sc_object())
    req(sc_fingerprint_scores())

    obj_sc <- souris_sc_object()[[1]]

    p_list <- Seurat::FeaturePlot(
      obj_sc,
      features = "SOURIS_fingerprint_score",
      split.by = "SOURIS_reference_group"
    ) 
    
    fix.sc <- scale_color_gradient2(
      low = "blue2",
      mid = "lavenderblush",
      high = "firebrick1",
      midpoint = 0.5,
      limits = c(0,1)
    )
    p_list <- lapply(
      p_list,
      function(x) {
        x + fix.sc +
          theme(
            axis.title.y = element_blank(),
            axis.title.x = element_blank()
          )
      }
    )
    
    # legend from first panel
    p_legend_src <- p_list[[1]] + theme(legend.position = "right")
    leg <- cowplot::get_legend(p_legend_src)
    
    # panels without legends
    p_list_nl <- lapply(p_list, function(x) x + theme(legend.position = "none"))
    p_panels  <- cowplot::plot_grid(plotlist = p_list_nl)
    
    cowplot::plot_grid(
      p_panels,
      leg,
      ncol = 2,
      rel_widths = c(1, 0.08)
    )
    
  })
  
  output$download_souris_sc_featureplot <- downloadHandler(
    filename = function() {
      paste0("souris_sc_featureplot_", Sys.Date(), ".pdf")
    },
    content = function(file) {
      # same plot logic, but drawn into an SVG device
      req(souris_sc_object())
      req(sc_fingerprint_scores())
      
      obj_sc <- souris_sc_object()[[1]]
      
      plist <- Seurat::FeaturePlot(
        obj_sc,
        features = "SOURIS_fingerprint_score",
        split.by = "SOURIS_reference_group"
      )
      
      fix.sc <- scale_color_gradient2(
        low = "blue2",
        mid = "lavenderblush",
        high = "firebrick1",
        midpoint = 0.5,
        limits = c(0, 1)
      )
      plist <- lapply(plist, function(x) x + fix.sc)
      p <- cowplot::plot_grid(plotlist = plist)
      
      pdf(file, width = 8, height = 6)
      print(p)
      dev.off()
    }
  )
  
  output$souris_sc_featureplot_subsystem_title <- renderUI({
    subsel <- clicked_subsystem_sc()
    req(subsel)
    titleText <- paste0("FeaturePlot · subsystem: ", subsel, " (SC)")
    tags$span(titleText)
  })
  
  output$souris_sc_featureplot_subsystem <- renderPlot({
    req(souris_sc_object())
    req(df_all_predictions())

    subsel <- clicked_subsystem_sc()
    req(subsel)

    obj_sc <- souris_sc_object()[[1]]
    df_all <- df_all_predictions()

    # filter to selected subsystem
    df_sub <- df_all[df_all$Subsystem == subsel, ]
    validate(
      need(
        nrow(df_sub) > 0,
        "No predictions found for the selected subsystem."
      )
    )

    # df_sub should have columns: Sample (cell ID), Reference, Prediction
    # align to Seurat cells

    idx <- match(colnames(obj_sc), df_sub$Sample)
    validate(
      need(
        any(!is.na(idx)),
        "No matching cells between subsystem predictions and single-cell object."
      )
    )

    scores <- rep(NA_real_, ncol(obj_sc))
    scores[!is.na(idx)] <- df_sub$Prediction[idx[!is.na(idx)]]

    # scale to [0,1] for a full gradient
    rng <- range(scores, na.rm = TRUE)
    if (diff(rng) == 0 || !is.finite(diff(rng))) {
      scores_scaled <- rep(0.5, length(scores))
    } else {
      scores_scaled <- (scores - rng[1]) / diff(rng)
    }

    obj_sc$SOURIS_subsystem_score <- scores_scaled
    
    p_list <- Seurat::FeaturePlot(
      obj_sc,
      features = "SOURIS_subsystem_score",
      split.by = "SOURIS_reference_group"
    ) 
    
    fix.sc <- scale_color_gradient2(
      low = "blue2",
      mid = "lavenderblush",
      high = "firebrick1",
      midpoint = 0.5,
      limits = c(0,1),
      name = "Subsystem score"
    )
    p_list <- lapply(
      p_list,
      function(x) {
        x + fix.sc +
          theme(
            axis.title.y = element_blank(),
            axis.title.x = element_blank()
          )
      }
    )
    
    # use first panel as legend source
    p_legend_src <- p_list[[1]] + theme(legend.position = "right")
    
    # extract legend grob
    leg <- cowplot::get_legend(p_legend_src)
    
    # remove legends from actual panels
    p_list_nl <- lapply(p_list, function(x) x + theme(legend.position = "none"))
    
    # combine panels (e.g. default grid layout)
    p_panels <- cowplot::plot_grid(plotlist = p_list_nl)
    
    title_grob <- ggplot2::ggplot() +
      ggplot2::annotate(
        "text",
        x = 0.5, y = 0.5,
        label = subsel,
        fontface = "bold",
        size = 4
      ) +
      ggplot2::theme_void()
    
    cowplot::plot_grid(
      title_grob,
      cowplot::plot_grid(
        p_panels,
        leg,
        ncol = 2,
        rel_widths = c(1, 0.08)
      ),
      ncol = 1,
      rel_heights = c(0.08, 0.92)
    )
  })

  observeEvent(
    plotly::event_data("plotly_click", source = "souris_auc_rxn_sc"),
    {
      d <- plotly::event_data("plotly_click", source = "souris_auc_rxn_sc")
      if (is.null(d)) {
        return()
      }

      rxn_label <- as.character(d$customdata) # or d$customdata if you used that
      clicked_reaction_sc(rxn_label)

      if (!sc_reaction_feature_panel_inserted()) {
        insertUI(
          selector = "#sc_reaction_featureplot_placeholder",
          where    = "afterEnd",
          ui = fluidRow(
            class = "fade-in-up",
            column(
              width = 12,
              box(
                width = 12, status = "primary", solidHeader = TRUE,
                title = NULL,
                div(
                  style = "display:flex; flex-direction:column; height:500px;",
                  div(
                    id = "souris-sc-feature-reaction-header",
                    style = "
                position: sticky;
                top: 0;
                z-index: 5;
                padding: 6px 10px;
                background: linear-gradient(135deg, #002B5C, #0D47A1);
                color: #E3F2FD;
                display: flex;
                align-items: center;
                justify-content: space-between;
                gap: 8px;
                box-shadow: 0 2px 6px rgba(0,0,0,0.15);
              ",
                    tags$div(
                      style = 'display:flex; flex-direction:column; gap:2px;',
                      tags$div(
                        style = 'font-size:12px; font-weight:600; letter-spacing:0.06em; text-transform:uppercase;',
                        uiOutput("souris_sc_featureplot_reaction_title")
                      ),
                      tags$div(
                        style = 'font-size:11px; color:#BBDEFB;',
                        "Reaction activity per reference group · single-cell"
                      )
                    )
                  ),
                  div(
                    style = "flex:1 1 auto; overflow-y:auto; padding:4px 0 0 0;",
                    plotOutput("souris_sc_featureplot_reaction",
                               height = "100%", width = "100%")
                  )
                )
              )
            )
          ),
          immediate = TRUE
        )
        sc_reaction_feature_panel_inserted(TRUE)
      }
    }
  )

  output$souris_sc_featureplot_reaction_title <- renderUI({
    rxnlab <- clicked_reaction_sc()
    req(rxnlab)
    titleText <- paste0("FeaturePlot · reaction: ", rxnlab, " (SC)")
    tags$span(titleText)
  })
  
  output$souris_sc_featureplot_reaction <- renderPlot({
    req(souris_sc_object())
    obj_sc <- souris_sc_object()[[1]]

    subsel <- clicked_subsystem_sc()
    rxnlab <- clicked_reaction_sc()
    req(subsel, rxnlab)

    rxn_id <- rownames(subset(react_meta_filt, RECON3D == rxnlab))[1]

    refs <- souris_reference_list()
    validate(
      need(
        !is.null(refs) && length(refs) > 0,
        "Reference assays must be loaded to retrieve reaction values."
      )
    )

    vals <- lapply(refs, function(df) {
      if (!rxn_id %in% colnames(df)) {
        return(NULL)
      }
      data.frame(
        Cell = df[[1]],
        Reaction = df[[rxn_id]],
        stringsAsFactors = FALSE
      )
    })
    vals <- vals[!vapply(vals, is.null, logical(1))]
    validate(
      need(
        length(vals) > 0,
        "Selected reaction not found in the reference assays."
      )
    )
    dfr <- dplyr::bind_rows(vals)

    idx <- match(colnames(obj_sc), dfr$Cell)
    validate(
      need(
        any(!is.na(idx)),
        "No matching cells between reaction values and single-cell object."
      )
    )

    rxn_scores <- rep(NA_real_, ncol(obj_sc))
    rxn_scores[!is.na(idx)] <- dfr$Reaction[idx[!is.na(idx)]]

    obj_sc$SOURIS_reaction_score <- rxn_scores

    rng <- range(obj_sc$SOURIS_reaction_score, na.rm = TRUE)

    p_list <- Seurat::FeaturePlot(
      obj_sc,
      features = "SOURIS_reaction_score",
      split.by = "SOURIS_reference_group"
    ) 
    
    fix.sc <- scale_color_gradient2(
      low = "blue2",
      mid = "lavenderblush",
      high = "firebrick1",
      midpoint = median(rxn_scores, na.rm = TRUE),
      limits = rng,
      name = rxnlab
    )
    p_list <- lapply(
      p_list,
      function(x) {
        x + fix.sc +
          theme(
            axis.title.y = element_blank(),  # remove vertical SOURIS_reaction_score
            axis.title.x = element_blank()   # optional: also remove x if desired
          )
      }
    )
    
    # Use the first panel as legend source
    p_legend_src <- p_list[[1]] + theme(legend.position = "right")
    
    # extract legend grob
    leg <- cowplot::get_legend(p_legend_src)
    
    # remove legends from actual panels
    p_list_nl <- lapply(p_list, function(x) x + theme(legend.position = "none"))
    
    # combine panels horizontally (or however you arrange them)
    p_panels <- cowplot::plot_grid(plotlist = p_list_nl)
    
    cowplot::plot_grid(
      p_panels,
      leg,
      ncol = 2,
      rel_widths = c(1, 0.08)  # adjust legend column width as needed
    )
  })
}

shinyApp(ui, server)
