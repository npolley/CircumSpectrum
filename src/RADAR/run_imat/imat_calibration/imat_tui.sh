#!/bin/bash

# Terminal-based GUI for imat_launcher.sh using whiptail.
# If whiptail is not available, you can swap it with dialog and adjust options.

# -------- ASSAY (free text) --------
ASSAY=$(whiptail --inputbox "Enter assay name (ASSAY):" 10 60 3>&1 1>&2 2>&3)
exitstatus=$?
if [ $exitstatus -ne 0 ]; then
  echo "Cancelled."
  exit 1
fi

# -------- MEDIUM (dropdown / menu) --------
MEDIUM=$(whiptail --title "Select medium" --menu "Choose MEDIUM:" 15 60 5 \
  "mediaPlasma"       "Plasma-like medium" \
  "mediaAMEM"         "AMEM medium" \
  "mediaDMEM_highGLUC" "DMEM high glucose" \
  "mediaIMDM"         "IMDM medium" \
  "mediaRPMI"         "RPMI medium" 3>&1 1>&2 2>&3)
exitstatus=$?
if [ $exitstatus -ne 0 ]; then
  echo "Cancelled."
  exit 1
fi

# -------- MODEL (dropdown / menu) --------
MODEL=$(whiptail --title "Select model" --menu "Choose MODEL:" 12 60 2 \
  "humanGEM" "Human model" \
  "mouseGEM" "Mouse model" 3>&1 1>&2 2>&3)
exitstatus=$?
if [ $exitstatus -ne 0 ]; then
  echo "Cancelled."
  exit 1
fi

# -------- N_SIMS (numeric input, default 100) --------
N_SIMS=$(whiptail --inputbox "Enter N_SIMS (default 100):" 10 60 "100" 3>&1 1>&2 2>&3)
exitstatus=$?
if [ $exitstatus -ne 0 ]; then
  echo "Cancelled."
  exit 1
fi

# -------- CORES (numeric input, default 10) --------
CORES=$(whiptail --inputbox "Enter number of CORES (default 10):" 10 60 "10" 3>&1 1>&2 2>&3)
exitstatus=$?
if [ $exitstatus -ne 0 ]; then
  echo "Cancelled."
  exit 1
fi

# -------- Confirm and run --------
SUMMARY="ASSAY:  $ASSAY
MEDIUM: $MEDIUM
MODEL:  $MODEL
N_SIMS: $N_SIMS
CORES:  $CORES

Run imat_launcher.sh with these parameters?"

whiptail --yesno "$SUMMARY" 15 60
exitstatus=$?
if [ $exitstatus -ne 0 ]; then
  echo "Cancelled."
  exit 1
fi

# Call your existing launcher
echo "Running: bash imat_launcher.sh \"$ASSAY\" \"$MEDIUM\" \"$MODEL\" \"$N_SIMS\" \"$CORES\""
bash imat_launcher.sh "$ASSAY" "$MEDIUM" "$MODEL" "$N_SIMS" "$CORES"
