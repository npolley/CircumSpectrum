#!/bin/bash

# Terminal-based GUI for expr_2_imat_launcher.sh using whiptail.
# If whiptail is not available, you can swap it with dialog and adjust options.


ASSAY_DIR="../../../../../data/sc_objects"

# Get unique filenames matching *_norm.csv
mapfile -t files < <(find "$ASSAY_DIR" -maxdepth 1 -type f -name "*.rds" -printf "%f\n" | sort -u)

if [ ${#files[@]} -eq 0 ]; then
  echo "No *.rds files found in $ASSAY_DIR"
  exit 1
fi

MENU_ITEMS=()
for f in "${files[@]}"; do
  base="${f%.rds}"      # strip suffix
  MENU_ITEMS+=("$base" "")  # only show the base; empty description
done

# -------- ASSAY (dropdown / menu) --------
ASSAY=$(whiptail --title "Select ASSAY" \
  --menu "Choose ASSAY:" \
  20 78 12 \
  "${MENU_ITEMS[@]}" \
  3>&1 1>&2 2>&3)

if [ $? -ne 0 ]; then
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

# -------- RESOLUTION (numeric input, default 40) --------
RES=$(whiptail --inputbox "Enter microclustering resolution (default 40):" 10 60 "40" 3>&1 1>&2 2>&3)
exitstatus=$?
if [ $exitstatus -ne 0 ]; then
  echo "Cancelled."
  exit 1
fi

# -------- Confirm and run --------
SUMMARY="ASSAY:  $ASSAY
MODEL:  $MODEL
RESOLUTION:  $RES


Run expr_2_imat_launcher.sh with these parameters?"

whiptail --yesno "$SUMMARY" 15 60
exitstatus=$?
if [ $exitstatus -ne 0 ]; then
  echo "Cancelled."
  exit 1
fi

# Call your existing launcher
echo "Running: bash expr_2_imat_sc_launcher.sh \"$ASSAY\" \"$MODEL\" \"$RES\""
Rscript expr_2_imat_sc.R "$ASSAY" "$MODEL" "$RES"
