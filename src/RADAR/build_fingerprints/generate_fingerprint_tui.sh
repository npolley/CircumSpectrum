#!/bin/bash

# Terminal-based GUI for generate_fingerprint.sh using whiptail.
# If whiptail is not available, you can swap it with dialog and adjust options.


OBJ_ASSAY_DIR="../../../data/fingerprint_prep_objects/RADAR_objects/"
COHORT_DIR="../../../data/RADAR_xCheck_cohort/"
ASSAY_DIR="../../../data/RADAR_xCheck_experimental_assay/"

# Get unique filenames matching *.rds
mapfile -t files < <(find "$OBJ_ASSAY_DIR" -maxdepth 1 -type f -name "*.rds" -printf "%f\n" | sort -u)

if [ ${#files[@]} -eq 0 ]; then
  echo "No *.rds files found in $OBJ_ASSAY_DIR"
  exit 1
fi

MENU_ITEMS=()
for f in "${files[@]}"; do
  base="${f%.rds}"      # strip suffix
  MENU_ITEMS+=("$base" "")  # only show the base; empty description
done

# -------- OBJ (dropdown / menu) --------
OBJ=$(whiptail --title "Select FINGERPRINT PREPARATION OBJECT" \
  --menu "Choose FINGERPRINT OBJECT:" \
  20 78 12 \
  "${MENU_ITEMS[@]}" \
  3>&1 1>&2 2>&3)

if [ $? -ne 0 ]; then
  echo "Cancelled."
  exit 1
fi

# -------- FLUX (combined folder names from two dirs) --------
# Collect directory basenames from both locations
mapfile -t flux_dirs < <(
  {
    find "$COHORT_DIR" -maxdepth 1 -mindepth 1 -type d -printf "%f\n"
    find "$ASSAY_DIR"  -maxdepth 1 -mindepth 1 -type d -printf "%f\n"
  } 2>/dev/null | sort -u
)

if [ ${#flux_dirs[@]} -eq 0 ]; then
  echo "No flux directories found in $COHORT_DIR or $ASSAY_DIR"
  exit 1
fi

FLUX_MENU_ITEMS=()
for d in "${flux_dirs[@]}"; do
  FLUX_MENU_ITEMS+=("$d" "")
done

FLUX=$(whiptail --title "Select FLUX ASSAY" \
  --menu "Choose FLUX assay:" \
  20 78 12 \
  "${FLUX_MENU_ITEMS[@]}" \
  3>&1 1>&2 2>&3)

if [ $? -ne 0 ]; then
  echo "Cancelled."
  exit 1
fi

# -------- Confirm and run --------
SUMMARY="FINGERPRINT OBJECT:  $OBJ

FLUX:               $FLUX


Run generate_fingerprint.sh with these parameters?"

if ! whiptail --yesno "$SUMMARY" 15 60; then
  echo "Cancelled."
  exit 1
fi

# Call your existing launcher
echo "Running: generate_fingerprint.sh \"$OBJ\" \"$FLUX\""
bash generate_fingerprint.sh "$OBJ" "$FLUX"
