#!/bin/bash

ASSAY="$1"
MODEL="$2"
RES="$3"

Rscript expr_2_imat_sc.R "$ASSAY" "$MODEL" "$RES"
