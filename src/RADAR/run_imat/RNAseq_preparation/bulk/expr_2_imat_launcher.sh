#!/bin/bash

ASSAY="$1"
MODEL="$2"

Rscript expr_2_imat.R "$ASSAY" "$MODEL"
