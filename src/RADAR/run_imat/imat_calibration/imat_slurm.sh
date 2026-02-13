#!/bin/bash
#SBATCH -J launchRscript
#SBATCH -t 0-6:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8        # can be overridden by sbatch
#SBATCH --mem-per-cpu=5G
#SBATCH --array=1-1              # placeholder; overridden by sbatch
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

module load statistics/R/4.3.0

ASSAY="$1"
MEDIUM="$2"
MODEL="$3"
N_SIMS="$4"
CORES="$5"

# ---------------- Retry logic ----------------
MAX_RETRIES=3
RETRY_FILE="logs/${ASSAY}_task${SLURM_ARRAY_TASK_ID}_retries.txt"

if [ -f "$RETRY_FILE" ]; then
  RETRIES=$(cat "$RETRY_FILE")
else
  RETRIES=0
fi

# Random seed for this attempt
SEED=$RANDOM

echo "Task $SLURM_ARRAY_TASK_ID, attempt $((RETRIES+1))/$MAX_RETRIES, seed=$SEED"

# Run R with the seed as last argument
Rscript imat.R "$SLURM_ARRAY_TASK_ID" "$ASSAY" "$MEDIUM" "$MODEL" "$SEED" "$N_SIMS" "$CORES" 
STATUS=$?

if [ $STATUS -ne 0 ]; then
  echo "Task $SLURM_ARRAY_TASK_ID failed with status $STATUS"

  if [ "$RETRIES" -lt "$MAX_RETRIES" ]; then
    RETRIES=$((RETRIES + 1))
    echo "$RETRIES" > "$RETRY_FILE"

    echo "Requeuing task $SLURM_ARRAY_TASK_ID (retry $RETRIES/$MAX_RETRIES)..."
    scontrol requeue "$SLURM_JOB_ID"   # requeue this array job/task [web:63][web:64][web:73]
    exit 0
  else
    echo "Max retries ($MAX_RETRIES) reached for task $SLURM_ARRAY_TASK_ID. Not requeuing."
    exit $STATUS
  fi
else
  rm -f "$RETRY_FILE"
  echo "Task $SLURM_ARRAY_TASK_ID completed successfully."
  exit 0
fi
