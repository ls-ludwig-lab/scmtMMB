#!/bin/bash
#SBATCH --export=ALL
#SBATCH --job-name=mgatk
#SBATCH --output=logs/%x-%j.log
#SBATCH --ntasks=8
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=16G
#SBATCH --partition=medium
#SBATCH --time=3-00:00:00

set -e
export TMPDIR=${HOME}/scratch/tmp
mkdir -p ${TMPDIR}

# Barcode files for each sample
declare -A BARCODE_TSV
BARCODE_TSV[CTRL]="CTRL_CB_mgatk.tsv"
BARCODE_TSV[KI36]="KI36_CB_mgatk.tsv"
BARCODE_TSV[KIA2]="KIA2_CB_mgatk.tsv"

# Loop over samples
for SAMPLE in CTRL KI36 KIA2; do
  bcfile=${BARCODE_TSV[$SAMPLE]}

  for BAM in ${SAMPLE}_*.bam; do
    [ -f "$BAM" ] || continue
    prefix=${BAM%.bam}

    echo "Submitting $BAM ..."
    sbatch mgatk.sh "$BAM" "$prefix" "$bcfile"
  done
done

echo "done"
>&2 date

