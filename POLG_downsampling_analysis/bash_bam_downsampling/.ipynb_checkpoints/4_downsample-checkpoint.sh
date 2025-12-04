#!/bin/bash
#SBATCH --export=ALL
#SBATCH --job-name=downsample
#SBATCH --output=logs/%x-%j.log
#SBATCH --ntasks=8
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=16G
#SBATCH --partition=medium
#SBATCH --time=1-00:00:00

set -e
export TMPDIR=${HOME}/scratch/tmp
mkdir -p ${TMPDIR}

# Compute fractions (use median coverage)
export CTRL=56
export KIA2=255
export KI36=154

frac_KIA2=$(echo "$CTRL / $KIA2" | bc -l)
frac_KI36=$(echo "$CTRL / $KI36" | bc -l)

echo "KIA2 fraction: $frac_KIA2"
echo "KI36 fraction: $frac_KI36"

# Subsample
# KIA2
samtools view -bs 42${frac_KIA2} KIA2_filtered.bam > KIA2_downsampled.bam
# KI36
samtools view -bs 42${frac_KI36} KI36_filtered.bam > KI36_downsampled.bam

# Count
KIA2_COUNT=$(samtools view -c KIA2_downsampled.bam)
echo "KIA2 has $KIA2_COUNT reads"

KI36_COUNT=$(samtools view -c KI36_downsampled.bam)
echo "KI36 has $KI36_COUNT reads"

echo "done"
>&2 date
