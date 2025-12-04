#!/bin/bash
#SBATCH --export=ALL
#SBATCH --job-name=multiDS
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
export CTRL=42701668
export KIA2=258414081
export KI36=135496613

seed=42
outdir=~/scratch/bam_DS
mkdir -p "$outdir"

subsample() {
  local inbam=$1 total=$2 target_reads=$3 outbam=$4

  # sanity checks
  if (( target_reads <= 0 )); then
    echo "Target (${target_reads}) <= 0 → skipping ${inbam}"
    return
  fi
  if (( target_reads >= total )); then
    echo "Target (${target_reads}) >= total (${total}) → copying ${inbam} to ${outbam}"
    samtools view -b "$inbam" > "$outbam"
    return
  fi

  # fraction (6 decimals), ensure 0 < frac < 1
  frac=$(awk -v t="$target_reads" -v n="$total" 'BEGIN{printf("%.6f", t/n)}')
  frac_no0=${frac#0}  # make SEED.FRAC

  echo "Subsampling ${inbam} to ~${target_reads} reads (fraction ${frac}) → ${outbam}"
  samtools view -bs ${seed}${frac_no0} "$inbam" > "$outbam"
}

# ------ 10–50M for all three ------
for DEPTH_M in 10 20 30 40; do
  echo "Downsampling to ${DEPTH_M}M reads…"
  TARGET=$(( DEPTH_M * 1000000 ))

  subsample CTRL_filtered.bam  "$CTRL"  "$TARGET"  "${outdir}/CTRL_${DEPTH_M}M.bam"
  subsample KIA2_filtered.bam  "$KIA2"  "$TARGET"  "${outdir}/KIA2_${DEPTH_M}M.bam"
  subsample KI36_filtered.bam  "$KI36"  "$TARGET"  "${outdir}/KI36_${DEPTH_M}M.bam"
done

for DEPTH_M in 50 100; do
  echo "Downsampling to ${DEPTH_M}M reads…"
  TARGET=$(( DEPTH_M * 1000000 ))
  subsample KI36_filtered.bam "$KI36" "$TARGET" "${outdir}/KI36_${DEPTH_M}M.bam"
  subsample KIA2_filtered.bam "$KIA2" "$TARGET" "${outdir}/KIA2_${DEPTH_M}M.bam"
done

# ------ 150/200M for KIA2 only ------
for DEPTH_M in 150 200; do
  echo "Downsampling KIA2 to ${DEPTH_M}M reads…"
  TARGET=$(( DEPTH_M * 1000000 ))
  subsample KIA2_filtered.bam "$KIA2" "$TARGET" "${outdir}/KIA2_${DEPTH_M}M.bam"
done

# Optional: index outputs and verify counts
for f in "${outdir}"/*.bam; do
  samtools index -@ "$threads" "$f"
  c=$(samtools view -c "$f")
  echo "$(basename "$f"): $c"
done

echo "done"
>&2 date