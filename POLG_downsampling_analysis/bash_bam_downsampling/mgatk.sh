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

source activate mgatk

export BAM=$1
export OUT=$2
export BCF=$3

mgatk tenx -i ${BAM} \
 -o ${OUT} \
 -bt CB -b ${BCF}

echo "done"
>&2 date
