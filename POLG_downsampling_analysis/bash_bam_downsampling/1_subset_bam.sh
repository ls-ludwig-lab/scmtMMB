#!/bin/bash
#SBATCH --export=ALL
#SBATCH --job-name=subam
#SBATCH --output=logs/%x-%j.log
#SBATCH --ntasks=8
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=16G
#SBATCH --partition=medium
#SBATCH --time=1-00:00:00

set -e
export TMPDIR=${HOME}/scratch/tmp
mkdir -p ${TMPDIR}

export OUTDIR=/data/cephfs-1/home/users/yhsieh_m/work/POLG_fastq
export SAMPLE=HEK_POLG_ATAC_hg38_v20-mtMask
export BAM_FILE="${OUTDIR}/${SAMPLE}/outs/possorted_bam.bam"

# extract reads map to chrM
echo "extract chrM reads ..."
samtools view -b $BAM_FILE chrM > chrM.bam

# Save the header lines
samtools view -H chrM.bam > SAM_header

# loop through each filter txt file
for FILTER in CTRL_CB.txt KIA2_CB.txt KI36_CB.txt
do
    # get prefix without .txt
    PREFIX=$(basename "$FILTER" _CB.txt)

    echo "Processing $PREFIX ..."

    # filter alignments
    samtools view chrM.bam | LC_ALL=C grep -F -f "$FILTER" > "${PREFIX}_filtered_SAM_body"

    # combine header and filtered body
    cat SAM_header "${PREFIX}_filtered_SAM_body" > "${PREFIX}_filtered.sam"

    # convert to BAM
    samtools view -b "${PREFIX}_filtered.sam" > "${PREFIX}_filtered.bam"

    # remove sam file 
    rm -rf "${PREFIX}_filtered.sam"
    rm -rf *SAM_body

    echo "Output: ${PREFIX}_filtered.bam"
    
done

echo "done"
>&2 date
