#!/bin/bash

process() {
    RUN_DIR=$1
    REF_DATA=$2
    RAW_READS=$3
    NUM_CELLS=$4

    mkdir -p "$RUN_DIR"/{inputs,outputs}

    for file in "$RAW_READS"/*_R1_raw.fastq.gz; do
        sample=$(basename "${file%_R1_raw.fastq.gz}")
        ln -s "$RAW_READS/${sample}_R1_raw.fastq.gz" "$RUN_DIR/inputs/${sample}_S1_L001_R1_001.fastq.gz"
        ln -s "$RAW_READS/${sample}_R2_raw.fastq.gz" "$RUN_DIR/inputs/${sample}_S1_L001_R2_001.fastq.gz"
    done

    cd "$RUN_DIR/outputs" || exit
    for file in "$RUN_DIR/inputs"/*_S1_L001_R1_001.fastq.gz; do
        sample=$(basename "${file%_S1_L001_R1_001.fastq.gz}")
        cellranger count --id="$sample" \
                        --fastqs="$RUN_DIR/inputs" \
                        --sample="$sample" \
                        --transcriptome="$REF_DATA" \
                        --expect-cells="$NUM_CELLS" \
                        --create-bam=true \
                        --jobmode=local \
                        --localcores=16
    done
}

process "../data/osd402" "../data/refdata-gex-mm10-2020-A" "../data/osd402/raw" 8000
process "../data/osd403" "../data/refdata-gex-mm10-2020-A" "../data/osd403/raw" 8000
process "../data/osd404" "../data/refdata-gex-mm10-2020-A" "../data/osd404/raw" 8000
process "../data/osd405" "../data/refdata-gex-mm10-2020-A" "../data/osd405/raw" 8000
