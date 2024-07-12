---
title: "UMI Deduplication to FastQ files"
author: "Stelios Gkionis"
date: "`r Sys.Date()`"
output: 
    rmdformats::robobook:
      fig_width: 7
      fig_height: 7
      fig_caption: true
      highlight: kate
      lightbox: true
      cards: true
      use_bookdown: false
      mathjax: null
      self_contained: true
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(
  echo = TRUE,
  eval = FALSE
  )
```

### Usage note

Replace `HOME_PATH=<PROJECT_DIR>` with your working directory.

## UMI Deduplication

Filter aligned `.bam` file to only keep primary alignments and properly-paired reads. This is required for correctly reconstructing the read pairs in the final `.fastq` files. Index the filtered `.bam` file and deduplicate using `umi_tools dedup`.

See [this section](https://umi-tools.readthedocs.io/en/latest/reference/dedup.html#unmapped-reads) in the documentation for more info, as the behavior of the algorithm and of the `--unmapped-reads --unpaired-reads --chimeric-pairs` flags is not very intuitive.

```{bash}
#!/bin/bash

HOME_PATH=<PROJECT_DIR>
BAM_PATH=$HOME_PATH/bam_umi
BAM_FILTERED=$HOME_PATH/bam_umi_filtered
BAM_OUTPATH=$HOME_PATH/bam_umi_filtered_dedup

if [ ! -d $BAM_FILTERED ]
then
    mkdir -p $BAM_FILTERED
fi
if [ ! -d $BAM_OUTPATH ]
then
    mkdir -p $BAM_OUTPATH
fi

UMI_TOOLS_COMMAND=`which umi_tools`
SAMTOOLS_COMMAND=`which samtools`

for FILE in `ls $BAM_PATH/*.bam`
do
    SAMPLE=`basename $FILE | sed s/\.bam//`
    echo "===== Processing $SAMPLE..."
    
    $SAMTOOLS_COMMAND view \
      -hb -f2 -F2304 \
      $BAM_PATH/$SAMPLE".bam" \
      > $BAM_FILTERED/$SAMPLE".bam"
    $SAMTOOLS_COMMAND index $BAM_FILTERED/$SAMPLE".bam"
    $UMI_TOOLS_COMMAND dedup \
      --output-stats $BAM_OUTPATH/$SAMPLE \
      --stdin $BAM_FILTERED/$SAMPLE".bam" \
      --stdout $BAM_OUTPATH/$SAMPLE".bam" \
      --unmapped-reads=discard --unpaired-reads=discard --chimeric-pairs=discard \
      --log $BAM_OUTPATH/$SAMPLE"_dedup.log" \
      --paired &
done

wait

for FILE in `ls $BAM_OUTPATH/*.bam`
do
    $SAMTOOLS_COMMAND index $FILE &
done
```

## Generate valid .fastq files

Re-order `.bam` file to group reads by name, then generate the `.fastq` files. Optionally, also generate a singletons file using `-0` to capture any errors.

```{bash}
#!/bin/bash

HOME_PATH=<PROJECT_DIR>
BAM_PATH=$HOME_PATH/bam_umi_filtered_dedup
FASTQ_OUTPATH=$HOME_PATH/fastq_umi_filtered_dedup

if [ ! -d $FASTQ_OUTPATH ]
then
    mkdir -p $FASTQ_OUTPATH
fi

SAMTOOLS_COMMAND=`which samtools`

for FILE in `ls $BAM_PATH/*.bam`
do
    SAMPLE=`basename $FILE | sed s/\.bam//`
    echo "===== Processing $SAMPLE..."

    $SAMTOOLS_COMMAND collate -f -O $BAM_PATH/$SAMPLE".bam" | \
    $SAMTOOLS_COMMAND fastq -n \
          -1 $FASTQ_OUTPATH/$SAMPLE"_1.fastq.gz" \
          -2 $FASTQ_OUTPATH/$SAMPLE"_2.fastq.gz" \
          -0 $FASTQ_OUTPATH/$SAMPLE"_0.fastq.gz" &
done
```
