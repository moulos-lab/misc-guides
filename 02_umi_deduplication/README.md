UMI Deduplication from raw .fastq files to deduplicated .fastq (or .bam) files
================================================================================

<!-- TOC start (generated with https://github.com/derlin/bitdowntoc) -->

- [Prerequisites](#prerequisites)
   * [Software](#software)
   * [Data](#data)
- [Process Outline](#process-outline)
   * [UMI Extraction](#umi-extraction)
   * [Preprocessing and Quality Control](#preprocessing-and-quality-control)
   * [Mapping](#mapping)
   * [UMI Deduplication](#umi-deduplication)
   * [Generation of FastQ Files (Optional)](#generation-of-fastq-files-optional)
- [Detailed Usage Example with Scripts](#detailed-usage-example-with-scripts)
      + [Usage Note](#usage-note)
   * [UMI Extraction](#umi-extraction-1)
   * [Preprocessing and Quality Control](#preprocessing-and-quality-control-1)
   * [Mapping](#mapping-1)
   * [UMI Deduplication](#umi-deduplication-1)
   * [Generate valid .fastq files (Optional)](#generate-valid-fastq-files-optional)

<!-- TOC end -->

<!-- TOC --><a name="prerequisites"></a>
# Prerequisites

<!-- TOC --><a name="software"></a>
## Software

* [UMI-tools](https://github.com/CGATOxford/UMI-tools?tab=readme-ov-file#installation) ([manual](https://umi-tools.readthedocs.io/en/latest/))
* [TrimGalore](https://github.com/FelixKrueger/TrimGalore/archive/refs/tags/0.6.10.tar.gz) ([manual](https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md))
* [MultiQC](https://github.com/MultiQC/MultiQC#installation) ([manual](https://multiqc.info/docs/))
* [samtools](https://github.com/samtools/samtools/releases/download/1.19.2/samtools-1.19.2.tar.bz2) ([manual](https://www.htslib.org/doc/samtools.html))
* [HISAT2](https://cloud.biohpc.swmed.edu/index.php/s/oTtGWbWjaxsQ2Ho/download) ([manual](https://daehwankimlab.github.io/hisat2/manual/))
* [Bowtie 2](https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.4.2/bowtie2-2.4.2-sra-linux-x86_64.zip/download) ([manual](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml))

<!-- TOC --><a name="data"></a>
## Data

* Paired-end `.fastq` files obtained from UMI (Unique Molecular Identifier) sequencing.

UMI sequences are assumed to be the first 12 nucleotides on the forward read of each pair. The protocol will need to be modified slightly if this is not the case for your data.

* Hisat2 and Bowtie2 indices for your genome to deduplicate reads using UMI and mapping coordinates.

In this protocol we will use the UCSC GRCm38/mm10 mouse genome.

To get the `hisat2` index for mm10:

```
wget https://genome-idx.s3.amazonaws.com/hisat/mm10_genome.tar.gz
tar –xvzf mm10_genome.tar.gz
```

To get the `bowtie2` index for mm10:

```
wget https://genome-idx.s3.amazonaws.com/bt/mm10.zip
unzip mm10.zip
```

<!-- TOC --><a name="process-outline"></a>
# Process Outline

<!-- TOC --><a name="umi-extraction"></a>
## UMI Extraction

1. Move the UMI sequence from the read sequence to the read header using `umi_tools extract`.

<!-- TOC --><a name="preprocessing-and-quality-control"></a>
## Preprocessing and Quality Control

2. Preprocess and run quality control with `trim_galore`.
3. Rename `trim_galore` output files.
4. Aggregate quality control reports for all samples with `multiqc`.
5. Clean up intermediate files.

<!-- TOC --><a name="mapping"></a>
## Mapping

6. Map the reads to the reference genome with `hisat2`, retaining the unmapped read pairs separately in `.fastq` files for the next step.
7. Remap the unmapped read pairs with `bowtie2`.
8. Convert the Hisat2 `.sam` output file to `.bam` format with `samtools view -bhS`.
9. Merge the Hisat2 and Bowtie2 `.bam` output files into a single file with `samtools merge`.
10. Sort the merged output file by coordinates for indexing with `samtools sort`.
11. Index final `.bam` file with `samtools index`.
12. Clean up intermediate files.

<!-- TOC --><a name="umi-deduplication"></a>
## UMI Deduplication

13. Filter `.bam` file to only keep primary alignments and properly-paired reads with `samtools view -hb -f2 -F2304`.
14. Index the filtered `.bam` file with `samtools index`.
15. Deduplicate reads with `umi_tools dedup`.

<!-- TOC --><a name="generation-of-fastq-files-optional"></a>
## Generation of FastQ Files (Optional)

16. Re-order `.bam` file to group reads by name with `samtools collate`.
17. Generate the `.fastq` files with `samtools fastq`.

<!-- TOC --><a name="detailed-usage-example-with-scripts"></a>
# Detailed Usage Example with Scripts

<!-- TOC --><a name="usage-note"></a>
### Usage Note

This pipeline starts with raw paired-end `.fastq` files from UMI sequencing.

Replace `HOME_PATH=<PROJECT_DIR>`, `HISAT2_INDEX=<HISAT2_INDEX_DIR>`, `BOWTIE2_INDEX=<BOWTIE2_INDEX_DIR>` in the scripts to run.

Also, pay attention to `.fastq` file naming patterns.

<!-- TOC --><a name="umi-extraction-1"></a>
## UMI Extraction

Extract the UMI sequence from the read sequence to the read header. Modify `BCPAT` regex if needed.

```{bash}
#!/bin/bash

HOME_PATH=<PROJECT_DIR>
FASTQ_PATH=$HOME_PATH/fastq
FASTQ_OUTPATH=$HOME_PATH/fastq_umi

if [ ! -d $FASTQ_OUTPATH ]
then
    mkdir -p $FASTQ_OUTPATH
fi

UMI_TOOLS_COMMAND=`which umi_tools`

BCPAT="NNNNNNNNNNNN"

for FILE in `ls $FASTQ_PATH/*_1.fastq.gz`
do
    SAMPLE=`basename $FILE | sed s/_1\.fastq\.gz//`
    echo "===== Processing $SAMPLE..."
    
    $UMI_TOOLS_COMMAND extract \
      --bc-pattern=$BCPAT \
      --stdin $FASTQ_PATH/$SAMPLE"_1.fastq.gz" \
      --stdout $FASTQ_OUTPATH/$SAMPLE"_1.fastq.gz" \
      --read2-in $FASTQ_PATH/$SAMPLE"_2.fastq.gz" \
      --read2-out=$FASTQ_OUTPATH/$SAMPLE"_2.fastq.gz" \
      --log=$FASTQ_OUTPATH/$SAMPLE".log" \
      --ignore-read-pair-suffixes &
done
```

<!-- TOC --><a name="preprocessing-and-quality-control-1"></a>
## Preprocessing and Quality Control

Trim adapter sequences, Ns and low-quality base calls, filter short reads and run quality control with `trim_galore`. Aggregate quality control reports with `multiqc`.

```{bash}
#!/bin/bash

HOME_PATH=<PROJECT_DIR>

TRIMGALORE_COMMAND=`which trim_galore`
CUTADAPT_COMMAND=`which cutadapt`
MULTIQC_COMMAND=`which multiqc`

FASTQ_PATTERN=*.fastq.gz

CORES=4

FASTQ_PATH=$HOME_PATH/fastq_umi

TRIMGALORE_OUTPUT=$HOME_PATH/fastq_umi_qual
if [ ! -d $TRIMGALORE_OUTPUT ]
then
    mkdir -p $TRIMGALORE_OUTPUT
fi

for FILE in $FASTQ_PATH/*_1.fastq.gz
do
    BASE=`basename $FILE | sed s/_1\.fastq\.gz//`
    F1=$FASTQ_PATH/$BASE"_1.fastq.gz"
    F2=$FASTQ_PATH/$BASE"_2.fastq.gz"
    $TRIMGALORE_COMMAND \
        --quality 20 \
        --length 50 \
        --output_dir $TRIMGALORE_OUTPUT \
        --path_to_cutadapt $CUTADAPT_COMMAND \
        --paired \
        --fastqc \
        --trim-n $F1 $F2 &
done

wait

for FILE in $FASTQ_PATH/*_1.fastq.gz
do
    BASE=`basename $FILE | sed s/_1\.fastq\.gz//`
    mv $TRIMGALORE_OUTPUT/$BASE"_1_val_1.fq.gz" \
        $TRIMGALORE_OUTPUT/$BASE"_1.fastq.gz"
    mv $TRIMGALORE_OUTPUT/$BASE"_2_val_2.fq.gz" \
        $TRIMGALORE_OUTPUT/$BASE"_2.fastq.gz"
    mv $TRIMGALORE_OUTPUT/$BASE"_1_val_1_fastqc.html" \
        $TRIMGALORE_OUTPUT/$BASE"_1_fastqc.html"
    mv $TRIMGALORE_OUTPUT/$BASE"_1_val_1_fastqc.zip" \
        $TRIMGALORE_OUTPUT/$BASE"_1_fastqc.zip"
    mv $TRIMGALORE_OUTPUT/$BASE"_2_val_2_fastqc.html" \
        $TRIMGALORE_OUTPUT/$BASE"_2_fastqc.html"
    mv $TRIMGALORE_OUTPUT/$BASE"_2_val_2_fastqc.zip" \
        $TRIMGALORE_OUTPUT/$BASE"_2_fastqc.zip"
done

$MULTIQC_COMMAND $TRIMGALORE_OUTPUT -o $TRIMGALORE_OUTPUT
rm -r $TRIMGALORE_OUTPUT/multiqc_data
```

<!-- TOC --><a name="mapping-1"></a>
## Mapping

First map with HISAT2, keeping all the read pairs that do not map properly in separate `.fastq` files, so they can be remapped with Bowtie2. Then merge the results, index the final `.bam` file and cleanup intermediates.

```{bash}
#!/bin/bash

HOME_PATH=<PROJECT_DIR>

HISAT2_COMMAND=`which hisat2`
BOWTIE2_COMMAND=`which bowtie2`
SAMTOOLS_COMMAND=`which samtools`

HISAT2_INDEX=<HISAT2_INDEX_DIR> 
BOWTIE2_INDEX=<BOWTIE2_INDEX_DIR>

CORES=32

FASTQ_PATH=$HOME_PATH/fastq_umi_qual
BAM_PATH=$HOME_PATH/bam_umi

if [ ! -d $BAM_PATH ]
then
    mkdir -p $BAM_PATH
fi

for FILE in `ls $FASTQ_PATH/*_1.fastq.gz`
do
    SAMPLE=`basename $FILE | sed s/_1\.fastq\.gz//`
    echo "===== Mapping with hisat2 for $SAMPLE..."
    $HISAT2_COMMAND -p $CORES \
        --un-conc $BAM_PATH/unmapped.fastq \
        -x $HISAT2_INDEX \
        -1 $FASTQ_PATH/$SAMPLE"_1.fastq.gz" \
        -2 $FASTQ_PATH/$SAMPLE"_2.fastq.gz" \
        -S $BAM_PATH/hisat2.sam \
        --no-unal \
        --no-mixed \
        --no-discordant
    echo " "
    
    echo "===== Trying to map unmapped reads with bowtie2 for $SAMPLE..."
    $BOWTIE2_COMMAND --local \
        --very-sensitive-local --dovetail \
        -p $CORES -x $BOWTIE2_INDEX \
        -1 $BAM_PATH/unmapped.1.fastq \
        -2 $BAM_PATH/unmapped.2.fastq | \
        $SAMTOOLS_COMMAND view -bhS \
        -o $BAM_PATH/unmapped_remap.uns -
    
    echo "===== Merging all reads for $SAMPLE..." 
    $SAMTOOLS_COMMAND view \
        -bhS $BAM_PATH/hisat2.sam > \
        $BAM_PATH/hisat2.bam
    $SAMTOOLS_COMMAND merge \
        -f $BAM_PATH/$SAMPLE".tmb" \
        $BAM_PATH/hisat2.bam \
        $BAM_PATH/unmapped_remap.uns
    
    echo "===== Coordinate sorting all reads for $SAMPLE..."
    $SAMTOOLS_COMMAND sort $BAM_PATH/$SAMPLE".tmb" > $BAM_PATH/$SAMPLE".bam"
    $SAMTOOLS_COMMAND index $BAM_PATH/$SAMPLE".bam"
    
    echo "===== Removing intermediate garbage for $SAMPLE..."
    rm $BAM_PATH/hisat2.bam \
        $BAM_PATH/hisat2.sam \
        $BAM_PATH/unmapped*.fastq \
        $BAM_PATH/unmapped_remap.uns \
        $BAM_PATH/$SAMPLE".tmb"
    echo " "
done
```

<!-- TOC --><a name="umi-deduplication-1"></a>
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

<!-- TOC --><a name="generate-valid-fastq-files-optional"></a>
## Generate valid .fastq files (Optional)

Re-order `.bam` file to group reads by name, then generate the `.fastq` files. Optionally, also generate a singletons file using `-0` to capture any errors as this file should be empty. Alternatively run `samtools flagstat` on the `.bam` file to see if `read1` and `read2` flags are equal.

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
