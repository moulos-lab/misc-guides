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
   * [Generation of FastQ Files (Optional)](#generate-valid-fastq-files-optional)
   * [Stranded genome tracks (Optional)](#stranded-genome-tracks)
- [Detailed Usage Example with Scripts](#detailed-usage-example-with-scripts)
      + [Usage Note](#usage-note)
   * [UMI Extraction](#umi-extraction-1)
   * [Preprocessing and Quality Control](#preprocessing-and-quality-control-1)
   * [Mapping](#mapping-1)
   * [UMI Deduplication](#umi-deduplication-1)
   * [Generate valid .fastq files (Optional)](#generate-valid-fastq-files-optional-1)
   * [Stranded genome tracks (Optional)](#stranded-genome-tracks-1)


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
* [Kent tools](https://github.com/ucscGenomeBrowser/kent/archive/refs/tags/v492_branch.1.tar.gz)([manual](https://github.com/ucscGenomeBrowser/kent))
* [bedtools](https://github.com/arq5x/bedtools2/releases/download/v2.31.1/bedtools-2.31.1.tar.gz)([manual](https://bedtools.readthedocs.io/en/latest/))

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

<!-- TOC --><a name="generate-valid-fastq-files-optional"></a>
## Generation of FastQ Files (Optional)

16. Re-order `.bam` file to group reads by name with `samtools collate`.
17. Generate the `.fastq` files with `samtools fastq`.

<!-- TOC --><a name="stranded-genome-tracks"></a>
## Stranded genome tracks (Optional)

18. Generation of `.bigWig` files for UCSC Genome Browser (Adapted from Alex Galaras and Panagiotis Moulos)

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
set -euo pipefail

HOME_PATH=<PROJECT_DIR>
UMI_TOOLS_COMMAND=$(which umi_tools)

BCPAT="NNNNNNNNNNNN"

# Loop over sample directories
for SAMPLE_DIR in "$HOME_PATH"/*/; do
    FASTQ_PATH="$SAMPLE_DIR/fastq"
    FASTQ_OUTPATH="$SAMPLE_DIR/fastq_umi"

    # Skip if fastq directory does not exist
    if [ ! -d "$FASTQ_PATH" ]; then
        echo "Skipping $SAMPLE (no fastq directory)"
        continue
    fi
    if [ ! -d $FASTQ_OUTPATH ]
    then
        mkdir -p $FASTQ_OUTPATH
    fi

    SAMPLE=$(basename "$SAMPLE_DIR")

    # Loop over read 1 files
    for FILE in "$FASTQ_PATH"/*_1.fq.gz; do
        [ -e "$FILE" ] || continue

        BASENAME=$(basename "$FILE" _1.fq.gz)

        echo "===== Processing $SAMPLE / $BASENAME..."

        "$UMI_TOOLS_COMMAND" extract \
            --bc-pattern="$BCPAT" \
            --stdin "$FASTQ_PATH/${BASENAME}_1.fq.gz" \
            --stdout "$FASTQ_OUTPATH/${BASENAME}_1.fq.gz" \
            --read2-in "$FASTQ_PATH/${BASENAME}_2.fq.gz" \
            --read2-out="$FASTQ_OUTPATH/${BASENAME}_2.fq.gz" \
            --log="$FASTQ_OUTPATH/${BASENAME}.log" \
            --ignore-read-pair-suffixes &
    done
done

wait
```

<!-- TOC --><a name="preprocessing-and-quality-control-1"></a>
## Preprocessing and Quality Control

Trim adapter sequences, Ns and low-quality base calls, filter short reads and run quality control with `trim_galore`. Aggregate quality control reports with `multiqc`.

```{bash}
#!/bin/bash
set -euo pipefail

HOME_PATH=<PROJECT_DIR>

TRIMGALORE_COMMAND=`which trim_galore`
CUTADAPT_COMMAND=`which cutadapt`
MULTIQC_COMMAND=`which multiqc`

CORES=1

# Loop over sample directories
for SAMPLE_DIR in "$HOME_PATH"/*/; do
    FASTQ_PATH="$SAMPLE_DIR/fastq_umi"
    TRIMGALORE_OUTPUT="$SAMPLE_DIR/fastq_umi_qual"

    # Skip if fastq_umi directory does not exist
    if [ ! -d "$FASTQ_PATH" ]; then
        echo "Skipping $SAMPLE (no fastq_umi directory)"
        continue
    fi
    if [ ! -d $TRIMGALORE_OUTPUT ]
    then
        mkdir -p $TRIMGALORE_OUTPUT
    fi

    SAMPLE=$(basename "$SAMPLE_DIR")

    for FILE in $FASTQ_PATH/*_1.fq.gz; do
        BASE=`basename $FILE | sed s/_1\.fq\.gz//`
        F1=$FASTQ_PATH/$BASE"_1.fq.gz"
        F2=$FASTQ_PATH/$BASE"_2.fq.gz"
        $TRIMGALORE_COMMAND \
            --quality 20 \
            --length 50 \
            --output_dir $TRIMGALORE_OUTPUT \
            --path_to_cutadapt $CUTADAPT_COMMAND \
            --paired \
            --fastqc \
            --cores $CORES \
            --trim-n $F1 $F2 &
    done
done

wait

for SAMPLE_DIR in "$HOME_PATH"/*/; do
    FASTQ_PATH="$SAMPLE_DIR/fastq_umi"
    TRIMGALORE_OUTPUT="$SAMPLE_DIR/fastq_umi_qual"

    # Skip if fastq_umi directory does not exist
    if [ ! -d "$FASTQ_PATH" ]; then
        echo "Skipping $SAMPLE (no fastq_umi directory)"
        continue
    fi

    SAMPLE=$(basename "$SAMPLE_DIR")

    for FILE in $FASTQ_PATH/*_1.fastq.gz; do
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
done

wait
```

<!-- TOC --><a name="mapping-1"></a>
## Mapping

First map with HISAT2, keeping all the read pairs that do not map properly in separate `.fastq` files, so they can be remapped with Bowtie2. Then merge the results, index the final `.bam` file and cleanup intermediates.

```{bash}
#!/bin/bash
set -euo pipefail

HOME_PATH=<PROJECT_DIR>

HISAT2_COMMAND=`which hisat2`
BOWTIE2_COMMAND=`which bowtie2`
SAMTOOLS_COMMAND=`which samtools`

HISAT2_INDEX=<HISAT2_INDEX_DIR> 
BOWTIE2_INDEX=<BOWTIE2_INDEX_DIR>

CORES=32

# Loop over sample directories
for SAMPLE_DIR in "$HOME_PATH"/*/; do
    FASTQ_PATH="$SAMPLE_DIR/fastq_umi_qual"
    BAM_PATH="$SAMPLE_DIR/bam_umi"

    # Skip if fastq_umi directory does not exist
    if [ ! -d "$FASTQ_PATH" ]; then
        echo "Skipping $SAMPLE (no fastq_umi directory)"
        continue
    fi
    if [ ! -d $BAM_PATH ]
    then
        mkdir -p $BAM_PATH
    fi

    SAMPLE=$(basename "$SAMPLE_DIR")

    for FILE in `ls $FASTQ_PATH/*_1.fastq.gz`; do
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
        $SAMTOOLS_COMMAND view -@ $CORES \
            -bhS $BAM_PATH/hisat2.sam > \
            $BAM_PATH/hisat2.bam
        $SAMTOOLS_COMMAND merge -@ $CORES \
            -f $BAM_PATH/$SAMPLE".tmb" \
            $BAM_PATH/hisat2.bam \
            $BAM_PATH/unmapped_remap.uns

        echo "===== Coordinate sorting all reads for $SAMPLE..."
        $SAMTOOLS_COMMAND sort -@ $CORES $BAM_PATH/$SAMPLE".tmb" > $BAM_PATH/$SAMPLE".bam"
        $SAMTOOLS_COMMAND index -@ $CORES $BAM_PATH/$SAMPLE".bam"
      
        echo "===== Removing intermediate garbage for $SAMPLE..."
        rm $BAM_PATH/hisat2.bam \
            $BAM_PATH/hisat2.sam \
            $BAM_PATH/unmapped*.fastq \
            $BAM_PATH/unmapped_remap.uns \
            $BAM_PATH/$SAMPLE".tmb"
        echo " "
    done
done

wait
```

<!-- TOC --><a name="umi-deduplication-1"></a>
## UMI Deduplication

Filter aligned `.bam` file to only keep primary alignments and properly-paired reads. This is required for correctly reconstructing the read pairs in the final `.fastq` files. Index the filtered `.bam` file and deduplicate using `umi_tools dedup`.

See [this section](https://umi-tools.readthedocs.io/en/latest/reference/dedup.html#unmapped-reads) in the documentation for more info, as the behavior of the algorithm and of the `--unmapped-reads --unpaired-reads --chimeric-pairs` flags is not very intuitive.

```{bash}
#!/bin/bash
set -euo pipefail

HOME_PATH=<PROJECT_DIR>

UMI_TOOLS_COMMAND=`which umi_tools`
SAMTOOLS_COMMAND=`which samtools`

# Loop over sample directories
for SAMPLE_DIR in "$HOME_PATH"/*/; do
    BAM_PATH=$SAMPLE_DIR/bam_umi
    BAM_FILTERED=$SAMPLE_DIR/bam_umi_filtered
    BAM_OUTPATH=$SAMPLE_DIR/bam_umi_filtered_dedup

    # Skip if fastq_umi directory does not exist
    if [ ! -d "$BAM_PATH" ]; then
        echo "Skipping $SAMPLE (no bam_umi directory)"
        continue
    fi
    if [ ! -d $BAM_FILTERED ]
    then
        mkdir -p $BAM_FILTERED
    fi
    if [ ! -d $BAM_OUTPATH ]
    then
        mkdir -p $BAM_OUTPATH
    fi

    SAMPLE=$(basename "$SAMPLE_DIR")

    for FILE in `ls $BAM_PATH/*.bam`; do
        FILENAME=`basename $FILE | sed s/\.bam//`
        echo "===== Processing $SAMPLE..."
      
        $SAMTOOLS_COMMAND view -@ 48 \
            -hb -f2 -F2304 \
            $BAM_PATH/$FILENAME".bam" \
            > $BAM_FILTERED/$FILENAME".bam"
        $SAMTOOLS_COMMAND index -@ 48 $BAM_FILTERED/$FILENAME".bam"
        $UMI_TOOLS_COMMAND dedup \
            --output-stats $BAM_OUTPATH/$FILENAME \
            --stdin $BAM_FILTERED/$FILENAME".bam" \
            --stdout $BAM_OUTPATH/$FILENAME".bam" \
            --unmapped-reads=discard --unpaired-reads=discard --chimeric-pairs=discard \
            --log $BAM_OUTPATH/$FILENAME"_dedup.log" \
            --paired &
    done
done

wait

for SAMPLE_DIR in "$HOME_PATH"/*/; do
    BAM_PATH=$SAMPLE_DIR/bam_umi
    BAM_FILTERED=$SAMPLE_DIR/bam_umi_filtered
    BAM_OUTPATH=$SAMPLE_DIR/bam_umi_filtered_dedup

    # Skip if fastq_umi directory does not exist
    if [ ! -d "$BAM_PATH" ]; then
        echo "Skipping $SAMPLE (no bam_umi directory)"
        continue
    fi
  
    SAMPLE=$(basename "$SAMPLE_DIR")

    for FILE in `ls $BAM_OUTPATH/*.bam`; do
        $SAMTOOLS_COMMAND index $FILE &
    done
done

wait
```

<!-- TOC --><a name="generate-valid-fastq-files-optional-1"></a>
## Generate valid .fastq files (Optional)

Re-order `.bam` file to group reads by name, then generate the `.fastq` files. Optionally, also generate a singletons file using `-0` to capture any errors as this file should be empty. Alternatively run `samtools flagstat` on the `.bam` file to see if `read1` and `read2` flags are equal.

```{bash}
#!/bin/bash
set -euo pipefail

HOME_PATH=<PROJECT_DIR>

SAMTOOLS_COMMAND=`which samtools`

# Loop over sample directories
for SAMPLE_DIR in "$HOME_PATH"/*/; do
    BAM_PATH=$SAMPLE_DIR/bam_umi_filtered_dedup
    FASTQ_OUTPATH=$SAMPLE_DIR/fastq_umi_filtered_dedup

    # Skip if fastq_umi directory does not exist
    if [ ! -d "$BAM_PATH" ]; then
        echo "Skipping $SAMPLE (no bam_umi directory)"
        continue
    fi
    if [ ! -d $FASTQ_OUTPATH ]
    then
        mkdir -p $FASTQ_OUTPATH
    fi

    SAMPLE=$(basename "$SAMPLE_DIR")

    for FILE in `ls $BAM_PATH/*.bam`; do
        SAMPLE=`basename $FILE | sed s/\.bam//`
        echo "===== Processing $SAMPLE..."

        $SAMTOOLS_COMMAND collate -f -O $BAM_PATH/$SAMPLE".bam" | \
        $SAMTOOLS_COMMAND fastq -n \
            -1 $FASTQ_OUTPATH/$SAMPLE"_1.fastq.gz" \
            -2 $FASTQ_OUTPATH/$SAMPLE"_2.fastq.gz" \
            -0 $FASTQ_OUTPATH/$SAMPLE"_0.fastq.gz" &
    done
done

wait
```

<!-- TOC --><a name="stranded-genome-tracks-1"></a>

## Stranded genome tracks (Optional)

Create stranded genome tracks for use with the UCSC Genome Browser in a Genome Hub hosted publically on the web. Use the chromosome length file for the genome used in the mapping step or the default chromInfo.txt file with any modifications necessary. Then create a hub directory with the required `hub.txt`, `trackDb.txt` and `genomes.txt` files. The `trackDb.txt` file is generated automatically from this script. Follow the documentation found online [here](https://genome.ucsc.edu/goldenpath/help/trackDb/trackDbHub.html) if needed.

```{bash}
#!/bin/bash
set -euo pipefail

HOME_PATH=<PROJECT_DIR>

LINK=<PUBLIC_URL>
# Required for Samples starting with numbers as is common in our Facility
PREFIX=Raw
GENOME=<chromInfo_FILE>
NORMALIZE=/media/samba/hatzis_lab/general_scripts/normalize_bedgraph.pl

CORES=32

SAMTOOLS=$(command -v samtools)
BEDTOOLS=$(command -v bedtools)
KENTTOOLS=$(command -v bedGraphToBigWig)

for SAMPLE_DIR in "$HOME_PATH"/*/; do
    BAM_PATH=$SAMPLE_DIR/bam_umi

    # Skip if bam_umi directory does not exist
    if [ ! -d "$BAM_PATH" ]; then
        echo "Skipping $SAMPLE (no bam_umi directory)"
        continue
    fi

    SAMPLE=$(basename "$SAMPLE_DIR")
    echo "========== Processing $SAMPLE"

    for FILE in $BAM_PATH/*.bam; do
        FILENAME=`basename $FILE | sed s/\.bam//`

        #### First mate
        # First mate - forward strand
        # Flags: -f 64 (first read in pair), -F 16 (exclude reverse strand reads)
        $SAMTOOLS view -@ $CORES -bh -f64 -F16 $FILE > $BAM_PATH/$FILENAME"_mate1_forward.bam"    # Create BAM for first read in pair mapped to forward strand

        # First mate - reverse strand
        # Flags: -f 80 (first read in pair & reverse strand)
        $SAMTOOLS view -bh -@ $CORES -f80 $FILE > $BAM_PATH/$FILENAME"_mate1_reverse.bam"   # Create BAM for second read in pair mapped to forward strand

        #### Second mate
        # Second mate - forward strand
        # Flags: -f 128 (second read in pair), -F 16 (exclude reverse strand reads)
        $SAMTOOLS view -bh -@ $CORES -f128 -F16 $FILE > $BAM_PATH/$FILENAME"_mate2_forward.bam"           # Create BAM for second read in pair mapped to forward strand

        # Second mate - reverse strand
        # Flags: -f 144 (second read in pair, mapped to reverse strand)
        $SAMTOOLS view -bh -@ $CORES -f144 $FILE > $BAM_PATH/$FILENAME"_mate2_reverse.bam"               # Create BAM for second read in pair mapped to reverse strand

        #### Merge
        $SAMTOOLS merge -@ $CORES $BAM_PATH/$FILENAME"_for.bam" $BAM_PATH/$FILENAME"_mate1_forward.bam" $BAM_PATH/$FILENAME"_mate2_reverse.bam"
        $SAMTOOLS merge -@ $CORES $BAM_PATH/$FILENAME"_rev.bam" $BAM_PATH/$FILENAME"_mate1_reverse.bam" $BAM_PATH/$FILENAME"_mate2_forward.bam"

        #### Tracks
        $BEDTOOLS genomecov -bg -split -ibam $BAM_PATH/$FILENAME"_for.bam" | grep -vP 'chrJH|chrMT|chrGL|_GL|_JH'  | sort -k1,1 -k2g,2 > $BAM_PATH/$FILENAME".plus.bedGraph" 

        $BEDTOOLS genomecov -bg -split -ibam $BAM_PATH/$FILENAME"_rev.bam" | grep -vP 'chrJH|chrMT|chrGL|_GL|_JH' | sort -k1,1 -k2g,2 | awk '{ print $1"\t"$2"\t"$3"\t-"$4 }' > $BAM_PATH/$FILENAME".minus.bedGraph" 

        #### Normalize
        perl $NORMALIZE --input $BAM_PATH/$FILENAME".plus.bedGraph" --sumto 500000000 --exportfactors $BAM_PATH/$FILENAME"_p_normfactors.txt" --ncores 8
        perl $NORMALIZE --input $BAM_PATH/$FILENAME".minus.bedGraph" --sumto -500000000 --exportfactors $BAM_PATH/$FILENAME"_m_normfactors.txt" --ncores 8

        #### Create bigwig files for every bedGraph
        for FILE in $BAM_PATH/${FILENAME}.*_norm.bedGraph; do
            BASE=`basename $FILE | sed s/\.bedGraph//`
            $KENTTOOLS $FILE $GENOME $BAM_PATH/$BASE".bigWig"
        done

        for FILE in $BAM_PATH/${FILENAME}.plus_norm.bigWig; do
            BASE=`basename $FILE | sed s/\.plus_norm.bigWig//`
            echo "
track $PREFIX"_"$SAMPLE
container multiWig
aggregate transparentOverlay
showSubtrackColorOnUi on
shortLabel $PREFIX"_"$SAMPLE
longLabel $PREFIX"_"$SAMPLE RNA-sequencing
boxedCfg on
autoScale on
maxHeightPixels 128:64:16
visibility full
type bigWig

track $PREFIX"_"$SAMPLE"_plus"
parent $PREFIX"_"$SAMPLE
type bigWig
bigDataUrl $LINK/$BASE".plus_norm.bigWig"
color 255,0,0

track $PREFIX"_"$SAMPLE"_minus"
parent $PREFIX"_"$SAMPLE
type bigWig
bigDataUrl $LINK/$BASE".minus_norm.bigWig"
color 153,0,0
" >> $HOME_PATH/trackdb.txt
        done

        rm $BAM_PATH/*_for* $BAM_PATH/*_rev*
    done
done

wait
```
