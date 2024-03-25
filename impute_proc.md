Imputation process for studies with existing PLINK files
================================================================================

# Prerequisites

## Software

* [bcftools](https://github.com/samtools/bcftools/releases/download/1.19/bcftools-1.19.tar.bz2) ([manual](https://samtools.github.io/bcftools/bcftools.html))
* [samtools](https://github.com/samtools/samtools/releases/download/1.19.2/samtools-1.19.2.tar.bz2) ([manual](https://www.htslib.org/doc/samtools.html))
* [PLINK](https://www.cog-genomics.org/plink/1.9/) (at least 1.90)
* The prisma R package or another way to modify PLINK files in terms of 
major/minor allele order to comply with the VCF format

## Data

* Genotypic data in PLINK or **proper** [VCF](https://samtools.github.io/hts-specs/VCFv4.2.pdf) format.
* A proper reference genome FASTA file and a `.fai` index to post-process VCF
files created with PLINK.

Regarding the 2nd point above, there are two options:

- The hg19 (GRCh37) human genome version
- The hg38 (GRCh38) human genome version

To create a `.fai` file for hg19:

```
wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
gunzip hs37d5.fa.gz
samtools faidx hs37d5.fa
rm hs37d5.fa
```

To create a `.fai` file for hg38:

```
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
gunzip GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
mv GCA_000001405.15_GRCh38_no_alt_analysis_set.fna hg38_no_alt.fa
samtools faidx hg38_no_alt.fa
rm hg38_no_alt.fa
```

# Process outline

## Pre-imputation

1. Read PLINK into a `GWASExperiment` object in prisma.
2. Run `prisma::harmonizeWithReference()` on the imported `GWASExperiment`.
3. Write PLINK files with the output of `prisma::harmonizeWithReference()` with
`prisma::writePlink()`.
4. Convert the output PLINK files to VCF with `plink`.
5. The output VCF does not contain a proper VCF header. Correct this with 
`bcftools reheader` and the proper human reference genome FASTA.
6. The output of (5) may not be properly sorted. Sort it with `bcftools sort`
and output compressed and indexed.
7. Split the sorted output per chromosome.
8. Upload the resulting chromosome-wise VCFs to TOPMED server for imputation.

## Post-imputation linear processing

9. Concatenate the resulting files (one large VCF per chromosome) with 
`bcftools concat`.
10. The resulting file does not have an ID for each SNP. Fix this with 
`bcftools annotate`.
11. Use `bcftools filter` to remove the `TYPED` variants from the imputation 
results as this will later cause problems with duplicates.
12. Use `bcftools filter` to filter variants with low imputation score less than
a threshold (e.g. 0.2).
13. Some deconvoluted duplicate complex indels remain (e.g. rs893063787,
rs538643934) remain because of their existence in the reference panel. Fix this
with `bcftools norm`.

## Post-imputation parallel processing

9. The resulting files do not have an ID for each SNP. Fix this with 
`bcftools annotate` in parallel for each chromosome.
10. Use `bcftools filter` to remove the `TYPED` variants from the imputation 
results as this will later cause problems with duplicates in parallel for each
chromosome.
11. Use `bcftools filter` to filter variants with low imputation score less than
a threshold (e.g. 0.2) in parallel for each chromosome.
12. Some deconvoluted duplicate complex indels remain (e.g. rs893063787,
rs538643934) remain because of their existence in the reference panel. Fix this
with `bcftools norm`.
13. Concatenate the resulting files (one large VCF per chromosome) with 
`bcftools concat`.

## Final conversion to PLINK format

14. After concatenation, convert to PLINK bid/bed/fam with `plink --vcf`.


# For the brave

The following bash script covers in almost one step, Steps 9-13 of the post
imputation process. It is also faster, as it firstly filters `TYPED` and 
'R2<0.2' variants and then assigning IDs to the remaining. The ID assignment
process of `bcftools` is slower than filtering, so if we assign IDs 
post-filtering, fewer IDs are assigned resulting in a faster completion. The 
final output is most probable suitable for PLINK operations as it is most
probably fully annotated and deduplicated.

It is advised that you study the whole process in detail before running the
following script. However, because of piping, this will run blazingly faster 
than individual steps but you must make sure you fully understand what is going
on.

```
#!/bin/bash

echo "Processing chromosomes..."
for CHR in `seq 1 22`
do
  bcftools filter -e 'INFO/TYPED=1 | INFO/R2<0.2' chr${CHR}.dose.vcf.gz | \
    bcftools norm --rm-dup indels - | \
    bcftools annotate --set-id +'%CHROM\:%POS\:%REF\:%FIRST_ALT' \
      -Oz -o chr${CHR}_proc.vcf.gz - &
done

wait

echo "Concatenating chromosomes..."
bcftools concat -Oz -o imputed_cleaned.vcf.gz chr{1..22}_proc.vcf.gz

echo "Cleaning..."
rm chr*proc.vcf.gz

echo "Done!"
```

The above script is a fast way to prepare PLINK files from TOPMED imputation
server results but it does not offer fine control over intermediate steps like
filtering threshold etc.

The final command to convert to PLINK files:

```
plink --vcf imputed_cleaned.vcf.gz --make-bed --out imputed_cleaned
```


# Detailed steps

We assume that we start from a bim/bed/fam triplet of PLINK files named
`data.bim`, `data.bed`, `data.fam` and the reference genome is hg19 (or GRCh37).
The other option is hg38 (or GRCh38).

## Steps 1-3

**1. Read PLINK into a `GWASExperiment` object in prisma.**
**2. Run `prisma::harmonizeWithReference()` on the imported `GWASExperiment`.**
**3. Write PLINK files with the output of `prisma::harmonizeWithReference()` 
with `prisma::writePlink()`.**

Within R:

```
library(prisma)

indata <- importGWAS(
   input="data",
   backend="snpStats",
   genome="hg19"
)

outdata <- harmonizeWithReference(indata)

writePlink(outdata,outBase="data_prevcf")
```

## Step 4

**4.Convert the output PLINK files to VCF with `plink`.**

```
plink --bfile data_prevcf --keep-allele-order --recode vcf --out data_unheaded
```

## Step 5-6

**5. The output VCF does not contain a proper VCF header. Correct this with 
`bcftools reheader` and the proper human reference genome FASTA.**
**6. The output of (5) may not be properly sorted. Sort it with `bcftools sort`
and output compressed and indexed.**

```
bcftools reheader --fai hs37d5.fa.fai -o data_unsorted.vcf data_unheaded.vcf
bcftools sort -Oz -o data.vcf.gz data_unsorted.vcf
tabix data.vcf.gz
```

## Step 7

**7. Split the sorted output per chromosome.**

```
bcftools index -s data.vcf.gz | cut -f1 | while read C; do bcftools view -Oz -o data_chr${C}.vcf.gz data.vcf.gz "${C}"; done
```

## Step 8

Create an account in [TOPMED imputation server](https://imputation.biodatacatalyst.nhlbi.nih.gov/).
Follow the instructions to upload the `data_chr....vcf.gz` files and start and
imputation job. Select "Skip" in the "Population" option and be **sure** about
the "Array Build" option (hg19/hg38). Wait for the results. Upon finished,
follow the instruction to download the results. These are 22 VCF files and 22
zipped text files with information for each imputed variant, contained in
22 zip encrypted zip files. The password to decrypt is emailed.

## Step 9

**9. Concatenate the resulting files (one large VCF per chromosome) with 
`bcftools concat`.**

```
bcftools concat -Oz -o imputed.vcf.gz chr{1..22}.dose.vcf.gz
```

**TIP**: After imputation, the resulting VCF files are huge. Therefore, commands
like the above will take time to finish. Consider detaching it from the terminal
so you can log out safely without interrupting your work. You can do this with
`nohup`. For example:

```
nohup bcftools concat -Oz -o imputed.vcf.gz chr{1..22}.dose.vcf.gz &
```

## Step 10

**10. The resulting file does not have an ID for each SNP. Fix this with 
`bcftools annotate`.**

The imputation resulting VCF files do not have IDs (usually [dbSNP](https://www.ncbi.nlm.nih.gov/snp/) ids) for each variant. That is 
because the TOPMED reference panel contains also novel variants, not yet 
verified and submitted to dbSNP. We need an ID for each variant. We correct this
with `bcftools annotate`. The ID will be of the format:

*chromosome:position:reference_genome_allele:alternative_allele*

```
bcftools annotate --set-id +'%CHROM\:%POS\:%REF\:%FIRST_ALT' -Oz -o imputed_ids.vcf.gz imputed.vcf.gz
```

## Step 11

**11. Use `bcftools filter` to remove the `TYPED` variants from the imputation 
results as this will later cause problems with duplicates.**

```
bcftools filter -e 'INFO/TYPED=1' -Oz -o imputed_ids_only.vcf.gz imputed_ids.vcf.gz
```

## Step 12

**12. Use `bcftools filter` to filter variants with low imputation score less 
than a threshold (e.g. 0.2).**

```
bcftools filter -e 'INFO/R2<0.2' -Oz -o imputed_ids_only_filtered.vcf.gz imputed_ids_only.vcf.gz
```

**HINT** We do not apply both filters in steps 11 and 12 together because we may 
want to experiment with different R<sup>2</sup> thresholds. However, we can 
apply them at the same time by changing the filter expression to 
`INFO/TYPED=1 | INFO/R2<0.2`.

## Step 13

**13. Some deconvoluted duplicate complex indels remain (e.g. rs893063787,
rs538643934) remain because of their existence in the reference panel. Fix this
with `bcftools norm`.**

This will be accomplished by specifycing the `--rm-dup` option of `bcftools` to
remove the second instance of `indels` regardless their ref/alt alleles. See 
also `[bcftools norm](https://samtools.github.io/bcftools/bcftools.html#norm)`
and `--collapse` in [Common Options](https://samtools.github.io/bcftools/bcftools.html#norm).

```
bcftools norm --rm-dup indels -Oz -o imputed_ids_only_filtered_rmdup.vcf.gz imputed_ids_only_filtered.vcf.gz
```

## Step 14

We are ready to create the final PLINK files. Most probably, they are ready for
downstream operations such as LD-pruning, PCA etc. If not, only minor processing
should be required.

```
plink --vcf imputed_ids_only_filtered_rmdup.vcf.gz --make-bed --out imputed_ids_only_filtered_rmdup
```

### Clean up

```
rm data_prevcf* data_unheaded* data_unsorted* imputed.vcf.gz imputed_ids.vcf.gz imputed_ids_only.vcf.gz imputed_ids_only_filtered.vcf.gz
```


# Alternative post-imputation processing

As the VCF files after imputation are huge and even `bcftools` appear slow. We 
will attempt to speed up the process in a multi-processor system with enough 
memory by re-arranging step 9 to the end of the process and do all id-setting 
and filtering per chromosome in parallel. This process may require also more 
storage for intermediate steps.

The following start after downloading the imputation results.

## Step 9

**9. The resulting files do not have an ID for each SNP. Fix this with 
`bcftools annotate` in parallel for each chromosome.**

The imputation resulting VCF files do not have IDs (usually [dbSNP](https://www.ncbi.nlm.nih.gov/snp/) ids) for each variant. That is 
because the TOPMED reference panel contains also novel variants, not yet 
verified and submitted to dbSNP. We need an ID for each variant. We correct this
with `bcftools annotate`. The ID will be of the format:

*chromosome:position:reference_genome_allele:alternative_allele*

```
for CHR in `seq 1 22`
do
  bcftools annotate --set-id +'%CHROM\:%POS\:%REF\:%FIRST_ALT' -Oz -o chr${CHR}_ids.dose.vcf.gz chr${CHR}.dose.vcf.gz &
done
```

Wait for this to finish before moving to the next step.

**HINT**: The above cannot run with `nohup`. Consider placing in a small shell 
script like `id.sh` and theh `nohup sh id.sh &`.

## Step 10

**10. Use `bcftools filter` to remove the `TYPED` variants from the imputation 
results as this will later cause problems with duplicates in parallel for each
chromosome.**

```
for CHR in `seq 1 22`
do
  bcftools filter -e 'INFO/TYPED=1' -Oz -o chr${CHR}_ids_only.dose.vcf.gz  chr${CHR}_ids.dose.vcf.gz &
done
```

Wait for this to finish before moving to the next step.

## Step 11

**11. Use `bcftools filter` to filter variants with low imputation score less 
than a threshold (e.g. 0.2) in parallel for each chromosome.**

```
for CHR in `seq 1 22`
do
  bcftools filter -e 'INFO/R2<0.2' -Oz -o chr${CHR}_ids_only_filtered.dose.vcf.gz  chr${CHR}_ids_only.dose.vcf.gz &
done
```

**HINT** We do not apply both filters in steps 11 and 12 together because we may 
want to experiment with different R<sup>2</sup> thresholds. However, we can 
apply them at the same time by changing the filter expression to 
`INFO/TYPED=1 | INFO/R2<0.2`.

## Step 12

**12. Some deconvoluted duplicate complex indels remain (e.g. rs893063787,
rs538643934) remain because of their existence in the reference panel. Fix this
with `bcftools norm`.**

This will be accomplished by specifycing the `--rm-dup` option of `bcftools` to
remove the second instance of `indels` regardless their ref/alt alleles. See 
also `[bcftools norm](https://samtools.github.io/bcftools/bcftools.html#norm)`
and `--collapse` in [Common Options](https://samtools.github.io/bcftools/bcftools.html#norm).

```
for CHR in `seq 1 22`
do
  bcftools norm --rm-dup indels -Oz -o chr${CHR}_ids_only_filtered_rmdup.dose.vcf.gz  chr${CHR}_ids_only_filtered.dose.vcf.gz &
done
```

Wait for this to finish.

## Step 13

**12. Concatenate the resulting files (one large VCF per chromosome) with 
`bcftools concat`.**

```
bcftools concat -Oz -o imputed_cleaned.vcf.gz chr{1..22}_ids_only_filtered_rmdup.dose.vcf.gz
```

## Step 14

We are ready to create the final PLINK files. Most probably, they are ready for
downstream operations such as LD-pruning, PCA etc. If not, only minor processing
should be required.

```
plink --vcf imputed_cleaned.vcf.gz --make-bed --out imputed_cleaned
```

### Clean up

```
rm data_prevcf* data_unheaded* data_unsorted* chr*_ids.dose.vcf.gz chr*_ids_only.dose.vcf.gz chr*_ids_only_filtered.dose.vcf.gz
```

# Other steps

After post-imputation processing of multiple datasets, it should also be easier
to merge them with `bcftools merge`. Then, the merged files can be converted
to PLINK format for population-statistics-based processing of the merged 
dataset.


### TODO

* Consider firstly removing the `TYPED` and low imputation score variants and
then adding ids as the id-addition step is much slower.

### Tests

1. Duplicate detection

```
plink --bfile imputed_cleaned --list-duplicate-vars

# Output
# CHR     POS     ALLELES IDS
# No duplicates! - OK
```

2. LD-pruning

```
plink --bfile imputed_cleaned --indep-pairwise 50 5 0.2 --out test 
```
