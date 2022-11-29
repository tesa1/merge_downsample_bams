# Zwart lab standard workflow to merge ChIP-seq bam files for use in comparison profile/tornado heatmap plots -- instructions for the RHPC server

Provided are the Zwart lab approved basic steps to merge bam files from a ChIP-seq experiment to make comparison profile plots when running on the RHPC server.

## Vignette Info

Here we show the basic steps for merging filtered bam files in two different categories so that they can be compared in an ESeq/deepTools profile or heatmap for specific regions. This is a real-world example using real-world data **(thanks, Yanyun Zhu!)**. The data come from FOXA1 ChIP-seq experiments of normal (healthy) and primary prostate cancer samples. Please do not disseminate or interpret these data, they are for instructional purposes only. 

Requirements for this tutorial are familiarity with:

- ChIP-seq 
- pre-/post-GCF course from Sebastian Gregoricchio and Tesa Severson
- Zwartlab snakePipes ChIP-seq pipeline written by former lab member Joe Siefert
- Upstream steps to merging including alignment of fastq.gz files to genome, filtering aligned files with MQ20

  
This vignette assumes you have aligned and filtered (MQ20) bam files for your samples of interest. This can be done by running peakcalling on your samples using the Zwartlab snakePipes ChIP-seq pipeline.

We also assume you have looked at the snakePipes ChIP-seq pipeline /QC_report/QC_report_all.tsv file for your experiments and have determined the read counts are enough (at least 20M reads/sample) and that QC metrics like Fraction of Reads in Peaks (FRiP) are acceptable and roughly similar between experiments (within ChIP-seq factors, eg. H3K27ac) and not significantly different between comparisons. If not or you do not know, ask for bioinformatic help, you may need to do additional sequencing.

In this tutorial for 2 categories of samples (Category A and B) we will merge the files of interest, index them, check their mapped reads and downsample accordingly.

**Zwart lab reproducibility means you observe these approved methods and do not deviate from them unless you have very specific scientific reasons**. In this case, ask for bioinformatic help.


 ## Merging and indexing of filtered bam files from Category A ##
Use samtools merge to create a new file `foxa1_healthy.bam` from the other filtered.bam files listed. 
Additionally, we use a flag to run this on 8 cores to speed up the process. Note, this will create a very 
big file as you are merging 10 bam files together. 

 ```bash
samtools merge -@8 foxa1_healthy.bam wz2086.filtered.bam wz2088.filtered.bam wz2090.filtered.bam 

samtools index foxa1_healthy.bam
```

 ## Check the mapped reads of the newly merged Category A file ##
Now use samtools flagstat to get the number of mapped reads in your new file. 
The new file `foxa1_healthy.bam` has 55524941 mapped reads. 

 ```bash
samtools flagstat foxa1_healthy.bam > foxa1_healthy.flag

cat foxa1_healthy.flag
```

![Screenshot](cat_foxa1_healthy_flagstat.png)

## Downsample Category A file to roughly 20M reads ##
Use samtools view to downample a file (-b) by a fraction (-s) to obtain roughly 20M reads.

```bash
samtools view -s 0.36 -b foxa1_healthy.bam > foxa1_healthy_ds.bam

samtools index foxa1_healthy_ds.bam
```

## Double-check you now have rougly 20M reads in Category A downsampled file

```bash
samtools flagstat foxa1_healthy_ds.bam > foxa1_healthy_ds.flag

cat foxa1_healthy_ds.flag
```

![Screenshot](cat_foxa1_healthy_ds_flagstat.png)

## Do the same for Category B files ## 

Based on the top line in the foxa1_primary.flag file you can calculate the fraction (N) to input into the dowsample command to obtain close to 20M reads.

For example 20000000/N = fraction_n

```bash
samtools merge -@8 foxa1_primary.bam file1.bam file2.bam file3.bam file4.bam 

samtools index foxa1_primary.bam

samtools flagstat foxa1_primary.bam > foxa1_primary.flag

samtools view -s fraction_n -b foxa1_primary.bam > foxa1_primary_ds.bam

samtools index foxa1_primary_ds.bam

samtools flagstat foxa1_primary_ds.bam > foxa1_primary_ds.flag
```

