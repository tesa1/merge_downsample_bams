# Zwart lab standard workflow to merge ChIP-seq bam files for use in comparison profile/tornado heatmap plots -- instructions for the RHPC server

Provided is a Zwart lab approved method to merge bam files from a ChIP-seq experiment to make comparison profile plots when running on the RHPC server.

## Vignette Info

With comparison plots, you can show the differences in signal strength for different groups of ChIP-seq data across a set (or sets) or regions of interest:
![Screenshot](profile_example.png)


Some things to consider when comparing signal across datasets: you need to make sure that you don't see a difference in signal because one set of data has more samples than the other dataset(s). To avoid that, for each data set individually, we merge the filtered aligned bam files within the set, check the number of reads in the merged file and then downsample to a target number of reads (~20M). This way, when you have done this with all the sets you want to compare, you have datasets of roughly equivalent read counts. The resulting files can be used in an EaSeq/deepTools profile or tornado/heatmap for your regions of interest. This is a realistic example of FOXA1 ChIP-seq data from normal (healthy) and primary prostate cancer samples. Please do not disseminate or interpret these data, they are for instructional purposes only. 

Requirements for this tutorial are familiarity with:

- ChIP-seq
- pre-/post-GCF course from Sebastian Gregoricchio and Tesa Severson
- Zwartlab snakePipes ChIP-seq pipeline written by former lab member Joe Siefert
- Upstream steps to merging including alignment of fastq.gz files to genome, filtering aligned files with MQ20
- samtools 

  
This vignette assumes you have aligned and filtered (MQ20) bam files for your samples of interest. This can be done by running peakcalling on your samples using the Zwartlab snakePipes ChIP-seq pipeline.

We also assume you have looked at the snakePipes ChIP-seq pipeline /QC_report/QC_report_all.tsv file for your experiments and have determined the read counts are enough (at least 20M reads/sample) and that QC metrics like Fraction of Reads in Peaks (FRiP) are acceptable and roughly similar between experiments (within ChIP-seq factors, eg. H3K27ac) and not significantly different between comparisons. If not or you do not know, ask for bioinformatic help, you may need to do additional sequencing.

In this tutorial we want to compare 2 datasets of samples (datasetA and datasetB). First we will merge the files of interest, index them, check their mapped reads and downsample accordingly.

**This is just one way to compare datasets in a normalized fashion. There are other (also correct) ways to do this, but this manner requires very few steps and is relatively convenient.**

## Merging and indexing of filtered bam (MQ20) files from datasetA (3 healthy tissue FOXA1 ChIP-seq samples) ##
Use samtools merge to create a new file `foxa1_healthy.bam` from the other filtered.bam files listed. 
Additionally, we use a flag to run this on 8 cores to speed up the process. Note, this will create a very 
big file as you are merging many bam files together. 

```bash
# merge files
samtools merge -@8 foxa1_healthy.bam file1.filtered.bam file2.filtered.bam file3.filtered.bam 

# index
samtools index foxa1_healthy.bam
```

## Check the mapped reads of the newly merged datasetA file ##
Now use samtools flagstat to get the number of mapped reads in your new file. 
The new file `foxa1_healthy.bam` has 55524941 high-quality mapped reads. You can see this from the top line of the .flag file  

```bash
# make a flagstat file to see the number of mapped reads
samtools flagstat foxa1_healthy.bam > foxa1_healthy.flag

# view the flagstat file on the terminal
cat foxa1_healthy.flag
```

![Screenshot](cat_foxa1_healthy_flagstat.png)

## Downsample datasetA file to roughly 20M reads ##
Note, you can never get exactly 20M reads with this downsample command because it downsamples reads randomly and the fraction you input is not precise. As long as it???s very close to 20M it???s fine. 

Use samtools view to downample a file (-b) by a fraction (-s) to obtain roughly 20M reads.

```bash
# downsample the merged file
samtools view -s 0.36 -b foxa1_healthy.bam > foxa1_healthy_ds.bam

# index the downsampled file
samtools index foxa1_healthy_ds.bam
```

## Double-check you now have rougly 20M reads in datasetA downsampled file

```bash
# make a flagstat file to see the number of mapped reads
samtools flagstat foxa1_healthy_ds.bam > foxa1_healthy_ds.flag

# view the flagstat file on the terminal
cat foxa1_healthy_ds.flag
```

![Screenshot](cat_foxa1_healthy_ds_flagstat.png)

## Clean up -- remove non-downsampled files because they're very big ##

```bash
# remove the big bam file and it's index
rm foxa1_healthy.bam
rm foxa1_healthy.bam.bai
```

## Do the same for datasetB files (4 primary tumor tissue FOXA1 ChIP-seq samples) ## 

Based on the top line in the foxa1_primary.flag (N) file you can calculate the fraction (fraction_n) to input into the dowsample command to obtain close to 20M reads.

For example 20000000/N = fraction_n

```bash
# merge files
samtools merge -@8 foxa1_primary.bam file1.filtered.bam file2.filtered.bam file3.filtered.bam file4.filtered.bam 

# index merged file
samtools index foxa1_primary.bam

# get number of reads mapped
samtools flagstat foxa1_primary.bam > foxa1_primary.flag

# check number of reads mapped by viewing on the terminal
cat foxa1_primary.flag

# downsample according to number of reads mapped and target number of reads
samtools view -s fraction_n -b foxa1_primary.bam > foxa1_primary_ds.bam

# index 
samtools index foxa1_primary_ds.bam

# double-check your dowsample
samtools flagstat foxa1_primary_ds.bam > foxa1_primary_ds.flag

# cleanup
rm foxa1_primary.bam
rm foxa1_primary.bam.bai
```

## You can now use these downsampled bam files to look at comparison plots ##
In EaSeq you can use the ds.bam files directly with a bed file of interest to make your profile or tornado plots.  For deeptools, you need to first convert to bigwig.

## Deeptools commands to make a comparison profile plot with bed file of interest
Assumes you have 
```bash
# activate deeptools
conda activate deeptools

# make bigwigs from your bam file (default data bin settings)
bamCoverage -b foxa1_healthy_ds.bam -o foxa1_healthy_ds.bw
bamCoverage -b foxa1_primary_ds.bam -o foxa1_primary_ds.bw

# make a matrix of your files and the coverage at your regions
computeMatrix reference-point -S foxa1_healthy_ds.bw foxa1_primary_ds.bw -R foxa1_peaks.bed  -b 2000 -a 2000 --numberOfProcessors 8 --outFileName foxa1_healthy_primary_vs_foxa1_peaks_matrix.gz --referencePoint center

# plot the data in a profile with standard-error shading and labeling
plotProfile -m foxa1_healthy_primary_vs_foxa1_peaks_matrix.gz --plotType se --perGroup -out foxa1_healthy_primary_vs_foxa1_peaks_profile.pdf --regionsLabel FOXA1_sites

# plot a basic heatmap/tornado of the data
plotHeatmap -m foxa1_healthy_primary_vs_foxa1_peaks_matrix.gz  -out foxa1_healthy_primary_vs_foxa1_peaks_heatmap.pdf --regionsLabel FOXA1_sites --samplesLabel healthy_FOXA1 primary_FOXA1

# close your environment
conda deactivate
```


