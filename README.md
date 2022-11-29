# Zwart lab standard workflow to merge ChIP-seq bam files for use in comparison profile plots -- instructions for the RHPC server

Provided are the Zwart lab approved basic steps to merge bam files from a ChIP-seq experiment to make comparison profile plots when running on the RHPC server.

## Vignette Info

Here we show the basic steps for merging filtered bam files in two different categories so that they can be compared in an ESeq/deepTools profile or heatmap for specific regions. This is a real-world example using real-world data **(thanks, Yanyun Zhu!)**. The data come from a ChIP-seq experiments of normal (healthy) and primary prostate cancer samples. Please do not disseminate or interpret these data, they are for instructional purposes only. 

Requirements for this tutorial are familiarity with:

- ChIP-seq 
- pre-/post-GCF course from Sebastian Gregoricchio and Tesa Severson
- Zwartlab snakePipes ChIP-seq pipeline written by former lab member Joe Siefert
- Upstream steps to merging including alignment of fastq.gz files to genome, filtering aligned files with MQ20

  
This vignette assumes you have aligned and filtered (MQ20) bam files for your samples of interest. This can be done by running peakcalling on your samples using the Zwartlab snakePipes ChIP-seq pipeline.

We also assume you have looked at the snakePipes ChIP-seq pipeline /QC_report/QC_report_all.tsv file for your experiments and have determined the read counts are enough (at least 20M reads/sample) and that QC metrics like Fraction of Reads in Peaks (FRiP) are acceptable and roughly similar between experiments (within ChIP-seq factors, eg. H3K27ac) and not significantly different between comparisons. If not or you do not know, ask for bioinformatic help, you may need to do additional sequencing.


 ## Installation and use of HiC-Pro ##
Hi-C Pro can be installed here from https://github.com/nservant/HiC-Pro

On darwin, using downloaded version of HiC-Pro.2.11 from https://github.com/nservant/HiC-Pro/releases and installed in /home/t.severson/tools/HiC-Pro-2.11.1.install/bin/. Config file was changed to suit the necessary requirements (bowtie2 and samtools 1.9 locations).

Move fastq sample folders into folder h3k27ac_rawdata. 

 ```bash
/home/t.severson/tools/HiC-Pro-2.11.1.install/bin/HiC-Pro -i h3k27ac_rawdata/ -o h3k27ac_outputs -c config-hicpro.txt
```

This will create h3k27ac_outputs folder with mapping and HiC analysis results


## Installation and use of hichipper for HiC-Pro results ##

Steps for analyzing of Hi-ChIP data from HiC Pro using hichipper using a sample manifest file (.yaml). Hichipper uses macs2 to call peaks. I generated peaks using the Mubach method (i.e. EACH,ALL).
Using called peaks and restriction fragment locations as input and the output from HiC-Pro output that can be used to:
  - Determine library quality
  - Visualize loops with strength and confidence metrics

Install hichipper using instructions from: https://hichipper.readthedocs.io/en/latest/

```bash
virtualenv -p /usr/bin/python2.7 venv
source venv/bin/activate
pip install hichipper
hichipper --version
```
sss
Run hichipper on the h3k27ac_outputs folder from HiC-Pro.
```bash
hichipper --out h3k27ac_hichipper --make-ucsc  hichipper.yaml
```


 ## Running a downsample experiment on the data ##
 
Need to use HTSeq to downsample the paired-end fastq files 

```bash
# install htseq with bioconda
conda create -c bioconda --name htseq htseq
conda activate htseq
```

For downsampling to 25% of reads
```bash
gunzip sample_R1.fastq.gz
gunzip sample_R1.fastq.gz
python subsample.py 0.25 sample_R1.fastq sample_R2.fastq sample_R1_ds25.fastq sample_R2_ds25.fastq
gzip -c sample_R1_ds25.fastq > sample_R1_ds25.fastq.gz
gzip -c sample_R2_ds25.fastq > sample_R2_ds25.fastq.gz
```

The subsample.py script is from the developer of HTSeq 
http://seqanswers.com/forums/archive/index.php/t-12070.html: 
It will throw a 'StopIteration' error but the files will be fine.

Next, run HiC-Pro on the downsampled samples

```bash
/home/t.severson/tools/HiC-Pro-2.11.1.install/bin/HiC-Pro -i ds_25/ -o ds_25_outputs -c config-hicpro.txt
````
Footer
