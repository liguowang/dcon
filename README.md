# Dcon: Estimate DNA contamination

Next-generation sequencing (NGS) of DNA is routinely used in clinical testing. However,
the complicated, multistep DNA extraction protocol and the unknown pre-analytic specimen
handling often lead to sample-to-sample DNA contaminations. DNA contamination poses
substantial challenges in NGS assays and affects diagnosis and treatment if contaminant
DNA harbors clinically significant variants.

**Dcon** is a program to detect and quantify sample-to-sample germline DNA contaminations using the
skewed (or shifted) distribution of variant allele frequency (VAF, VAF refers to the
fraction of sequencing reads overlapping a genomic coordinate of a DNA variant such as SNP).

**Unlike other tools, Dcon does NOT require any prior knowledge such as "genotype" and/or "allele frequencies in human population" calculated from 1000 Genome project. Dcon works on single BAM file and does not require tumor-normal pairs. Therefore, Dcon is also able to detect contamination from DNA sequencing data generated from non-human samples, where allele frequencies and genotype are generally not available.**

## Getting Started

### How Dcon works

![workflow](https://github.com/liguowang/dcon/blob/master/img/worflow.png?raw=true)


### Prerequisites
Dcon requires [python2.7](https://www.python.org/download/releases/2.7/) and these python packages (These packages will be automatically installed, please see instructions below). 

1. [**numpy**](http://www.numpy.org/)
2. [**scipy**](https://www.scipy.org/)
3. [**pysam**](https://pypi.python.org/pypi/pysam)
4. [**scikit-learn**](https://pypi.python.org/pypi/scikit-learn)

### Installing

Following instructions [here](https://pip.pypa.io/en/stable/installing/) to install **pip**. Please note: **pip** is already installed if you're using Python 2 >=2.7.9 or Python 3 >=3.4 binaries downloaded from [python.org](https://www.python.org/), but you'll need to upgrade **pip**. 

Users have two ways to install Dcon:

1. Install Dcon from [pypi](https://pypi.python.org/pypi):
    * Open a terminal and type command: ```pip install Dcon```

2. Install Dcon from local directory
    * Download Dcon directly from: [this link](https://github.com/liguowang/dcon/archive/0.1.5.zip)
    * Go to the directory where **dcon-version.zip** was saved (for example, dcon-0.1.5.zip).
    * Open a terminal and type command:`pip install dcon-0.1.5.zip`

### General usage

User provides a [BAM](https://genome.ucsc.edu/FAQ/FAQformat.html#format5.1) file and a [VCF](https://genome.ucsc.edu/FAQ/FAQformat.html#format10.1) file. Variants in VCF file can be called from this particular BAM file or extracted from [dbSNP](https://www.ncbi.nlm.nih.gov/SNP/) (i.e. extract all SNPs that falling into captured genomic regions, Dcon automatically checks the VAF for each SNP and skips irrelevant SNPs).

   `python2.7 Dcon.py -b input.bam -v input.vcf -o output`
    
**Tips: How do I get a sub-section of a VCF file?**
1. Download and install [tabix](http://sourceforge.net/projects/samtools/files/tabix/)
2. Using tabix to get a sub-section of a VCF file.

`tabix -h ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20100804/ALL.2of4intersection.20100804.genotypes.vcf.gz 2:39967768-39967768 > output.vcf`

General usage instruction is listed below:   
```
Usage: Dcon.py [options]

Dcon: Estimate DNA Contamination from BAM file. 

Options:
  --version             show program's version number and exit
  -h, --help            show this help message and exit
  -b BAM_FILE, --bam=BAM_FILE
                        Alignment file in BAM format. BAM file must be indexed
                        and sorted using samTools
                        (http://samtools.sourceforge.net/).
  -v VCF_FILE, --variant=VCF_FILE
                        VCF format file defining SNPs. VCF file can be a plain
                        text, compressed (.gz, .z, .Z, .bz, .bz2 and .gzip2)
                        or remote file (http://, https://, ftp://). Exclusive
                        with '-r'.
  -o OUTPUT_FILE, --output=OUTPUT_FILE
                        Prefix of output file. Will generate three files:
                        "prefix.SNP.tsv", "prefix.PI.tsv" "prefix.dcon.R".
                        "prefix.dcon.pdf" will be generated if R exists.
  -q MIN_MAPQ, --mapq=MIN_MAPQ
                        Minimum mapping quality
                        (http://maq.sourceforge.net/qual.shtml). Mapping
                        quality is Phred-scaled probability of an alignment
                        being wrong. default=30
  -m MIN_SEQQ, --seqq=MIN_SEQQ
                        Minimum base phred quality
                        (http://maq.sourceforge.net/qual.shtml). Base quality
                        is Phread-scaled probability of a base calling being
                        wrong. default=30
  -c MIN_COVERAGE, --cvg=MIN_COVERAGE
                        Minimum number of reads supporting variant. default=30
  -a MIN_ALLELE_READ_COUNT, --allele=MIN_ALLELE_READ_COUNT
                        Minimum number of reads supporting one allele.
                        default=5
  -p PROCESSOR_NUM, --processor=PROCESSOR_NUM
                        Number of processes used to call SNPs. default=1
  -s, --skipXY          Skip SNPs on X and Y chromosomes. default=False
  -y PROBABILITY_CUT, --prob-cutoff=PROBABILITY_CUT
                        Probability cutoff. default=0.5


```

## Create testing datasets
Targeted-capture sequencing data of 10 chronic myelomonocytic leukemia patients were downloaded from [NCBI Sequence Read Archive (SRA)](https://www.ncbi.nlm.nih.gov/sra).

* [SRR6756023](https://www.ncbi.nlm.nih.gov/sra/?term=SRR6756023)
* [SRR6756025](https://www.ncbi.nlm.nih.gov/sra/?term=SRR6756025)
* [SRR6756028](https://www.ncbi.nlm.nih.gov/sra/?term=SRR6756028)
* [SRR6756033](https://www.ncbi.nlm.nih.gov/sra/?term=SRR6756033)
* [SRR6756036](https://www.ncbi.nlm.nih.gov/sra/?term=SRR6756036)
* [SRR6756040](https://www.ncbi.nlm.nih.gov/sra/?term=SRR6756040)
* [SRR6756042](https://www.ncbi.nlm.nih.gov/sra/?term=SRR6756042)
* [SRR6756045](https://www.ncbi.nlm.nih.gov/sra/?term=SRR6756045)
* [SRR6756048](https://www.ncbi.nlm.nih.gov/sra/?term=SRR6756048)  
* [SRR6756053](https://www.ncbi.nlm.nih.gov/sra/?term=SRR6756053)

A total of 45 synthetic datasets were created by mixing the above 10 original datasets in pairwise at 5%, 10%, ..., 45%, 50%. Mixed reads were then aligned to human genome (hg19) using BWA. Below are mixed BAM files from "SRR6756025" and "SRR6756028":

| File Name                   |  Mix percentage  |   Reads from SRR6756025 |  Reads from SRR6756028  |  Total Reads  | Download Links                                                                                                                                      |  MD5Sum of BAM file            |
|:----------------------------|:----------------:|:-----------------------:|:-----------------------:|:-------------:|:---------------------------------------------------------------------------------------------------------------------------------------------------:|:------------------------------:|
|SRR6756025_SRR6756028_P05.bam|      5%          |1,900,360                |100,122                  |2,000,482      | [bam](https://drive.google.com/open?id=1ueMZYaw927Wa1V-QJ2ir-L-59nhVezdr), [bai](https://drive.google.com/open?id=1uQZ_ZMG9JzRb7bcKFRV0JtAYeUGAMzHK)|0fa86559117ee9f3070986593ff3b667|
|SRR6756025_SRR6756028_P10.bam|      10%         |1,800,190                |199,734                  |1,999,924      | [bam](https://drive.google.com/open?id=1r1OHZ04SQnGygyslEkRurCk_2jTv7CVR), [bai](https://drive.google.com/open?id=1fV7TGq9zb0xQ981F0ZI_45GYcgHTmyrp)|3e8705858900060836eccf984e4fe778|
|SRR6756025_SRR6756028_P15.bam|      15%         |1,700,516                |298,978                  |1,999,494      | [bam](https://drive.google.com/open?id=1WMkmwRlA6OEYq_fOn8jsBa8sGoKVwDKv), [bai](https://drive.google.com/open?id=1P91gOw31QcBl-wGqNiLFezuxmDjlg3Sw)|9ea836d52e7b78ec2e21d7fe8a547bdf|
|SRR6756025_SRR6756028_P20.bam|      20%         |1,600,992                |399,158                  |2,000,150      | [bam](https://drive.google.com/open?id=16pUsoZ97kw_nw8p4xA65Pf9NJ_h-b-Mq), [bai](https://drive.google.com/open?id=1wHolvdxHa3srDjPssOaABf1X-jUVRSfJ)|f2ed2b1cfbf89b6331f2a86065248767|
|SRR6756025_SRR6756028_P25.bam|      25%         |1,501,164                |499,074                  |2,000,238      | [bam](https://drive.google.com/open?id=1IptltVhzRO7UxdkGqLqJgkLPFrvJFvFI), [bai](https://drive.google.com/open?id=13hYujxRbD3z6rVRn5gGcU6FzrjC8r2e1)|3f838b28fa8898d12c6313de82473694|
|SRR6756025_SRR6756028_P30.bam|      30%         |1,401,194                |600,108                  |2,001,302      | [bam](https://drive.google.com/open?id=1wgGSXZacF4GHOylOEY9lYAOJvGluE4vm), [bai](https://drive.google.com/open?id=1cNhU1l80YDrkgDOQ0pbQI-eTEmF3y2yg)|2b9392ef983e1e789b5d818c9f604b19|
|SRR6756025_SRR6756028_P35.bam|      35%         |1,300,672                |700,448                  |2,001,120      | [bam](https://drive.google.com/open?id=1eJysW6twpyykkKRcA0vD4Y-2--yJmPVP), [bai](https://drive.google.com/open?id=1FxLrk6PVrGWV1HUGiQIKEt_Vy12fB9vl)|1d3a0fe4990bfa4da599df44d8f4c25a|
|SRR6756025_SRR6756028_P40.bam|      40%         |1,200,038                |799,920                  |1,999,958      | [bam](https://drive.google.com/open?id=1nTbFVqPJHQEZDo54lB4dY0F9RV-azxrR), [bai](https://drive.google.com/open?id=1ACX5f_xprOV-x2iIcAd-shEYzLieL6Wr)|e52bbeec296c2f933019a49e59639d46|
|SRR6756025_SRR6756028_P45.bam|      45%         |1,099,964                |899,424                  |1,999,388      | [bam](https://drive.google.com/open?id=1rIx5XpVx5epf58QAIm73ZcJ-D4XKUyIc), [bai](https://drive.google.com/open?id=1s8acjTvUi-ZSZtwCdwDQJT2qgJceQ5at)|72c8e3fca913735d1270875a3be59b17|
|SRR6756025_SRR6756028_P50.bam|      50%         | 999,998                 |998,942                  |1,998,940      | [bam](https://drive.google.com/open?id=1ejAisBka-CX7J_3sJz5lJbGJiISd3hMy), [bai](https://drive.google.com/open?id=15roTtddOUjC64zfcGJD3tQVDMSvkqzeL)|75562d5c2ae4db37d2c97fe09c368222|


## Running test

1. Download BAM files, BAM index files from the table above.
2. Download the VCF file from [here](https://drive.google.com/open?id=1kKIebM7iIGBHJpZdtNztJ3Ot5Na3VCd2).
3. Running Dcon to estimate contamination level.

```

#####################################
# Run Dcon on BAM file mixed at 5%  #
#####################################

$ Dcon.py -a 4 -b SRR6756025_SRR6756028_P05.bam -v gene83exon.snp.vcf.gz -o P05

@ 2018-03-19 21:52:08: Running Dcon v0.1.5
@ 2018-03-19 21:52:08: Reading VCF file "gene83exon.snp.vcf.gz" ...
@ 2018-03-19 21:52:08: Using 9256 SNPs in file "gene83exon.snp.vcf.gz"
@ 2018-03-19 21:52:09: Counting reads from "SRR6756025_SRR6756028_P05.bam"
@ 2018-03-19 22:00:56: Writing variants to "P05.SNP.tsv"
@ 2018-03-19 22:00:56: Filtering SNPs ...
@ 2018-03-19 22:00:56: Estimating contamination for each variant. Saved to "P05.SNP.PI.tsv"
@ 2018-03-19 22:00:56: Building Bayesian Gaussian Mixture model (BGMM) ...
	Converge status: True
	Iterations: 4

@ 2018-03-19 22:00:56: Summerzie BGMM model ...
	Means of homozygous component: 0.047493
	Means of heterozygous component: 0.745989
	Weight of homozygous component: 0.374071
	Weight of heterozygous component: 0.625929

@ 2018-03-19 22:00:56: Classify variants ...

@ 2018-03-19 22:00:56: Writing to P05.PI.tsv ...
@ 2018-03-19 22:00:56: Estimating overall contamination using maximum likelihood estimation (MLE) ...

	Overall contamination: 0.056


#####################################
# Run Dcon on BAM file mixed at 10% #
#####################################

$ Dcon.py -a 4 -b SRR6756025_SRR6756028_P10.bam -v gene83exon.snp.vcf.gz -o P10

@ 2018-03-19 21:41:28: Running Dcon v0.1.5
@ 2018-03-19 21:41:28: Reading VCF file "gene83exon.snp.vcf.gz" ...
@ 2018-03-19 21:41:28: Using 9256 SNPs in file "gene83exon.snp.vcf.gz"
@ 2018-03-19 21:41:28: Counting reads from "SRR6756025_SRR6756028_P10.bam"
@ 2018-03-19 21:50:19: Writing variants to "P10.SNP.tsv"
@ 2018-03-19 21:50:19: Filtering SNPs ...
@ 2018-03-19 21:50:19: Estimating contamination for each variant. Saved to "P10.SNP.PI.tsv"
@ 2018-03-19 21:50:19: Building Bayesian Gaussian Mixture model (BGMM) ...
	Converge status: True
	Iterations: 10

@ 2018-03-19 21:50:19: Summerzie BGMM model ...
	Means of homozygous component: 0.078710
	Means of heterozygous component: 0.703040
	Weight of homozygous component: 0.388537
	Weight of heterozygous component: 0.611463

@ 2018-03-19 21:50:19: Classify variants ...

@ 2018-03-19 21:50:19: Writing to P10.PI.tsv ...
@ 2018-03-19 21:50:20: Estimating overall contamination using maximum likelihood estimation (MLE) ...

	Overall contamination: 0.100


#####################################
# Run Dcon on BAM file mixed at 15% #
#####################################

$ Dcon.py -a 4 -b SRR6756025_SRR6756028_P15.bam -v gene83exon.snp.vcf.gz -o P15

@ 2018-03-19 21:55:51: Running Dcon v0.1.5
@ 2018-03-19 21:55:51: Reading VCF file "gene83exon.snp.vcf.gz" ...
@ 2018-03-19 21:55:51: Using 9256 SNPs in file "gene83exon.snp.vcf.gz"
@ 2018-03-19 21:55:51: Counting reads from "SRR6756025_SRR6756028_P15.bam"
@ 2018-03-19 22:05:02: Writing variants to "P15.SNP.tsv"
@ 2018-03-19 22:05:02: Filtering SNPs ...
@ 2018-03-19 22:05:02: Estimating contamination for each variant. Saved to "P15.SNP.PI.tsv"
@ 2018-03-19 22:05:02: Building Bayesian Gaussian Mixture model (BGMM) ...
	Converge status: True
	Iterations: 17

@ 2018-03-19 22:05:02: Summerzie BGMM model ...
	Means of homozygous component: 0.121970
	Means of heterozygous component: 0.687535
	Weight of homozygous component: 0.387907
	Weight of heterozygous component: 0.612093

@ 2018-03-19 22:05:02: Classify variants ...

@ 2018-03-19 22:05:02: Writing to P15.PI.tsv ...
@ 2018-03-19 22:05:02: Estimating overall contamination using maximum likelihood estimation (MLE) ...

	Overall contamination: 0.166


#####################################
# Run Dcon on BAM file mixed at 20% #
#####################################

$ Dcon.py -a 4 -v gene83exon.snp.vcf.gz -b SRR6756025_SRR6756028_P20.bam -o P20

@ 2018-03-19 22:04:26: Running Dcon v0.1.5
@ 2018-03-19 22:04:26: Reading VCF file "gene83exon.snp.vcf.gz" ...
@ 2018-03-19 22:04:26: Using 9256 SNPs in file "gene83exon.snp.vcf.gz"
@ 2018-03-19 22:04:26: Counting reads from "SRR6756025_SRR6756028_P20.bam"
@ 2018-03-19 22:13:34: Writing variants to "P20.SNP.tsv"
@ 2018-03-19 22:13:34: Filtering SNPs ...
@ 2018-03-19 22:13:34: Estimating contamination for each variant. Saved to "P20.SNP.PI.tsv"
@ 2018-03-19 22:13:34: Building Bayesian Gaussian Mixture model (BGMM) ...
	Converge status: True
	Iterations: 17

@ 2018-03-19 22:13:34: Summerzie BGMM model ...
	Means of homozygous component: 0.151335
	Means of heterozygous component: 0.658398
	Weight of homozygous component: 0.386997
	Weight of heterozygous component: 0.613003

@ 2018-03-19 22:13:34: Classify variants ...

@ 2018-03-19 22:13:34: Writing to P20.PI.tsv ...
@ 2018-03-19 22:13:34: Estimating overall contamination using maximum likelihood estimation (MLE) ...

	Overall contamination: 0.207


#####################################
# Run Dcon on BAM file mixed at 25% #
#####################################

$ Dcon.py -a 4 -v gene83exon.snp.vcf.gz  -b SRR6756025_SRR6756028_P25.bam -o  P25

@ 2018-03-19 22:08:07: Running Dcon v0.1.5
@ 2018-03-19 22:08:07: Reading VCF file "gene83exon.snp.vcf.gz" ...
@ 2018-03-19 22:08:07: Using 9256 SNPs in file "gene83exon.snp.vcf.gz"
@ 2018-03-19 22:08:07: Counting reads from "SRR6756025_SRR6756028_P25.bam"
@ 2018-03-19 22:17:12: Writing variants to "P25.SNP.tsv"
@ 2018-03-19 22:17:12: Filtering SNPs ...
@ 2018-03-19 22:17:12: Estimating contamination for each variant. Saved to "P25.SNP.PI.tsv"
@ 2018-03-19 22:17:12: Building Bayesian Gaussian Mixture model (BGMM) ...
	Converge status: True
	Iterations: 24

@ 2018-03-19 22:17:12: Summerzie BGMM model ...
	Means of homozygous component: 0.171741
	Means of heterozygous component: 0.632972
	Weight of homozygous component: 0.377965
	Weight of heterozygous component: 0.622035

@ 2018-03-19 22:17:12: Classify variants ...

@ 2018-03-19 22:17:12: Writing to P25.PI.tsv ...
@ 2018-03-19 22:17:12: Estimating overall contamination using maximum likelihood estimation (MLE) ...

	Overall contamination: 0.232


#####################################
# Run Dcon on BAM file mixed at 30% #
#####################################

$ Dcon.py -a 4 -v gene83exon.snp.vcf.gz  -b SRR6756025_SRR6756028_P30.bam -o  P30

@ 2018-03-19 22:34:46: Running Dcon v0.1.5
@ 2018-03-19 22:34:46: Reading VCF file "gene83exon.snp.vcf.gz" ...
@ 2018-03-19 22:34:46: Using 9256 SNPs in file "gene83exon.snp.vcf.gz"
@ 2018-03-19 22:34:47: Counting reads from "SRR6756025_SRR6756028_P30.bam"
@ 2018-03-19 22:43:41: Writing variants to "P30.SNP.tsv"
@ 2018-03-19 22:43:41: Filtering SNPs ...
@ 2018-03-19 22:43:41: Estimating contamination for each variant. Saved to "P30.SNP.PI.tsv"
@ 2018-03-19 22:43:41: Building Bayesian Gaussian Mixture model (BGMM) ...
	Converge status: True
	Iterations: 35

@ 2018-03-19 22:43:42: Summerzie BGMM model ...
	Means of homozygous component: 0.205044
	Means of heterozygous component: 0.618505
	Weight of homozygous component: 0.397182
	Weight of heterozygous component: 0.602818

@ 2018-03-19 22:43:42: Classify variants ...

@ 2018-03-19 22:43:42: Writing to P30.PI.tsv ...
@ 2018-03-19 22:43:42: Estimating overall contamination using maximum likelihood estimation (MLE) ...

	Overall contamination: 0.278
	
	
#####################################
# Run Dcon on BAM file mixed at 35% #
#####################################

$ Dcon.py -a 4 -v gene83exon.snp.vcf.gz  -b SRR6756025_SRR6756028_P35.bam -o  P35

@ 2018-03-19 22:22:12: Running Dcon v0.1.5
@ 2018-03-19 22:22:12: Reading VCF file "gene83exon.snp.vcf.gz" ...
@ 2018-03-19 22:22:12: Using 9256 SNPs in file "gene83exon.snp.vcf.gz"
@ 2018-03-19 22:22:12: Counting reads from "SRR6756025_SRR6756028_P35.bam"
@ 2018-03-19 22:31:24: Writing variants to "P35.SNP.tsv"
@ 2018-03-19 22:31:24: Filtering SNPs ...
@ 2018-03-19 22:31:24: Estimating contamination for each variant. Saved to "P35.SNP.PI.tsv"
@ 2018-03-19 22:31:24: Building Bayesian Gaussian Mixture model (BGMM) ...
	Converge status: True
	Iterations: 95

@ 2018-03-19 22:31:24: Summerzie BGMM model ...
	Means of homozygous component: 0.242770
	Means of heterozygous component: 0.629033
	Weight of homozygous component: 0.463618
	Weight of heterozygous component: 0.536382

@ 2018-03-19 22:31:24: Classify variants ...

@ 2018-03-19 22:31:24: Writing to P35.PI.tsv ...
@ 2018-03-19 22:31:24: Estimating overall contamination using maximum likelihood estimation (MLE) ...

	Overall contamination: 0.354


#####################################
# Run Dcon on BAM file mixed at 40% #
#####################################

$ Dcon.py -a 4 -v gene83exon.snp.vcf.gz  -b SRR6756025_SRR6756028_P40.bam -o  P40

@ 2018-03-19 22:22:18: Running Dcon v0.1.5
@ 2018-03-19 22:22:18: Reading VCF file "gene83exon.snp.vcf.gz" ...
@ 2018-03-19 22:22:18: Using 9256 SNPs in file "gene83exon.snp.vcf.gz"
@ 2018-03-19 22:22:18: Counting reads from "SRR6756025_SRR6756028_P40.bam"
@ 2018-03-19 22:31:29: Writing variants to "P40.SNP.tsv"
@ 2018-03-19 22:31:29: Filtering SNPs ...
@ 2018-03-19 22:31:29: Estimating contamination for each variant. Saved to "P40.SNP.PI.tsv"
@ 2018-03-19 22:31:29: Building Bayesian Gaussian Mixture model (BGMM) ...
	Converge status: True
	Iterations: 49

@ 2018-03-19 22:31:29: Summerzie BGMM model ...
	Means of homozygous component: 0.286721
	Means of heterozygous component: 0.707425
	Weight of homozygous component: 0.373812
	Weight of heterozygous component: 0.626188

@ 2018-03-19 22:31:29: Classify variants ...

@ 2018-03-19 22:31:29: Writing to P40.PI.tsv ...
@ 2018-03-19 22:31:29: Estimating overall contamination using maximum likelihood estimation (MLE) ...

	Overall contamination: 0.418


#####################################
# Run Dcon on BAM file mixed at 45% #
#####################################

$ Dcon.py -a 4 -v gene83exon.snp.vcf.gz  -b SRR6756025_SRR6756028_P45.bam -o  P45

@ 2018-03-19 22:35:03: Running Dcon v0.1.5
@ 2018-03-19 22:35:03: Reading VCF file "gene83exon.snp.vcf.gz" ...
@ 2018-03-19 22:35:04: Using 9256 SNPs in file "gene83exon.snp.vcf.gz"
@ 2018-03-19 22:35:04: Counting reads from "SRR6756025_SRR6756028_P45.bam"
@ 2018-03-19 22:43:58: Writing variants to "P45.SNP.tsv"
@ 2018-03-19 22:43:58: Filtering SNPs ...
@ 2018-03-19 22:43:58: Estimating contamination for each variant. Saved to "P45.SNP.PI.tsv"
@ 2018-03-19 22:43:58: Building Bayesian Gaussian Mixture model (BGMM) ...
	Converge status: True
	Iterations: 22

@ 2018-03-19 22:43:58: Summerzie BGMM model ...
	Means of homozygous component: 0.303108
	Means of heterozygous component: 0.766116
	Weight of homozygous component: 0.287331
	Weight of heterozygous component: 0.712669

@ 2018-03-19 22:43:58: Classify variants ...

@ 2018-03-19 22:43:58: Writing to P45.PI.tsv ...
@ 2018-03-19 22:43:58: Estimating overall contamination using maximum likelihood estimation (MLE) ...

	Overall contamination: 0.435


#####################################
# Run Dcon on BAM file mixed at 50% #
#####################################

$ Dcon.py -a 4 -v gene83exon.snp.vcf.gz  -b SRR6756025_SRR6756028_P50.bam -o  P50

@ 2018-03-19 22:50:13: Running Dcon v0.1.5
@ 2018-03-19 22:50:13: Reading VCF file "gene83exon.snp.vcf.gz" ...
@ 2018-03-19 22:50:13: Using 9256 SNPs in file "gene83exon.snp.vcf.gz"
@ 2018-03-19 22:50:14: Counting reads from "SRR6756025_SRR6756028_P50.bam"
@ 2018-03-19 22:59:12: Writing variants to "P50.SNP.tsv"
@ 2018-03-19 22:59:12: Filtering SNPs ...
@ 2018-03-19 22:59:12: Estimating contamination for each variant. Saved to "P50.SNP.PI.tsv"
@ 2018-03-19 22:59:12: Building Bayesian Gaussian Mixture model (BGMM) ...
	Converge status: True
	Iterations: 11

@ 2018-03-19 22:59:12: Summerzie BGMM model ...
	Means of homozygous component: 0.292286
	Means of heterozygous component: 0.761119
	Weight of homozygous component: 0.299293
	Weight of heterozygous component: 0.700707

@ 2018-03-19 22:59:12: Classify variants ...

@ 2018-03-19 22:59:12: Writing to P50.PI.tsv ...
@ 2018-03-19 22:59:12: Estimating overall contamination using maximum likelihood estimation (MLE) ...

	Overall contamination: 0.428
				
```

## Output files

Dcon will generate 3 output files
1. **prefix.SNP.tsv**
  * column-1: Chromosome ID
  * column-2: SNP Position 
  * column-3: Allele 1
  * column-4: Number of reads supporting Allele 1
  * column-5: Allele 2
  * column-6: Number of reads supporting Allele 1

2. **prefix.PI.tsv**  
  * column-1: Chromosome ID
  * column-2: SNP Position 
  * column-3: Allele 1
  * column-4: Number of reads supporting Allele 1
  * column-5: Allele 2
  * column-6: Number of reads supporting Allele 2
  * column-7: Ratio. {Number of reads supporting less frequent allele}/{Total reads covering this SNP}  
  * column-8: Probability of homozygous SNP
  * column-9: Probability of heterozygous SNP
  * column-10: Label. "Hom" or "Het"
  * column-11: Contamination percentage calculated from this individual SNP.

3. **prefix.overall_contamination.R**
  * R script to generate **prefix.overall_contamination.pdf**
  


```
$ head P05.SNP.tsv

1	36937059	A	126	G	117
1	36937065	A	127	G	124
2	25469502	T	116	C	95
3	105439026	G	477	A	11
3	128199380	G	220	A	14
3	128199662	G	170	A	11
3	128204951	C	97	T	90
3	128205860	C	30	G	15
3	128206618	C	36	A	7
3	136056033	G	671	A	5
...

$ head P05.PI.tsv

Chrom	Ref_pos	Allele_1	Allele_1_count	Allele_2	Allele_2_count	Ratio	Prob_of_Het	Prob_of_Hom	Label	Contamination_p
1	36937059	A	126	G	117	0.928571428571	1.0	1.17505373058e-32	Het	0.037037037037
1	36937065	A	127	G	124	0.976377952756	1.0	4.08974298614e-36	Het	0.0119521912351
2	25469502	T	116	C	95	0.818965517241	1.0	2.52236024344e-25	Het	0.0995260663507
3	105439026	G	477	A	11	0.0230607966457	0.000202183057686	0.999797816942	Hom	0.0450819672131
3	128199380	G	220	A	14	0.0636363636364	0.000486510687758	0.999513489312	Hom	0.119658119658
3	128199662	G	170	A	11	0.0647058823529	0.000499670278562	0.999500329721	Hom	0.121546961326
3	128204951	C	97	T	90	0.927835051546	1.0	1.32461743765e-32	Het	0.0374331550802
3	128205860	C	30	G	15	0.5	0.999999989858	1.01421802114e-08	Het	0.333333333333
3	128206618	C	36	A	7	0.194444444444	0.0470112000123	0.952988799988	Hom	0.325581395349
...

```  

![workflow](https://github.com/liguowang/dcon/blob/master/img/P05_overall_contamination.png?raw=true)


## speed benchmark

**CPU model**

Intel(R) Core(TM) i7-3720QM CPU @ 2.60GHz

**Speed**:

When VCF (9265 candidate SNPs within captured regions) file was provided , Dcon took **6 minutes and 5 seconds** to estimate DNA contamination from a 177 Mb BAM file containing 2 million alignments. Using multiple threads (-p) will increase speed significantly. 

## Performance

In below figure, each boxplot is generated from 100 BAM files (Capture sequencing target 37 genes). And each BAM file is purposely mixed from two different samples at pre-defined percentages (10%, 20%, 30%, 40%, 50%).
For example, in the leftmost boxplot(green), 10% of reads were from one sample, and 90% of reads were from another sample.
 
![](https://github.com/liguowang/dcon/blob/master/img/strip_chart.png?raw=true)

## Contributing Authors

1. Liguo Wang <Wang.Liguo@mayo.edu>
2. Tao Ma <Ma.Tao@mayo.edu>
3. Rohan Gnanaolivu <gnanaolivu.rohandavid@mayo.edu>
4. Jagadheshwar Balan <Balan.Jagadheshwar@mayo.edu>
5. Eric Klee Klee <Klee.Eric@mayo.edu>
6. Jean-Pierre Kocher <kocher.jeanpierre@mayo.edu>

## License
This project is licensed under the GPL License - see the [LICENSE](LICENSE) file for details

## Acknowledgments

This project is partially supported by the Clinical Genome Sequencing Laboratory (CGSL), Mayo Clinic.