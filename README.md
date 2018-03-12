# Dcon: Estimate DNA contamination

Next-generation sequencing (NGS) of DNA is routinely used in clinical testing. However,
the complicated, multistep DNA extraction protocol and the unknown pre-analytic specimen
handling often lead to sample-to-sample DNA contaminations. DNA contamination poses
substantial challenges in NGS assays and affects diagnosis and treatment if contaminant
DNA harbors clinically significant variants.

**Dcon** is a program to detect and quantify sample-to-sample germline DNA contaminations using the
skewed (or shifted) distribution of variant allele frequency (VAF, VAF refers to the
fraction of sequencing reads overlapping a genomic coordinate of a DNA variant such as SNP).

**Unlike other tools, Dcon does NOT require any prior knowledge such as "genotype" and/or "allele frequencies" calculated from 1000 Genome project. Dcon works on single BAM file and does not require tumor-normal pairs.**

## Getting Started

### How Dcon works

![workflow](https://github.com/liguowang/dcon/blob/master/img/worflow.png?raw=true)


### Prerequisites
Dcon requires [python2.7](https://www.python.org/download/releases/2.7/) and these python packages. 

1. [**numpy**](http://www.numpy.org/)
2. [**scipy**](https://www.scipy.org/)
3. [**pysam**](https://pypi.python.org/pypi/pysam)
4. [**scikit-learn**](https://pypi.python.org/pypi/scikit-learn)
5. [**pyBigWig**](https://pypi.python.org/pypi/pyBigWig)

### Installing

Following instructions [here](https://pip.pypa.io/en/stable/installing/) to install **pip**. Please note: **pip** is already installed if you're using Python 2 >=2.7.9 or Python 3 >=3.4 binaries downloaded from [python.org](https://www.python.org/), but you'll need to upgrade **pip**. 

Users have three ways to install Dcon:

1. Install Dcon from pypi:
  * Open a terminal and type command: ```pip install Dcon```

2. Install Dcon from local directory
  * Download Dcon directly from: [this link](https://github.com/liguowang/dcon/archive/master.zip)
  * Go to the directory where **dcon-master.zip** was saved.
  * Open a terminal and type command:`pip install dcon-master.zip`

3. Install Dcon from local directory
  * Clone Dcon. Open a terminal and type command: `git clone https://github.com/liguowang/dcon.git`
  * Type command: `cd ./dcon`
  * Type command: `pip install .` 

### General usage

Dcon has two run modes:

1. Use provides a [BAM](https://genome.ucsc.edu/FAQ/FAQformat.html#format5.1) file and a [BED](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) file. For targeted sequencing, BED file defines the captured genomic regions. For whole exome sequencing, BED file defines the whole exons. 

   `python2.7 Dcon.py -b input.bam -r input.bed -o output`
    
   This mode is not only slow but also produces less accurate estimation. Only use this mode if no known SNPs were found in the captured regions. 
    
2. Use provides a [BAM](https://genome.ucsc.edu/FAQ/FAQformat.html#format5.1) file and a [VCF](https://genome.ucsc.edu/FAQ/FAQformat.html#format10.1) file. VCF file can be called from this particular BAM file or extracted from [dbSNP](https://www.ncbi.nlm.nih.gov/SNP/) (i.e. extract all SNPs that falling into captured genomic regions, Dcon automatically checks the VAF for each SNP and skips irrelevant SNPs).

   `python2.7 Dcon.py -b input.bam -v input.vcf -o output`

   This mode is fast and gives accurate estimation [recommended]. 
   
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
  -r BED_FILE, --region=BED_FILE
                        BED format file defining genomic regions from which
                        SNPs will be called. The first 3 columns ('chrom_ID',
                        'start', 'end') are required. BED file can be a plain
                        text, compressed (.gz, .z, .Z, .bz, .bz2 and .gzip2)
                        or remote file (http://, https://, ftp://). Note: if
                        genomic regions are overlapped, some SNPs might be
                        reported multiple times. Exclusive with '-v'.
  -v VCF_FILE, --variant=VCF_FILE
                        VCF format file defining SNPs. VCF file can be a plain
                        text, compressed (.gz, .z, .Z, .bz, .bz2 and .gzip2)
                        or remote file (http://, https://, ftp://). Exclusive
                        with '-r'.
  -o OUTPUT_FILE, --output=OUTPUT_FILE
                        Prefix of output file. Will generate three files:
                        "prefix.SNP.tsv", "prefix.PI.tsv" "prefix.dcon.R".
                        "prefix.dcon.pdf" will be generated if R exists.
  -u BIGWIG_FILE, --alignability=BIGWIG_FILE
                        Alignability file in BigWig format. ENCODE
                        alignability score S (0 < S <= 1), measures how often
                        a read found at the particular location will align
                        within the whole genome. S = 1/(number of matches
                        found in the genome): S=1 means one match in the
                        genome, S=0.5 is two matches in the genome, and so on.
                        Note: BAM file, BED (or VCF) file, and Alignability
                        file should base on the same reference genome.
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
  -a MIN_ALLLE_READ_COUNT, --allele=MIN_ALLLE_READ_COUNT
                        Minimum number of reads supporting one allele.
                        default=5
  -p PROCESSOR_NUM, --processor=PROCESSOR_NUM
                        Number of processes used to call SNPs. default=1
  -s, --skipXY          Skip SNPs on X and Y chromosomes. default=False
  -x SNP_NUM, --snp-num=SNP_NUM
                        Number of SNPs used to calcualte contamination. Only
                        applied to targeted sequencing. default=10
  -w ALIGNABILITY_SCORE, --align-cutoff=ALIGNABILITY_SCORE
                        Alignability score cutoff. An SNP will be filtered out
                        if the average alignability score (of the 20-bp window
                        centering on this SNP) is less than this cutoff.
                        default=1
  -y PROBABILITY_CUT, --prob-cutoff=PROBABILITY_CUT
                        Probability cutoff. default=0.5

```

## Running the tests

## speed benchmark

**CPU model**

Intel(R) Core(TM) i7-3720QM CPU @ 2.60GHz

**Speed**:

1. When VCF file was provided (940 candidate SNPs), Dcon took **1 minutes and 51 seconds** to estimate DNA contamination from a 318Mb BAM file containing 4,436,375 alignments.
2. When BED file was provided, Dcon took **32 minutes and 2 seconds** to call variants and estimate DNA contamination from the same BAM file mentioned above.

## Performance

In below figure, each boxplot is generated from 100 BAM files. And each BAM file is purposely mixed from two different samples at pre-defined percentages (10%, 20%, 30%, 40%, 50%).
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