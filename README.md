# Dcon: Estimate DNA contamination

Next-generation sequencing (NGS) of DNA is routinely used in clinical testing. However,
the complicated, multistep DNA extraction protocol and the unknown pre-analytic specimen
handling often lead to sample-to-sample DNA contaminations. DNA contamination poses
substantial challenges in NGS assays and affects diagnosis and treatment if contaminant
DNA harbors clinically significant variants.

**Dcon** is a program to detect and quantify sample-to-sample germline DNA contaminations using the
skewed (or shifted) distribution of variant allele frequency (VAF, VAF refers to the
fraction of sequencing reads overlapping a genomic coordinate of a DNA variant such as SNP).  

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

1. Following instructions [here](https://pip.pypa.io/en/stable/installing/) to install **pip**. Please note: **pip** is already installed if you're using Python 2 >=2.7.9 or Python 3 >=3.4 binaries downloaded from [python.org](https://www.python.org/), but you'll need to upgrade **pip**. 
2. Use `pip install Dcon` to install.

## Running the tests

## Deployment

## Contributing

## Authors

## License
This project is licensed under the GPL License - see the [LICENSE](LICENSE) file for details

## Acknowledgments