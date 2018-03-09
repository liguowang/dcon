# Dcon: Estimate DNA contamination

Next-generation sequencing (NGS) of DNA is routinely used in clinical testing. However,
the complicated, multistep DNA extraction protocol and the unknown pre-analytic specimen
handling often lead to human-to-human DNA contaminations. DNA contamination poses
substantial challenges in NGS assays and affects diagnosis and treatment if contaminant
DNA harbors clinically significant variants.

**Dcon** is a program to detect and quantify sample-to-sample DNA contaminations using the
skewed (or shifted) distribution of Variant allele frequency (VAF, VAF refers to the
fraction of sequencing reads overlapping a genomic coordinate of a DNA variant such as SNP).  
