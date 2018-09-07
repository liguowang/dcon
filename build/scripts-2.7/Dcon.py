#!/usr/bin/env python

"""
Dcon: Estimate DNA Contamination from BAM file. 
"""

import os,sys
if sys.version_info[0] != 2 or sys.version_info[1] != 7:
        print >>sys.stderr, "\nYou are using python" + str(sys.version_info[0]) + '.' + \
        str(sys.version_info[1]) + " This program needs python2.7!\n"
        sys.exit()

try:
	import pysam
except ImportError:
	print >>sys.stderr, "Import Error: Please install pysam (https://pypi.python.org/pypi/pysam)!"
	
from optparse import OptionParser
from time import strftime
import subprocess
import collections
import operator
import pysam
from numpy import mean
from multiprocessing import Process, Manager, current_process
from DconModule import ireader
from DconModule.utils import *
from DconModule.bgmm import *


__author__ = "Liguo Wang"
__contributor__="Liguo Wang"
__copyright__ = "Copyright 2017-2018, Mayo Clinic"
__credits__ = []
__license__ = "GPL"
__version__="0.1.8"
__maintainer__ = "Liguo Wang"
__email__ = "wangliguo78@gmail.com"
__status__ = "4 - Beta"
	
def bamsnv(called_SNPs, sam_file, g_chr, g_st, g_end, min_mapq, min_seqq):
	'''
	
	Calculate variant allele frequency from BAM file.
	
	* called_SNPs: shared variable. Send back results to main()
	* g_chr: chromsome ID
	* g_start: start postion (0-based, not included)
	* g_end: end position (1-based, included)
	* min_mapq: minimum mapping quality score in phred sclale
	* min_seqq: minimum sequencing quality score in phred sclale
	'''
	samfile = pysam.Samfile(sam_file,'rb')
	name = current_process().name
	
	for pileupcolumn in samfile.pileup(g_chr, g_st, g_end, truncate=True):
				
		ref_pos = pileupcolumn.reference_pos + 1
		
		# all read bases at a particular position
		read_alleles = []							
				
		for pileupread in pileupcolumn.pileups:
			read = pileupread.alignment
			if read.is_qcfail: continue
			if read.is_secondary: continue
			if read.is_unmapped: continue
			if read.is_duplicate: continue
			
			# filter mapping quality
			#if read.mapping_quality == 255: continue	#255 indicates mapping quality is not available
			if read.mapping_quality < min_mapq: continue	

			# filter base quality
			if pileupread.query_position is None:
				continue
			base_qual = ord(pileupread.alignment.qual[pileupread.query_position])-33
			if base_qual < min_seqq:
				continue
			
			# base must be A/C/G/T
			base = pileupread.alignment.seq[pileupread.query_position].upper()
			if base not in ('A','C','G','T'): continue
						
			read_alleles.append(base)
							
		#bfreq is a dict with base as key, base frequency (count) as value
		# Counter({'T': 53, 'A': 12})
		bfreq = collections.Counter(read_alleles)	
		
		# only consider bi-allelic sites
		if len(bfreq) != 2:
			continue		
		
		output = [g_chr, ref_pos]
		#sorted alleles by frequency (in decreasing order)
		for i,j in sorted(bfreq.items(), key = lambda x:x[1], reverse=True):
			output.append(i)
			output.append(j)
						
		called_SNPs.append(output)	#[g_chr, ref_pos, base1, base1_count, base2, base2_count]	


def main():
	usage="%prog [options]" + '\n' + __doc__ + "\n"
	parser = OptionParser(usage,version="%prog " + __version__)
	parser.add_option("-b","--bam",action="store",type="string",dest="bam_file", help="Alignment file in BAM format. BAM file must be indexed and sorted using samTools (http://samtools.sourceforge.net/).")	
	parser.add_option("-v","--variant",action="store",type="string",dest="vcf_file",help="VCF format file containing *candiate* SNPs (Note: Not every SNP in this file has to be the real SNP in your BAM file. They could be SNPs extracted from dbSNP/1000Genome database, as long as they are located in the captured region. However, too many spurious SNPs will reduce the accuracy.). VCF file can be a plain text, compressed (.gz, .z, .Z, .bz, .bz2 and .gzip2) or remote file (http://, https://, ftp://). Exclusive with '-r'. ")
	parser.add_option("-o","--output",action="store",type="string",dest="output_file", help="Prefix of output file. Will generate four files: \"prefix.SNP.tsv\", \"prefix.PI.tsv\",\"prefix.overall_contamination.R\".  \"prefix.overall_contamination.pdf\" will also be generated if R exists.")
	parser.add_option("-m","--mapq",action="store",type="int",dest="min_mapq",default=30,help="Minimum mapping quality (http://maq.sourceforge.net/qual.shtml). Mapping quality is Phred-scaled probability of an alignment being wrong. default=%default")	
	parser.add_option("-q","--seqq",action="store",type="int",dest="min_seqq",default=30,help="Minimum base phred quality (http://maq.sourceforge.net/qual.shtml). Base quality is Phread-scaled probability of a base calling being wrong. default=%default")	
	parser.add_option("-c","--cvg",action="store",type="int",dest="min_coverage",default=30,help="Minimum number of reads supporting variant. default=%default")	
	parser.add_option("-n","--processor",action="store",type="int",dest="processor_num",default=1,help="Number of processes. default=%default")	
	parser.add_option("-p","--prob",action="store",type="float",dest="probability_cut",default=0.5,help="Cutoff of probability.  Probability is calcualted from Bayesian Gaussian Mixture model to decicde if a SNP is homozygous or heterozygous. default=%default")
	parser.add_option("-z","--zscore",action="store",type="float",dest="zscore_cut",default=2.5,help="Cutoff of Z-score. The modified Z-score is calculated from median absolute deviation. SNPs with Z-score greater than this cutoff will be considered as outlier and not used to estimate overall contamination level. default=%default")
	
	(options,args)=parser.parse_args()
	
	if options.probability_cut >= 1 or options.probability_cut <= 0:
		print >>sys.stderr, "\nProbability cutoff (-p) must be in (0,1).\n\n"
		sys.exit(0)
	if options.zscore_cut <= 0:
		print >>sys.stderr, "\nZ-score cutoff (-z) must be a positiv value.\n\n"		
		sys.exit(0)
	if options.min_mapq <= 0:
		print >>sys.stderr, "\nMinimum mapping quality score (-m) must be a positiv value.\n\n"	
		sys.exit(0)
	if options.min_seqq <= 0:
		print >>sys.stderr, "\Minimum sequencing quality score (-q) must be a positiv value.\n\n"	
		sys.exit(0)

	# check input BAM and bai
	if not (options.bam_file):
		print >>sys.stderr, "\nBAM file (-b) is required\n\n"
		parser.print_help()
		sys.exit(0)
	else:
		if not os.path.isfile(options.bam_file):
			print >>sys.stderr, "\nError: BAM file: %s does not exist!" % (options.bam_file)
			parser.print_help()
			sys.exit(0)
		if not os.path.isfile(options.bam_file + '.bai'):
			print >>sys.stderr, "\nError: BAM index file: %s does not exist!" % (options.bam_file+ '.bai')
			parser.print_help()
			sys.exit(0)

	# check output file
	if not (options.output_file):
		print >>sys.stderr, "\nThe prefix of output file (-o) is required\n\n"
		parser.print_help()
		sys.exit(0)			
	
	# check input VCF file
	if options.vcf_file:
		if not os.path.isfile(options.vcf_file):
			print >>sys.stderr, "\nError: VCF file: %s does not exist!" % (options.vcf_file)
			parser.print_help()
			sys.exit(0)


	print >>sys.stderr, "\n"		
	printlog(["Running Dcon", 'v'+ __version__])
	
	
	# targeted SNPs
	all_regions = []
	if options.vcf_file:
		printlog(['Reading VCF file \"%s\" ...' % options.vcf_file])	
		all_regions = read_vcf(options.bam_file,options.vcf_file)	#generator [chr,st,end,name],[chr2,st2,end2,name2],...
		
	manager = Manager()
	#list of list. shared variable between main() and bamcnv()
	called_SNPs = manager.list()	
	
	printlog(['Calculate variant allele frequencies from \"%s\"' % options.bam_file])
	count = 0
	for (chr, st, end, name) in all_regions:
		count += 1
		jobs = []
		job_name = '\t'.join([str(i) for i in (chr,st,end,name)])
		p = Process(name = job_name,target = bamsnv, args = (called_SNPs, options.bam_file, chr, st, end, options.min_mapq, options.min_seqq))
		jobs.append(p)
		p.start()
		
		if count == options.processor_num:
			for proc in jobs:
				proc.join()
				count = 0	
	
	printlog(['Writing %d variants to \"%s\"' % (len(called_SNPs), options.output_file + '.SNP.tsv')])
	OUT_TAB1 = open(options.output_file + '.SNP.tsv', 'w')
	seen = set()
	for snp in called_SNPs:
		snp_coord = snp[0] + ':' + str(snp[1])
		if snp_coord not in seen:
			print >>OUT_TAB1, '\t'.join([str(i) for i in snp])
			seen.add(snp_coord)
	OUT_TAB1.close()	
	
	
	printlog(['Read SNPs ...'])	
	(snp_list, ratio_list) = read_SNPs(infile = options.output_file + '.SNP.tsv',min_mvaf = 0.01, min_cvg = options.min_coverage, min_alleleRC=3)
	
	if len(snp_list) <= 3:
		print >>sys.stderr, "No SNPs passed selection criteria. Try to increase '-z' to get more candidate SNPs"
		print >>sys.stderr, "Aborted"
		sys.exit(0)
	
	printlog(['Estimating contamination for each variant. Saved to \"%s\"' % (options.output_file + '.PI.tsv')])	
	build_GMM(outfile = options.output_file + '.PI.tsv', d = ratio_list, names=snp_list, n=2, iter=1000, tol=0.0001, rnd=99, prob_cut = options.probability_cut)
	
	printlog(['Filter outlier SNPs ...'])	
		
	filter_outlier(options.output_file + '.PI.tsv', options.output_file + '.PI.filtered.tsv', zcut = options.zscore_cut)
	
	if os.path.exists(options.output_file + '.PI.filtered.tsv') and os.stat(options.output_file + '.PI.filtered.tsv').st_size > 0:
		printlog(['Delete \"%s\"' % (options.output_file + '.PI.tsv')])
		os.remove(options.output_file + '.PI.tsv')
		
	
	printlog(['Estimating overall contamination using maximum likelihood estimation (MLE) ...'])
	
	pi_hom = mle(infile = options.output_file + '.PI.filtered.tsv')
	
	print "\n\nOverall contamination level of %s is %.3f\n"  % (options.bam_file.replace('.bam',''), pi_hom)
	
if __name__=='__main__':
	main()

