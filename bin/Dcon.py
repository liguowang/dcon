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
import pysam
import operator
import collections
from numpy import mean
#from scipy.stats import binom_test
from multiprocessing import Process, Manager, current_process
from DconModule import ireader
from DconModule.utils import *
from DconModule.bgmm import build_GMM


__author__ = "Liguo Wang"
__contributor__="Liguo Wang"
__copyright__ = "Copyright 2017-2018, Mayo Clinic"
__credits__ = []
__license__ = "GPL"
__version__="0.1.5"
__maintainer__ = "Liguo Wang"
__email__ = "wang.liguo@mayo.edu; wangliguo78@gmail.com"
__status__ = "4 - Beta"


		
def bamsnv(called_SNPs, sam_file, g_chr, g_st, g_end, min_mapq, min_seqq, min_cvg, min_alleleRC):
	'''
	calling SNPs from BAM file. (g_st, g_end]
	* called_SNPs: shared variable. Send back results to main()
	* g_chr: chromsome ID
	* g_start: start postion (0-based, not included)
	* g_end: end position (1-based, included)
	* min_mapq: minimum mapping quality score in phred sclale
	* min_seqq: minimum sequencing quality score in phred sclale
	* min_cvg: minimum read coverage.
	* min_alleleRC: allele's minimum read count. Alleles with read count than this will not report.
	* min_allelePercent: allele's minimum read percent. Alleles with read percent less than this will not report. 
	'''
	samfile = pysam.Samfile(sam_file,'rb')
	name = current_process().name
	
	for pileupcolumn in samfile.pileup(g_chr, g_st, g_end, truncate=True):
		#not enough read coverage
		if pileupcolumn.nsegments < min_cvg:
			continue
			
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
			if read.mapping_quality == 255: continue
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
		
		# filter overall coverage
		tmp = bfreq.values()
		if sum(tmp) < min_cvg:
			continue
		
		# filter allele's minimum read count
		if min(tmp) < min_alleleRC:
			continue

		
		output = [g_chr, ref_pos]
		#sorted alleles by frequency (in decreasing order)
		for i,j in sorted(bfreq.items(), key = lambda x:x[1], reverse=True):
			output.append(i)
			output.append(j)
						
		called_SNPs.append(output)	#[g_chr, ref_pos, base1, base1_count, base2, base2_count]
	#print >>sys.stderr, "@ " + strftime("%Y-%m-%d %H:%M:%S") + ": " + name + " -- fisinshed"
	
	
	
def main():
	usage="%prog [options]" + '\n' + __doc__ + "\n"
	parser = OptionParser(usage,version="%prog " + __version__)
	parser.add_option("-b","--bam",action="store",type="string",dest="bam_file", help="Alignment file in BAM format. BAM file must be indexed and sorted using samTools (http://samtools.sourceforge.net/).")	
	parser.add_option("-v","--variant",action="store",type="string",dest="vcf_file",help="VCF format file defining SNPs. VCF file can be a plain text, compressed (.gz, .z, .Z, .bz, .bz2 and .gzip2) or remote file (http://, https://, ftp://). Exclusive with '-r'. ")
	parser.add_option("-o","--output",action="store",type="string",dest="output_file", help="Prefix of output file. Will generate three files: \"prefix.SNP.tsv\", \"prefix.PI.tsv\" \"prefix.dcon.R\". \"prefix.dcon.pdf\" will be generated if R exists.")
	parser.add_option("-q","--mapq",action="store",type="int",dest="min_mapq",default=30,help="Minimum mapping quality (http://maq.sourceforge.net/qual.shtml). Mapping quality is Phred-scaled probability of an alignment being wrong. default=%default")	
	parser.add_option("-m","--seqq",action="store",type="int",dest="min_seqq",default=30,help="Minimum base phred quality (http://maq.sourceforge.net/qual.shtml). Base quality is Phread-scaled probability of a base calling being wrong. default=%default")	
	parser.add_option("-c","--cvg",action="store",type="int",dest="min_coverage",default=30,help="Minimum number of reads supporting variant. default=%default")	
	parser.add_option("-a","--allele",action="store",type="int",dest="min_allele_read_count",default=5,help="Minimum number of reads supporting one allele. default=%default")	
	parser.add_option("-p","--processor",action="store",type="int",dest="processor_num",default=1,help="Number of processes used to call SNPs. default=%default")	
	parser.add_option("-s","--skipXY",action="store_true",dest="skip_XY",default=False,help="Skip SNPs on X and Y chromosomes. default=%default")	
	parser.add_option("-y","--prob-cutoff",action="store",type="float",dest="probability_cut",default=0.5,help="Probability cutoff. default=%default")
	
	(options,args)=parser.parse_args()
	
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
		#[[chr,st,end,name],[chr2,st2,end2,name2],...]
		all_regions = read_vcf(options.bam_file,options.vcf_file, skipxy = options.skip_XY)	
		printlog(['Using %d SNPs in file \"%s\"' % (len(all_regions), options.vcf_file)])
	
	if len(all_regions) == 0:
		printlog(['Error: no genomic regins or SNPs. Exit!'])
		sys.exit(0)
		
		
	manager = Manager()
	#list of list. shared variable between main() and bamcnv()
	called_SNPs = manager.list()	
	
	if options.vcf_file:
		printlog(['Counting reads from \"%s\"' % options.bam_file])
		
	for chunk in chunks(l = all_regions, n = options.processor_num):
		jobs = []
		for chr, st, end, name in chunk:
			job_name = '\t'.join([str(i) for i in (chr,st,end,name)])
			p = Process(name = job_name,target = bamsnv, args = (called_SNPs, options.bam_file, chr, st, end, options.min_mapq, options.min_seqq, options.min_coverage, options.min_allele_read_count))
			jobs.append(p)
			p.start()

		for proc in jobs:
			proc.join()	
	
	printlog(['Writing variants to \"%s\"' % (options.output_file + '.SNP.tsv')])
	OUT_TAB1 = open(options.output_file + '.SNP.tsv', 'w')
	seen = set()
	for snp in called_SNPs:
		snp_coord = snp[0] + ':' + str(snp[1])
		if snp_coord not in seen:
			print >>OUT_TAB1, '\t'.join([str(i) for i in snp])
			seen.add(snp_coord)
	OUT_TAB1.close()	
	
	
	printlog(['Filtering SNPs ...'])	
	(snp_list, ratio_list) = filter_SNPs(infile = options.output_file + '.SNP.tsv')
	
	printlog(['Estimating contamination for each variant. Saved to \"%s\"' % (options.output_file + '.SNP.PI.tsv')])	
	build_GMM(outfile = options.output_file + '.PI.tsv', d = ratio_list, names=snp_list, n=2, iter=1000, tol=0.001, rnd=99, prob_cut = options.probability_cut)
	
	printlog(['Estimating overall contamination using maximum likelihood estimation (MLE) ...'])
	pi = estimate_overall_pi(infile = options.output_file + '.PI.tsv')
	#print pi
	print "\n\tOverall contamination: %f"  % pi


if __name__=='__main__':
	main()
	
