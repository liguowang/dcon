from time import strftime
import os,sys
import collections
import pysam
import ireader
from numpy import mean,std
from scipy import stats
from scipy.stats import binom
import subprocess
import operator


def get_reference_name(bamfile):
	'''
	get the name and length of all chromosome
	'''
	chroms = {}	#id to length
	samfile = pysam.Samfile(bamfile,'rb')
	a = samfile.header['SQ']
	for i in a:
		chroms[i['SN']] = i['LN']
	return chroms
	
def chunks(l, n):
	"""Yield successive n-sized chunks from l."""
	for i in range(0, len(l), n):
		yield l[i:i + n]
        
def printlog (mesg_lst):
	'''print message to STDERR'''
	if len(mesg_lst)==1:
		msg = "@ " + strftime("%Y-%m-%d %H:%M:%S") + ": " +  mesg_lst[0]
	else:
		msg = "@ " + strftime("%Y-%m-%d %H:%M:%S") + ": " + ' '.join([str(i) for i in mesg_lst])
	print >>sys.stderr,msg

def read_bed(bam_file, bed_file, skipxy):
	'''
	Read BED file, BED file could regular or compressed format or remote file
	Compressed files support: .gz, .z, .Z, .bz, .bz2, .gzip2
	Remote files support: http://, https://, ftp://
	skipxy <BOOLEAN>: skip regions on X and Y chromosomes
	'''
	
	chrom_infor = get_reference_name(bam_file)	#id to length
	regions = []
	for line in ireader.reader(bed_file):
		if line.startswith('#'):
			continue
		f = line.split()
		if len(f) < 3:
			continue
		chrom = f[0]
		if chrom not in chrom_infor:
			print >>sys.stderr, "Error: cannot find \"%s\" from the \"%s\" file. Please check your chromosome IDs in \"%s\"." % (chrom, bam_file, bed_file)
			sys.exit(0)
		if (skipxy is True) and (chrom.upper() in ['CHRX','CHRY','X','Y']):
			continue
		try:
			start = int(f[1])
			end = int(f[2])
		except ValueError:
			print >>sys.stderr, "Error in line %s, not a BED file" % line
			sys.exit(0)
		
		start = max(start,0)
		end = min(end, chrom_infor[chrom])

		try:
			name = '\t'.join(f[3:])
		except ValueError:
			name = ''

		
		regions.append((chrom, start, end, name))
	return regions

def read_vcf(bam_file, vcf_file, skipxy):
	'''
	Read VCF file. VCF file could regular or compressed format or remote file
	Compressed files support: .gz, .z, .Z, .bz, .bz2, .gzip2
	Remote files support: http://, https://, ftp://
	skipxy <BOOLEAN>: skip regions on X and Y chromosomes
	'''
	
	chrom_infor = get_reference_name(bam_file)	#id to length
	regions = []
	for line in ireader.reader(vcf_file):
		if line.startswith('#'):
			continue
		f = line.split()
		
		if f[0] in chrom_infor:
			chrom = f[0]
		else:
			chrom = 'chr' + f[0]
			if chrom not in chrom_infor:
				print >>sys.stderr, "Error: cannot find \"%s\"  (chromosome ID) from the \"%s\" file. Please check your chromosome IDs in \"%s\"." % (f[0], bam_file, bed_file)
				sys.exit(0)
			
		
		if (skipxy is True) and (chrom.upper() in ['CHRX','CHRY', 'X', 'Y']):
			continue
		try:
			end = int(f[1])
			start = end - 1
		except ValueError:
			continue
		try:
			name = f[2]
		except ValueError:
			name = ''
		
		start = max(start,0)
		end = min(end, chrom_infor[chrom])
		
		regions.append((chrom, start, end, name))
	return regions

#=======

def filter_SNPs(infile, bigwig_file = None, cutoff = 0):
	'''
	cutoff: Alignability score (between 0 and 1) cutoff. 
	'''
	
	ratios = []
	SNPs = []
	
	for l in open(infile,'r'):
		l = l.strip()
		f = l.split()
		
		#filter SNPs whose alignability score less than cutoff. 
		if bigwig_file:
			bw = pyBigWig.open(bigwig_file)
			chrom = f[0]
			SNP_pos = int(f[1])
			start = SNP_pos - 11
			end = SNP_pos + 10
			if start < 0:
				start = 0
			if end > bw.chroms(chrom):
				end = bw.chroms(chrom)
			avg_alignability = bw.stats(chrom, start, end)
			if avg_alignability[0] < cutoff:
				continue
		
		base1_count = float(f[3])
		base2_count = float(f[5])
		
		#Minor variant allele fraction (or percentage)
		MVAF = min(base1_count, base2_count)/(base1_count + base2_count)
		
		if MVAF < 0.01:
			continue
		# Minor variant allele fraction/Major variant allele fraction ratio. 
		MVAR = min(base1_count, base2_count)/max(base1_count, base2_count)	

		ratios.append(MVAR)
		SNPs.append('\t'.join(f[0:6]))
	return(SNPs, ratios)

#=======	
			
def estimate_overall_pi(infile):
	'''
	Estimate the overall contamination percentage using maximum likelihood estimation (MLE)
	'''
	candidate_PIs = [i/1000.0 for i in range(0,501)]
	snp_hom = []
	snp_het = []
	for l in open(infile,'r'):
		l = l.strip()
		if l.startswith('Chrom'):continue
		f = l.split()
		allele_1_count = int(f[3])
		allele_2_count = int(f[5])
		if allele_1_count < allele_2_count:
			(allele_1_count, allele_2_count) = (allele_2_count,allele_1_count)
		
		if f[9] == 'Het':
			continue
			#snp_het.append([allele_1_count + allele_2_count, allele_2_count, 'Het']) #n,k
		elif f[9] == 'Hom':
			snp_hom.append([allele_1_count + allele_2_count, allele_2_count, 'Hom']) #n,k
		else:
			continue
	
	#key is pi, value is the product of pmfs across all snps
	
	
	# hom SNPs
	prob = -float("inf")
	pi_of_max_prob1 = 0
	for pi in candidate_PIs:
		# Homozygous SNP contaminated by homozygous SNP
		#p1 = pi	
		# Homozygous SNP contaminated by heterozygous SNP
		p2 = pi/2.0
		
		joint_prob = 0
		for n,k,t in snp_hom:
			#pmf_1 = binom.logpmf(k,n,p1)
			pmf_2 = binom.logpmf(k,n,p2)
			#joint_prob += max(pmf_1, pmf_2)
			joint_prob += pmf_2
			
		if joint_prob > prob:
			prob = joint_prob
			pi_of_max_prob1 = pi
	
	
	return pi_of_max_prob1
	#print >>sys.stderr, "Contamination estimated from Homozygous SNPs: %f" % 	pi_of_max_prob1		
	
	"""
	# het SNPs
	prob = -float("inf")
	pi_of_max_prob2 = 0	
	for pi in candidate_PIs:
		p = (1.0 - pi)/2.0
		joint_prob = 0
		for n,k,t in snp_het:
			joint_prob += binom.logpmf(k,n,p)
		
		if joint_prob > prob:
			prob = joint_prob
			pi_of_max_prob2 = pi
	#print >>sys.stderr, "Contamination estimated from Heterozygous SNPs: %f" % 	pi_of_max_prob2		
	
	# all SNPs
	prob = -float("inf")
	pi_of_max_prob3 = 0
	for pi in candidate_PIs:
		joint_prob = 0
		for n,k,t in (snp_hom + snp_het):
			if t == 'Hom':
				p1 = pi	
				p2 = pi/2.0
				joint_prob += max(binom.logpmf(k,n,p1), binom.logpmf(k,n,p2))
			elif t == 'Het':
				p = (1.0 - pi)/2.0
				joint_prob += binom.logpmf(k,n,p)
		if joint_prob > prob:
			prob = joint_prob
			pi_of_max_prob3 = pi
	
	#print >>sys.stderr, "Contamination estimated from All SNPs: %f" % 	pi_of_max_prob3	
	"""
	#print str(infile) + '\t' + str(pi_of_max_prob1)