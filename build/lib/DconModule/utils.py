from time import strftime
import os,sys
import collections
import pysam
import numpy as np
from numpy import mean,std
from scipy import stats
from scipy.stats import binom
import subprocess
import operator
from DconModule import ireader




def get_reference_name(bamfile):
	'''
	get the id and length of all chromosomes from bam file.
	'''
	chroms = {}	#id to length
	samfile = pysam.Samfile(bamfile,'rb')
	a = samfile.header['SQ']
	for i in a:
		chroms[i['SN']] = i['LN']
	return chroms
        
def printlog (mesg_lst):
	'''
	print message to STDERR.
	'''
	if len(mesg_lst)==1:
		msg = "@ " + strftime("%Y-%m-%d %H:%M:%S") + ": " +  mesg_lst[0]
	else:
		msg = "@ " + strftime("%Y-%m-%d %H:%M:%S") + ": " + ' '.join([str(i) for i in mesg_lst])
	print >>sys.stderr,msg

def read_vcf(bam_file, vcf_file):
	'''
	Read VCF file. VCF file could regular or compressed format or remote file
	Compressed files support: .gz, .z, .Z, .bz, .bz2, .gzip2
	Remote files support: http://, https://, ftp://
	'''
	
	chrom_infor = get_reference_name(bam_file)	#id to length
	for line in ireader.reader(vcf_file):
		if line.startswith('#'):
			continue
		f = line.replace('\r','').split()
		
		if f[0] in chrom_infor:
			chrom = f[0]
		else:
			chrom = 'chr' + f[0]
			if chrom not in chrom_infor:
				print >>sys.stderr, "Error: cannot find \"%s\"  (chromosome ID) from the \"%s\" file. Please check your chromosome IDs in \"%s\"." % (f[0], bam_file, bed_file)
				sys.exit(0)
			
		try:
			end = int(f[1])
			start = end - 1
		except ValueError:
			continue
		
		try:
			name = f[2]
		except ValueError:
			name = ''
		
		if end > chrom_infor[chrom]:
			print >>sys.stderr, "%d is larger than the size of %s." % (end, chrom)
			sys.exit(0)
		
		yield [chrom, start, end, name]	

def read_SNPs(infile, min_mvaf = 0.01, min_cvg=30, min_alleleRC=3):
	'''
	remove SNPs with Minor variant allele fraction (or percentage) < 0.01
	'''
	
	ratios = []
	SNPs = []
	
	for l in open(infile,'r'):
		l = l.strip()
		f = l.split()
				
		base1_count = float(f[3])
		base2_count = float(f[5])
		
		#Total coverage
		if base1_count + base2_count < min_cvg:
			continue
		#Coverage over minor allele
		if min(base1_count, base2_count) < min_alleleRC:
			continue
			
		#Minor variant allele fraction (or percentage)
		MVAF = min(base1_count, base2_count)/(base1_count + base2_count)
		if MVAF < min_mvaf:
			continue
		
		# Minor variant allele fraction/Major variant allele fraction ratio. 
		MVAR = min(base1_count, base2_count)/max(base1_count, base2_count)	

		ratios.append(MVAR)
		SNPs.append('\t'.join(f[0:6]))
	return(SNPs, ratios)

#=======
def outlier(point_keys, points, thresh):
	'''
	Detect outlier using median-absolute-deviation (MAD)
	'''
	points = np.array(points)
	if len(points.shape) == 1:
		points = points[:,None]
	imedian = np.median(points, axis=0)
	diff = np.sum((points - imedian)**2, axis=-1)
	diff = np.sqrt(diff)
	med_abs_deviation = np.median(diff)
	
	modified_z_score = 0.6745 * diff / med_abs_deviation
	filtered_keys = []
	for k,z in zip(point_keys, list(modified_z_score)):
		if z > thresh:
			continue
		filtered_keys.append(k)
	return filtered_keys


def filter_outlier(infile, outfile, zcut):
	'''
	filter out heterozygous and outlier SNPs. 
	'''
	hom_snp_wise_pi = {}
	het_snp_wise_pi = {}
	FIN = open(infile,'r')
	for l in FIN:
		l = l.strip()
		f = l.split()
		k = f[0] + '_' + f[1]
				
		if l.startswith('Chrom'):
			continue
		if f[9] == 'Hom':
			hom_snp_wise_pi[k] = float(f[11])
		elif f[9] == 'Het':
			het_snp_wise_pi[k] = float(f[11])
			
	FIN.close()
	
	
	FOUT = open(outfile,'w')
	hom_kept = outlier(hom_snp_wise_pi.keys(), hom_snp_wise_pi.values(), zcut)	
	het_kept = outlier(het_snp_wise_pi.keys(), het_snp_wise_pi.values(), zcut)	
	FIN = open(infile,'r')
	for l in FIN:
		l = l.strip()
		if l.startswith('Chrom'):
			print >>FOUT, l + '\t' + 'Outlier'
			continue
		f = l.split()
		k = f[0] + '_' + f[1]
		
		if f[9] == 'Hom':
			if  k in hom_kept:
				print >>FOUT, l + '\t' + 'Pass'
			else:
				print >>FOUT, l + '\t' + 'Fail'
		elif f[9] == 'Het':
			if  k in het_kept:
				print >>FOUT, l + '\t' + 'Pass'
			else:
				print >>FOUT, l + '\t' + 'Fail'
			
	FIN.close()
	FOUT.close()
	
	#return (len(kept))			
	


#=======	
			
def mle(infile):
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
		if f[12] == 'Fail':continue
		
		allele_1_count = int(f[3])
		allele_2_count = int(f[5])

		if f[9] == 'Hom':
			snp_hom.append([allele_1_count + allele_2_count, allele_2_count, "Hom"]) #n,k
		elif f[9] == 'Het':
			snp_het.append([allele_1_count + allele_2_count, allele_2_count, "Het"]) #n,k
		else:
			continue
			
	print >>sys.stderr, '@ ' + strftime("%Y-%m-%d %H:%M:%S") + ": Estimating contamination from homozygous SNPs ..." 
	prob = -float("inf")
	pi_of_max_prob_hom = 0.0
	for pi in candidate_PIs:
		p2 = pi/2.0
			
		joint_prob = 0
		for n,k,t in snp_hom:			
			pmf_2 = binom.logpmf(k,n,p2)
			joint_prob += pmf_2
							
		if joint_prob > prob:
			prob = joint_prob
			pi_of_max_prob_hom = pi
	return pi_of_max_prob_hom
	
	#print >>sys.stderr, '@ ' + strftime("%Y-%m-%d %H:%M:%S") + ": Estimating contamination from heterozygous SNPs ..." 
	#prob = -float("inf")
	#pi_of_max_prob_het = 0.0
	#for pi in candidate_PIs:
	#	joint_prob = 0
	#	for n,k,t in snp_het:
	#		p = (1.0 - pi)/2.0
	#		joint_prob += binom.logpmf(k,n,p)
	#	if joint_prob > prob:
	#		prob = joint_prob
	#		pi_of_max_prob_het = pi
	#
	#print >>sys.stderr, '@ ' + strftime("%Y-%m-%d %H:%M:%S") + ": Estimating contamination from all SNPs ..." 
	#prob = -float("inf")
	#pi_of_max_prob_both = 0.0
	#for pi in candidate_PIs:
	#	joint_prob = 0
	#	for n,k,t in (snp_hom + snp_het):
	#		if t == 'Hom':
	#			p2 = pi/2.0
	#			pmf_2 = binom.logpmf(k,n,p2)
	#			joint_prob += pmf_2
	#		elif t == 'Het':
	#			p = (1.0 - pi)/2.0
	#			joint_prob += binom.logpmf(k,n,p)
	#	if joint_prob > prob:
	#		prob = joint_prob
	#		pi_of_max_prob_both = pi
	#
	#return (pi_of_max_prob_hom, pi_of_max_prob_het, pi_of_max_prob_both)
			
def gradient_chart(R_outfile, pdf_outfile, arrow_pos):	
	'''
	arrow_pos: indicates the overall contamination
	'''
	ROUT = open(R_outfile,'w')
	
	#print >>ROUT, "library(RColorBrewer)"

	print >>ROUT, "#Overall estimated contamination level is: %s%%" % str(arrow_pos*100)
	print >>ROUT, "cc=c('#a50026','#d73027','#f46d43','#fdae61','#fee08b','#d9ef8b','#a6d96a','#66bd63','#1a9850','#006837')"
	print >>ROUT, "newcol <- colorRampPalette(rev(cc))"
	print >>ROUT, "cols2 <- newcol(501)"
	print >>ROUT, "pdf(\"%s\", width=10, height=4)" % pdf_outfile
	print >>ROUT, "barplot(0:500, col = cols2, border = NA,space=0,yaxt='n',xlab='Contamination percentage',main=\"%s\")" % (pdf_outfile.replace('.pdf',''))

	print >>ROUT, "arrows(%f, -1, %f, 0, xpd = TRUE,col='red')" % (arrow_pos*1000,arrow_pos*1000)
	print >>ROUT, "text(%.2f, -70, labels=c('%.1f%%'),xpd=TRUE)" % (arrow_pos*1000,arrow_pos*100)	
	print >>ROUT, "dev.off()"