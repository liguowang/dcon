from time import strftime
import os,sys
import collections
import pysam
import ireader
import numpy as np
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
def outlier(point_keys, points, thresh):
	'''
	Detect outlier using median-absolute-deviation (MAD)
	'''
	points = np.array(points)
	if len(points.shape) == 1:
		points = points[:,None]
	median = np.median(points, axis=0)
	diff = np.sum((points - median)**2, axis=-1)
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
	filter out outlier SNPs 
	'''
	hom_pi = {}
	het_pi = {}
	FIN = open(infile,'r')
	FOUT = open(outfile,'w')
	for l in FIN:
		l = l.strip()
		f = l.split()
		k = f[0] + '_' + f[1]
				
		if l.startswith('Chrom'):
			print >>FOUT, l + '\t' + 'Outlier'
			continue
		else:	
			if f[9] == 'Hom':
				hom_pi[k] = float(f[11]) 
			elif f[9] == 'Het':
				het_pi[k] = float(f[11]) 
	FIN.close()
	
	hom_filtered_keys = outlier(hom_pi.keys(), hom_pi.values(), zcut)
	het_filtered_keys = outlier(het_pi.keys(), het_pi.values(), zcut)
	
	FIN = open(infile,'r')
	for l in FIN:
		l = l.strip()
		f = l.split()
		k = f[0] + '_' + f[1]
		if f[9] == 'Hom':
			if  k in hom_filtered_keys:
				print >>FOUT, l + '\t' + 'Pass'
			else:
				print >>FOUT, l + '\t' + 'Failed'
				
		if f[9] == 'Het':
			if  k in het_filtered_keys:
				print >>FOUT, l + '\t' + 'Pass'
			else:
				print >>FOUT, l + '\t' + 'Failed'
		else:
			pass
	FIN.close()
	FOUT.close()
	
	return (len(hom_filtered_keys), len(het_filtered_keys))			
	
	
def read_SNPs(infile, cutoff = 0.01):
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
		
		#Minor variant allele fraction (or percentage)
		MVAF = min(base1_count, base2_count)/(base1_count + base2_count)
		
		if MVAF < cutoff:
			continue
		# Minor variant allele fraction/Major variant allele fraction ratio. 
		MVAR = min(base1_count, base2_count)/max(base1_count, base2_count)	

		ratios.append(MVAR)
		SNPs.append('\t'.join(f[0:6]))
	return(SNPs, ratios)

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
		allele_1_count = int(f[3])
		allele_2_count = int(f[5])
		if f[12] == 'Failed':
			continue
		if allele_1_count < allele_2_count:
			(allele_1_count, allele_2_count) = (allele_2_count,allele_1_count)
		
		if allele_2_count*1.0/(allele_2_count + allele_1_count) < 0.01:
			continue
		if f[9] == 'Het':
			snp_het.append([allele_1_count + allele_2_count, allele_2_count, 'Het']) #n,k
		elif f[9] == 'Hom':
			snp_hom.append([allele_1_count + allele_2_count, allele_2_count, 'Hom']) #n,k
		else:
			continue
	
	#key is pi, value is the product of pmfs across all snps
	
	# number of SNPs used for MLE
	size_cut = 10
	
	if len(snp_hom) >= size_cut:
		prob = -float("inf")
		pi_of_max_prob = 0.0
		#tmp = {}
		for pi in candidate_PIs:
			joint_prob = 0
			# Homozygous SNP contaminated by homozygous SNP
			#p1 = pi	
			# Homozygous SNP contaminated by heterozygous SNP
			p2 = pi/2.0
		
			for n,k,t in snp_hom:
				#pmf_1 = binom.logpmf(k,n,p1)
				pmf_2 = binom.logpmf(k,n,p2)
				#joint_prob += max(pmf_1, pmf_2)
				joint_prob += pmf_2
			#tmp[pi] = joint_prob
			if joint_prob > prob:
				prob = joint_prob
				pi_of_max_prob = pi
		#for i,j in sorted(tmp.items(), key=lambda x:x[1]):
		#	print i,j
		return pi_of_max_prob
	else:
		prob = -float("inf")
		pi_of_max_prob = 0.0
		for pi in candidate_PIs:
			joint_prob = 0
			for n,k,t in (snp_hom + snp_het):
				if t == 'Hom':
					#p1 = pi	
					#pmf_1 = binom.logpmf(k,n,p1)
					
					p2 = pi/2.0
					pmf_2 = binom.logpmf(k,n,p2)
					joint_prob += pmf_2
				elif t == 'Het':
					p = (1.0 - pi)/2.0
					joint_prob += binom.logpmf(k,n,p)
			if joint_prob > prob:
				prob = joint_prob
				pi_of_max_prob = pi
	
		return pi_of_max_prob

def estimate_overall_pi(infile, R_outfile):
	'''
	Estimate the overall contamination percentage using mode.
	'''
	
	all = []
	homo = []
	hete = []
	overall_pi = 0.0
	with open(infile,'r') as L:
		for line in L:
			line = line.strip()
			if line.startswith('Chrom'):continue
			f = line.split()
			label = f[9]
			distance = float(f[10])
			pi = float(f[11])
			if pi < 0.01:
				continue
				
			filter = f[12]
			if filter == 'Failed':
				continue
				
			all.append((distance, pi))
			if label == 'Het':
				hete.append((distance, pi))
			elif label == 'Hom':
				homo.append((distance, pi))
			
	all_sorted = [i[1] for i in sorted(all,key=lambda x:x[0])]
	homo_sorted = [i[1] for i in sorted(homo,key=lambda x:x[0])]
	hete_sorted = [i[1] for i in sorted(hete,key=lambda x:x[0])]
	
	
	if len(all_sorted) == 0:
		 return overall_pi
	elif len(all_sorted) <= 2:
		overall_pi = mean(all_sorted)
		return overall_pi
	else:
	
		ROUT = open(R_outfile,'w')
		print >>ROUT, "ALL_SNP_PI <- c(%s)" % ','.join([str(i) for i in all_sorted])
		print >>ROUT, "ALL_D <- density(ALL_SNP_PI,na.rm=T)"
		print >>ROUT, "m1 <- ALL_D$x[which.max(ALL_D$y)]"
		print >>ROUT, "cat( round(m1, 4) )"		
		ROUT.close()
		
		command = ['Rscript', R_outfile, 'shell=True']
		try:
			overall_pi = subprocess.check_output(command, universal_newlines=True)
		except:
			print >>sys.stderr, 'Cannot run ' + '"' + R_outfile + '". ' 
			pass
		return float(overall_pi)
			
def gradient_chart(R_outfile, pdf_outfile, arrow_pos):	
	'''
	arrow_pos: indicates the overall contamination
	'''
	ROUT = open(R_outfile,'w')
	
	#print >>ROUT, "library(RColorBrewer)"

	print >>ROUT, "cc=c('#a50026','#d73027','#f46d43','#fdae61','#fee08b','#d9ef8b','#a6d96a','#66bd63','#1a9850','#006837')"
	print >>ROUT, "newcol <- colorRampPalette(rev(cc))"
	print >>ROUT, "cols2 <- newcol(501)"
	print >>ROUT, "pdf(\"%s\", width=10, height=4)" % pdf_outfile
	print >>ROUT, "barplot(0:500, col = cols2, border = NA,space=0,yaxt='n',xlab='Contamination percentage',main=\"%s\")" % (pdf_outfile.replace('.pdf',''))

	print >>ROUT, "arrows(%f, -1, %f, 0, xpd = TRUE,col='red')" % (arrow_pos*1000,arrow_pos*1000)
	print >>ROUT, "text(%.2f, -70, labels=c('%.1f%%'),xpd=TRUE)" % (arrow_pos*1000,arrow_pos*100)	
	print >>ROUT, "dev.off()"