from time import strftime
import os,sys
import string
import collections
import pysam
import ireader
from numpy import mean,std
from scipy import stats
import subprocess
import pyBigWig



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
		# Minor variant allele fraction/Major variant allele fraction ratio. 
		MVAR = min(base1_count, base2_count)/max(base1_count, base2_count)	
		if MVAF < 0.01:
			continue
		ratios.append(MVAR)
		SNPs.append('\t'.join(f[0:6]))
	return(SNPs, ratios)

#=======	
def estimate_overall_pi(infile, outfile, size_cut):
	'''
	estimate overall contaminant percentage
	size_cut:
		* when the number of usable SNP is less then size_cut. Use their mean to estimate
		* when the number of usable SNP is larger then size_cut. Use mode 
	'''
	all = []
	homo = []
	hete = []
	container = collections.defaultdict(list)
	
	dist_cut = 30
	
	with open(infile,'r') as L:
		for line in L:
			line = line.strip()
			if line.startswith('Chrom'):continue
			f = line.split()
			coord = f[0] + ':' + f[1]
			
			chrom = f[0]
			end = int(f[1])
			start = end -1 
			
			if chrom not in container:
				container[chrom].append(start)
			else:
				if min([abs(i - start) for i in container[chrom]]) > dist_cut: #at least 30 nt away
					container[chrom].append(start)
				else:
					continue
			
			
			label = f[9]
			distance = float(f[10])
			pi = float(f[11])
			if pi < 0.01:
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
		 print  "\n Cannot find SNPs to estimated contamination level" 
		 sys.exit(0)
	
	# when the number of usable SNP is less than 3. Use their mean to estimate	
	elif (len(all_sorted) > 0) and (len(all_sorted) < 3):
		pi_mean = mean(all_sorted)
		print  "\nEstimated contamination level (All SNP):\t%.3f" % (pi_mean)
		sys.exit(0)
		
	# when the number of usable SNP is larger than 3. Use mode 
	elif len(all_sorted) >= 3:
	
		ROUT = open(outfile + '.dcon.R','w')
		
		print >>ROUT, "pdf(\"%s\",width=6,height=8)" % (outfile + '.dcon.pdf')
		print >>ROUT, 'par(mfrow=c(2,2))'
	
		print >>ROUT, "ALL_SNP_PI <- c(%s)" % ','.join([str(i*100.0) for i in all_sorted])
		print >>ROUT, "ALL_D <- density(ALL_SNP_PI[1:%d],na.rm=T)" % size_cut
		print >>ROUT, "m1 <- ALL_D$x[which.max(ALL_D$y)]"
		print >>ROUT, "plot(ALL_D,xlab='Contamination percentage',col='blue',main=paste(\"ALL SNP, Estimated = \", round(m1,2), \"%\"))"
		print >>ROUT, "abline(v=m1,col='red',lty='dashed')"
		
		if len(homo_sorted) >= 3:
			print >>ROUT, "\n"
			print >>ROUT, "HOMO_SNP_PI <- c(%s)" % ','.join([str(i*100.0) for i in homo_sorted])
			print >>ROUT, "HOMO_D <- density(HOMO_SNP_PI[1:%d],na.rm=T)" % size_cut
			print >>ROUT, "m2 <- HOMO_D$x[which.max(HOMO_D$y)]"
			print >>ROUT, "plot(HOMO_D,xlab='Contamination percentage',col='blue',main=paste(\"HOMO SNP, Estimated = \", round(m2,2), \"%\"))"
			print >>ROUT, "abline(v=m2,col='red',lty='dashed')"
		
		if len(hete_sorted) >= 3:
			print >>ROUT, "\n"
			print >>ROUT, "HETE_SNP_PI <- c(%s)" % ','.join([str(i*100.0) for i in hete_sorted])
			print >>ROUT, "HETE_D <- density(HETE_SNP_PI[1:%d],na.rm=T)" % size_cut
			print >>ROUT, "m3 <- HETE_D$x[which.max(HETE_D$y)]"
			print >>ROUT, "plot(HETE_D,xlab='Contamination percentage',col='blue',main=paste(\"HETE SNP, Estimated = \", round(m3,2), \"%\"))"
			print >>ROUT, "abline(v=m3,col='red',lty='dashed')"
		
		
		print >>ROUT, 'dev.off()'
		
		if len(all_sorted) >= 3:
			print >>ROUT, "cat('Estimated contamination level (All SNPs): ', round(m1/100.0, 4), '\\n', sep=' ')"
			#print >>ROUT, 'print ("\n")'
		elif len(all_sorted) >0:
			print >>ROUT, "cat('Estimated contamination level (All SNPs): ', round(%s,4), '\\n', sep=' ')" % mean(all_sorted)
			#print >>ROUT, 'print ("\n")'
		else:
			print >>ROUT, "cat('Estimated contamination level (ALL SNPs): unknown', '\\n', sep=' ')"
			#print >>ROUT, 'print ("\n")'
		
		if len(homo_sorted) >= 3:
			print >>ROUT, "cat('Estimated contamination level (HOMO SNPs only): ', round(m2/100.0, 4),'\\n', sep=' ')"
		elif len(homo_sorted) > 0:
			print >>ROUT, "cat('Estimated contamination level (HOMO SNPs): ', round(%s,4), '\\n', sep=' ')" % mean(homo_sorted)
		else:
			print >>ROUT, "cat('Estimated contamination level (HOMO SNPs): unknown', '\\n', sep=' ')"
		
		if len(hete_sorted) >= 3:
			print >>ROUT, "cat('Estimated contamination level (HETE SNPs only): ', round(m3/100.0, 4),'\\n', sep=' ')"
		elif len(hete_sorted) > 0 :
			print >>ROUT, "cat('Estimated contamination level (HETE SNPs only): ', round(%s,4), '\\n', sep=' ')" % mean(hete_sorted)
		else:
			print >>ROUT, "cat('Estimated contamination level (HETE SNPs): unknown', '\\n', sep=' ')"
		ROUT.close()
		
		
		try:
			subprocess.call('Rscript ' + outfile + '.dcon.R', shell=True)
		except:
			print >>sys.stderr, 'Cannot generate pdf file from ' + '"' + outfile + '.dcon.R' + '". ' 
			pass
			
				
