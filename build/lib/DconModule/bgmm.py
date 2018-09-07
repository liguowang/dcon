import sys
from time import strftime
import numpy as np
from sklearn import mixture

def build_GMM(outfile, d, names, n, iter, tol, rnd, prob_cut):
	"""
	Return means of components of Gaussian Mixture Model.
	d is a list of ratios (minor variant allele fraction / major variant allele fraction). 
	names is a list of SNP IDs
	n is the number of components
	rnd is a random number. You get exactly the same results when running multiple times using the same random number. Must be integer. 
	"""
	
	print >>sys.stderr, '@ ' + strftime("%Y-%m-%d %H:%M:%S") + ": Building Bayesian Gaussian Mixture model (BGMM) ..."
	bgmm = mixture.BayesianGaussianMixture(n_components = n, covariance_type = 'full', max_iter = iter, tol = tol, random_state=rnd)
	bgmm_model = bgmm.fit(np.array(d).reshape(-1,1))
	
	print >>sys.stderr, "\tConverge status: " + str(bgmm_model.converged_)
	print >>sys.stderr, "\tIterations: " + str(bgmm_model.n_iter_)
	
	print >>sys.stderr,'\n'
	print >>sys.stderr, '@ ' + strftime("%Y-%m-%d %H:%M:%S") + ": Summerzie BGMM model ..." 
	
	component_means = bgmm_model.means_[:,0]
	#print component_means
	hom_mean = sorted(component_means)[0]
	het_mean = sorted(component_means)[1]
	print >>sys.stderr, "\tMeans of homozygous component: %f" % hom_mean
	print >>sys.stderr, "\tMeans of heterozygous component: %f" % het_mean
	
	component_weights = sorted(bgmm_model.weights_)
	#print component_weights
	print >>sys.stderr, "\tWeight of homozygous component: %f" % component_weights[0]
	print >>sys.stderr, "\tWeight of heterozygous component: %f" %component_weights[1]
	
	
	print >>sys.stderr, '\n'
	print >>sys.stderr, '@ ' + strftime("%Y-%m-%d %H:%M:%S") + ": Classify variants ...\n"
	labels = {}
	for idx,val in enumerate(component_means):
		if val == max(component_means):
			labels[idx] = "Het"
		elif val == min(component_means):
			labels[idx] = "Hom"
			
	probs = bgmm_model.predict_proba(np.array(d).reshape(-1,1))
	
	
	print >>sys.stderr, '@ ' + strftime("%Y-%m-%d %H:%M:%S") + ": Writing to %s ..." % outfile
	OUT_FILE = open(outfile,'w')
	header = ['Chrom', 'Ref_pos', 'Allele_1', 'Allele_1_count', 'Allele_2', 'Allele_2_count', 'Ratio', 'Prob_of_' + labels[0], 'Prob_of_' + labels[1], 'Label', 'Distance', 'Contamination_pi']
	print >>OUT_FILE, '\t'.join(header)
	
	for nm,ratio, p in zip(names, d, probs):
		f = nm.split('\t')
		count1 = float(f[3])	
		count2 = float(f[5])
		minor = min(count1, count2)	#minor variant allele count
		count_total = count1 + count2
		p_list = list(p)	#list of two probabilities
		if max(p_list) < prob_cut:
			continue
		index_of_max_p = p_list.index(max(p_list))	
		lab = labels[index_of_max_p]	#Assign variant to the group with maximum probability
		if lab == 'Het':			
			distance = abs(ratio - het_mean)
			contamination = 1.0 - (minor*2.0)/count_total
		elif lab == 'Hom':		
			distance = abs(ratio - hom_mean)
			contamination = (minor*2.0)/count_total
			#if contamination > 0.5:
				#contamination = contamination/2.0
		
		#remove contamination estimation (pi) less than 0.01
		#when contamination  (pi) is larger than 0.5. Use 1 - pi instead. 
		#if contamination < 0.01:
		#	continue
		if contamination > 0.5:
			contamination = 1.0 - contamination
			
		print >>OUT_FILE, nm + '\t' + str(ratio) + '\t' + '\t'.join([str(i) for i in p_list]) + '\t' + lab + '\t' + str(distance) + '\t' + str(contamination)
		#print >>OUT_FILE, nm + '\t' + str(ratio) + '\t' + '\t'.join([str(i) for i in p_list]) + '\t' + lab + '\t' + str(contamination)
		
	OUT_FILE.close()
	
