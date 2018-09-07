#!/usr/bin/env python

#import built-in modules
import os,sys
import re
import string
from optparse import OptionParser
import warnings
import string
import collections
import math

#import third-party modules

#changes to the paths

#changing history to this module
#05/26/2011: suppport multiple spliced mapped reads

__author__ = "Liguo Wang"
__copyright__ = "Copyright 2010, Wei Li's Lab"
__credits__ = []
__license__ = "GPL"
__version__ = "1.1.0"
__maintainer__ = "Liguo Wang"
__email__ = "liguow@bcm.edu"
__status__ = "Development" #Prototype or Production


def RSS(arg):
	'''calculate Square root of sum of square. Input is ',' separated numbers'''
	lst=arg.split(',')
	lst_sum=0
	for i in [ int(i)**2 for i in lst]:
		lst_sum += i 
	#nsr=10*math.log10((1+noi_sum**0.5)/(1+sig_sum**0.5))
	return lst_sum**0.5
	
def H_mean(arg):
	'''calculate harmornic mean. Input is ',' separated numbers'''
	lst=[1/float(i) for i in arg.split(',') if float(i) !=0]
	if len(lst) == 0:
		return "NA"
	else:
		return len(lst)/(sum(lst))

def shannon_entropy(arg):
	'''calculate shannon's entropy (or Shannon-Wiener index). Input is list of numbers'''
	lst=[float(i) for i in arg if float(i)>0]
	if sum(lst)<=0 or min(lst)<0:return "NA"
	entropy=0.0
	for i in lst:
		entropy += (i/sum(lst)) * math.log((i/sum(lst)))
	return -entropy

	
def shannon_entropy_es(arg):
	'''calculate estimator of shannon's entropy (Chao & Shen, 2003)'''
	
	lst=[float(i) for i in arg if float(i)>0]
	if sum(lst)<=0 or min(lst)<0:return "NA"	#if there is no fragmental splicing
	if (len(lst)==1): return 0					#if there is only 1 fragmental splicing
	lst.append(2)
	
	#estimate C_bar
	singleton=0
	entropy=0.0
	for i in lst:
		if i ==1:singleton +=1
	
	C_bar = 1- (singleton/sum(lst))
	for i in lst:entropy += ( (C_bar*i/sum(lst)) * math.log((C_bar*i/sum(lst))) )/(1-(1-C_bar*i/sum(lst))**sum(lst))
	return -entropy

def shannon_entropy_ht(arg):
	'''calculate estimator of shannon's entropy based on Horzitz-Thompson'''
	lst=arg.split(',')
	lst=[float(i) for i in lst if float(i)>0]
	if sum(lst)<=0 or min(lst)<0:return "NA"	#if there is no fragmental splicing
	if (len(lst)==1): return 0					#if there is only 1 fragmental splicing
	
	#estimate C_bar
	entropy=0.0
	for i in lst:
		entropy += ( (i/sum(lst)) * math.log((i/sum(lst))) )/(1-(1-i/sum(lst))**sum(lst))
	return -entropy
	
def simpson_index(arg):
	'''calculate Gini-Simpson's index. Input is ',' separated numbers'''
	lst=arg.split(',')
	lst=[float(i) for i in lst if float(i)>0]
	simpson=0.0
	
	try:
		for i in lst:
			simpson = simpson + (i/sum(lst))**2
		return 1-simpson
	except: return 0
	
def simpson_index_es(arg):
	'''calculate estimator Gini-Simpson's index. Input is ',' separated numbers'''
	lst=arg.split(',')
	lst=[float(i) for i in lst if float(i)>0]
	simpson=0.0
	
	try:
		for i in lst:
			simpson = simpson + i*(i-1)
		return 1- (simpson/(sum(lst)*(sum(lst)-1)))
	except: return 0
	
def Hill_number(arg,qvalue=1):
	'''Calculate real diversity (Hill's number). Input is ',' separated numbers. qvalue is the only
	parameter for Hill's function. When q=1, it return exp(H) which is the effective number of junctions
	calculated by Shannon's entropy. When q<1, Hill's function was favors low frequency junctions. 
	When q>1, Hill's function was favors high frequency junctions (common junctions). Simpon's Index
	is particular case of Hill's function as q=2'''
	
	lst=arg.split(',')
	lst=[float(i) for i in lst if float(i)>0]
	freq=[(i/sum(lst))**qvalue for i in lst]
	try:
		return (sum(freq))**(1/(1-qvalue))
	except:
		return math.exp(shannon_entropy(arg))
def fisher_test(a,alt='two-sided'):
	'''
	Fisher Exact Test
	alternative : {'two-sided', 'less', 'greater'}, optional
	'''
	if len(a) !=4:
		return None
	p = stats.fisher_exact([[ int(a[0]),int(a[1]) ],[int(a[2]),int(a[3]) ]],alternative=alt)
	return p[1]

def bh_adjust(pvalues, correction_type = "Benjamini-Hochberg"):                
    """                                                                                                   
    consistent with R - print correct_pvalues_for_multiple_testing
    ([0.0, 0.01, 0.029, 0.03, 0.031, 0.05, 0.069, 0.07, 0.071, 0.09, 0.1]) 
    """
    from numpy import array, empty                                                                        
    pvalues = array(pvalues) 
    n = float(pvalues.shape[0])                                                                           
    new_pvalues = empty(n)
    if correction_type == "Bonferroni":                                                                   
        new_pvalues = n * pvalues
    elif correction_type == "Bonferroni-Holm":                                                            
        values = [ (pvalue, i) for i, pvalue in enumerate(pvalues) ]                                      
        values.sort()
        for rank, vals in enumerate(values):                                                              
            pvalue, i = vals
            new_pvalues[i] = (n-rank) * pvalue                                                            
    elif correction_type == "Benjamini-Hochberg":                                                         
        values = [ (pvalue, i) for i, pvalue in enumerate(pvalues) ]                                      
        values.sort()
        values.reverse()                                                                                  
        new_values = []
        for i, vals in enumerate(values):                                                                 
            rank = n - i
            pvalue, index = vals                                                                          
            new_values.append((n/rank) * pvalue)                                                          
        for i in xrange(0, int(n)-1):  
            if new_values[i] < new_values[i+1]:                                                           
                new_values[i+1] = new_values[i]                                                           
        for i, vals in enumerate(values):
            pvalue, index = vals
            new_pvalues[index] = new_values[i]                                                                                                                  
    return new_pvalues
