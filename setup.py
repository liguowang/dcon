import sys, os, platform, glob
from distutils.core import setup
from setuptools import *

"""
Setup script for Dcon  -- A program to estimate DNA contamination from BAM file.
"""

if sys.version_info[0] != 2 or sys.version_info[1] < 7:
	print >> sys.stderr, "ERROR: Dcon requires Python 2.7"
	sys.exit()
	

def main():
	setup(  name = "Dcon",
			version = "0.1.8",
			py_modules = [ 'psyco_full' ],
			packages = find_packages( 'lib' ),
			package_dir = { '': 'lib' },
			package_data = { '': ['*.ps'] },
			scripts = glob.glob( "bin/*.py"),
			ext_modules = [],
			test_suite = 'nose.collector',
			setup_requires = ['nose>=0.10.4','cython>=0.12'],
			author = "Liguo Wang",
			author_email ="wangliguo78@gmail.com",
			platforms = ['Linux','MacOS'],
			install_requires = ['numpy>=1.13.3','scipy','pysam>=0.13','sklearn'], 
			description = " A program to assess DNA contamination level",
			url = "",
			zip_safe = False,
			dependency_links = [],
			license='MIT',
			python_requires='>=2.6, !=3.0.*, !=3.1.*, !=3.2.*, <4',
			classifiers=[
				'Development Status :: 5 - Production/Stable',
				'Environment :: Console',
				'Intended Audience :: Science/Research',
				'License :: OSI Approved :: GNU General Public License v2 (GPLv2)',
				'Operating System :: MacOS :: MacOS X',
				'Operating System :: POSIX',
				'Programming Language :: Python',
				'Topic :: Scientific/Engineering :: Bio-Informatics',
			],
			
			keywords='DNA contamination, Capture gene panel sequencing, Whole exome sequencing, whole genome sequencing',
             )


if __name__ == "__main__":
	main()
