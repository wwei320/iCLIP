#===============================================================================
#         FILE: iclip_seq step1
#
#  DESCRIPTION: Fully automatic pipline downloading iclip-seq data SRA files and converting to fastq files
#
#      OPTIONS:  --help;                                       show the help message and exit
#                --project(required) <project_name>; define the project name;
#                --rootdir(required) <path_to_your_root_dir>; define the rootdir;  
#                --sradir(required) <path_to_your_sra_dir>; define the sra dir;
#                --refdir(required) <path_to_your_ref_dir>; define the reference dir;#                
#                --threads(required) <# of threads>; define # of threads;
#                --genome(Optional) <mm9(default)/hg19/rn5>; define genome version;
#
# REQUIREMENTS: Python_Modules as 'import section';"STAR" as RNAseq mapper; SRA toolkit; RSeQ;
#
#       AUTHOR: Wei Wang  wwei320@gmail.com
#       
#      VERSION: 2.0
#      CREATED: 2019-June-05
#===============================================================================
####import section
import matplotlib; matplotlib.use('agg')
import itertools
import numpy as np
import pandas as pd
from pylab import *
import os
import sys
import ntpath
import re
import glob
import scipy
import collections
import subprocess
import time
# import multiprocessing
import argparse
import unittest

def strandinfo_SINGLE(logfile,outfile):
	f=open( logfile ).readlines()
	F = 1
	R = 1
	for line in f:
		if "++,--" in line:
			inputline=line.split( ": " )[1]
			F=float(inputline.replace('\n',''))
		elif "+-,-+" in line:
			inputline=line.split( ": " )[1]
			R=float(inputline.replace('\n',''))
	with open(outfile, 'w') as f:	
		if F/R>=3:
			f.write('forward\n')
		elif R/F>=3:
			f.write('invert\n')
		else:
			f.write('NONE\n')
		f.close()	
	os.system('rm '+logfile)

def strandinfo_PAIR(logfile):
	f=open( logfile ).readlines()
	STRINFO = ''
	F = 1
	R = 1
	for line in f:
		if "1++,1--,2+-,2-+" in line:
			inputline=line.split( ": " )[1]
			F=float(inputline.replace('\n',''))
		elif "1+-,1-+,2++,2--" in line:
			inputline=line.split( ": " )[1]
			R=float(inputline.replace('\n',''))	
		if F/R>=3:
			STRINFO='FR' 
		elif R/F>=3:
			STRINFO='RF'		
		else:
			STRINFO='P_NONE'
	return STRINFO
	
def handel_PAIR(strKEY,fastqfile1,fastqfile2,bamfile):	
	if strKEY=='RF':		
		fastqfile=fastqfile1		
	else:
		fastqfile=fastqfile2
	return fastqfile

def downloadcheck(fqmdir,sample,seqtype,tail1,tail2):
	if seqtype=='SINGLE':
		exists = os.path.isfile(fqmdir+sample+'.fastq')
		if exists:
			print sample+"Downloading finished"
			DNCHECK='YES'
		else:
			print sample+"Downloading error, try re-download once"
			DNCHECK='NO'

	elif seqtype=='PAIRED':
		exists1 = os.path.isfile(fqmdir+sample+tail1+'.fastq')
		exists2 = os.path.isfile(fqmdir+sample+tail2+'.fastq')
		if exists1 and exists2:
			print sample+"Downloading finished"
			DNCHECK='YES'
		else:
			print sample+"Downloading error, try to re-download once"
			DNCHECK='NO'
	return 	DNCHECK
		

##arg setting parts###
parser = argparse.ArgumentParser(description="Use this to run RNAseq-APA automatic pipeline.")
parser.add_argument("--project", action="store", dest='project',default='', metavar='<name_of_project>', help="define the project name")
parser.add_argument("--rootdir", action="store", dest='rootdir',default='', metavar='<rootdir>', help="define the rootdir")
parser.add_argument("--sradir", action="store", dest='sradir',default='', metavar='<sradir>', help="define the sradir")
parser.add_argument("--refdir", action="store", dest='refdir',default='', metavar='<refdir>', help="define the refdir")
parser.add_argument("--genodir", action="store", dest='genodir',default='', metavar='<genodir>', help="define the genome dir for mapping")
parser.add_argument("--genome", action="store", dest='genome',default='mm9', metavar='<mm9(default)/hg19/rn4/...>', help="define the genome version")
parser.add_argument("--threads", action="store", dest='threads',default='', metavar='<threads>', help="define the number of threads")

args=parser.parse_args()

#####python setting #####
rootdir=args.rootdir
project=args.project
sradir=args.sradir
refdir=args.refdir
geno=args.genome
genoDir=args.genodir
CPUS=args.threads

CPUS=str(CPUS)

####automatic setting section####
scrdir=os.path.dirname(os.path.abspath(__file__))
mapper='_star'
tail1='_1'
tail2='_2'
UTRdir=refdir
fqmdir=os.path.join(rootdir, project+'/rawfastq/')
# genoDir=os.path.join(rootdir, 'data/ucsc/genomes/'+geno+mapper)
##toggled by ww: samoutdir=os.path.join(rootdir, project+'/rawsam/')
samoutdir=os.path.join(rootdir, project+'/rawsam')
rawoutdir=os.path.join(rootdir, project+'/rawout/')
testdir=os.path.join(rootdir, project+'/test/')
plotdir=os.path.join(rootdir, project+'/plot/')
reportdir=os.path.join(rootdir, project+'/report/')
samplefile=rootdir+project+'/sample_list.txt'

if not os.path.exists(fqmdir):
    os.makedirs(fqmdir)	
if not os.path.exists(samoutdir):
    os.makedirs(samoutdir)	
if not os.path.exists(rawoutdir):
    os.makedirs(rawoutdir)	
if not os.path.exists(testdir):
    os.makedirs(testdir)	
if not os.path.exists(plotdir):
    os.makedirs(plotdir)	
if not os.path.exists(reportdir):
    os.makedirs(reportdir)	
	
dfsample = pd.read_table(samplefile)
samples=dfsample.Run.unique().tolist()


# os.chdir(dndir)

#for sample in samples:
for sample in samples:
###pre-step1 Download sra and convert to Fastq files####
	print "\n--------------Pre-step1: Download sra and convert to Fastq files, Sample: "+ sample +"----------------"	
	srafile=sample+'.sra'
	cmd1='prefetch -q '+ sample
	os.system(cmd1)
	
	cmd2='mv '+sradir+srafile+' ' +fqmdir
	os.system(cmd2)
	
	seqtype=dfsample[dfsample.Run==sample].LibraryLayout.values[0]
	samplename=dfsample[dfsample.Run==sample].samplename.values[0]
	
	if seqtype=='SINGLE':	
		cmd3='fastq-dump --outdir ' + fqmdir +' '+ fqmdir + srafile
	elif seqtype=='PAIRED':
		cmd3='fastq-dump --outdir ' +fqmdir +' --split-files --split-3 ' + fqmdir+srafile
	os.system(cmd3)
	
	cmd4 = 'rm '+ fqmdir+srafile
	os.system(cmd4)	

	DNCHECK=downloadcheck(fqmdir,sample,seqtype,tail1,tail2)
	if DNCHECK=='NO':
		os.system(cmd1)
		os.system(cmd2)
		os.system(cmd3)
		os.system(cmd4)
		
