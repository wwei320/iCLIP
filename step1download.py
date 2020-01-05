#!/bin/bash
#===============================================================================
#         FILE: all_steps.RNAseq_APA.pip.py
#
#        USAGE: python all_steps.RNAseq_APA.pip.py --project --<OPTIONS> 
#
#  DESCRIPTION: Fully automatic pipline profiling APA,IPA and gene expression using RNAseq data from GEO
#
#      OPTIONS:  --help;                                       show the help message and exit
#                --project(required) <project_name>; define the project name;
#                --rootdir(required) <path_to_your_root_dir>; define the rootdir;  
#                --sradir(required) <path_to_your_sra_dir>; define the sra dir;
#                --refdir(required) <path_to_your_ref_dir>; define the reference dir;#                
#                --threads(required) <# of threads>; define # of threads;
#                --genome(Optional) <mm9(default)/hg19/rn5>; define genome version;
#
# REQUIREMENTS: Pythod_Modules as 'import section'; "APAanalizer" in bioconductor; "STAR" as RNAseq mapper; SRA toolkit; RSeQ;
#
#       AUTHOR: Ruijia Wang  rjwang.bioinfo@gmail.com
#      ADVISOR: Bin Tian     btian@njms.rutgers.edu
# ORGANIZATION: 
#      VERSION: 2.0
#      CREATED: 2018-June-05
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
		
#### step1 Mapping ####	
# 	bamfile=samoutdir+'/'+sample+'.Aligned.sortedByCoord.out.bam'
# 	print "\n--------------Step1: Mapping to the Genome("+seqtype+"),Sample:"+sample+"----------------"	
	
# 	if seqtype=='SINGLE':
# 		fastqfile=fqmdir+sample+'.fastq'
# 		cmd5="STAR --runThreadN "+CPUS+" --genomeDir "+genoDir+" --outSAMtype BAM SortedByCoordinate --readFilesIn "+fastqfile+ " --outFileNamePrefix "+ samoutdir+'/' +sample+'.'
# 	elif seqtype=='PAIRED':
# 		fastqfile1=fqmdir+sample+tail1+'.fastq'	
# 		fastqfile2=fqmdir+sample+tail2+'.fastq'	
# 		# fastqfile=fqmdir+sample+tail2+'.fastq'		
# 		cmd5="STAR --runThreadN "+CPUS+" --genomeDir "+genoDir+" --outSAMtype BAM SortedByCoordinate --readFilesIn "+fastqfile1 +" "+fastqfile2+ " --outFileNamePrefix "+ samoutdir +'/'+sample+'.'
		
# 	# cmd5="STAR --runThreadN "+CPUS+" --genomeDir "+genoDir+" --outSAMtype BAM SortedByCoordinate --readFilesIn "+fastqfile+ " --outFileNamePrefix "+ samoutdir +'/'+sample+'.'
# 	print cmd5
# 	os.system(cmd5)

# #### step2 checking strand ####	
# 	print "\n--------------Step2: Determing strand info("+seqtype+"),Sample:"+sample+"----------------"	
# 	logfile=testdir+'check.txt'
# 	outfile=testdir+'strand.txt'
# 	cmd6 = 'infer_experiment.py -r '+refdir+'/'+geno+'.RefSeq.bed -i '+bamfile+' > '+logfile
# 	print cmd6
# 	os.system(cmd6)		
# 	if seqtype=='SINGLE': 
# 		strandinfo_SINGLE(logfile,outfile)
# 	elif seqtype=='PAIRED':
# 		STRINFO=strandinfo_PAIR(logfile)
# 		# os.system('rm '+bamfile)
# 		# os.system('rm '+logfile)
# 		# os.system('rm '+outfile)
# 		fastqfile=handel_PAIR(STRINFO,fastqfile1,fastqfile2,bamfile)
# 		cmd5="STAR --runThreadN "+CPUS+" --genomeDir "+genoDir+" --outSAMtype BAM SortedByCoordinate --readFilesIn "+fastqfile+ " --outFileNamePrefix "+ samoutdir +'/'+sample+'.'
# 		os.system(cmd5)
# 		os.system(cmd6)
# 		strandinfo_SINGLE(logfile,outfile)		
# 	else:
# 		print " Library type ERROR, check your sample table!!! "
	
# ####step3 Count reads using R####
# 	newbamfile=samoutdir+'/'+samplename+'.bam'
# 	os.system('mv '+bamfile +' '+ newbamfile)
	
# 	print "\n--------------Step3: Counting reads("+seqtype+"),Sample:"+sample+"----------------"		
# 	prjdir=os.path.join(rootdir, project+'/')
# 	arg1=prjdir
# 	arg2=geno
# 	arg3=samoutdir+'/'
# 	arg4=UTRdir
# 	arg5=CPUS
# 	###Rscriptfile0=scrdir+'/RNAseq.APAprofile.R'
# 	Rscriptfile0=scrdir+'RNAseq.APAprofile.R'

# 	####  RUN Rscript ###
# 	subprocess.call(['Rscript', Rscriptfile0, arg1, arg2, arg3,arg4,arg5])	
# 	os.system('rm '+newbamfile)
# 	if seqtype=='SINGLE':
# 		os.system('rm '+fqmdir+sample+'.*')
# 	elif seqtype=='PAIRED':
# 		os.system('rm '+fqmdir+sample+'_*')
		
# ####step4 Combine all sample together ####
# print "\n--------------Step4:  Combine all sample together ----------------"
# statdir=rawoutdir
# outputdir=reportdir

# statfiles=sorted(glob.glob(os.path.join(statdir, '*.CDS_reads_TPM.conserved.txt')))
# for statfile in statfiles:
# 	df1 = pd.read_table(statfile)
# 	if statfile==statfiles[0]:
# 		dfcds=df1.copy()
# 	else:
# 		dfcds=pd.merge(dfcds,df1,on=['gene_symbol'],how='left')
		
# statfiles=sorted(glob.glob(os.path.join(statdir, '*.UTR_reads_RPM.conserved.txt')))		
# for statfile in statfiles:
# 	df1 = pd.read_table(statfile)
# 	if statfile==statfiles[0]:
# 		dfUTR=df1.copy()
# 	else:
# 		dfUTR=pd.merge(dfUTR,df1,on=['gene_symbol'],how='left')
		
# statfiles=sorted(glob.glob(os.path.join(statdir, '*.IPA_LE_reads_RPK.txt')))
# for statfile in statfiles:
# 	df1 = pd.read_table(statfile)
# 	if statfile==statfiles[0]:
# 		dfIPA=df1.copy()
# 	else:
# 		dfIPA=pd.merge(dfIPA,df1,on=['gene_symbol', 'IPAID', 'PASid'],how='left')
		
# dfcds.to_csv(outputdir+project+'.CDS.allsample.txt', index=None, mode='w', sep='\t')
# dfUTR.to_csv(outputdir+project+'.UTR.allsample.txt', index=None, mode='w', sep='\t')		
# dfIPA.to_csv(outputdir+project+'.IPA_LE.allsample.txt', index=None, mode='w', sep='\t')

