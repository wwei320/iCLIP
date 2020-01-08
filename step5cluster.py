#===============================================================================
#         FILE: iclip_seq step5
#
#  DESCRIPTION: cluster uniquely mapped reads in Stau1 iCLIP data by modifying a few funcs from current PASS pipline
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
#      CREATED: 2019-Oct-05
#===============================================================================
 




def split_sam(sam_file, min_mapq = 10, direction = 'reverse', spike_in = None):
    '''add FM tag to show first aligement position for both plus and minus strand in unique mapped reads sam file
    '''
    pass_file = open(sam_file.replace('.sam', '.FM.sam'), 'w')
    
    
    lap = 0

    with open(sam_file, 'r') as fin:
        for line in fin:
            # Skip header
            if line[0] == '@':
                pass_file.write(line)
                # nonpass_file.write(line)
                # if spike_in: spike_in_file.write(line)
                continue
            # Process each line
            (readname, flag, chromosome, position, mapq, cigar) = line.split()[:6]
            position = int(position)
            mapq = int(mapq)

            if (direction.lower() == 'reverse' and flag == '16') or \
                    (direction.lower() == 'forward' and flag == '0'):
                strand = '+'
               
                lap=position
            elif (direction.lower() == 'reverse' and flag == '0') or \
                    (direction.lower() == 'forward' and flag == '16'):
                strand = '-'

                
                ###M: all mismatch and match ; N: skipped region on ref; I: insertion; D: deletion
                nums = re.split('[MND]', cigar)[:-1]
                covered = sum(int(re.search(r'(\d+)\D+$', x + 'M').group(1)) for x in nums)
              
                # downstream_seq = ps.get_seq(chromosome, strand,
                #                          start = position - t_stretch_len, 
                #                          end = position - 1,
                #                          genome = genome)
                # lap = position
                lap = position + covered - 1
              
            else:
                continue

            pass_file.write(line.strip() + '\tFM:i:%d\n' % (lap))
            

    pass_file.close()


def find_nearby_indexes(v, max_distance):
    """ 
    Yield the indexes of lists containing neighboring positions
    
    Arguments:
    v is like [[100395423, [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 56]],...],
    sorted by the position (100395423, etc). 
    Output:
    """
    # A container to hold list of liss
    index_list = []
    for i in range(len(v))[:-1]:
        if v[i + 1][0] - v[i][0] <= max_distance:
            index_list.append(i)
            continue
        if len(index_list) > 0:
            index_list.append(i)  
            yield index_list
            index_list = []


def merge_nearby_clusters(cs_cluster, max_distance):
    """
    Recursively merge neighboring clusters using a devide and conqurer algorithm
    Neighboring positions located within max_distance from the peak cluster with 
    max RPM merged into the peak cluster. 
    Arguments:
    cs_cluster is like:[position, [num in pass_files]] (Ordered by position)
    [[155603441, [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0]],
     [155603444, [3, 1, 3, 14, 13, 8, 6, 1, 1, 10, 1, 6]],
     [155603445, [8, 1, 4, 13, 20, 8, 6, 0, 5, 19, 5, 19]]]
    Output:
    """
    sample_number = len(cs_cluster[0][1]) // 2  
    # Check base case
    clustered = True
    for i in range(len(cs_cluster))[:-1]:
        if cs_cluster[i + 1][0] - cs_cluster[i][0] <= max_distance:
            clustered = False

    if clustered:
        return True
    # When the clustering is not completed:
    else:
        # Get max of normalized read numbers from all samples and the index for max
        _, i = max((v, i) for i, v in enumerate((sum(pos_num[1][sample_number:])
                                                   for pos_num in cs_cluster)))
        # Merge cs within max_distance from the cs with max read number
        left_index = i  # leftmost position merged
        right_index = i  # rightmost position merged
        for j, (pos, nums) in enumerate(cs_cluster):
            if abs(pos - cs_cluster[i][0]) <= max_distance and \
                    abs(pos - cs_cluster[i][0]) > 0:
                # Combine the read numbers for two CS, sample by sample
                # cs_cluster[i] is like
                #              [155603471, [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0]]
                # nums is like [0, 0, 0, 0, 0, 5, 1, 0, 0, 0, 0, 0]
                cs_cluster[i][1] = [sum(x)
                                    for x in zip(cs_cluster[i][1], nums)]
                cs_cluster[j][1] = 0  
                left_index = min(left_index, j)
                right_index = max(right_index, j)
        # Devide and conqure
        if left_index > 0:
            merge_nearby_clusters(cs_cluster[:left_index],
                                               max_distance)
        if right_index < len(cs_cluster) - 1:
            merge_nearby_clusters(cs_cluster[right_index + 1:],
                                               max_distance)



def cluster_dict(readcounts, max_distance):
    '''Converts a dict of readcount into a dict of clustered readcount.
    
    This function is only meant to be called by cluster_pass_reads(),
    rescue_nonpass_reads(), cluster_cleavage_sites().
    Arguments:
    readcounts, a dict containing something like:
    'chr11:+:31270274': [6, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    '''
    print('Clustering reads ...')
    # Calculate the normalized read numbers and attach to the read numbers
    df = pd.DataFrame(readcounts).T
    normalized = df / df.sum(0)
    df = pd.concat([df, normalized], axis = 1, ignore_index = True)
    readcounts = df.T.to_dict('list')
    # Separate positions based on chrmosome & strand combination
    cpcounts = collections.defaultdict(list)
    while len(readcounts) > 0:
        k, v = readcounts.popitem()  
        k = k.split(':')
        cpcounts[':'.join(k[:2])].append([int(k[2]), v]) 
        # cpcounts.popitem() is like
        # {'chr9:-':[[100395423, [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 56]],...]}
    del readcounts
    print('start cluster')
    # Sort the list of lists for each chromosome & strand combination
    for k, v in cpcounts.items():
        # Sort in place to save memory
        v.sort(key=lambda val: val[0])
        # Get the indexes of lists containing neighboring positions
        for indexes in find_nearby_indexes(v, max_distance):
            
            cs_cluster = v[indexes[0]:(indexes[-1] + 1)]
            # The original list in the dict will be edited in place:
            merge_nearby_clusters(cs_cluster, max_distance)
        # Delete positions with 0 read number after clustering

        print('finish cluster')
        cpcounts[k] = [posi_num for posi_num in v if not posi_num[1] == 0]
    return cpcounts


def cluster_pass_reads(pass_files,
                       output = 'clusters.csv',
                       direction = 'reverse',
                       max_distance = 24):
    """
    Cluster PASS reads within max_distance.
    First generate a table with the following columns: chromosome, strand, pos, 
    num. Then recursively combine reads with LAP within max_distance nt (from 
    position with the highest read number to the position with next highest 
    read number).
    Build a dict of readcounts like {'chr9:-:100395423':[0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 56]}. The [list of int] saves the number of reads from each 
    pass_file. 
    """
    print('Reading sam files containing PASS reads ...')
    file_num = len(pass_files)
    # Pattern for finding first mapped position (FM) in read names
    re_pattern = re.compile('FM:i:(\d+)')

    print('\nre.compile Done!') 
    # Readcounts is a read id counter
    readcounts = {}
    # Go through sam files and get read numbers for each calculated read id
    for i, pass_file in enumerate(pass_files):
       # Read sam file, calculate read id, and count number of reads
        fin = open(pass_file, 'r')
        for line in fin:
            if line[0] == '@':
                continue
            elements = line.split()
            chromosome = elements[2]
            flag = elements[1]

            if (direction == 'reverse' and flag == '16') or \
                    (direction == 'forward' and flag == '0'):
                strand = '+'
            elif (direction == 'reverse' and flag == '0') or \
                    (direction == 'forward' and flag == '16'):
                strand = '-'

            position = re.search(re_pattern, elements[-1]).group(1)

            read_id = ':'.join([chromosome, strand, position])
            # Initialize the list for holding number of reads for each file
            readcounts.setdefault(read_id, [0] * file_num)[i] += 1
            # readcount.popitem() may return something like:
            # ('chr11:+:31270274', [6, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
        fin.close()

    # Converts a dict of readcount into a dict of clustered readcount.
    cpcounts = cluster_dict(readcounts, max_distance)

    # Write the result to disk in csv format
    print('Writing to output file ...')
    with open(output, 'w') as fout:
        # Write header
        sample_string = ','.join([Path(pass_file).stem.split('.')[0].split('/')[-1]
                                  for pass_file in pass_files])
        fout.write(f'chromosome,strand,position,{sample_string}\n')
        # Write read counts for each cluster
        for k in cpcounts:
            chromosome, strand = k.split(':')
            for record in cpcounts[k]:
                position = str(record[0])
                counts = ','.join([str(int(count)) for count
                                   in record[1][:file_num]])
                fout.write(f'{chromosome},{strand},{position},{counts}\n')
    print('CLUSTER Done!')



## 1.1. Load packages
import os
import re
from pathlib import Path
from subprocess import check_output
import collections
import itertools
import numpy as np
import pandas as pd
import seaborn as sns
import multiprocessing as mp
from time import sleep
import sys

sam_dir = Path('/scratch/ww346/analyze/Str_data/Stau1_iCLIP/rawsam')

 # Sequencing direction relative to mRNA sequence. 
SEQ_DIRECTION = 'forward' # Either 'reverse' (antisense) or 'forward' (sense)

# Maximum distance between clustered reads in the same cluster
MAX_DISTANCE = 15

#### step1 add first alignment (FM) tag to sam files 

sam_files = sorted(str(sam_file) for sam_file in sam_dir.glob('*out.unique.sam'))
l = len(sam_files)

print('\nadding FM tags to unique mapped reads from the following files:\n') 
print('\n'.join([sam_file.split('/')[-1] for sam_file in sam_files]))



for sam_file in sam_files:
    print("Processing ", sam_file)
    split_sam(sam_file = sam_file)


print('\nCIGAR Done!') 

upass_files = sorted(str(pass_file) for pass_file in sam_dir.glob('*.FM.sam'))

upass_reads_output = sam_dir/'clusters.unique.reads.csv'

cluster_pass_reads(upass_files, output = upass_reads_output,
                      direction = SEQ_DIRECTION, max_distance = MAX_DISTANCE)




