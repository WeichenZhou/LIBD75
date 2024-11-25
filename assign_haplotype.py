#!/usr/bin/env python
# coding: utf-8

import numpy as np
import random
import pybedtools
import sys
import os
from IPython.utils import io
import pysam
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import seaborn as sns
import math
from scipy.stats import zscore, bernoulli
from sklearn.preprocessing import MinMaxScaler, Normalizer, PowerTransformer
from sklearn.preprocessing import Normalizer
from sklearn.preprocessing import PowerTransformer
from sklearn.preprocessing import StandardScaler
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from itertools import islice, product
from bisect import bisect, bisect_left, bisect_right
from collections import defaultdict
import pickle
import argparse

def InRange(mylist, pos):
    for idx, val in enumerate(mylist):
        if(pos in range(val[0],val[1]+1)):
            return True
    return False

def WhichRange(mylist, pos):
    diff = np.zeros(len(mylist))
    for idx, val in enumerate(mylist):
        diff[idx] = val[1] - val[0] + 1
    for idx, val in enumerate(mylist):
        if(pos in range(val[0],val[1]+1)):
            return [idx, sum(islice(diff, idx))]
    return -1

def LargestSmaller(myList, pos):
    i = bisect(myList, pos)
    if i:
        return myList[i-1]
    else:
        return -1

def SmallestLarger(myList, pos):
    i = bisect(myList, pos)
    if i >= len(myList):
        return np.inf
    elif i:
        return myList[i]
    else:
        return -1

def WindowPercent(cigarList, seq, window_size):
    length = len(seq)
    lastpos = cigarList[len(cigarList)-1][1]
    firstpos = cigarList[0][0]
    lastwindow = lastpos // window_size
    firstwindow = firstpos // window_size

    windows = []
    for i in range(firstwindow,lastwindow+1):
        window_start = i * window_size
        window_end = window_start + window_size - 1

        read_coverage = 0
        for pair in cigarList:
            if pair[1] <= window_start:
                continue
            if pair[0] >= window_end:
                break

            overlap_start = max(pair[0], window_start)
            overlap_end = min(pair[1], window_end)
            read_coverage += overlap_end - overlap_start + 1

        coverage_percentage = read_coverage / length
        windows.append((window_start, coverage_percentage))

    return windows

def GenerateCIGARList(read):
    cigarLine = read.cigar;
    cigarList = []
    cur_pos = read.reference_start + 1
    seq = read.seq
    if len(cigarLine) == 1 and cigarLine[0][0] == 0:
        pass
    else:
        for (cigarType,cigarLength) in cigarLine:
            try:
                if(cigarType == 0 or cigarType == 7 or cigarType == 8): # Match
                    next_pos = cur_pos + cigarLength - 1
                    cigarList.append([cur_pos, next_pos])
                    cur_pos = next_pos + 1
                elif(cigarType == 1): # Insertions # Should I remove it?
                    ins_pos = cur_pos - (read.reference_start + 1)
                    seq = seq[:ins_pos] + seq[ins_pos+cigarLength:]
                elif(cigarType == 2): # Deletion
                    cur_pos = cur_pos + cigarLength
                elif(cigarType == 3): # Skip
                    cur_pos = cur_pos + cigarLength
                elif(cigarType == 4): # Soft clipping
                    sc_pos = cur_pos - (read.reference_start + 1)
                    seq = seq[:sc_pos] + seq[sc_pos+cigarLength:]
                elif(cigarType == 5): # Hard clipping
                    continue
                elif(cigarType == 6): # Padding
                    continue
                else: # Catch and expect later?
                    print("Wrong CIGAR number")
                    sys.exit(1);
            except:
                print("Problem")
                sys.exit(1);
    return cigarList, seq

def AssignHaplotype(read, hetsnp, cigarList, seq):
    haplotype = 0
    for snp in hetsnp.fetch(read.reference_name, read.reference_start, read.reference_end):
        # If SNP is in gap of read (N or D)
        if len(cigarList) != 0 and InRange(cigarList, snp.pos) == False:
            continue
        # If only matches or in range of read
        else:
            gts = [s['GT'] for s in snp.samples.values()]
            # Only matches
            if len(cigarList) == 0:
                startpos = read.reference_start + 1
                diffsum = 0
            else:
                # Start of other ranges
                idx, diffsum = WhichRange(cigarList, snp.pos)
                startpos = cigarList[idx][0]
            # Get the base in read at the SNP position
            #print(read, cigarList, seq, snp.pos, startpos, diffsum, '\n')
            base = seq[snp.pos - startpos + int(diffsum)]

            # If read have not been visited (first SNP overlapping this read)
            if haplotype == 0:
                if gts == [(0, 1)]:
                    if base == snp.alleles[0]:
                        haplotype = 1
                    elif base == snp.alleles[1]:
                        haplotype = 2
                    else: # No match
                        pass
                elif gts == [(1, 0)]:
                    if base == snp.alleles[1]:
                        haplotype = 1
                    elif base == snp.alleles[0]:
                        haplotype = 2
                    else: # No match
                        pass
            else: # If multiple SNPs overlapping with this read
                if gts == [(0, 1)]:
                    # If same haplotype then do nothing
                    if base == snp.alleles[0] and haplotype == 1:
                        pass
                    elif base == snp.alleles[1] and haplotype == 2:
                        pass
                    elif base != snp.alleles[0] or base != snp.alleles[1]:
                        pass
                    else: # If a read overlaps with two different haplotypes, ignore this read
                        haplotype = 3
                        break
                elif gts == [(1, 0)]:
                    if base == snp.alleles[1] and haplotype == 1:
                        pass
                    elif base == snp.alleles[0] and haplotype == 2:
                        pass
                    elif base != snp.alleles[0] or base != snp.alleles[1]: # if miss match, ignore this snp
                        pass
                    else: # If a read overlaps with two different haplotypes, ignore this read
                        haplotype = 3
                        break
    return haplotype

def ExtractNth(myList, i):
    return [item[i] for item in myList]

def GenerateResult(cellname, readfile, outfile, hapcount, hetsnp, snp_pos):
    waiting_reads = {}
    hapcount[cellname] = {'Haplotype1':0, 'Haplotype2':0, 'Multisnp':0, 'PairDisagree':0, 'Total':0, 'UnmappedPair':0}

    for read in readfile.fetch():  # parallel each chromosome
        hapcount[cellname]['Total'] += 1
        if (not (read.is_unmapped) or not (read.mate_is_unmapped)) and read.is_proper_pair:
            seqname = read.query_name
            startpos = read.reference_start
            endpos = read.reference_end
            chromosome = read.reference_name

            closestsnp = SmallestLarger(snp_pos[chromosome], startpos)# Smallest SNP pos larger than startpos

            ###### non-informative #######
            if (
                                closestsnp > endpos
                        ):  # When the read is not overlapping with a SNP but its pair might have haplotype
                if seqname not in waiting_reads:  # Seqname has not been recorde
                    waiting_reads[seqname] = [read, -1]
                else:  # Seqname has been recorded
                    pairs = waiting_reads[seqname]
                    visited_read = pairs[0]
                    visited_haplotype = pairs[1]

                    # Haplotype
                    if visited_haplotype not in [
                        0,
                        -1,
                        3
                    ]:  # If pair has haplotype
                        read.set_tag('HP', visited_haplotype, replace = False)
                        visited_read.set_tag('HP', visited_haplotype, replace = False)
                        outfile.write(read)
                        outfile.write(visited_read)
                    else:  # If pair doesn't have haplotype or 'unphased'
                        read.set_tag('HP', 'Unphased', replace = False)
                        visited_read.set_tag('HP', 'Unphased', replace = False)
                        outfile.write(read)
                        outfile.write(visited_read)
                        #pass

                    waiting_reads.pop(seqname)

            ###### Informative #########
            else:  # When this read overlaps with a SNP
                cigarList, seq = GenerateCIGARList(read)
                haplotype = AssignHaplotype(read, hetsnp, cigarList, seq)
                if haplotype == 3:
                    hapcount[cellname]['Multisnp'] += 1

                if (seqname not in waiting_reads):  # Seqname has not been recorded
                    waiting_reads[seqname] = [read, haplotype]
                else:  # Seqname has been recorded
                    pairs = waiting_reads[seqname]
                    #print(pairs)
                    visited_read = pairs[0]
                    #print(visited_read)
                    visited_haplotype = pairs[1]

                    # Haplotype
                    if visited_haplotype not in [
                        -1,
                        0,
                        3
                    ]:  # Pair has haplotype
                        if (
                                visited_haplotype == haplotype
                        ):  # If same haplotype and none false
                            read.set_tag('HP', visited_haplotype, replace = False)
                            visited_read.set_tag('HP', visited_haplotype, replace = False)
                            outfile.write(read)
                            outfile.write(visited_read)
                            hapcount[cellname]['Haplotype' + str(visited_haplotype)] += 2
                        elif haplotype in [-1, 0, 3]:  # If current false
                            read.set_tag('HP', visited_haplotype, replace = False)
                            visited_read.set_tag('HP', visited_haplotype, replace = False)
                            outfile.write(read)
                            outfile.write(visited_read)
                            hapcount[cellname]['Haplotype' + str(visited_haplotype)] += 2
                        elif ( # Non match
                                visited_haplotype != haplotype
                                and haplotype in [1, 2]
                        ):
                            hapcount[cellname]['PairDisagree'] += 2
                            read.set_tag('HP', 'PairDisagree', replace = False)
                            visited_read.set_tag('HP', 'PairDisagree', replace = False)
                            outfile.write(read)
                            outfile.write(visited_read)
                            #pass
                    else:  # Pair doesn't have haplotype
                        if haplotype not in [-1, 0, 3]:
                            read.set_tag('HP', haplotype, replace = False)
                            visited_read.set_tag('HP', haplotype, replace = False)
                            outfile.write(read)
                            outfile.write(visited_read)
                            hapcount[cellname]['Haplotype' + str(haplotype)] += 2
                        else:
                            pass
                    waiting_reads.pop(seqname)
        else:
            read.set_tag('HP', 'Unphased', replace = False)
            outfile.write(read)

    if len(waiting_reads) != 0:
        for k,v in waiting_reads.items():
            visited_read = v[0]
            visited_haplotype = v[1]
            if visited_haplotype in [1,2]:
                visited_read.set_tag('HP', 'Unphased', replace = False)
                hapcount[cellname]['UnmappedPair'] += 1
                outfile.write(visited_read)
            else:
                visited_read.set_tag('HP', 'Unphased', replace = False)
                outfile.write(visited_read)


    outfile.close()
    return hapcount

# Set up argument parsing
parser = argparse.ArgumentParser(description='Process BAM files and generate haplotype counts.')
parser.add_argument('--input_bam', required=True, help='Path to the input BAM file')
parser.add_argument('--output_bam', required=True, help='Path to the output BAM file')
parser.add_argument('--output_hapcount', required=True, help='Path to the output hapcount pickle file')
parser.add_argument('--hetsnp_vcf', required=True, help='Path to the hetsnp VCF file')

args = parser.parse_args()

# VCF
hetsnp = pysam.VariantFile(args.hetsnp_vcf)

hapcount = {}

# SNP positions
snp_pos = defaultdict(list)
for snp in hetsnp.fetch():
    snp_pos[snp.chrom].append(snp.pos)

# Input BAM file
readfile = pysam.AlignmentFile(args.input_bam, mode='rb')

# Output BAM file
outfile = pysam.AlignmentFile(args.output_bam, "wb", template=readfile)

hapcount = GenerateResult('WGS', readfile, outfile, hapcount, hetsnp, snp_pos)

# Write hapcount to a pickle file
with open(args.output_hapcount, 'wb') as f:
    pickle.dump(hapcount, f)
