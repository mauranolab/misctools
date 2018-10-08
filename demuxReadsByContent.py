#!/bin/env python3

#Split out reads from a pairs of fastq files based on regexps matching R1 and R1
#This can be useful for amplicon sequencing or other deterministic libraries
#...especially if you forget to program the index reads into basespace

#Usage: demuxReadsByContent.py Undetermined_S0_R1_001.fastq.gz Undetermined_S0_R2_001.fastq.gz

from sys import argv
import sys
import re

import gzip

verbose = False


read1filename = argv[1]
#read1filename = "Undetermined_S0_R1_001.fastq.gz"
read1file = gzip.open(read1filename, 'rb')
read2filename = argv[2]
#read2filename = "Undetermined_S0_R2_001.fastq.gz"
read2file = gzip.open(read2filename, 'rb')


def initializeParser(sample, r1re, r2re, bc1re, bc2re):
    read1filename = sample + "_R1_001.fastq.gz"
    read2filename = sample + "_R2_001.fastq.gz"
    reads1 = gzip.open(read1filename, 'wt')
    reads2 = gzip.open(read2filename, 'wt')
    return (sample, r1re, r2re, bc1re, bc2re, reads1, reads2)


def shutdownParser(parsespec):
    try:
        reads1 = parsespec[5]
        reads1.close()
    except:
        print("Failed to close R1 for " + parsespec[0] + "!!", file=sys.stderr)
    
    try:
        reads2 = parsespec[6]
        reads2.close()
    except:
        print("Failed to close R2 for " + parsespec[0] + "!!", file=sys.stderr)



try:
    #Add your own regexps here, they will be matched greedily against read pairs
    #Format: Sample ID, R1, R2, BC1, BC2
    parsespecs = (initializeParser("BS01572B", "", "", "^CGTACTAG$", ""),
    initializeParser("BS01579A", "", "", "^TCCTGAGC$", ""),
    initializeParser("BS01580A", "", "", "^GGACTCCT$", ""),
    initializeParser("Undetermined_S0", "", "", "", "")
    )
    
    numskipped=0
    numread=0
    #BUGBUG never ends
    while(True or numread<500000000):
        reads1 = [None] * 4
        reads1[0] = read1file.readline().decode().rstrip('\n')
        reads1[1] = read1file.readline().decode().rstrip('\n')
        reads1[2] = read1file.readline().decode().rstrip('\n')
        reads1[3] = read1file.readline().decode().rstrip('\n')
        
        reads2 = [None] * 4
        reads2[0] = read2file.readline().decode().rstrip('\n')
        reads2[1] = read2file.readline().decode().rstrip('\n')
        reads2[2] = read2file.readline().decode().rstrip('\n')
        reads2[3] = read2file.readline().decode().rstrip('\n')
        
        bcs = reads1[0].split(' ')[1].split(":")[3].split("+")
        
        numread+=1
        foundMatch=False
        
        for parsespec in parsespecs:
            if parsespec[1] is None:
                r1match = None
            else:
                r1match = re.search(parsespec[1], reads1[1], flags=re.IGNORECASE) is not None
            
            if parsespec[2] is None:
                r2match = None
            else:
                r2match = re.search(parsespec[2], bcs[0], flags=re.IGNORECASE) is not None
            
            if parsespec[3] is None:
                bc1match = None
            else:
                bc1match = re.search(parsespec[3], bcs[0], flags=re.IGNORECASE) is not None
            
            if parsespec[4] is None:
                bc2match = None
            else:
                bc2match = re.search(parsespec[2], bcs[1], flags=re.IGNORECASE) is not None
            
            if (r1match is not None or r1match is not None) and (r1match is None or r1match) and (r2match is None or r2match) and (bc1match is None or bc1match) and (bc2match is None or bc2match):
                if verbose:
                    print(reads1[0], "\t", reads1[1], "\t", reads2[1], "\t", bcs[0], "\t", bcs[1], "\tr1=", r1match, "; r2=", r2match, "; bc1=", bc1match, "; bc2=", bc2match, "; assigned " + parsespec[0], sep="")
                parsespec[5].write('\n'.join(reads1) + '\n')
                parsespec[6].write('\n'.join(reads2) + '\n')
                
                foundMatch=True
                
                break
        
        if not foundMatch:
            numskipped+=1
            if verbose:
                print(reads1[0], "\t", reads1[1], "\t", reads2[1], "\t", bcs[0], "\t", bcs[1], "\tr1=", r1match, "; r2=", r2match, "; bc1=", bc1match, "; bc2=", bc2match, "; skipped", sep="")

finally:
    read1file.close()
    read2file.close()
    shutdownParser(parsespec)


print("Successfully parsed", numread, ". Skipped", numskipped, "reads.", file=sys.stderr)
