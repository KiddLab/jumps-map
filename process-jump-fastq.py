#!/usr/bin/env python
# Jeff Kidd
# process reads from miseq from a jump library, prepare for mapping
# added read length test 07-10-2015


import fastqstats
import subprocess
import sys
from Bio import pairwise2
import os
from optparse import OptionParser


###############################################################################
# Helper function to run commands, handle return values
def runCMD(cmd):
    val = subprocess.Popen(cmd, shell=True).wait()
    if val == 0:
        pass
    else:
        print 'command failed'
        print cmd
        sys.exit(1)
###############################################################################
#####################################################################
# some utility functions
def open_gzip_read(fileName):
    gc = 'gunzip -c ' + fileName
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)  # To deal with fact that might close file before reading all
    try:
        inFile = os.popen(gc, 'r')
    except:
        print "ERROR!! Couldn't open the file " + fileName + " using gunzip -c\n"
        sys.exit(1)
    return inFile
#####################################################################
def open_gzip_write(fileName):
    try:
        gc = 'gzip > ' + fileName
        outFile = os.popen(gc, 'w')
    except:
        print "ERROR!! Couldn't open the output file " + fileName+ " (with gzip)\n"
        sys.exit(1)
    return outFile
#####################################################################
##############################################################################
# Returns complement of a bp.  If not ACGT then return same char
def complement(c):
    if c == 'A':
        return 'T'
    if c == 'T':
        return 'A'
    if c == 'C':
        return 'G'
    if c == 'G':
        return 'C'
    if c == 'a':
        return 't'
    if c == 't':
        return 'a'
    if c == 'c':
        return 'g'
    if c == 'g':
        return 'c'
    # If not ACTGactg simply return same character
    return c
##############################################################################
# Returns the reverse compliment of sequence 
def revcomp(seq):
    c = ''
    seq = seq[::-1] #reverse
    # Note, this good be greatly sped up using list operations
    seq = [complement(i) for i in seq]
    c = ''.join(seq)
    return c
##############################################################################
#############################################################################
# Takes an open fastq file and returns a dictionary of information about
# the next sequence record.  Returns None if the end of the file
# has been reached.  Assumes that PHRED quality is ord(i)-33, can use different offset if needed
# The structure of the dictionary is:
# seqName, seq, qualString, qualList
def get_next_seq_record(myFile, int qoffSet = 33):
    record = {}
    record['readName'] = myFile.readline()
    if record['readName'] == '':
        return None
    record['readName'] = record['readName'].rstrip()
    if record['readName'][0] != '@':
        print "ERROR!! In get_next_seq_record expect first char of '"+ record['readName'] +"' to be '@'"
        sys.exit(1)
    record['readName'] = record['readName'][1:]
    record['seq'] = myFile.readline()
    if record['seq'] == '':
        print "ERROR!! In get_next_seq_record expect sequence line from fastq file"
        sys.exit(1)
    record['seq'] = record['seq'].rstrip()
    plusLine = myFile.readline()
    if plusLine == '' or plusLine[0] != '+':
        print "ERROR!! In get_next_seq_record expect line from fastq file to begin with '+'"
        sys.exit(1)
    record['qualString'] = myFile.readline()
    if record['qualString'] == '':
        print "ERROR!! In get_next_seq_record expect qual line from fastq file"
        sys.exit(1)
    record['qualString'] = record['qualString'].rstrip()
    lord = ord #local function pointer for speedup    
    record['qualList'] = [lord(i) - qoffSet for i in record['qualString'] ]
    record['qual33Str'] = ''.join([chr(i+33) for i in record['qualList']])
    return record
###############################################################################
def run_pear(myData):
    # PEAR aligns/merges overlapping read pairs, which is the case that we have here
    myData['pearBase'] = myData['outDir'] + myData['sampleName'] + '.pear'

    cmd = 'pear --nbase -f %s -r %s -o %s' % (myData['r1fq'],myData['r2fq'],myData['pearBase'])
    
    myData['assembledFQ'] = myData['pearBase'] + '.assembled.fastq'
    myData['discardedFQ'] = myData['pearBase'] + '.discarded.fastq'
    myData['notAssemF'] = myData['pearBase'] + '.unassembled.forward.fastq'
    myData['notAssemR'] = myData['pearBase'] + '.unassembled.reverse.fastq'

    # check to see if should run
    outgz = myData['assembledFQ'] + '.gz'
    if os.path.isfile(outgz) is True:
        print 'found gzip output already, will not rerun'
        myData['assembledFQ'] += '.gz'
        myData['discardedFQ'] += '.gz'
        myData['notAssemF'] += '.gz'
        myData['notAssemR'] += '.gz'
    else:
        print cmd
        runCMD(cmd)
        cmd = 'gzip ' + myData['assembledFQ']
        print cmd
        runCMD(cmd)
        myData['assembledFQ'] += '.gz'

        cmd = 'gzip ' + myData['discardedFQ']
        print cmd
        runCMD(cmd)
        myData['discardedFQ'] += '.gz'

        cmd = 'gzip ' + myData['notAssemF']
        print cmd
        runCMD(cmd)
        myData['notAssemF'] += '.gz'
        
        cmd = 'gzip ' + myData['notAssemR']
        print cmd
        runCMD(cmd)
        myData['notAssemR'] += '.gz'        
###############################################################################
def count_num_not_assembled(myData):
    myData['numNotAssem'] = 0
    fqFile = open_gzip_read(myData['notAssemF'])
    while True:
        R1 = get_next_seq_record(fqFile)
        if R1 is None: break    
        myData['numNotAssem'] += 1
    fqFile.close()
###############################################################################
def count_num_discarded(myData):
    myData['numDiscarded'] = 0
    fqFile = open_gzip_read(myData['discardedFQ'])
    while True:
        R1 = get_next_seq_record(fqFile)
        if R1 is None: break    
        myData['numDiscarded'] += 1
    fqFile.close()
###############################################################################
def process_assembled(myData):     
    # take each read contig (from PEAR alignment) and search for the expected
    # linker sequence.  If we find it, extract the r1 and r2 and process
    # into separate files.  Takes of reverse complement of r1, and imposes
    # min score on linker alignment.  Tests if reads are longer than 22 bp.  If not, they are
    # written in the lenFail fastq files

    myData['newR1FileName'] = myData['outDir'] + myData['sampleName'] + '.processed.R1.fq.gz'
    myData['newR2FileName'] = myData['outDir'] + myData['sampleName'] + '.processed.R2.fq.gz'    
    outR1 = open_gzip_write(myData['newR1FileName'])
    outR2 = open_gzip_write(myData['newR2FileName'])

    myData['lenFail_R1FileName'] = myData['outDir'] + myData['sampleName'] + '.lenFail.R1.fq.gz'
    myData['lenFail_R2FileName'] = myData['outDir'] + myData['sampleName'] + '.lenFail.R2.fq.gz'    
    lenFailR1 = open_gzip_write(myData['lenFail_R1FileName'])
    lenFailR2 = open_gzip_write(myData['lenFail_R2FileName'])
    
    myData['numAssembled'] = 0
    myData['numOK'] = 0
    myData['numFail'] = 0
    myData['lenFail'] = 0

    fqFile = open_gzip_read(myData['assembledFQ'])
    while True:
        R1 = get_next_seq_record(fqFile)
        if R1 is None: break    
        myData['numAssembled'] += 1
        res = check_seq(R1,myData)
        if res['passChecks'] is True:
            passLength = read_len_test(res)
            if passLength is True:  
                myData['numOK'] += 1 
                name = R1['readName']
                name = name.split()[0]
                name1 = name + ' 1'
                name2 = name + ' 2'
                outR1.write('@%s\n%s\n+\n%s\n' % (name1,res['seq1'],res['seq1Qual']))
                outR2.write('@%s\n%s\n+\n%s\n' % (name2,res['seq2'],res['seq2Qual']))
        else:
            myData['numFail']  += 1

        if myData['numAssembled']  % 25000 == 0:
            print '\tProcesssed %i assembled seqs...' % (myData['numAssembled'])
        
        if res['passChecks'] is True and passLength is False:
            myData['lenFail'] += 1
            name = R1['readName']
            name = name.split()[0]
            name1 = name + ' 1'
            name2 = name + ' 2'
            lenFailR1.write('@%s\n%s\n+\n%s\n' % (name1,res['seq1'],res['seq1Qual']))           
            lenFailR2.write('@%s\n%s\n+\n%s\n' % (name2,res['seq2'],res['seq2Qual']))   

#        if myData['numAssembled']  >= 1000:
#            break
    fqFile.close()    
    myData['totReads'] = myData['numAssembled'] + myData['numDiscarded'] + myData['numNotAssem']
    outR1.close()
    outR2.close()    
###############################################################################
def check_seq(fq,myData):
    minScore = 49.0
    result = {}
    result['passChecks'] = False    

    # do not know the alignment orientation, so check both and take best score
    alignRes = pairwise2.align.globalms(myData['linkerSeq'], fq['seq'], 2, -1, -.5, -.2,penalize_end_gaps=False)
    alignResLinkerRC = pairwise2.align.globalms(myData['linkerSeqRC'], fq['seq'], 2, -1, -.5, -.2,penalize_end_gaps=False)    
    # compare scores
    if alignRes[0][2] >= alignResLinkerRC[0][2]:
        ls = myData['linkerSeq']
    else:
        ls = myData['linkerSeqRC']
        alignRes = alignResLinkerRC        
    result['align'] = alignRes

    # should only be one alignment.  Otherwise 
    if len(alignRes) != 1:
#        print 'have mulitple potential alignments'
#        print fq['seq']
        result['passChecks'] = False
        return result
        
 
    # check score
    if alignRes[0][2] < minScore:
        result['passChecks'] = False
        return result

    #figure out coordinates
    # go through keeping track of pos to go form column number to bp number
    # do the sequence in 1 based coordinates
    seq1ColToPos = []
    current = 0
    for i in range(len(alignRes[0][0])): #[1,2,3,4], looking at alignment
        if alignRes[0][0][i] != '-':
            current += 1
        seq1ColToPos.append(current)
    seq2ColToPos = []
    current = 0
    for i in range(len(alignRes[0][1])):
        if alignRes[0][1][i] != '-':
            current += 1
        seq2ColToPos.append(current) 

    linkerColStart = -1
    linkerColEnd = -1
    for i in range(len(seq1ColToPos)):
        if seq1ColToPos[i] == 1 and linkerColStart == -1:
            linkerColStart = i
        if seq1ColToPos[i] == len(myData['linkerSeq']) and linkerColEnd == -1:
            linkerColEnd = i
    

    # extract sequences -1 because python 0 based, colToSeq 1 based, 
    leftSeq = fq['seq'][0:seq2ColToPos[linkerColStart]-1]
    linkerSeq = fq['seq'][seq2ColToPos[linkerColStart]-1:seq2ColToPos[linkerColEnd]]
    rightSeq = fq['seq'][seq2ColToPos[linkerColEnd]:]    

    
    

    result['passChecks'] = True    
    # passess, so take out the sequence and quals
    leftSeqQual = fq['qual33Str'][0:seq2ColToPos[linkerColStart]-1]
    rightSeqQual = fq['qual33Str'][seq2ColToPos[linkerColEnd]:]
    
    # need to reverse comp R1 due to structure of the library
    leftSeq = revcomp(leftSeq)
    leftSeqQual = leftSeqQual[::-1]
    
    result['seq1'] = leftSeq
    result['seq1Qual'] = leftSeqQual
    result['seq2'] = rightSeq
    result['seq2Qual'] = rightSeqQual

    

    
    return result
###############################################################################

#Counts the length of both reads, if either is shorter than 22 it fails this test

def read_len_test(res):
    readLen = len(res['seq1']) > 22 and len(res['seq2']) > 22 and len(res['seq1']) < 30 and len(res['seq2']) < 30 
    #nCount = res['seq1'].count('N') < 22 and res['seq2'].count('N') < 22
    return readLen

###############################################################################

USAGE = """
process-jump-fastq.py  --r1fq <read 1 fq.gz>  --r2fq <read 2 fq.gz>  --sample <name of sample>  --outdir <dir of output>

    This script sets up reads for processing of Talkowski style jumping libraries.
    Starts with paired end fastq.  Merges together reads,  searches for adapter sequence,
    then splits out to new paired-end files suitable for mapping.
                       

"""
parser = OptionParser(USAGE)
parser.add_option('--r1fq',dest='r1fq', help = 'name of f1 fq.gz')
parser.add_option('--r2fq',dest='r2fq', help = 'name of f2 fq.gz')
parser.add_option('--sample',dest='sampleName', help = 'name of sample')
parser.add_option('--outdir',dest='outDir', help = 'name of output dir')



(options, args) = parser.parse_args()

if options.r1fq is None:
    parser.error('r1fq name not given')
if options.r2fq is None:
    parser.error('r2fq not given')
if options.sampleName is None:
    parser.error('sampleName not given')
if options.outDir is None:
    parser.error('out put dir not given')

###############################################################################
if options.outDir[-1] != '/':
    options.outDir += '/'

# setup file location info
myData = {}
myData['filesToDelete'] = []
myData['filesToGzip'] = []
myData['r1fq'] = options.r1fq
myData['r2fq'] = options.r2fq
myData['sampleName'] = options.sampleName
myData['outDir'] = options.outDir
myData['linkerSeq'] = 'CTGCTGTACCGTTCTCCGTACAGCAG'
# rev of linker is also possible, since do not know orietnation
myData['linkerSeqRC'] = revcomp(myData['linkerSeq'])


print 'Processing %s' % myData['sampleName']
#run pear to join together reads that overlap
run_pear(myData)
count_num_not_assembled(myData)
print '%i reads were not assembled' % myData['numNotAssem']
count_num_discarded(myData)
print '%i reads were discarded' % myData['numDiscarded']

process_assembled(myData)
print '%i reads were assembled' % myData['numAssembled']
print '%i assembled reads failed the check' % myData['numFail']
print '%i assembled reads that passed the check but failed the length test (< 23bp or > 29)' % myData['lenFail']
print '%i total reads in original fastq' % myData['totReads']


