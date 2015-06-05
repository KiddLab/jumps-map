#!/usr/bin/env python
# Jeff Kidd
# process reads from miseq from a jump library, prepare for mapping
import genutils
import fastqstats
import sys
from Bio import pairwise2
import os
from optparse import OptionParser


###############################################################################
def run_pear(myData):
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
        genutils.runCMD(cmd)
        cmd = 'gzip ' + myData['assembledFQ']
        print cmd
        genutils.runCMD(cmd)
        myData['assembledFQ'] += '.gz'

        cmd = 'gzip ' + myData['discardedFQ']
        print cmd
        genutils.runCMD(cmd)
        myData['discardedFQ'] += '.gz'

        cmd = 'gzip ' + myData['notAssemF']
        print cmd
        genutils.runCMD(cmd)
        myData['notAssemF'] += '.gz'
        
        cmd = 'gzip ' + myData['notAssemR']
        print cmd
        genutils.runCMD(cmd)
        myData['notAssemR'] += '.gz'        
###############################################################################
def count_num_not_assembled(myData):
    myData['numNotAssem'] = 0
    fqFile = genutils.open_gzip_read(myData['notAssemF'])
    while True:
        R1 = fastqstats.get_next_seq_record(fqFile)
        if R1 is None: break    
        myData['numNotAssem'] += 1
    fqFile.close()
###############################################################################
def count_num_discarded(myData):
    myData['numDiscarded'] = 0
    fqFile = genutils.open_gzip_read(myData['discardedFQ'])
    while True:
        R1 = fastqstats.get_next_seq_record(fqFile)
        if R1 is None: break    
        myData['numDiscarded'] += 1
    fqFile.close()
###############################################################################
def process_assembled(myData):     
    
    myData['newR1FileName'] = myData['outDir'] + myData['sampleName'] + '.processed.R1.fq.gz'
    myData['newR2FileName'] = myData['outDir'] + myData['sampleName'] + '.processed.R2.fq.gz'
    
    outR1 = genutils.open_gzip_write(myData['newR1FileName'])
    outR2 = genutils.open_gzip_write(myData['newR2FileName'])
    
    
    myData['numAssembled'] = 0
    myData['numOK'] = 0
    myData['numFail'] = 0
    fqFile = genutils.open_gzip_read(myData['assembledFQ'])
    while True:
        R1 = fastqstats.get_next_seq_record(fqFile)
        if R1 is None: break    
        myData['numAssembled'] += 1
        res = check_seq(R1,myData)
        if res['passChecks'] is True:
            myData['numOK'] += 1
            name = R1['readName']
            name = name.split()[0]
            name1 = name + ' 1'
            name2 = name + ' 2'
            outR1.write('@%s\n%s\n+\n%s\n' % (name1,res['seq1'],res['seq1Qual']))
            outR2.write('@%s\n%s\n+\n%s\n' % (name2,res['seq2'],res['seq2Qual']))


#            print myData['numAssembled']
#            print '>%i_f' % myData['numAssembled']
#            print res['seq1']
#            print '>%i_r' % myData['numAssembled']
#            print res['seq2']
        else:
            myData['numFail']  += 1
#            print res['align']
        
            
            
        
        if myData['numAssembled']  % 10000 == 0:
            print '\tProcesssed %i assembled seqs...' % (myData['numAssembled'])
        
        if myData['numAssembled']  >= 50:
            break
    fqFile.close()    
    myData['totReads'] = myData['numAssembled'] + myData['numDiscarded'] + myData['numNotAssem']
    outR1.close()
    outR2.close()    
###############################################################################
def check_seq(fq,myData):
    minScore = 49.0
    result = {}
    result['passChecks'] = False    

#    alignRes = pairwise2.align.globalms(myData['linkerSeq'], fq['seq'], 2, -1, -.5, -.1)
    alignRes = pairwise2.align.globalms(myData['linkerSeq'], fq['seq'], 2, -1, -.5, -.2,penalize_end_gaps=False)
    alignResLinkerRC = pairwise2.align.globalms(myData['linkerSeqRC'], fq['seq'], 2, -1, -.5, -.2,penalize_end_gaps=False)
    
    # compare scores
    if alignRes[0][2] >= alignResLinkerRC[0][2]:
        ls = myData['linkerSeq']
    else:
        ls = myData['linkerSeqRC']
        alignRes = alignResLinkerRC
        

    result['align'] = alignRes

    if len(alignRes) != 1:
        # just take the first one
        print 'have mulitple potential alignments'
        print fq['seq']
        result['passChecks'] = False
        return result
        
 
    # check score
    if alignRes[0][2] < minScore:
        result['passChecks'] = False
        return result

    #figure out coordinates
    seq1ColToPos = []
    current = 0
    for i in range(len(alignRes[0][0])):
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
    

    
    leftSeq = fq['seq'][0:seq2ColToPos[linkerColStart]-1]
    linkerSeq = fq['seq'][seq2ColToPos[linkerColStart]-1:seq2ColToPos[linkerColEnd]]
    rightSeq = fq['seq'][seq2ColToPos[linkerColEnd]:]

    # no longer needed -- since we are going with the score
    # check linker sequence length
#    if len(linkerSeq) != len(myData['linkerSeq']):
#        result['passChecks'] = False
#        return result
#    mismatch = 0
#    for i in range(len(linkerSeq)):
#        if linkerSeq[i] != myData['linkerSeq'][i]:
#            mismatch += 1
#    if mismatch > 1:
#        result['passChecks'] = False
#        return result
    
    result['passChecks'] = True
    
    # passess, so take out the sequence and quals
    leftSeqQual = fq['qual33Str'][0:seq2ColToPos[linkerColStart]-1]
    rightSeqQual = fq['qual33Str'][seq2ColToPos[linkerColEnd]:]
    
    # need to reverse r1
    leftSeq = genutils.revcomp(leftSeq)
    leftSeqQual = leftSeqQual[::-1]
    
    result['seq1'] = leftSeq
    result['seq1Qual'] = leftSeqQual
    result['seq2'] = rightSeq
    result['seq2Qual'] = rightSeqQual
    
    return result
###############################################################################


###############################################################################

USAGE = """
process-jump-fastq.py  --r1fq <read 1 fq.gz>  --r2fq <read 2 fq.gz>  --sample <name of sample>  --outdir <dir of output>
                       


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
# rev of linker is also possible...
myData['linkerSeqRC'] = genutils.revcomp(myData['linkerSeq'])


print 'Processing %s' % myData['sampleName']
#run pear to join together reads that overlap
run_pear(myData)
count_num_not_assembled(myData)
print '%i reads were not assembled' % myData['numNotAssem']
count_num_discarded(myData)
print '%i reads were discarded' % myData['numDiscarded']


process_assembled(myData)
print '%i reads were assembled' % myData['numAssembled']
print '%i total reads in original fastq' % myData['totReads']

