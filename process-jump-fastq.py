#!/usr/bin/env python
# Jeff Kidd
# process reads from miseq from a jump library, prepare for mapping
import genutils
import sys
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

#run pear to joint together reads that overlap
run_pear(myData)

