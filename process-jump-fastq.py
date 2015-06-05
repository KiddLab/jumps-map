
# Jeff Kidd
# process reads from miseq from a jump library, prepare for mapping
import genutils
import sys
from optparse import OptionParser



###############################################################################

USAGE = """
process-jump-fastq.py  --r1fq <read 1 fq.gz>  --r2fq <read 2 fq.gz>  
                       


"""
parser = OptionParser(USAGE)
parser.add_option('--r1fq',dest='r1fq', help = 'name of f1 fq.gz')
parser.add_option('--r2fq',dest='r2fq', help = 'name of f2 fq.gz')

(options, args) = parser.parse_args()

if options.r1fq is None:
    parser.error('r1fq name not given')
if options.r2fq is None:
    parser.error('r2fq not given')

###############################################################################


print 'hello'