import sys
import gff_parser
import gmap_pygr
import pprint
from pygr import seqdb, cnestedlist
from itertools import izip

genome = '/scratch/titus/lamprey/supercontigs.fa'
ests = '/scratch/titus/lamprey/Pma200805_collapsed_ESTs.fasta'

genome_db = seqdb.SequenceFileDB(genome)
ests_db = seqdb.SequenceFileDB(ests)

def counter(g, when):
    for n, x in enumerate(g):
        if n % when == 0:
            print '...', n
        yield x

gmap_fp = open(sys.argv[1])
al_name = sys.argv[2]

al = cnestedlist.NLMSA(al_name, 'w', pairwiseMode=True)
gmap_pygr.build_alignment(al, counter(gmap_fp, 1000), genome_db, ests_db)
al.build(saveSeqDict=True)
