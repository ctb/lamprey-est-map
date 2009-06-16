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

fp = open(sys.argv[1])
paths = gmap_pygr.build_ivals(fp, genome_db, ests_db)

al = cnestedlist.NLMSA('foo', 'memory', pairwiseMode=True)
for n, (path, ivals) in enumerate(paths):
    if n > 100:
        break
    al.add_aligned_intervals(ivals)

al.build()
