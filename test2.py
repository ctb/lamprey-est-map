import sys
import gff_parser
import gmap_pygr
import pprint
from pygr import seqdb
from itertools import izip

genome = '/scratch/titus/lamprey/supercontigs.fa'
ests = '/scratch/titus/lamprey/Pma200805_collapsed_ESTs.fasta'

genome_db = seqdb.SequenceFileDB(genome)
ests_db = seqdb.SequenceFileDB(ests)

fp = open(sys.argv[1])
ivals = gmap_pygr.build_ivals(fp, genome_db, ests_db)

first100 = list(izip(range(100), ivals))

pprint.pprint(first100)
