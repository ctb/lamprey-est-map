import sys
import bowtie_parser
import bowtie_pygr
from pygr import seqdb, cnestedlist

genome_db = seqdb.SequenceFileDB('/scratch/titus/lamprey/supercontigs.fa')
reads_db = seqdb.SequenceFileDB('/scratch/titus/lamprey/brain-rnaseq/lamprey30M2H/s_1_sequence.fasta')

def counter(g, when):
    for n, x in enumerate(g):
        if n % when == 0:
            print '...', n
        yield x

al = cnestedlist.NLMSA(sys.argv[1], mode='w', pairwiseMode=True)
fp = open(sys.argv[2])
bowtie_pygr.build_alignment(al, counter(fp, 1000), genome_db, reads_db)
al.build(saveSeqDict=True)
