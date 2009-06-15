import sys
import gff_parser
from pygr import seqdb

genome = '/scratch/titus/lamprey/supercontigs.fa'
ests = '/scratch/titus/lamprey/Pma200805_collapsed_ESTs.fasta'

genome_db = seqdb.SequenceFileDB(genome)
ests_db = seqdb.SequenceFileDB(ests)

fp = open(sys.argv[1])

for n, row in enumerate(gff_parser.read(fp)):
    src_seq = genome_db[row.seqid][row.start - 1:row.end]
    if row.strand == '-':
        src_seq = -src_seq

    ##

    target_name, target_start, target_stop, target_ori = \
                 gff_parser.parse_target(row.attributes.Target)

    target_seq = ests_db[target_name][target_start - 1:target_stop]
    if target_ori == '-':
        target_seq = -target_seq
    
    ##

    ivals = []
    gaps = row.attributes.Gap.split(' ')

    src_i = 0
    dst_i = 0
    for x in gaps:
        typ, m = x[0], int(x[1:])

        if typ == 'M':
            src_ival = src_seq[src_i:src_i + m]
            dst_ival = target_seq[dst_i:dst_i + m]
            ivals.append((src_ival, dst_ival))

            src_i += m
            dst_i += m
        elif typ == 'I':
            dst_i += m
        elif typ == 'D':
            src_i += m
        else:
            raise Exception, "unknown char in Gap attr: %s%d" % (typ, m)

    for (src, dst) in ivals:
        print row.seqid, row.start, row.end, ' == ', row.attributes.Target
        print src
        print dst

        s = str(src)
        d = str(dst)
        ident = 0
        for i in range(len(s)):
            if s[i] == d[i]:
                ident += 1

        print row.score, int(ident / float(len(s)) * 100)
        print row.strand, row.attributes.Target
        print ''

    print '*' * 20

        
    if n > 100:
        break
