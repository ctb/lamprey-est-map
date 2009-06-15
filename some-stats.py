import sys
import gff_parser
import bincount

###

scores = bincount.BinCount()
paths_by_contig = {}
cdna_by_contig = {}
paths_by_cdna = {}

for row in gff_parser.read(open(sys.argv[1])):
    scores.add(row.score)

    contig = row.seqid
    path = row.attributes.ID
    cdna = row.attributes.Name

    s = paths_by_contig.setdefault(contig, set())
    s.add(path)

    s = cdna_by_contig.setdefault(contig, set())
    s.add(cdna)
    
    s = paths_by_cdna.setdefault(cdna, set())
    s.add(path)

###

scores.bin(1).write(open('scores.bin', 'w'), center=False)

paths_by_contig_count = bincount.BinCount()
for v in paths_by_contig.itervalues():
    paths_by_contig_count.add(len(v))
paths_by_contig_count.bin(1).write(open('paths_by_contig.bin', 'w'), center=False)

cdna_by_contig_count = bincount.BinCount()
for v in cdna_by_contig.itervalues():
    cdna_by_contig_count.add(len(v))
cdna_by_contig_count.bin(1).write(open('cdna_by_contig.bin', 'w'), center=False)

paths_by_cdna_count = bincount.BinCount()
for v in paths_by_cdna.itervalues():
    paths_by_cdna_count.add(len(v))
paths_by_cdna_count.bin(1).write(open('paths_by_cdna.bin', 'w'), center=False)
