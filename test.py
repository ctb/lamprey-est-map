import sys
import gff_parser

target_collector_fn = lambda row: row.attributes.Name

fp = open(sys.argv[1])

n = 0
for collect in gff_parser.read_groups(fp, target_collector_fn):
    n += 1
    print set([ c.attributes.Name for c in collect ])
    if n > 100:
        break
