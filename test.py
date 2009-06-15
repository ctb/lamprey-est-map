import sys
import gff_parser

fp = open(sys.argv[1])

n = 0
for row in gff_parser.read(fp):
    n += 1
    print row.attributes
    if n > 100:
        break
