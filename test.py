import sys
import gff_parser

fp = open(sys.argv[1])

for n, row in enumerate(gff_parser.read(fp)):
    print row.attributes.Gap
    if n > 100:
        break
