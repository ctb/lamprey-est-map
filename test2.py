import sys
from pygr import seqdb, cnestedlist

al = cnestedlist.NLMSA(sys.argv[1])

contig1 = al.seqDict['supercontigs.Contig0']
slice = al[contig1]

for (src, dest, edge) in slice.edges():
    print repr(src), repr(dest), edge.pIdentity()
