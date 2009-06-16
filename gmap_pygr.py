"""
Parse GFF3 output from GMAP, and build pygr NLMSAs with it.

For GMAP, see:

  http://www.gene.com/share/gmap/
  http://www.gene.com/share/gmap/src/README
  http://bioinformatics.oxfordjournals.org/cgi/content/full/21/9/1859

For pygr, see:
  http://code.google.com/p/pygr/

Author: C. Titus Brown, <ctb@msu.edu>
"""

import gff_parser

def get_src_sequence(db, row):
    """Extract the source sequence from the given GFF row and seqdb."""
    seq = db[row.seqid][row.start - 1:row.end]
    if row.strand == '-':
        seq = -seq

    return seq

def get_dest_sequence(db, row):
    """Parse the Target attribute and extract the associated sequence."""
    name, start, stop, ori = gff_parser.parse_target(row)

    seq = db[name][start - 1:stop]
    if ori == '-':
        seq = -seq

    return seq

def gaps(row):
    """Parse the gmap 'Tags' attribute."""
    g = row.attributes.Gap
    g = g.split(' ')
    for x in g:
        typ, n = x[0], int(x[1:])
        yield typ, n

def build_ivals(fp, genome_db, ests_db):
    """Generate aligned intervals based on the GMAP/GFF3 data in 'fp'.

    Yields (pathname, ivals) where ivals is a list of interval pairs.
    """

    for n, row in enumerate(gff_parser.read(fp)):
        ivals = []
        src_seq = get_src_sequence(genome_db, row)
        dest_seq = get_dest_sequence(ests_db, row)
        
        src_i = 0
        dst_i = 0
        for typ, m in gaps(row):
            if typ == 'M':              # full match
                src_ival = src_seq[src_i:src_i + m]
                dst_ival = dest_seq[dst_i:dst_i + m]
                ivals.append((src_ival, dst_ival))

                src_i += m
                dst_i += m
            elif typ == 'I':            # insertion on dest seq
                dst_i += m
            elif typ == 'D':            # insertion on src seq
                src_i += m
            else:
                raise Exception, "unknown char in Gap attr: %s%d" % (typ, m)

        yield row.attributes.ID, ivals

if __name__ == '__main__':
    pass
