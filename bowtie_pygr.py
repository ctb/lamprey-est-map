import bowtie_parser

def get_src_sequence(db, row):
    seq = db[row.seqid]
    start = row.start
    stop = start + len(row.read)
    return seq[start:stop]

def get_read_sequence(db, row):
    seq = db[row.readname]
    if row.strand == '-':
        seq = -seq
        
    return seq

def build_ivals(fp, genome_db, reads_db):
    for row in bowtie_parser.read(fp):
        src_seq = get_src_sequence(genome_db, row)
        read_seq = get_read_sequence(reads_db, row)

        yield src_seq, read_seq
        
def build_alignment(al, fp, genome_db, ests_db):
    """Build an NLMSA out of the bowtie mapping in fp."""
    al.add_aligned_intervals(build_ivals(fp, genome_db, ests_db))
