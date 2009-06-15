"""See http://www.sequenceontology.org/gff3.shtml.

Note, for gaps, see

http://may2005.archive.ensembl.org/Docs/wiki/html/EnsemblDocs/CigarFormat.html
"""

import csv
import urllib

DELIMITER='\t'
FIELDNAMES = ['seqid', 'source', 'type', 'start', 'end', 'score',
              'strand', 'phase', 'attributes' ]

class Bag(dict):
    """dict-like class that supports attribute access as well as getitem.

    >>> x = Bag()
    >>> x['foo'] = 'bar'
    >>> x.foo
    'bar'
    
    """
    def __init__(self, *args, **kw):
        dict.__init__(self, *args, **kw)
        self.__dict__.update(kw)

    def __setitem__(self, k, v):
        dict.__setitem__(self, k, v)
        self.__dict__[k] = v

def read(fp):
    """Parse the given fp as GFF3, yielding Bag objects containing row info.
    """
    
    reader = csv.DictReader(fp, delimiter=DELIMITER, fieldnames=FIELDNAMES)

    for row in reader:
        if row['seqid'].startswith('#'):
            continue

        # pick out the attributes and deal with them specially:
        #  - only unquote after splitting on ;
        #  - make them into a dictionary
        
        attr = row['attributes']
        attr = attr.split(';')
        attr_d = Bag()
        for kv in attr:
            k, v = kv.split('=', 1)
            k, v = urllib.unquote(k), urllib.unquote(v)
            attr_d[k] = v

        # unquote all of the other values & store them in a new dict
        row_d = Bag([ (k, urllib.unquote(v)) for (k, v) in row.items() \
                       if k != 'attributes' ])
        
        # save attributes back into it
        row_d['attributes'] = attr_d
        
        # done!
        yield row_d

if __name__ == '__main__':
    import doctest
    doctest.testmod()
    
