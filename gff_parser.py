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
        for k in self.keys():
            self.__dict__[k] = self.__getitem__(k)

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

parent_collector_fn = lambda x: x.attributes.Parent
    
def read_groups(fp, collector_fn=None):
    """Read collections, grouping by provided criterion.  Default is Parent.

    'collector_fn' is an optional callable that returns a cookie; this
    cookie is used to decide whether or not to start a new collection.  It
    should return something with reasonable '==' behavior (e.g. a string!)
    and should never return None.

    Note, groupings must be contiguous within file.
    
    Returns lists of Bag objects representing each row.
    """
    if collector_fn is None:
        collector_fn = parent_collector_fn
        
    parent = None
    collect = []
    for row in read(fp):
        if parent != collector_fn(row):
            if collect:
                yield collect
            collect = []
            parent = collector_fn(row)
            
        collect.append(row)

    if collect:
        yield collect

if __name__ == '__main__':
    import doctest
    doctest.testmod()
    
