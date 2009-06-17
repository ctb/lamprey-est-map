from pygr import seqdb, cnestedlist
import pygr_draw
from pygr_draw.annotation import SequenceWrapperFactory, SequenceWrapper

# load in the NLMSAs
est_map = '/u/t/dev/lamprey-est-map/data/52/aln-52'
al1 = cnestedlist.NLMSA(est_map)

# plot the ESTs with names == sequence.id, color == 'red'
wrapper_est = SequenceWrapperFactory(color='red')

reads_filename = '/u/t/dev/lamprey-est-map/data/brain-rnaseq/lamprey30M2H/s_1_al-trim3-20'
al2 = cnestedlist.NLMSA(reads_filename)
wrapper_read = SequenceWrapperFactory(name='')

# because we're not using pygr.Data, the NLMSAs have different seqdbs
# for the supercontigs.  We have to munge them to make them the same
# so that we can do a single draw.  EVIL.  Ignore Me.
sc = al1.seqDict.prefixDict['supercontigs']
al2.seqDict.prefixDict['supercontigs'] = al1.seqDict.prefixDict['supercontigs']
al2.seqDict.dicts[sc] = 'supercontigs'

# choose the sequence to plot
contig1 = al1.seqDict['supercontigs.Contig0'][:5000]

# build a png
picture_class = pygr_draw.BitmapSequencePicture
p = pygr_draw.draw_annotation_maps(contig1, (al1, al2),
                                   picture_class=picture_class,
                                   wrappers=(wrapper_est, wrapper_read))

# write out
image = p.finalize()
filename = 'draw.png'
open(filename, 'w').write(image)
