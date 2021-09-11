import numpy as np
import sys
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from collections import Counter
from os.path import isfile
import RNA
import parameters as param
from scipy.special import comb


def Hamming_dist(seq1, seq2):
   assert len(seq1) == len(seq2)
   return len([i for i in range(len(seq1)) if seq1[i] != seq2[i]])




L = int(sys.argv[1])

###############################################################################################
print 'phenotypic robustness - for order of subplots' 
###############################################################################################
df_psample = pd.read_csv('./GPmapdata/Psample_seq_struct_data'+param.L_vsstring_parametersPsample[L]+'.csv')
shape_list_sample = df_psample['shape'].tolist()
sequence_list_sample = df_psample['sequence'].tolist()
shape_vs_seq_list = {shape: [sequence_list_sample[i] for i, s in enumerate(shape_list_sample) if s == shape] for shape in set(shape_list_sample) if shape.count('[')  > 0}
rho_sample = df_psample['robustness substitutions'].tolist()
shape_vs_rhoS = {shape: np.mean([rho_sample[i] for i, s in enumerate(shape_list_sample) if s == shape]) for shape in set(shape_list_sample)}
shapes_sorted = sorted(shape_vs_rhoS.keys(), key=shape_vs_rhoS.get, reverse=True)
print len(shapes_sorted)
if len(shapes_sorted)//20 > 1:
   shape_indices_sorted_forsubplots = [int(i) for i in np.arange(0, len(shapes_sorted)-1, step=len(shapes_sorted)//25)]
   while len(shape_indices_sorted_forsubplots) % 5 != 0: # want al columns filled
      shape_indices_sorted_forsubplots = [[i for i in range(len(shapes_sorted)) if i not in shape_indices_sorted_forsubplots][0],] + [s for s in shape_indices_sorted_forsubplots] #add first neglected shape
   shapes_sorted_forsubplots = [shapes_sorted[i] for i in sorted(shape_indices_sorted_forsubplots)]
   print 'subplots contains shape indices:', sorted(shape_indices_sorted_forsubplots)
   assert len(set(shapes_sorted_forsubplots)) == len(shapes_sorted_forsubplots)
else:
   shapes_sorted_forsubplots = shapes_sorted
nrows = len(shapes_sorted_forsubplots)//5 
if len(shapes_sorted_forsubplots)%5 != 0:
   nrows += 1
###############################################################################################
print 'Hamming distance distribution in p-sample' 
###############################################################################################
f, ax = plt.subplots(nrows=nrows, ncols=5, figsize=(12, 1.6*len(shapes_sorted_forsubplots)//5 + 1.6))
for plotindex, shape in enumerate(shapes_sorted_forsubplots):
   seq_list = shape_vs_seq_list[shape]
   Hamming_dist_list = [Hamming_dist(seq1, seq2) for seqindex, seq1 in enumerate(seq_list) for seq2 in seq_list[:seqindex]]
   if len(shapes_sorted_forsubplots)//5+1 > 1:
      axindex = tuple((plotindex//5, plotindex%5))
   else:
      axindex = plotindex
   Hamming_dist_range = [H - 0.5 for H in range(L + 2)]
   ax[axindex].hist(Hamming_dist_list, bins=Hamming_dist_range, density=False)
   ax[axindex].plot(list(range(L + 1)), [len(Hamming_dist_list)*comb(L, H)*3**H/float(4**L) for H in range(L + 1)], c='k')
   ax[axindex].set_xlabel('Hamming dist. in seq. sample')
   ax[axindex].set_ylabel('frequency')
   ax[axindex].set_title(shape)
f.tight_layout()
f.savefig('./Pplots/Psampling_Hamming_dist_test'+param.L_vsstring_parametersPsample[L]+'.png', bbox_inches='tight')
plt.close('all')


     
