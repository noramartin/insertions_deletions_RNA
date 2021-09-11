import numpy as np
from functools import partial
import sys
import pandas as pd
from functions.general_functions import get_shape_from_random_sequence
from collections import Counter
from multiprocessing import Pool
from os.path import isfile
import parameters as param


n_cpus = int(sys.argv[1])
L = int(sys.argv[2])
###################################################################################################
shape_funct_random_seq = partial(get_shape_from_random_sequence, L=L, GPfunction=param.GPfunction)
###################################################################################################
if not isfile('./GPmapdata/gsample_neutralsetdata'+param.L_vsstring_parametersGsample[L]+'.csv'): 
   ###############################################################################################
   print 'g-sample' 
   ###############################################################################################
   pool = Pool(processes = n_cpus)
   pool_result = pool.map(shape_funct_random_seq, np.arange(param.sample_sizeGsample))
   pool.close()
   pool.join()
   ###############################################################################################
   print 'extract data for neutral set size' 
   ###############################################################################################
   sequence_sample_list, shape_list_sample = zip(*pool_result)
   shape_vs_f = Counter(shape_list_sample)
   shapes_sorted_by_f = [s for s in sorted(shape_vs_f.keys(), key=shape_vs_f.get, reverse=True) if s not in ['|', '_']]
   frequency_list = [shape_vs_f[shape] for shape in shapes_sorted_by_f]
   df_sets = pd.DataFrame.from_dict({'structure': shapes_sorted_by_f, 
                                     'number of times in sample': frequency_list})
   df_sets.to_csv('./GPmapdata/gsample_neutralsetdata'+param.L_vsstring_parametersGsample[L]+'.csv')
   del sequence_sample_list, shape_list_sample
   del shape_vs_f, shapes_sorted_by_f, df_sets, pool_result


