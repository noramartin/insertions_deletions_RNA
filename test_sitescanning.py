import numpy as np
from functools import partial
import sys
import pandas as pd
from functions.general_functions import get_shape_from_random_sequence, decimal_to_scientific_notation
from collections import Counter
from multiprocessing import Pool
from os.path import isfile
import parameters as param
import functions.Psampling_shapes as inv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

def get_phi_pq_from_df(phi_list, norm):
   phi_pq = {}
   for phi in phi_list:
      for shape2 in phi.split('o'):
         if shape2.count('[') > 0:
            try:
               phi_pq[shape2] += 1/float(norm * len(phi_list))
            except KeyError:
               phi_pq[shape2] = 1/float(norm * len(phi_list))
   return phi_pq




n_cpus = int(sys.argv[1])
L = int(sys.argv[2])
###################################################################################################
shape_funct_random_seq = partial(get_shape_from_random_sequence, L=L, GPfunction=param.GPfunction)

###############################################################################################
print 'test site-scanning'
###############################################################################################
if not isfile('./GPmapdata/Gsample_seq_struct_data'+param.L_vsstring_parametersPsample[L]+'.csv'): 
   ###############################################################################################
   print 'g-sample' 
   ###############################################################################################
   pool = Pool(processes = n_cpus)
   pool_result = pool.map(shape_funct_random_seq, np.arange(param.sample_sizeGsample))
   pool.close()
   pool.join()
   sequence_sample_list, shape_list_sample = zip(*pool_result)
   shape_vs_f = Counter(shape_list_sample)
   shape_vs_sequence_sample = {shape: [seq[:] for i, seq in enumerate(sequence_sample_list) if shape_list_sample[i] == shape][:min(500, shape_vs_f[shape])] for shape, count in shape_vs_f.iteritems() if shape.count('[') > 0 and count >= 10}
   del sequence_sample_list, shape_list_sample, pool_result
   ###############################################################################################
   print 'extract data for p-robustness' 
   ###############################################################################################
   function_for_parallelised_calc = partial(inv.get_shape_robustness_from_sequence, GPfunction=param.GPfunction)
   sequence_sample_list = [s for seqlist in shape_vs_sequence_sample.values() for s in seqlist]
   pool = Pool(processes = n_cpus)
   pool_result = pool.map(function_for_parallelised_calc, sequence_sample_list)
   pool.close()
   pool.join()
   shape_list_sample, rho_listS, rho_listI, rho_listD, phi_listS, phi_listI, phi_listD = zip(*pool_result)

   df_G = pd.DataFrame.from_dict({'shape': shape_list_sample, 
                                'sequence': sequence_sample_list,
                                'robustness substitutions': rho_listS,
                                'robustness deletions': rho_listD,
                                'robustness insertions': rho_listI,
                                'phi substitutions': phi_listS,
                                'phi deletions': phi_listD,
                                'phi insertions': phi_listI})
   df_G.to_csv('./GPmapdata/Gsample_seq_struct_data'+param.L_vsstring_parametersPsample[L]+'.csv')
   del shape_list_sample, rho_listS, rho_listI, rho_listD, phi_listS, phi_listI, phi_listD
###############################################################################################
print 'plot - robustness'
###############################################################################################
df_G = pd.read_csv('./GPmapdata/Gsample_seq_struct_data'+param.L_vsstring_parametersPsample[L]+'.csv')
df_P = pd.read_csv('./GPmapdata/Psample_seq_struct_data'+param.L_vsstring_parametersPsample[L]+'.csv')
shape_list_all = [shape for shape in set(df_G['shape'].tolist()) if shape.count('[') > 0 and df_G['shape'].tolist().count(shape) > 10]
shape_list_allP = [shape for shape in set(df_P['shape'].tolist()) if shape.count('[') > 0]
print 'number of shapes g-sample', len(shape_list_all)
print 'number of shapes site-scanning', len(shape_list_allP)
for shape in shape_list_all:
   assert shape in shape_list_allP # nono of the shapes in the g-sample is not in the site-scanning (or p-sample)
f, ax = plt.subplots(ncols=3, figsize=(8,2.2))
for plotindex, type_mutation in enumerate(['substitutions', 'insertions', 'deletions']):
   type_sampling_vs_Probustness_list = {}
   for type_sampling, df in [['site-scanning', df_P], ['g-sampling', df_G]]:
      rho_list_df = df['robustness ' + type_mutation].tolist()
      shape_list_df = df['shape'].tolist()
      type_sampling_vs_Probustness_list[type_sampling] = [np.mean([rho_list_df[i] for i, s in enumerate(shape_list_df) if s == shape]) for shape in shape_list_all]
   ax[plotindex].scatter(type_sampling_vs_Probustness_list['g-sampling'], type_sampling_vs_Probustness_list['site-scanning'], s=5, zorder=0, c=['lime', 'b', 'r'][plotindex])
   ax[plotindex].plot([0, 1], [0, 1], c='k', zorder=-1)
   ax[plotindex].set_xlabel(r'$\rho_{P}$ ' + type_mutation + '\ngenotype-sampling')
   ax[plotindex].set_ylabel(r'$\rho_{P}$ ' + type_mutation + '\nsite-scanning')
   ax[plotindex].set_xlim(0, 1)
   ax[plotindex].set_ylim(0, 1)
   corr_coeff, pvalue = pearsonr(type_sampling_vs_Probustness_list['g-sampling'], type_sampling_vs_Probustness_list['site-scanning'])
   #ax[plotindex].set_title(r'correlation: '+str(round(corr_coeff,2)) + '\np-value='+ decimal_to_scientific_notation(pvalue),fontdict={'fontsize': 'small'})
   ax[plotindex].annotate('ABCDEF'[plotindex], xy=(0.04, 0.86), xycoords='axes fraction', fontweight='bold')   #, fontsize='large', size=14
f.tight_layout()
f.savefig('./Pplots/testrobustness.png')
###############################################################################################
print 'plot - phipq'
###############################################################################################
shape = '[[[_[]]_]_]'
dfG_shape = df_G.loc[df_G['shape'] == shape]
dfP_shape = df_P.loc[df_P['shape'] == shape]
f, ax = plt.subplots(ncols=3, figsize=(8,2.4))
for plotindex, type_mutation in enumerate(['substitutions', 'insertions', 'deletions']):
   type_sampling_vs_phi_pq_list = {}
   for type_sampling, df in [['site-scanning', dfP_shape], ['g-sampling', dfG_shape]]:
      print 'shape for phipq', shape, 'number of sequences', type_sampling, len(df['phi ' + type_mutation].tolist())
      phi_pq = get_phi_pq_from_df(df['phi ' + type_mutation].tolist(), norm ={'substitutions': 3*L, 'deletions': L, 'insertions': 4* (L+1)}[type_mutation])
      type_sampling_vs_phi_pq_list[type_sampling] = [phi_pq[shape2] if shape2 in phi_pq else np.nan for shape2 in shape_list_allP if shape2!= shape]
   ax[plotindex].scatter(type_sampling_vs_phi_pq_list['g-sampling'], type_sampling_vs_phi_pq_list['site-scanning'], s=5, zorder=0, c=['lime', 'b', 'r'][plotindex])
   
   ax[plotindex].set_xlabel([r'$\phi_{pq, S}$ ', r'$\phi_{pq, I}$ ', r'$\phi_{pq, D}$ '][plotindex] + '\ngenotype-sampling')
   ax[plotindex].set_ylabel([r'$\phi_{pq, S}$ ', r'$\phi_{pq, I}$ ', r'$\phi_{pq, D}$ '][plotindex] + '\nsite-scanning')
   ax[plotindex].set_xscale('log')
   ax[plotindex].set_yscale('log')
   min_x = min(np.nanmin(type_sampling_vs_phi_pq_list['g-sampling']), np.nanmin(type_sampling_vs_phi_pq_list['site-scanning']))
   ax[plotindex].plot([0.8*min_x, 2], [0.8*min_x, 2], c='k', zorder=-1)
   ax[plotindex].set_xlim(0.8*min_x, 2)
   ax[plotindex].set_ylim(0.8*min_x, 2)
   x_positive, y_positive = zip(*[(x, y) for (x, y) in zip(type_sampling_vs_phi_pq_list['g-sampling'], type_sampling_vs_phi_pq_list['site-scanning']) if not (np.isnan(x) or np.isnan(y))])
   corr_coeff, pvalue = pearsonr(np.log10(x_positive), np.log10(y_positive))
   #ax[plotindex].set_title(r'log-log correlation: '+str(round(corr_coeff, 2)) + '\np-value='+ decimal_to_scientific_notation(pvalue),fontdict={'fontsize': 'small'})
   ax[plotindex].annotate('ABCDEF'[plotindex], xy=(0.04, 0.86), xycoords='axes fraction', fontweight='bold')
f.tight_layout()
f.savefig('./Pplots/test_phipq.png')
###############################################################################################
print 'plot - robustness versus subsample'
###############################################################################################
subsample_size = 125
rho_list_df = df_P['robustness ' + type_mutation].tolist()
shape_list_df = df_P['shape'].tolist()
shape_list_allP = [shape for shape in set(shape_list_df) if shape.count('[') > 0]
shape_vs_index_list = {shape: [i for i, s in enumerate(shape_list_df) if s == shape] for shape in shape_list_allP}
f, ax = plt.subplots(ncols=3, figsize=(8,2.2))
for plotindex, type_mutation in enumerate(['substitutions', 'insertions', 'deletions']):   
   Probustness_listfull = [np.mean([rho_list_df[i] for i in shape_vs_index_list[shape]]) for shape in shape_list_allP]
   for repetition in range(10):
      Probustness_list_sample = []
      for shape in shape_list_allP:      
         index_list = np.random.choice(shape_vs_index_list[shape], size=subsample_size, replace=False)
         rho_g_list = [rho_list_df[i] for i in index_list]
         assert len(set(index_list)) == len(rho_g_list) == subsample_size
         Probustness_list_sample.append(np.mean(rho_g_list))
      ax[plotindex].scatter(Probustness_listfull, Probustness_list_sample, alpha=0.5, s=5, zorder=0, c=['lime', 'b', 'r'][plotindex])
   ax[plotindex].plot([0, 1], [0, 1], c='k', zorder=-1)
   ax[plotindex].set_xlabel(r'$\rho_{P}$ ' + type_mutation + '\nfull sample\n' + str(param.number_seq_per_str) + ' sequences')
   ax[plotindex].set_ylabel(r'$\rho_{P}$ ' + type_mutation + '\nsubsampling\n' + str(subsample_size) + ' sequences')
   ax[plotindex].set_xlim(0, 1)
   ax[plotindex].set_ylim(0, 1)
   corr_coeff, pvalue = pearsonr(type_sampling_vs_Probustness_list['g-sampling'], type_sampling_vs_Probustness_list['site-scanning'])
   #ax[plotindex].set_title(r'correlation: '+str(round(corr_coeff,2)) + '\np-value='+ decimal_to_scientific_notation(pvalue),fontdict={'fontsize': 'small'})
   ax[plotindex].annotate('ABCDEF'[plotindex], xy=(0.04, 0.86), xycoords='axes fraction', fontweight='bold')   #, fontsize='large', size=14
f.tight_layout()
f.savefig('./Pplots/testrobustness_subsampling.png', dpi=200)
###############################################################################################
print 'plot - phipq versus subsample'
###############################################################################################
shape = '[[[_[]]_]_]'
dfP_shape = df_P.loc[df_P['shape'] == shape]
f, ax = plt.subplots(ncols=3, figsize=(8,2.4))
for plotindex, type_mutation in enumerate(['substitutions', 'insertions', 'deletions']):
   type_sampling_vs_phi_pq_list = {}
   phi_pq_info_sequence_info = df['phi ' + type_mutation].tolist()
   ###
   phi_pq = get_phi_pq_from_df(phi_pq_info_sequence_info, norm ={'substitutions': 3*L, 'deletions': L, 'insertions': 4* (L+1)}[type_mutation])
   phi_pq_list = [phi_pq[shape2] if shape2 in phi_pq else np.nan for shape2 in shape_list_allP if shape2!= shape]
   ###
   phi_pq_info_sequence_info_subsample = np.random.choice(phi_pq_info_sequence_info, size=subsample_size, replace=False)
   assert len(phi_pq_info_sequence_info_subsample) == subsample_size
   phi_pq_subsample = get_phi_pq_from_df(phi_pq_info_sequence_info_subsample, norm ={'substitutions': 3*L, 'deletions': L, 'insertions': 4* (L+1)}[type_mutation])
   phi_pq_list_subsample = [phi_pq_subsample[shape2] if shape2 in phi_pq_subsample else np.nan for shape2 in shape_list_allP if shape2!= shape]
   ###
   ax[plotindex].scatter(phi_pq_list, phi_pq_list_subsample, s=5, zorder=0, c=['lime', 'b', 'r'][plotindex])
   ax[plotindex].set_xlabel([r'$\phi_{pq, S}$ ', r'$\phi_{pq, I}$ ', r'$\phi_{pq, D}$ '][plotindex] + '\nfull sample\n' + str(param.number_seq_per_str) + ' sequences')
   ax[plotindex].set_ylabel([r'$\phi_{pq, S}$ ', r'$\phi_{pq, I}$ ', r'$\phi_{pq, D}$ '][plotindex] + '\nsubsampling\n' + str(subsample_size) + ' sequences')
   ax[plotindex].set_xscale('log')
   ax[plotindex].set_yscale('log')
   min_x = min(np.nanmin(phi_pq_list), np.nanmin(phi_pq_list_subsample))
   ax[plotindex].plot([0.8*min_x, 2], [0.8*min_x, 2], c='k', zorder=-1)
   ax[plotindex].set_xlim(0.8*min_x, 2)
   ax[plotindex].set_ylim(0.8*min_x, 2)
   x_positive, y_positive = zip(*[(x, y) for (x, y) in zip(phi_pq_list, phi_pq_list_subsample) if not (np.isnan(x) or np.isnan(y))])
   corr_coeff, pvalue = pearsonr(np.log10(x_positive), np.log10(y_positive))
   #ax[plotindex].set_title(r'log-log correlation: '+str(round(corr_coeff, 2)) + '\np-value='+ decimal_to_scientific_notation(pvalue),fontdict={'fontsize': 'small'})
   ax[plotindex].annotate('ABCDEF'[plotindex], xy=(0.04, 0.86), xycoords='axes fraction', fontweight='bold')
f.tight_layout()
f.savefig('./Pplots/test_phipq_subsampling.png', dpi=200)





