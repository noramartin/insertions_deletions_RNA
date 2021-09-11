import numpy as np
from functools import partial
import sys
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from os.path import isfile
from multiprocessing import Pool
from scipy.stats import pearsonr, spearmanr
import parameters as param
import random
from functions.general_functions import decimal_to_scientific_notation
from functions.functions_phipqlocal_analysis import *
from copy import deepcopy



def update_dict_all_data(data_all, similarity_list, seq_list, shape, type_correlation, type_comparison):
   data_all['shape'] += deepcopy([shapes_sorted[shapeindex],] * len(similarity_list))
   data_all['similarity'] += similarity_list[:]
   data_all['sequence'] += seq_list[:]
   data_all['type of correlation'] += deepcopy([type_correlation,] * len(similarity_list))
   data_all['comparison'] += deepcopy([type_comparison,] * len(similarity_list))
   return data_all

###############################################################################################
###############################################################################################
n_cpus = int(sys.argv[1])
L = int(sys.argv[2])
min_nonneutral_neighbours = 3
min_nonneutral_Boltz_sum = 0
#############
df_psample = pd.read_csv('./GPmapdata/Psample_seq_struct_data'+param.L_vsstring_parametersPsample[L]+'.csv')
shape_list_sample = df_psample['shape'].tolist()
sequence_list_sample = df_psample['sequence'].tolist()
rho_sample = df_psample['robustness substitutions'].tolist()
rho_insertions_sample = df_psample['robustness insertions'].tolist()
rho_deletions_sample = df_psample['robustness deletions'].tolist()
shape_vs_rhoS = {shape: np.mean([rho_sample[i] for i, s in enumerate(shape_list_sample) if s == shape]) for shape in set(shape_list_sample)}
shapes_sorted = sorted(shape_vs_rhoS.keys(), key=shape_vs_rhoS.get, reverse=True)
###############################################################################################
if len(shapes_sorted)//20 > 1:
   shape_indices_forsubplots = [int(i) for i in np.arange(0, len(shapes_sorted)-1, step=len(shapes_sorted)//25)]
   while len(shape_indices_forsubplots) % 5 != 0: # want al columns filled
      shape_indices_forsubplots = [[i for i in range(len(shapes_sorted)) if i not in shape_indices_forsubplots][0],] + [s for s in shape_indices_forsubplots] #add first neglected shape
   shapes_sorted_forsubplots = [shapes_sorted[i] for i in sorted(shape_indices_forsubplots)]
   print 'subplots contains shape indices:', sorted(shape_indices_forsubplots)
   assert len(set(shapes_sorted_forsubplots)) == len(shapes_sorted_forsubplots)
else:
   shapes_sorted_forsubplots = shapes_sorted
nrows = len(shapes_sorted_forsubplots)//5
if len(shapes_sorted_forsubplots)%5 != 0:
   nrows += 1
###############################################################################################
phi_listS = df_psample['phi substitutions'].astype('str').tolist()
phi_listD = df_psample['phi deletions'].astype('str').tolist() 
phi_listI = df_psample['phi insertions'].astype('str').tolist()
types_corr = [ 'correlation', 'Bhattacharyya coeff.',
              'Jaccard index (for shapes with frequencies > 0.05)', 
              'Jaccard index (for shapes with frequencies > 0.01)'] 

##########################################################################################################################
##########################################################################################################################
print 'get statistics: generalised genetic correlations and plastogenetic congruence' 
##########################################################################################################################
##########################################################################################################################
all_param = param.L_vsstring_parametersPsample[L] + 'minN' +str(min_nonneutral_neighbours) + 'Boltzsum' +str(int(round(min_nonneutral_Boltz_sum * 100)))
localphipq_df = './GPmapdata/localphipqcomparison'+ all_param + '.csv'
if not isfile(localphipq_df):
   data_all = {'shape': [],
               'similarity': [],
               'type of correlation': [],
               'comparison': [],
               'sequence': []}
   typecorr_vs_locphipq_similaritySI_sameseq, typecorr_vs_locphipq_similaritySD_sameseq = {t: [] for t in types_corr}, {t: [] for t in types_corr}
   typecorr_vs_locphipq_similaritySI_diffseq, typecorr_vs_locphipq_similaritySD_diffseq = {t: [] for t in types_corr}, {t: [] for t in types_corr}
   typecorr_vs_Boltzphipq_similarityS_sameseq, typecorr_vs_Boltzphipq_similarityI_sameseq, typecorr_vs_Boltzphipq_similarityD_sameseq = [{t: [] for t in types_corr} for i in range(3)]
   typecorr_vs_Boltzphipq_similarityS_diffseq, typecorr_vs_Boltzphipq_similarityI_diffseq, typecorr_vs_Boltzphipq_similarityD_diffseq = [{t: [] for t in types_corr} for i in range(3)]
   for type_correlation in types_corr:
      print type_correlation
      function_for_phipqlocal_data = partial(read_in_shape_data_and_find_phiqplocal_correlation_list,
                                             min_nonneutral_neighbours=min_nonneutral_neighbours, min_nonneutral_Boltz_sum=min_nonneutral_Boltz_sum, 
                                             shape_list_sample = shape_list_sample, shapes_list_full= shapes_sorted, sequence_list_sample = sequence_list_sample,
                                             phi_listS=deepcopy(phi_listS), phi_listI=deepcopy(phi_listI), phi_listD = deepcopy(phi_listD), 
                                             type_correlation=type_correlation, Boltzmann_prob_function=param.shape_and_prob_function, range_kbT=param.range_kbT)
      pool = Pool(processes = n_cpus)
      pool_result = pool.map(function_for_phipqlocal_data, shapes_sorted)
      pool.close()
      pool.join()
      for shapeindex, shape in enumerate(shapes_sorted):
         similarity_data = deepcopy(pool_result[shapeindex])
         typecorr_vs_locphipq_similaritySI_sameseq[type_correlation] += similarity_data[0][:] 
         data_all = update_dict_all_data(data_all, similarity_data[0][:], shape=shape, type_correlation=type_correlation, seq_list=similarity_data[10][:], 
                                         type_comparison='phipq similarity same sequence substitutions/insertions')
         typecorr_vs_locphipq_similaritySD_sameseq[type_correlation] +=  similarity_data[1][:] 
         data_all = update_dict_all_data(data_all, similarity_data[1][:], shape=shape, type_correlation=type_correlation, seq_list=similarity_data[10][:],  
                                         type_comparison='phipq similarity same sequence substitutions/deletions')
         typecorr_vs_locphipq_similaritySI_diffseq[type_correlation] +=  similarity_data[2][:] 
         data_all = update_dict_all_data(data_all, similarity_data[2][:], shape=shape, type_correlation=type_correlation, seq_list=similarity_data[10][:], 
                                         type_comparison='phipq similarity different sequence substitutions/insertions')
         typecorr_vs_locphipq_similaritySD_diffseq[type_correlation] +=  similarity_data[3][:] 
         data_all = update_dict_all_data(data_all, similarity_data[3][:], shape=shape, type_correlation=type_correlation, seq_list=similarity_data[10][:],
                                         type_comparison='phipq similarity different sequence substitutions/deletions')
         typecorr_vs_Boltzphipq_similarityS_sameseq[type_correlation] +=  similarity_data[4][:] 
         data_all = update_dict_all_data(data_all, similarity_data[4][:], shape=shape, type_correlation=type_correlation, seq_list=similarity_data[10][:], 
                                         type_comparison='Boltzmann/phipq similarity same sequence substitutions')
         typecorr_vs_Boltzphipq_similarityI_sameseq[type_correlation] +=  similarity_data[5][:] 
         data_all = update_dict_all_data(data_all, similarity_data[5][:], shape=shape, type_correlation=type_correlation, seq_list=similarity_data[10][:], 
                                         type_comparison='Boltzmann/phipq similarity same sequence insertions')
         typecorr_vs_Boltzphipq_similarityD_sameseq[type_correlation] +=  similarity_data[6][:] 
         data_all = update_dict_all_data(data_all, similarity_data[6][:], shape=shape, type_correlation=type_correlation, seq_list=similarity_data[10][:], 
                                         type_comparison='Boltzmann/phipq similarity same sequence deletions')
         typecorr_vs_Boltzphipq_similarityS_diffseq[type_correlation] +=  similarity_data[7][:] 
         data_all = update_dict_all_data(data_all, similarity_data[7][:], shape=shape, type_correlation=type_correlation, seq_list=similarity_data[10][:],
                                         type_comparison='Boltzmann/phipq similarity different sequence substitutions')
         typecorr_vs_Boltzphipq_similarityI_diffseq[type_correlation] +=  similarity_data[8][:] 
         data_all = update_dict_all_data(data_all, similarity_data[8][:], shape=shape, type_correlation=type_correlation, seq_list=similarity_data[10][:],  
                                         type_comparison='Boltzmann/phipq similarity different sequence insertions')
         typecorr_vs_Boltzphipq_similarityD_diffseq[type_correlation] +=  similarity_data[9][:] 
         data_all = update_dict_all_data(data_all, similarity_data[9][:], shape=shape, type_correlation=type_correlation, seq_list=similarity_data[10][:],
                                          type_comparison='Boltzmann/phipq similarity different sequence deletions')
      del function_for_phipqlocal_data
   df_localphipq = pd.DataFrame.from_dict(data_all)
   df_localphipq.to_csv(localphipq_df)
df_localphipq = pd.read_csv(localphipq_df)
##########################################################################################################################
##########################################################################################################################
print 'plot data: generalised genetic correlations and plastogenetic congruence - general stats' 
##########################################################################################################################
##########################################################################################################################
norm_S, norm_I, norm_D = L*(param.K-1), param.K*(L+1), L
for layout in ['', 'two_rows']:
   for type_correlation in types_corr:
      df_specific_correlation = df_localphipq.loc[df_localphipq['type of correlation'] == type_correlation]
      f1, ax1 = plt.subplots(ncols = 2, figsize=(8.4, 2.8))
      ax_freq1, ax_freq2 = ax1
      if layout == '':
         f2, ax2 = plt.subplots(ncols = 3, figsize=(8.4, 2.5))
         ax_Boltz1, ax_Boltz2, ax_Boltz3 = ax2
      else:
         f2, ax2 = plt.subplots(ncols = 2, nrows=2, figsize=(6.5, 5))
         ax_Boltz1, ax_Boltz2, ax_Boltz3 = ax2[0, 1], ax2[1, 0], ax2[1, 1]
         ax2[0, 0].set_visible(False)  
      # first row: generalised genetic correlations
      locphipq_similaritySI_sameseq = df_specific_correlation.loc[df_specific_correlation['comparison'] == 'phipq similarity same sequence substitutions/insertions']['similarity'].tolist()
      locphipq_similaritySI_diffseq = df_specific_correlation.loc[df_specific_correlation['comparison'] == 'phipq similarity different sequence substitutions/insertions']['similarity'].tolist()
      locphipq_similaritySD_sameseq = df_specific_correlation.loc[df_specific_correlation['comparison'] == 'phipq similarity same sequence substitutions/deletions']['similarity'].tolist()
      locphipq_similaritySD_diffseq = df_specific_correlation.loc[df_specific_correlation['comparison'] == 'phipq similarity different sequence substitutions/deletions']['similarity'].tolist()
      print type_correlation
      list_nonan = [x for x in locphipq_similaritySI_sameseq if not np.isnan(x)]
      print 'SI_sameseq', min(list_nonan), np.percentile(list_nonan, 25), np.percentile(list_nonan, 50), np.percentile(list_nonan, 75)
      list_nonan = [x for x in locphipq_similaritySI_diffseq if not np.isnan(x)]
      print 'SI_diffseq', min(list_nonan), np.percentile(list_nonan, 25), np.percentile(list_nonan, 50), np.percentile(list_nonan, 75)
      list_nonan = [x for x in locphipq_similaritySD_sameseq if not np.isnan(x)]
      print 'SD_sameseq', min(list_nonan), np.percentile(list_nonan, 25), np.percentile(list_nonan, 50), np.percentile(list_nonan, 75)
      list_nonan = [x for x in locphipq_similaritySD_diffseq if not np.isnan(x)]
      print 'SD_diffseq', min(list_nonan), np.percentile(list_nonan, 25), np.percentile(list_nonan, 50), np.percentile(list_nonan, 75)
      min_value = min(locphipq_similaritySI_sameseq + locphipq_similaritySI_diffseq + locphipq_similaritySD_sameseq + locphipq_similaritySD_diffseq)
      if type_correlation.startswith ('Jaccard index'):
         min_freq = float(type_correlation.split('>')[1][:-1])
         if abs(min_freq - 0.05) < 0.0001:
            step = 0.1
         elif  abs(min_freq - 0.01) < 0.0001:
            step = 0.05
         bins = np.arange(0, 1 + 0.999 * step, step=step) # because there can be only up to 1/min_freq frequent shapes (if rho=0), the Jaccard index has a low resolution
      else:
         bins = np.linspace(min(0, min_value), 1, 20)
      ax_freq1.hist(locphipq_similaritySI_sameseq, facecolor='b', label='subst. vs ins.', alpha=0.3, density=False, bins=bins)
      ax_freq1.hist(locphipq_similaritySI_diffseq, facecolor='grey', label='baseline model', alpha=0.5, density=False, bins=bins)
      ax_freq2.hist(locphipq_similaritySD_sameseq, facecolor='r', label='subst. vs del.', alpha=0.3, density=False, bins=bins)
      ax_freq2.hist(locphipq_similaritySD_diffseq, facecolor='grey', label='baseline model', alpha=0.5, density=False, bins=bins)
      # second row: Boltzmann
      Boltzphipq_similarityS_sameseq = df_specific_correlation.loc[df_specific_correlation['comparison'] == 'Boltzmann/phipq similarity same sequence substitutions']['similarity'].tolist()
      Boltzphipq_similarityS_diffseq = df_specific_correlation.loc[df_specific_correlation['comparison'] == 'Boltzmann/phipq similarity different sequence substitutions']['similarity'].tolist()
      Boltzphipq_similarityI_sameseq = df_specific_correlation.loc[df_specific_correlation['comparison'] == 'Boltzmann/phipq similarity same sequence insertions']['similarity'].tolist()
      Boltzphipq_similarityI_diffseq = df_specific_correlation.loc[df_specific_correlation['comparison'] == 'Boltzmann/phipq similarity different sequence insertions']['similarity'].tolist()      
      Boltzphipq_similarityD_sameseq = df_specific_correlation.loc[df_specific_correlation['comparison'] == 'Boltzmann/phipq similarity same sequence deletions']['similarity'].tolist()
      Boltzphipq_similarityD_diffseq = df_specific_correlation.loc[df_specific_correlation['comparison'] == 'Boltzmann/phipq similarity different sequence deletions']['similarity'].tolist()
      min_value = min(Boltzphipq_similarityS_sameseq + Boltzphipq_similarityS_diffseq + Boltzphipq_similarityI_sameseq + Boltzphipq_similarityI_diffseq + Boltzphipq_similarityD_sameseq + Boltzphipq_similarityD_diffseq)
      if type_correlation.startswith ('Jaccard index'):
         min_freq = float(type_correlation.split('>')[1][:-1])
         bins = np.arange(0, 1 + 0.999 * step, step=step) # because there can be only up to 1/min_freq frequent shapes (if rho=0), the Jaccard index has a low resolution
      else:
         bins = np.linspace(min(0, min_value), 1, 20)
      ax_Boltz1.hist(Boltzphipq_similarityS_sameseq, facecolor='lime', label='Boltzmann vs subst.', alpha=0.3, density=False, bins=bins)
      ax_Boltz1.hist(Boltzphipq_similarityS_diffseq, facecolor='grey', label='baseline model', alpha=0.5, density=False, bins=bins)
      ax_Boltz2.hist(Boltzphipq_similarityI_sameseq, facecolor='b', label='Boltzmann vs ins.', alpha=0.3, density=False, bins=bins)
      ax_Boltz2.hist(Boltzphipq_similarityI_diffseq, facecolor='grey', label='baseline model', alpha=0.5, density=False, bins=bins)
      ax_Boltz3.hist(Boltzphipq_similarityD_sameseq, facecolor='r', label='Boltzmann vs del.', alpha=0.3, density=False, bins=bins)
      ax_Boltz3.hist(Boltzphipq_similarityD_diffseq, facecolor='grey', label='baseline model', alpha=0.5, density=False, bins=bins)
      for ax_hist in [ax_freq1, ax_freq2, ax_Boltz1, ax_Boltz2, ax_Boltz3]:
         ax_hist.set_ylabel('frequency')
         xmin, xmax = ax_hist.get_xlim()
         ax_hist.set_xlim(min(0, xmin), max(1, xmax))
         ax_hist.legend(frameon=False, loc='upper left', fontsize='small')
         y1, y2 = ax_hist.get_ylim()
         ax_hist.set_ylim(y1, y2 * 1.25)
      for plotindex, ax_hist in enumerate([ax_freq1, ax_freq2]):
         if len(type_correlation) > 15:
            ax_hist.set_xlabel(type_correlation+' for\n'+r'$\widetilde{\phi}^{(g)}_{q}$ (subst.)'+' vs. '+[r'$\widetilde{\phi}^{(g)}_{q}$ (ins.)', r'$\widetilde{\phi}^{(g)}_{q}$ (del.)'][plotindex])  
         else:
            ax_hist.set_xlabel(type_correlation+' for '+r'$\widetilde{\phi}^{(g)}_{q}$ (subst.)'+' vs. '+[r'$\widetilde{\phi}^{(g)}_{q}$ (ins.)', r'$\widetilde{\phi}^{(g)}_{q}$ (del.)'][plotindex])             
      for plotindex, ax_hist in enumerate([ax_Boltz1, ax_Boltz2, ax_Boltz3]):
         if type_correlation.startswith('Jaccard'):
            min_freq = float(type_correlation.split('>')[1][:-1])
            ax_hist.set_xlabel('Jaccard index (for shapes\nwith frequencies > ' + str(min_freq) + ')'+' for\n'+r'Boltzmann frequ.'+' vs.'+[r' $\widetilde{\phi}^{(g)}_{q}$ (subst.)', r' $\widetilde{\phi}^{(g)}_{q}}$ (ins.)', r' $\widetilde{\phi}^{(g)}_{q}}$ (del.)'][plotindex])        
         else:
            if type_correlation == 'correlation of non-zero frequencies':
               ax_hist.set_xlabel('correlation of non-zero values'+'\n'+r'for Boltzmann frequ.'+' vs.'+[r' $\widetilde{\phi}^{(g)}_{q}$ (subst.)', r' $\widetilde{\phi}^{(g)}_{q}$ (ins.)', r' $\widetilde{\phi}^{(g)}_{q}}$ (del.)'][plotindex])        
            else:
               ax_hist.set_xlabel(type_correlation+'\n'+r'for Boltzmann frequ.'+' vs.'+[r' $\widetilde{\phi}^{(g)}_{q}$ (subst.)', r' $\widetilde{\phi}^{(g)}_{q}$ (ins.)', r' $\widetilde{\phi}^{(g)}_{q}}$ (del.)'][plotindex])                    
      for axindex, ax_tolabel in enumerate([ ax_freq1, ax_freq2]):
         ax_tolabel.annotate('ABCDEFGH'[axindex], xy=(0.90, 0.86), xycoords='axes fraction', fontsize='large', fontweight='bold', size=15)
      for axindex, ax_tolabel in enumerate([ax_Boltz1, ax_Boltz2, ax_Boltz3]):
         ax_tolabel.annotate('ABCDEFGH'[axindex], xy=(0.90, 0.86), xycoords='axes fraction', fontsize='large', fontweight='bold', size=15)
      f1.tight_layout()
      f1.savefig('./Pplots/Psampling_phi_pq'+all_param+'_'.join(type_correlation.split())+'_1.png', bbox_inches='tight', dpi=300)
      f2.tight_layout()
      f2.savefig('./Pplots/Psampling_phi_pq'+all_param+'_'.join(type_correlation.split())+ layout +'_2.png', bbox_inches='tight', dpi=300)
      plt.close('all')
      del f1, f2, ax1, ax2

##########################################################################################################################
##########################################################################################################################
print 'plot data: generalised genetic correlations and plastogenetic congruence -  stats per shape' 
##########################################################################################################################
##########################################################################################################################
bins = np.linspace(-1, 1, 20)
for type_correlation in types_corr:
   df_specific_correlation = df_localphipq.loc[df_localphipq['type of correlation'] == type_correlation]
   f, ax = plt.subplots(nrows=5, figsize=(10, 10.3), sharex=True)
   for comparison_index, comparison in enumerate(['phipq similarity same sequence substitutions/insertions', 'phipq similarity same sequence substitutions/deletions',
                   'Boltzmann/phipq similarity same sequence substitutions', 'Boltzmann/phipq similarity same sequence insertions', 
                   'Boltzmann/phipq similarity same sequence deletions']):
      
      df_specific_comparison = df_specific_correlation.loc[df_specific_correlation['comparison'] == comparison]
      df_nullmodel = df_specific_correlation.loc[df_specific_correlation['comparison'] == comparison.replace('same', 'different')]
      fraction_data_higher_nullmodel_list = []
      for plotindex, shape in enumerate(shapes_sorted):
         data_for_shape = df_specific_comparison.loc[df_specific_comparison['shape'] == shape]
         data_for_shape_nullmodel = df_nullmodel.loc[df_nullmodel['shape'] == shape]
         seq_list1, similarity_list1 = data_for_shape['sequence'].tolist(), data_for_shape['similarity'].tolist()
         seq_list2, similarity_list2 = data_for_shape_nullmodel['sequence'].tolist(), data_for_shape_nullmodel['similarity'].tolist()
         data_higher_nullmodel_list = []
         for seq in set(seq_list1):
            data_for_seq = [similarity for sequence, similarity in zip(seq_list1, similarity_list1) if seq == sequence]
            data_for_seq_nullmodel = [similarity for sequence, similarity in zip(seq_list2, similarity_list2) if seq == sequence]
            assert len(data_for_seq) == 1 and len(data_for_seq_nullmodel) == 1
            if data_for_seq[0] > data_for_seq_nullmodel[0]:
               data_higher_nullmodel_list.append(1)
            elif not (np.isnan(data_for_seq[0]) or np.isnan(data_for_seq_nullmodel[0])):
               data_higher_nullmodel_list.append(0)
            if np.isnan(data_for_seq[0]) or np.isnan(data_for_seq_nullmodel[0]):
               print 'npnan in data', shape, type_correlation, comparison
         fraction_data_higher_nullmodel_list.append(np.mean(data_higher_nullmodel_list))
         if np.mean(data_higher_nullmodel_list) < 0.5:
            print type_correlation, comparison, 'fraction data higher than baseline model ', np.mean(data_higher_nullmodel_list), 'for', shape
      ax[comparison_index].plot(np.arange(len(shapes_sorted)), fraction_data_higher_nullmodel_list, ls='-', marker='o', c=['b', 'r', 'lime', 'b', 'r'][comparison_index], ms=2) 
      ax[comparison_index].plot([-1, len(shapes_sorted)], [0.5, 0.5], c='grey')
      if comparison_index < 2:
         annotation_text = r'$\widetilde{\phi}^{(g)}_{q}$ (subst.)'+' vs. '+[r'$\widetilde{\phi}^{(g)}_{q}$ (ins.)', r'$\widetilde{\phi}^{(g)}_{q}$ (del.)'][comparison_index]
      else:
         annotation_text = r'Boltzmann frequ.'+' vs.'+[r' $\widetilde{\phi}^{(g)}_{q}$ (subst.)', r' $\widetilde{\phi}^{(g)}_{q}$ (ins.)', r' $\widetilde{\phi}^{(g)}_{q}$ (del.)'][comparison_index - 2]

      ax[comparison_index].set_ylabel('fraction of sequences:\ncorrelation higher\nin sequence-specific data\nthan in baseline model\n'+ annotation_text, fontsize='small') 
      ax[comparison_index].set_ylim(0, 1)
      ax[comparison_index].set_xlim(-1, len(shapes_sorted))
      ax[comparison_index].annotate(annotation_text, xy=(0.95, 0.05), xycoords='axes fraction', ha='right', size=13) 
      ax[comparison_index].set_xlabel('shape rank\n(sorted in descending order '+r'by phenotype robustness to substitutions)')
   #f.text(-0.05, 0.5, r'fraction: correlation higher than in the baseline model', va='center', rotation=90, size=13)
   f.tight_layout()
   f.savefig('./Pplots/Psampling_phi_pq'+all_param+'summary'+'_'.join(type_correlation.split())+'_pershape.png', bbox_inches='tight', dpi=300)
   plt.close('all')

###############################################################################################
###############################################################################################
print 'plot data: phi_pq per neutral set with one subplot per p' 
###############################################################################################
###############################################################################################
gsampling_df = pd.read_csv('./GPmapdata/gsample_neutralsetdata'+param.L_vsstring_parametersGsample[L]+'.csv')
structure_vs_freq_estimate = {row['structure']: row['number of times in sample']/float(param.sample_sizeGsample) for rowindex, row in gsampling_df.iterrows()}
###############################################################################################
f, ax = plt.subplots(nrows=nrows, ncols=5, figsize=(12, 1.8 * len(shapes_sorted_forsubplots)//5 + 1.8))
for plotindex, shape in enumerate(shapes_sorted_forsubplots):
   phi_pqS, phi_pqI, phi_pqD = {s: 0.0 for s in shapes_sorted if s!=shape}, {s: 0.0 for s in shapes_sorted if s!=shape}, {s: 0.0 for s in shapes_sorted if s!=shape}
   for i, s in enumerate(shape_list_sample):
      if s == shape:
         list_neighboursS = phi_listS[i].split('o')
         list_neighboursD = phi_listD[i].split('o')
         list_neighboursI = phi_listI[i].split('o')
         for s2 in phi_pqS:
            phi_pqS[s2] += list_neighboursS.count(s2)/float((param.K-1) * L * param.number_seq_per_str)
            phi_pqD[s2] += list_neighboursD.count(s2)/float(L * param.number_seq_per_str)
            phi_pqI[s2] += list_neighboursI.count(s2)/float((L+1) * param.K * param.number_seq_per_str)
   if len(shapes_sorted_forsubplots)//5+1 > 1:
      axindex = tuple((plotindex//5, plotindex%5))
   else:
      axindex = plotindex
   if plotindex == 0: #for legend
      ax[axindex].scatter([phi_pqS[s] for s in phi_pqS.keys()], [phi_pqI[s] for s in phi_pqS.keys()], label='insertions', s=5, c='b', alpha=0.7, lw=0)
      ax[axindex].scatter([phi_pqS[s] for s in phi_pqS.keys()], [phi_pqD[s] for s in phi_pqS.keys()], label='deletions', s=5, c='r', alpha=0.7, lw=0)
   else:
      ax[axindex].scatter([phi_pqS[s] for s in phi_pqS.keys()], [phi_pqI[s] for s in phi_pqS.keys()], s=5, c='b', alpha=0.7, lw=0)
      ax[axindex].scatter([phi_pqS[s] for s in phi_pqS.keys()], [phi_pqD[s] for s in phi_pqS.keys()], s=5, c='r', alpha=0.7, lw=0)    
   if plotindex%5 == 0:
      ax[axindex].set_ylabel(r'$\phi_{pq, D}$ and $\phi_{pq, I}$')
   ax[axindex].set_xlabel(r'$\phi_{qp, S}$')
   ax[axindex].set_title(r'shape, $p$: ' + shape)
   ax[axindex].set_xscale('log')
   ax[axindex].set_yscale('log')
   min_axlim = 0.5 * min([min([phi for phi in phi_pqS.values() if phi >0]), min([phi for phi in phi_pqI.values() if phi >0]), min([phi for phi in phi_pqD.values() if phi >0]), 1.0/param.sample_sizeGsample])
   ax[axindex].set_xlim(min_axlim, 2)
   ax[axindex].set_ylim(min_axlim, 2)
   ax[axindex].plot([min_axlim, 1], [min_axlim, 1], c='k', lw=0.5)
   del phi_pqS, phi_pqI, phi_pqD
f.legend(bbox_to_anchor=(0.5, 1.02), ncol=2, loc='upper center')
f.tight_layout()
f.savefig('./Pplots/Psampling_phi_pq_allshapes'+param.L_vsstring_parametersPsample[L]+'.png', bbox_inches='tight', dpi=200)
plt.close('all')
del f, ax


####################################################################################################################
####################################################################################################################
print 'plot data: phi_pq per neutral set' 
####################################################################################################################
####################################################################################################################
for type_correlation, correlation_function in [('Pearson', logPearsonR_nonzero_values), ('Spearman', SpearmanR_nonzero_values)]:
   ####################################################################################################################
   ####################################################################################################################
   print 'plot data: phi_pq per neutral set - correlations with mean Boltzmann probabilities' 
   ####################################################################################################################
   ####################################################################################################################
   ###############################################################################################
   print 'get data: phi_pq vs mean Boltzmann probabilities  - correlations' 
   ###############################################################################################
   shape_vs_phi_pqS, shape_vs_phi_pqD, shape_vs_phi_pqI, shape_vs_mean_Boltzmann = {}, {}, {}, {}
   shape_vs_DFindex_list = {shape: [dfi for dfi, s in enumerate(shape_list_sample) if s == shape] for shape in set(shape_list_sample)}
   shape_vs_DFindex_listshuffled = {shape: random.sample(dfindex_list, k=len(dfindex_list)) for shape, dfindex_list in shape_vs_DFindex_list.iteritems()}
   num_seq_per_type_mutation = param.number_seq_per_str//4
   shape_vs_DFindex_listS = {shape: dfindex_list[:num_seq_per_type_mutation] for shape, dfindex_list in shape_vs_DFindex_listshuffled.iteritems()}
   shape_vs_DFindex_listI = {shape: dfindex_list[num_seq_per_type_mutation:2*num_seq_per_type_mutation] for shape, dfindex_list in shape_vs_DFindex_listshuffled.iteritems()}
   shape_vs_DFindex_listD = {shape: dfindex_list[2*num_seq_per_type_mutation:3*num_seq_per_type_mutation] for shape, dfindex_list in shape_vs_DFindex_listshuffled.iteritems()}
   shape_vs_DFindex_listBoltz = {shape: dfindex_list[3*num_seq_per_type_mutation:] for shape, dfindex_list in shape_vs_DFindex_listshuffled.iteritems()}
   corr_Boltz_phipqS, corr_Boltz_phipqD, corr_Boltz_phipqI, corr_phipqSD, corr_phipqSI = [], [], [], [], []
   for plotindex, shape in enumerate(shapes_sorted):
      phi_pqS, phi_pqI, phi_pqD, Boltz_f = {s: 0 for s in shapes_sorted if s!=shape}, {s: 0 for s in shapes_sorted if s!=shape}, {s: 0 for s in shapes_sorted if s!=shape}, {s: 0 for s in shapes_sorted if s!=shape}
      for dfindex in shape_vs_DFindex_listS[shape]:
         list_neighboursS = phi_listS[dfindex].split('o')
         for s2 in phi_pqS:
            if list_neighboursS.count(s2):
               phi_pqS[s2] += list_neighboursS.count(s2)/float((param.K-1) * L * num_seq_per_type_mutation)
      for dfindex in shape_vs_DFindex_listI[shape]:
         list_neighboursI = phi_listI[dfindex].split('o')
         for s2 in phi_pqI:
            if list_neighboursI.count(s2):
               phi_pqI[s2] += list_neighboursI.count(s2)/float((L + 1) * param.K * num_seq_per_type_mutation)
      for dfindex in shape_vs_DFindex_listD[shape]:
         list_neighboursD = phi_listD[dfindex].split('o')
         for s2 in phi_pqD:
            if list_neighboursD.count(s2):
               phi_pqD[s2] += list_neighboursD.count(s2)/float(L * num_seq_per_type_mutation)   
      for dfindex in shape_vs_DFindex_listBoltz[shape]:
         shape_vsP = param.shape_and_prob_function(sequence_list_sample[dfindex], range_kbT=param.range_kbT)
         for s2 in shape_vsP:
            if s2 in Boltz_f and s2 not in [shape, '_', '|']: 
               Boltz_f[s2] += shape_vsP[s2]/float(num_seq_per_type_mutation) 
            elif s2 not in [shape, '_', '|']:
               Boltz_f[s2] = shape_vsP[s2]/float(num_seq_per_type_mutation)           
      Boltz_list = [Boltz_f[s] for s in phi_pqS.keys()]
      phi_list_shapeI = [phi_pqI[s] for s in phi_pqS.keys()]
      phi_list_shapeD = [phi_pqD[s] for s in phi_pqS.keys()]
      phi_list_shapeS = [phi_pqS[s] for s in phi_pqS.keys()] 
      corr_phipqSD.append(correlation_function(phi_list_shapeS, phi_list_shapeD)[0])
      corr_phipqSI.append(correlation_function(phi_list_shapeS, phi_list_shapeI)[0])
      corr_Boltz_phipqS.append(correlation_function(Boltz_list, phi_list_shapeS)[0])
      corr_Boltz_phipqD.append(correlation_function(Boltz_list, phi_list_shapeD)[0])
      corr_Boltz_phipqI.append(correlation_function(Boltz_list, phi_list_shapeI)[0])
      shape_vs_phi_pqS[shape] = {s2: value for s2, value in phi_pqS.iteritems()}
      shape_vs_phi_pqD[shape] = {s2: value for s2, value in phi_pqD.iteritems()}
      shape_vs_phi_pqI[shape] = {s2: value for s2, value in phi_pqI.iteritems()}
      shape_vs_mean_Boltzmann[shape] = {s2: value for s2, value in Boltz_f.iteritems()}
   df_corr =  pd.DataFrame.from_dict({'shape': shapes_sorted, type_correlation + ' correlation phiqp: substitutions vs deletions': corr_phipqSD, 
                                      type_correlation + ' correlation phiqp: substitutions vs insertions': corr_phipqSI, 
                                      type_correlation + ' correlation mean Boltzmann frequ. vs phiqp substitutions': corr_Boltz_phipqS,
                                      type_correlation + ' correlation mean Boltzmann frequ. vs phiqp deletions': corr_Boltz_phipqD,
                                      type_correlation + ' correlation mean Boltzmann frequ. vs phiqp insertions': corr_Boltz_phipqI})
   df_corr.to_csv('./GPmapdata/Psampling_phenotype_phipq_correlations'+type_correlation+param.L_vsstring_parametersPsample[L]+'.csv')
   ###############################################################################################
   print 'plot data: phi_pq correlations - one example and stats (use seperate sequences to evaluate phipq for the three types of mutation)' 
   ###############################################################################################
   #indiced_shapes = list(set([0, 1, 2, 3, 4, 5, 6, 7] + [int(i) for i in np.arange(0, len(shapes_sorted)-1, step=len(shapes_sorted)//10)]))
   for shapeindex in sorted(shape_indices_forsubplots)[:12]:
      f, ax = plt.subplots(ncols=2, nrows=2, figsize=(8, 5.7))
      shape = shapes_sorted[shapeindex]   
      phi_pqS, phi_pqD, phi_pqI, Boltz_f  =  shape_vs_phi_pqS[shape], shape_vs_phi_pqD[shape], shape_vs_phi_pqI[shape], shape_vs_mean_Boltzmann[shape]  
      sc1 = ax[0,0].scatter([phi_pqS[s] for s in phi_pqS.keys()], [phi_pqI[s] for s in phi_pqS.keys()], label=r'$\phi_{qp, I}$', s=5, c='b')
      sc2 = ax[0,0].scatter([phi_pqS[s] for s in phi_pqS.keys()], [phi_pqD[s] for s in phi_pqS.keys()], label=r'$\phi_{qp, D}$', s=5, c='r')
      ax[0,0].set_ylabel(r'$\phi_{qp, I}$  and $\phi_{qp, D}$')
      ax[0,0].set_xlabel(r'$\phi_{qp, S}$')
      #ax[0].set_title(r'shape, $p$: ' + shape)
      ax[0,0].set_xscale('log')
      ax[0,0].set_yscale('log')
      min_axlim = 0.5 * min([min([phi for phi in phi_pqS.values() if phi >0]), min([phi for phi in phi_pqI.values() if phi >0]), min([phi for phi in phi_pqD.values() if phi >0])]) #, 1.0/param.sample_sizeGsample])
      ax[0,0].plot([min_axlim, 2], [min_axlim, 2], c='k', lw=0.5)
      ax[0,0].set_xlim(min_axlim, 2)
      ax[0,0].set_ylim(min_axlim, 2)
      #ax[0,0].legend(edgecolor='k', facecolor='darkgrey', loc='lower right', title=r'$p=$ '+shape, framealpha=0.5, fontsize='x-small')
      ax[0,1].legend(handles=[sc1, sc2], edgecolor='k', facecolor='darkgrey', loc='lower left', title=r'$p=$ '+shape, framealpha=0.5, fontsize='x-small')
      corr_coeffI, pvalueI = spearmanr([phi_pqS[s] for s in phi_pqS.keys()], [phi_pqI[s] for s in phi_pqS.keys()])
      corr_coeffD, pvalueD = spearmanr([phi_pqS[s] for s in phi_pqS.keys()], [phi_pqD[s] for s in phi_pqS.keys()])
      #ax[0,0].set_title(r'Insertions: '+str(round(corr_coeffI,2)) + ', p-value='+ str(decimal_to_scientific_notation(pvalueI)) + '\n' + r'Deletions: '+str(round(corr_coeffD,2)) + ', p-value='+ str(decimal_to_scientific_notation(pvalueD)),
      #                   fontdict={'fontsize': 'medium'})
      ###
      df_corr = pd.DataFrame.from_dict({type_correlation + ' correlation coeff.': corr_phipqSD + corr_phipqSI, 'mutations': ['subst. vs\ndel.',]*len(corr_phipqSD) + ['subst. vs\nins.',]*len(corr_phipqSI)})
      sns.boxplot(data=df_corr, x='mutations', y=type_correlation+ ' correlation coeff.', ax=ax[0,1], whis=(5,95), fliersize=2, palette={'subst. vs\ndel.':'r', 'subst. vs\nins.': 'b'})
      ax[0,1].set_xlabel('')
      ax[0,1].set_ylabel(type_correlation + ' correlation coeff.' + '\n' + r'between $\log\ \phi_{qp}$ for'+'\ndifferent mutation types')
      ax[0,1].set_ylim(0,1)
      ###
      sc3 = ax[1,0].scatter([Boltz_f[s] for s in phi_pqS.keys()], [phi_pqS[s] for s in phi_pqS.keys()], label=r'$\phi_{qp, S}$', s=5, c='g')
      sc4 = ax[1,0].scatter([Boltz_f[s] for s in phi_pqS.keys()], [phi_pqI[s] for s in phi_pqS.keys()], label=r'$\phi_{qp, I}$', s=5, c='b')
      sc5 = ax[1,0].scatter([Boltz_f[s] for s in phi_pqS.keys()], [phi_pqD[s] for s in phi_pqS.keys()], label=r'$\phi_{qp, D}$', s=5, c='r')
      ax[1,0].set_ylabel(r'$\phi_{qp}$ for $p=$ '+shape)
      ax[1,0].set_xlabel('mean Boltzmann freq of\n' + r'shape $q$' + '\nin neutral set of shape' + r' $p=$ '+shape)
      #ax[0].set_title(r'shape, $p$: ' + shape)
      ax[1,0].set_xscale('log')
      ax[1,0].set_yscale('log')
      min_axlim = 0.7 * min([min([0.5,]+[g for g in Boltz_f.values() if g >0]), min([0.5,]+[phi for phi in phi_pqS.values() if phi >0]), min([0.5,]+[phi for phi in phi_pqI.values() if phi >0]), min([0.5,]+[phi for phi in phi_pqD.values() if phi >0])]) #, 1.0/param.sample_sizeGsample])
      ax[1,0].plot([min_axlim, 2], [min_axlim, 2], c='k', lw=0.5)
      ax[1,0].set_xlim(min_axlim, 2)
      ax[1,0].set_ylim(min_axlim, 2)
      #ax[1,0].legend(edgecolor='k', facecolor='darkgrey', loc='lower right', framealpha=0.5, fontsize='x-small')
      ax[1,1].legend(handles=[sc3, sc4, sc5], edgecolor='k', facecolor='darkgrey', loc='lower left', framealpha=0.5, fontsize='x-small')
      corr_coeffS, pvalueS = spearmanr([Boltz_f[s] for s in phi_pqS.keys()], [phi_pqS[s] for s in phi_pqS.keys()])
      corr_coeffI, pvalueI = spearmanr([Boltz_f[s] for s in phi_pqS.keys()], [phi_pqI[s] for s in phi_pqS.keys()])
      corr_coeffD, pvalueD = spearmanr([Boltz_f[s] for s in phi_pqS.keys()], [phi_pqD[s] for s in phi_pqS.keys()])
      #ax[1,0].set_title(r'Substitutions: '+str(round(corr_coeffS,2)) + ', p-value='+ str(decimal_to_scientific_notation(pvalueS)) + '\n' +r'Insertions: '+str(round(corr_coeffI,2)) + ', p-value='+ str(decimal_to_scientific_notation(pvalueI)) + '\n' + r'Deletions: '+str(round(corr_coeffD,2)) + ', p-value='+ str(decimal_to_scientific_notation(pvalueD)),
      #                   fontdict={'fontsize': 'small'})
      ###
      df_corr2 = pd.DataFrame.from_dict({type_correlation+ ' correlation coeff': corr_Boltz_phipqS + corr_Boltz_phipqI + corr_Boltz_phipqD, 
                                          'mutations': ['subst.',]*len(corr_Boltz_phipqS) + ['ins.',]*len(corr_Boltz_phipqI) + ['del.',]*len(corr_Boltz_phipqD)})
      sns.boxplot(data=df_corr2, x='mutations', y=type_correlation + ' correlation coeff', ax=ax[1,1], whis=(5,95), 
                      fliersize=2, palette={'subst.':'lime', 'del.': 'r', 'ins.': 'b'})
      ax[1,1].set_ylabel(type_correlation + ' correlation coeff. between' + '\n' + r'mean Boltzmann freq. of $q$'+ '\n' + r'and $\log\ \phi_{qp}$')
      ax[1,1].set_xlabel('')
      ax[1,1].set_ylim(0,1)
      for i in range(2):
         x1, x2 = ax[i,1].get_xlim()
         ax[i,1].set_xlim(x1 - 0.1 * abs(x2-x1), x2)
      for axindex, axi in enumerate([ax[0,0], ax[0,1], ax[1,0], ax[1,1],]):
            axi.annotate('ABCDEF'[axindex], xy=(0.04, 0.86), xycoords='axes fraction', fontsize='large', fontweight='bold', size=13)       
      f.tight_layout()
      f.savefig('./Pplots/Psampling_phi_pq_Boltzmann'+type_correlation+param.L_vsstring_parametersPsample[L]+'example'+str(shapeindex)+'.png', bbox_inches='tight')
      plt.close('all')
      del f, ax
      df_plot =  pd.DataFrame.from_dict({'shape': phi_pqS.keys(), 'Boltzmann freq.': [Boltz_f[s] for s in phi_pqS.keys()], 
                                         'phiqp substitutions': [phi_pqS[s] for s in phi_pqS.keys()], 
                                          'phiqp insertions': [phi_pqI[s] for s in phi_pqS.keys()], 
                                          'phiqp deletions': [phi_pqD[s] for s in phi_pqS.keys()]})
      df_plot.to_csv('./data_plots/Psampling_phi_pq_Boltzmann'+type_correlation+param.L_vsstring_parametersPsample[L]+'example'+str(shapeindex)+'.csv')
      del phi_pqS, phi_pqI, phi_pqD, Boltz_f

