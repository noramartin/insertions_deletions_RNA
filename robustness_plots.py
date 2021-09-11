import numpy as np
import sys
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from os.path import isfile
from scipy.stats import pearsonr
import RNA
import parameters as param
from functions.general_functions import decimal_to_scientific_notation


L = int(sys.argv[1])
###############################################################################################
print 'load data' 
###############################################################################################
type_m_vs_no_pos = {'substitutions': L, 'deletions': L, 'insertions': L+1}
df_psample = pd.read_csv('./GPmapdata/Psample_seq_struct_data'+param.L_vsstring_parametersPsample[L]+'.csv')
###############################################################################################
print 'phenotype robustness' 
###############################################################################################
df_psample = pd.read_csv('./GPmapdata/Psample_seq_struct_data'+param.L_vsstring_parametersPsample[L]+'.csv')
shape_list_sample = df_psample['shape'].tolist()
seq_sample = df_psample['sequence'].tolist()
rhoS_list, rhoI_list, rhoD_list = df_psample['robustness substitutions'].tolist(), df_psample['robustness insertions'].tolist(), df_psample['robustness deletions'].tolist()
shape_vs_seq_list = {shape: [seq_sample[i] for i, s in enumerate(shape_list_sample) if s == shape] for shape in set(shape_list_sample) if shape.count('[')  > 0}
typem_vs_shape_vs_rho = {'substitutions': {shape: np.mean([rhoS_list[i] for i, s in enumerate(shape_list_sample) if s == shape]) for shape in set(shape_list_sample) if shape.count('[')  > 0},
                          'insertions': {shape: np.mean([rhoI_list[i] for i, s in enumerate(shape_list_sample) if s == shape]) for shape in set(shape_list_sample) if shape.count('[')  > 0},
                          'deletions': {shape: np.mean([rhoD_list[i] for i, s in enumerate(shape_list_sample) if s == shape]) for shape in set(shape_list_sample) if shape.count('[')  > 0}}
shapes_sorted = [s for s in sorted(typem_vs_shape_vs_rho['substitutions'].keys(), key=typem_vs_shape_vs_rho['substitutions'].get, reverse=True)]

###############################################################################################
print 'shape sample for subplot' 
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
print 'read out p-sample data' 
###############################################################################################
rho_sample = df_psample['robustness substitutions'].tolist()
rho_insertions_sample = df_psample['robustness insertions'].tolist()
rho_deletions_sample = df_psample['robustness deletions'].tolist()

###############################################################################################
###############################################################################################
print 'p-robustness' 
###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
print 'p-robustness with frequency changes plot' 
###############################################################################################
###############################################################################################
pheno_rhoS_list = [typem_vs_shape_vs_rho['substitutions'][shape] for shape in shapes_sorted]
pheno_rhoI_list = [typem_vs_shape_vs_rho['insertions'][shape]  for shape in shapes_sorted]
pheno_rhoD_list = [typem_vs_shape_vs_rho['deletions'][shape]  for shape in shapes_sorted]
f, ax = plt.subplots(ncols=2, figsize=(7, 2.7))
ax[0].scatter(pheno_rhoS_list, pheno_rhoI_list, c='b', s=8, label='ins.', alpha=0.7, lw=0)
ax[0].scatter(pheno_rhoS_list, pheno_rhoD_list, c='r', s=8, label='del.', alpha=0.7, lw=0)
ax[0].plot([0, 1], [0, 1], c='grey', alpha=0.7)
ax[0].legend(loc='lower right', fontsize='small')
ax[0].set_xlabel('phenotype robustness\nsubstitutions')
ax[0].set_ylabel('phenotype robustness\ninsertions/\ndeletions')
ax[0].set_ylim(-0.03, 1)
ax[0].set_xlim(-0.03, 1)
print 'number of shapes in phenotype robustness plot:', len(pheno_rhoS_list)
assert 0 <= min(pheno_rhoS_list) <= 1
assert 0 <= min(pheno_rhoI_list) <= 1
assert 0 <= min(pheno_rhoD_list) <= 1
assert 0 <= max(pheno_rhoS_list) <= 1
assert 0 <= max(pheno_rhoI_list) <= 1
assert 0 <= max(pheno_rhoD_list) <= 1
L_list = [L-1, L, L+1]
L_vs_shape_vs_freq = {}
for L_Gsample in L_list:
   df_sets = pd.read_csv('./GPmapdata/gsample_neutralsetdata'+param.L_vsstring_parametersGsample[L_Gsample]+'.csv')
   structure_list = df_sets['structure'].tolist()
   N_list = df_sets['number of times in sample'].tolist()
   L_vs_shape_vs_freq[L_Gsample] = {shape: N_list[i]/float(param.sample_sizeGsample) if  N_list[i] >= param.minimum_number_found else np.nan for i, shape in enumerate(structure_list)}
   if L_Gsample == L:
      print 'shapes in G-sample, but not in P-sample:', {shape: N_list[i]/float(param.sample_sizeGsample) for i, shape in enumerate(structure_list) if shape not in shape_list_sample}
pheno_freq_list_Lminusone = [L_vs_shape_vs_freq[L_list[0]][shape] if shape in L_vs_shape_vs_freq[L_list[0]] else np.nan for shape in shapes_sorted]
pheno_freq_list_Lplusone = [L_vs_shape_vs_freq[L_list[2]][shape] if shape in L_vs_shape_vs_freq[L_list[2]] else np.nan for shape in shapes_sorted]
ax[1].scatter(np.true_divide(pheno_freq_list_Lminusone, pheno_freq_list_Lplusone), np.true_divide(pheno_rhoD_list, pheno_rhoI_list), s=12, c='dimgrey', alpha=0.7, lw=0)
ax[1].set_xlabel('phenotype frequency ratio:\n'+r'$L=$'+str(L_list[0])+r' compared to $L=$'+str(L_list[2]))
ax[1].set_ylabel('phenotype robustness ratio:\ndeletions vs insertions')
x1, x2 = ax[1].get_xlim()
y1, y2 = ax[1].get_ylim()
y2 = max(y2, 1.1)
x2 = max(x2, 1.1)
ax[1].plot([x1, x2], [x1, x2], c='k', alpha=0.5)
#ax[1].plot([x1, x2], [1, 1], c='k', alpha=0.5)
#ax[1].plot([1, 1], [y1, y2], c='k', alpha=0.5)
ax[1].set_xlim(x1, x2)
ax[1].set_ylim(y1, y2)
for i in range(2):
   ax[i].annotate('ABCDEF'[i], xy=(0.04, 0.89), xycoords='axes fraction', fontsize='large', fontweight='bold', size=14)    
f.tight_layout()
f.savefig('./Pplots/Prho_robustness_ratio'+param.L_vsstring_parametersGsample[L]+'.png', bbox_inches='tight', dpi=300)
plt.close('all')
print 'phenotype robustness correlation: substitutions/insertions', pearsonr(pheno_rhoS_list, pheno_rhoI_list)
print 'phenotype robustness correlation: substitutions/deletions', pearsonr(pheno_rhoS_list, pheno_rhoD_list)
###############################################################################################
###############################################################################################
print 'frequency changes plot' 
###############################################################################################
###############################################################################################
f, ax = plt.subplots(ncols=3, figsize=(10, 3))
pheno_freq_list_Lminusone = [L_vs_shape_vs_freq[L_list[0]][shape] if shape in L_vs_shape_vs_freq[L_list[0]] else np.nan for shape in shapes_sorted]
pheno_freq_list_L = [L_vs_shape_vs_freq[L_list[1]][shape] if shape in L_vs_shape_vs_freq[L_list[1]] else np.nan for shape in shapes_sorted]
pheno_freq_list_Lplusone = [L_vs_shape_vs_freq[L_list[2]][shape] if shape in L_vs_shape_vs_freq[L_list[2]] else np.nan for shape in shapes_sorted]
for label, freq_list in [(str(L-1), pheno_freq_list_Lminusone), (str(L), pheno_freq_list_L), (str(L+1), pheno_freq_list_Lplusone)]:
   print 'number of phenotype values for L=' + label, ':', len([x for x in freq_list if not (np.isnan(x) or np.isinf(x))])
x_list, y_list = np.true_divide(pheno_freq_list_Lminusone, pheno_freq_list_L), np.true_divide(pheno_rhoD_list, pheno_rhoS_list)
ax[0].scatter(x_list, y_list, s=5, c='grey')
ax[0].set_xlabel('phenotype frequency ratio:\n'+r'$L=$'+str(L_list[0])+r' compared to $L=$'+str(L_list[1]))
ax[0].set_ylabel('phenotype robustness ratio:\ndeletions vs substitutions')
data_for_corr = [xy for xy in zip(x_list, y_list) if not (np.isnan(xy[0]) or np.isnan(xy[1]) or np.isinf(max(np.abs(xy))))]
print 'number of data points in phenotypic frequency change plot', len(data_for_corr)
#corr_coeff, pvalue = pearsonr(zip(*data_for_corr)[0], zip(*data_for_corr)[1])
#ax[0].set_title(r'correlation: '+str(round(corr_coeff,2)) + ', p-value='+ decimal_to_scientific_notation(pvalue), fontdict={'fontsize': 'small'})
print 'A)', 'phenotype frequency ratio:\n'+r'$L=$'+str(L_list[1])+r' compared to $L=$'+str(L_list[0]), 'vs', 'phenotype robustness ratio:\nsubstitutions vs deletions'
x_list, y_list = np.true_divide(pheno_freq_list_Lminusone, pheno_freq_list_Lplusone), np.true_divide(pheno_rhoD_list, pheno_rhoI_list)
ax[1].scatter(x_list, y_list, s=5, c='grey')
ax[1].set_xlabel('phenotype frequency ratio:\n'+r'$L=$'+str(L_list[0])+r' compared to $L=$'+str(L_list[2]))
ax[1].set_ylabel('phenotype robustness ratio:\ndeletions vs insertions')
data_for_corr = [xy for xy in zip(x_list, y_list) if not (np.isnan(xy[0]) or np.isnan(xy[1]) or np.isinf(max(np.abs(xy))))]
#corr_coeff, pvalue = pearsonr(zip(*data_for_corr)[0], zip(*data_for_corr)[1])
#ax[1].set_title(r'correlation: '+str(round(corr_coeff,2)) + ', p-value='+ decimal_to_scientific_notation(pvalue), fontdict={'fontsize': 'small'})
print 'B)', 'phenotype frequency ratio:\n'+r'$L=$'+str(L_list[2])+r' compared to $L=$'+str(L_list[0]), 'vs', 'phenotype robustness ratio:\ninsertions vs deletions'
x_list, y_list = np.true_divide(pheno_freq_list_L, pheno_freq_list_Lplusone), np.true_divide(pheno_rhoS_list, pheno_rhoI_list)
ax[2].scatter(x_list, y_list, s=5, c='grey')
ax[2].set_xlabel('phenotype frequency ratio:\n'+r'$L=$'+str(L_list[1])+r' compared to $L=$'+str(L_list[2]))
ax[2].set_ylabel('phenotype robustness ratio:\nsubstitutions vs insertions')
data_for_corr = [xy for xy in zip(x_list, y_list) if not (np.isnan(xy[0]) or np.isnan(xy[1]) or np.isinf(max(np.abs(xy))))]
#corr_coeff, pvalue = pearsonr(zip(*data_for_corr)[0], zip(*data_for_corr)[1])
#ax[2].set_title(r'correlation: '+str(round(corr_coeff,2)) + ', p-value='+ decimal_to_scientific_notation(pvalue), fontdict={'fontsize': 'small'})
print 'C)', 'phenotype frequency ratio:\n'+r'$L=$'+str(L_list[2])+r' compared to $L=$'+str(L_list[1]), 'vs', 'phenotype robustness ratio:\ninsertions vs substitutions'
for i in range(3):
   x1, x2 = ax[i].get_xlim()
   y1, y2 = ax[i].get_ylim()
   #ax[i].plot([x1, x2], [1, 1], c='k', alpha=0.5)
   #ax[i].plot([1, 1], [y1, y2], c='k', alpha=0.5)
   ax[i].plot([x1, x2], [x1, x2], c='k', alpha=0.5)
for i in range(3):
   ax[i].annotate('ABCDEF'[i], xy=(0.05, 0.84), xycoords='axes fraction', fontsize='large', fontweight='bold', size=14)    
f.tight_layout()
f.savefig('./Pplots/freq_change_robustness_ratio'+str(L)+'_'+param.L_vsstring_parametersGsample[L_Gsample]+'.png', bbox_inches='tight')
plt.close('all')
df_probustness = pd.DataFrame.from_dict({'shape': shapes_sorted,
                                          'phenotype robustness substitutions': pheno_rhoS_list,
                                          'phenotype robustness insertions': pheno_rhoI_list,
                                          'phenotype robustness deletions': pheno_rhoD_list,
                                          'phenotypic frequency L=' + str(L-1): pheno_freq_list_Lminusone,
                                          'phenotypic frequency L=' + str(L): pheno_freq_list_L,
                                          'phenotypic frequency L=' + str(L+1): pheno_freq_list_Lplusone})
df_probustness.to_csv('./data_plots/Prhofreq_change_robustness_ratio'+param.L_vsstring_parametersPsample[L]+'.csv')
###############################################################################################
###############################################################################################
print 'g-robustness' 
###############################################################################################
###############################################################################################
###############################################################################################
print 'plot data: g-robustness per phenotype' 
###############################################################################################
f, ax = plt.subplots(nrows=nrows, ncols=5, figsize=(15, 2.25 * len(shapes_sorted_forsubplots)//5 + 2.25))
for plotindex, shape in enumerate(shapes_sorted_forsubplots):
   geno_rho_listS = [rho_sample[i] for i, s in enumerate(shape_list_sample) if s == shape]
   if len(shapes_sorted_forsubplots)//5+1 > 1:
      axindex = tuple((plotindex//5, plotindex%5))
   else:
      axindex = plotindex
   for type_m in ['insertions', 'deletions']:      
      if type_m == 'insertions':
         geno_rho_list = [rho_insertions_sample[i] for i, s in enumerate(shape_list_sample) if s == shape]
      elif type_m == 'deletions':
         geno_rho_list = [rho_deletions_sample[i] for i, s in enumerate(shape_list_sample) if s == shape] 
      if plotindex == 0:
          ax[axindex].scatter(geno_rho_listS, geno_rho_list, label=type_m, c=param.type_m_vs_color[type_m], s=5, alpha=0.7, lw=0)
      else:
          ax[axindex].scatter(geno_rho_listS, geno_rho_list, c=param.type_m_vs_color[type_m], s=5, alpha=0.7, lw=0)
   if plotindex%5 == 0:
      ax[axindex].set_ylabel(r'$\widetilde{\rho}^{(g)}$'+'\ninsertions/deletions')
   ax[axindex].set_xlabel(r'$\widetilde{\rho}^{(g)}$'+'\nsubstitutions')
   ax[axindex].set_title('shape: ' + shape)
   del geno_rho_list, geno_rho_listS
   ax[axindex].plot([0,1], [0,1], c='k', zorder=-2)
   ax[axindex].set_xlim(-0.02, 1.02)
   ax[axindex].set_ylim(-0.02, 1.02)
f.tight_layout()
f.savefig('./Pplots/Psampling_Grobustness'+param.L_vsstring_parametersPsample[L]+'.png', bbox_inches='tight', dpi=300)
plt.close('all')
del f, ax
###############################################################################################
print 'get data: g-robustness correlations' 
###############################################################################################
zero_onlyD = 0
variance_rhoS, variance_rhoD, variance_rhoI = [], [], [] 
pearsonrPvalue_list_SI, pearsonrPvalue_list_SD = [], []
for shape in shapes_sorted:
   geno_rhoS_list = [rho_sample[i] for i, s in enumerate(shape_list_sample) if s == shape]
   geno_rhoI_list = [rho_insertions_sample[i] for i, s in enumerate(shape_list_sample) if s == shape]
   geno_rhoD_list = [rho_deletions_sample[i] for i, s in enumerate(shape_list_sample) if s == shape] 
   variance_rhoS.append(np.std(geno_rhoS_list)) 
   variance_rhoD.append(np.std(geno_rhoD_list)) 
   variance_rhoI.append(np.std(geno_rhoI_list))
   if max(geno_rhoS_list) - min(geno_rhoS_list) > 10**(-6) and max(geno_rhoI_list) - min(geno_rhoI_list) > 10**(-6):
      pearsonrPvalue_list_SI.append(tuple(pearsonr(geno_rhoS_list, geno_rhoI_list)))
   else:
      pearsonrPvalue_list_SI.append(tuple([np.nan, np.nan]))
   if max(geno_rhoS_list) - min(geno_rhoS_list) > 10**(-6) and max(geno_rhoD_list) - min(geno_rhoD_list) > 10**(-6):
      pearsonrPvalue_list_SD.append(tuple(pearsonr(geno_rhoS_list, geno_rhoD_list)))
   else:
      pearsonrPvalue_list_SD.append(tuple([np.nan, np.nan]))
   if max(geno_rhoD_list) < 0.2/float(L):
      zero_onlyD += 1 
      assert np.isnan(pearsonrPvalue_list_SD[-1][0])
df_corr_grobustess = pd.DataFrame.from_dict({'shape': shapes_sorted,
                                              'Pearson correlation: genotype robustness substitutions vs insertions': zip(*pearsonrPvalue_list_SI)[0],
                                              'p-value: genotype robustness substitutions vs insertions': zip(*pearsonrPvalue_list_SI)[1],
                                              'Pearson correlation: genotype robustness substitutions vs deletions': zip(*pearsonrPvalue_list_SD)[0],
                                              'p-value: genotype robustness substitutions vs deletions': zip(*pearsonrPvalue_list_SD)[1]})
for columnname in ['p-value: genotype robustness substitutions vs insertions',
                   'p-value: genotype robustness substitutions vs deletions']:
   number_p_vales_significant = len([p for p in df_corr_grobustess[columnname].tolist() if p < 0.05 and not np.isnan(p)])
   number_p_vales_total = len([p for p in df_corr_grobustess[columnname].tolist()])
   print columnname
   print 'fraction p-values < 0.05: ', float(number_p_vales_significant)/number_p_vales_total
   print 'or', number_p_vales_significant, 'out of', number_p_vales_total, '\n\n'
df_corr_grobustess.to_csv('./data_plots/genotype_robustness_correlations' + param.L_vsstring_parametersPsample[L] + '.csv')
print 'number of neutral sets with zero robustness to deletions', zero_onlyD
print 'min correlation SI', np.nanmin(zip(*pearsonrPvalue_list_SI)[0])
print 'min correlation SD', np.nanmin(zip(*pearsonrPvalue_list_SD)[0])
assert len(shapes_sorted) == 227
###############################################################################################
print 'correlations and variation in rho' 
###############################################################################################
f, ax = plt.subplots(ncols=2, figsize=(20,4))
for axindex, corr_list in enumerate([pearsonrPvalue_list_SI, pearsonrPvalue_list_SD]):
   ax[axindex].scatter(variance_rhoS, zip(*corr_list)[0], c='g', s=4, alpha=0.3)
   ax[axindex].scatter(variance_rhoI, zip(*corr_list)[0], c='b', s=4, alpha=0.3)
   ax[axindex].scatter(variance_rhoD, zip(*corr_list)[0], c='r', s=4, alpha=0.3)
   ax[axindex].set_xlabel('standard deviation: genotypic robustness distribution')
for axindex, ylabel in enumerate([r'correlation: $\widetilde{\rho}^{(g)}$ substitutions- insertions', r'correlation: $\widetilde{\rho}^{(g)}$ substitutions- deletions']):
   ax[axindex].set_ylabel(ylabel)
f.tight_layout()
f.savefig('./Pplots/grobustness_variations_and_correlations'+param.L_vsstring_parametersPsample[L]+'.png', bbox_inches='tight')
plt.close('all')
del f, ax   
###############################################################################################
print 'plot data: g-robustness per phenotype: one example and correlations' 
###############################################################################################
#indiced_shapes = list(set([0, 1, 2, 3, 4, 5] + [int(i) for i in np.arange(0, len(shapes_sorted)-1, step=len(shapes_sorted)//5)]))
for shapeindex in sorted(shape_indices_forsubplots):
   shape = shapes_sorted[shapeindex]
   f, ax = plt.subplots(ncols=3, figsize=(8.7, 3), gridspec_kw={'width_ratios': [1, 1, 0.3]})
   typem_vs_geno_rho_list = {'substitutions': [rho_sample[i] for i, s in enumerate(shape_list_sample) if s == shape]}
   typem_vs_geno_rho_list['insertions'] = [rho_insertions_sample[i] for i, s in enumerate(shape_list_sample) if s == shape]
   typem_vs_geno_rho_list['deletions'] = [rho_deletions_sample[i] for i, s in enumerate(shape_list_sample) if s == shape] 
   sc1 = ax[0].scatter(typem_vs_geno_rho_list['substitutions'], typem_vs_geno_rho_list['insertions'], label='insertions', s=3, c='b')
   sc2 = ax[0].scatter(typem_vs_geno_rho_list['substitutions'], typem_vs_geno_rho_list['deletions'], label='deletions', s=3, c='r')
   ax[0].set_ylabel(r'$\widetilde{\rho}^{(g)}$'+' insertions / deletions', fontsize='large')
   ax[0].set_xlim(-0.05, 1.02)
   ax[0].set_ylim(-0.05, 1.02)
   ax[0].set_title(shape)
   ax[0].plot([0, 1], [0, 1], c='grey', label='substitutions')
   ax[0].set_xlabel(r'$\widetilde{\rho}^{(g)}$'+' substitutions', fontsize='large')
   ax[0].annotate(shape, xy=(0.99, 0.05), xycoords='axes fraction', fontsize='medium', fontweight='bold', size=12, horizontalalignment='right')
   corr_coeffI, pvalueI = pearsonr(typem_vs_geno_rho_list['substitutions'], typem_vs_geno_rho_list['insertions'])
   corr_coeffD, pvalueD = pearsonr(typem_vs_geno_rho_list['substitutions'], typem_vs_geno_rho_list['deletions'])
   #ax[0].set_title(r'Insertions: '+str(round(corr_coeffI,2)) + ', p-value='+ decimal_to_scientific_notation(pvalueI) + '\n' + r'Deletions: '+str(round(corr_coeffD,2)) + ', p-value='+ decimal_to_scientific_notation(pvalueD),
   #                fontdict={'fontsize': 'small'})
   ax[0].set_title(r'correlation (insertions): '+str(round(corr_coeffI,2)) + '\n' + r'correlation (deletions): '+str(round(corr_coeffD,2)), fontdict={'fontsize': 'small'})
   hist1 = ax[1].hist([P for P in zip(*pearsonrPvalue_list_SI)[0] if not np.isnan(P)], color='b', alpha=0.5, label='insertions')
   hist2 = ax[1].hist([P for P in zip(*pearsonrPvalue_list_SD)[0] if not np.isnan(P)], color='r', alpha=0.5, label='deletions')
   ax[1].set_xlabel(r'Pearson correlation coefficient' + '\n'+r'$\widetilde{\rho}^{(g)}$ subst. vs ins./del.')
   ax[1].set_ylabel('number of neutral sets')
   minP = min([P for P in zip(*pearsonrPvalue_list_SI)[0] if not np.isnan(P)] + [P for P in zip(*pearsonrPvalue_list_SD)[0] if not np.isnan(P)])
   assert ax[1].get_xlim()[0] <= minP
   ax[1].set_xlim(min(0, ax[1].get_xlim()[0]), 1)
   ax[2].axis('off')
   ax[2].legend(handles=hist1[2] + hist2[2] + [sc1, sc2], loc='center')
   for axindex, axi in enumerate(ax[:-1]):
      axi.annotate('ABCDEF'[axindex], xy=(0.06, 0.86), xycoords='axes fraction', fontsize='large', fontweight='bold', size=15)
   f.tight_layout()
   f.savefig('./Pplots/Psampling_Grobustness'+param.L_vsstring_parametersPsample[L]+'example'+str(shapeindex)+'.png', bbox_inches='tight', dpi=200)
   plt.close('all')
   del f, ax


