import numpy as np
from scipy.stats import pearsonr, spearmanr
from copy import deepcopy




def read_in_shape_data_and_find_phiqplocal_correlation_list(shape, min_nonneutral_neighbours, min_nonneutral_Boltz_sum, shape_list_sample, shapes_list_full, sequence_list_sample, phi_listS, phi_listI, phi_listD, type_correlation, Boltzmann_prob_function, range_kbT):
   """perform analysis for phipq data"""
   g_vs_localphi_pqS, g_vs_localphi_pqI, g_vs_localphi_pqD, g_vs_BoltzmannP = {}, {}, {}, {}
   for i, s in enumerate(shape_list_sample): # get local phipq for each genotype in neutral network
      if s == shape:
         seq = sequence_list_sample[i][:]
         list_neighboursS = phi_listS[i].split('o')
         g_vs_localphi_pqS[seq] = [list_neighboursS.count(s2) for s2 in shapes_list_full if s2 != shape]
         list_neighboursD = phi_listD[i].split('o')
         g_vs_localphi_pqD[seq] = [list_neighboursD.count(s2) for s2 in shapes_list_full if s2 != shape]
         list_neighboursI = phi_listI[i].split('o')
         g_vs_localphi_pqI[seq] = [list_neighboursI.count(s2) for s2 in shapes_list_full if s2 != shape]
         shape_vsP = Boltzmann_prob_function(seq, range_kbT=2*range_kbT)
         g_vs_BoltzmannP[seq] = [shape_vsP[s2] if s2 in shape_vsP else 0 for s2 in shapes_list_full if s2 != shape]
         del list_neighboursS, list_neighboursD, list_neighboursI, shape_vsP
         if min([sum(g_vs_localphi_pqS[seq]), sum(g_vs_localphi_pqD[seq]), sum(g_vs_localphi_pqI[seq])]) < min_nonneutral_neighbours  or sum(g_vs_BoltzmannP[seq]) < min_nonneutral_Boltz_sum:
            del g_vs_localphi_pqS[seq] # do not consider sequences with no alternative phenotypes in neighbourhood
            del g_vs_localphi_pqI[seq] # do not consider sequences with no alternative phenotypes in neighbourhood
            del g_vs_localphi_pqD[seq] # do not consider sequences with no alternative phenotypes in neighbourhood
            del g_vs_BoltzmannP[seq] # do not consider sequences with no alternative phenotypes in neighbourhood   
         elif sum(g_vs_BoltzmannP[seq]) < 10.0**(-5):
            print 'no non-neutral subopt shapes', s, seq
   return [deepcopy(c[:]) for c in find_correlation_list_phiqplocal(type_correlation, g_vs_localphi_pqS, g_vs_localphi_pqI, g_vs_localphi_pqD, g_vs_BoltzmannP)]


def find_correlation_list_phiqplocal(type_correlation, g_vs_localphi_pqS, g_vs_localphi_pqI, g_vs_localphi_pqD, g_vs_BoltzmannP):
   """measure similarity of sequence-specific phi_pq"""
   seq_list = g_vs_localphi_pqS.keys()
   L = len(seq_list[0])
   normS, normI, normD = 1/float(3*L), 1/float(4*(L+1)), 1/float(L)
   seq1_to_seq2 = {seq1: np.random.choice([seq for seq in seq_list if seq != seq1]) for seq1 in seq_list}
   similaritySI_same_seq = [evaluate_type_correlation(g_vs_localphi_pqS[seq], g_vs_localphi_pqI[seq], type_correlation, normlist1=normS, normlist2=normI) for seq in seq_list]
   similaritySD_same_seq = [evaluate_type_correlation(g_vs_localphi_pqS[seq], g_vs_localphi_pqD[seq], type_correlation, normlist1=normS, normlist2=normD) for seq in seq_list]
   similaritySI_diff_seq = [evaluate_type_correlation(g_vs_localphi_pqS[seq1], g_vs_localphi_pqI[seq2], type_correlation, normlist1=normS, normlist2=normI) for seq1, seq2 in seq1_to_seq2.iteritems()]
   similaritySD_diff_seq = [evaluate_type_correlation(g_vs_localphi_pqS[seq1], g_vs_localphi_pqD[seq2], type_correlation, normlist1=normS, normlist2=normD)  for seq1, seq2 in seq1_to_seq2.iteritems()]
   # plastogenetic congruence
   similarityS_same_seq = [evaluate_type_correlation(g_vs_BoltzmannP[seq], g_vs_localphi_pqS[seq], type_correlation, normlist2=normS) for seq in seq_list]
   similarityI_same_seq = [evaluate_type_correlation(g_vs_BoltzmannP[seq], g_vs_localphi_pqI[seq], type_correlation, normlist2=normI) for seq in seq_list]
   similarityD_same_seq = [evaluate_type_correlation(g_vs_BoltzmannP[seq], g_vs_localphi_pqD[seq], type_correlation, normlist2=normD) for seq in seq_list]
   similarityS_diff_seq = [evaluate_type_correlation(g_vs_BoltzmannP[seq1], g_vs_localphi_pqS[seq2], type_correlation, normlist2=normS) for seq1, seq2 in seq1_to_seq2.iteritems()]
   similarityI_diff_seq = [evaluate_type_correlation(g_vs_BoltzmannP[seq1], g_vs_localphi_pqI[seq2], type_correlation, normlist2=normI) for seq1, seq2 in seq1_to_seq2.iteritems()]
   similarityD_diff_seq = [evaluate_type_correlation(g_vs_BoltzmannP[seq1], g_vs_localphi_pqD[seq2], type_correlation, normlist2=normD) for seq1, seq2 in seq1_to_seq2.iteritems()]
   return similaritySI_same_seq, similaritySD_same_seq, similaritySI_diff_seq, similaritySD_diff_seq, similarityS_same_seq, similarityI_same_seq, similarityD_same_seq, similarityS_diff_seq, similarityI_diff_seq, similarityD_diff_seq, seq_list


def evaluate_type_correlation(list1, list2, type_correlation, zero_threshold=10.0**(-5), normlist1=1, normlist2=1):
   assert len(list1) == len(list2)
   if type_correlation == 'correlation':
      if len(list1) < 3 or max(list1) - min(list1) == 0 or max(list2) - min(list2) == 0:
         print 'no correlation found for', list1, list2
         return np.nan
      return pearsonr(list1, list2)[0]
   elif type_correlation.startswith ('Jaccard index'):
      min_freq = float(type_correlation.split('>')[1][:-1])
      overlapping = len([j for j in range(len(list1)) if list1[j] * normlist1 > min_freq and list2[j] * normlist2 > min_freq])
      total = len([j for j in range(len(list1)) if list1[j] * normlist1 > min_freq or list2[j] * normlist2 > min_freq])
      if total > 0:
         return overlapping/float(total)
      else:
         assert overlapping == 0 and total == 0
         return 1.0
   elif type_correlation ==  'Bhattacharyya coeff. with no normalisation':
      list1_renormalised = [e * normlist1 for e in list1]
      list2_renormalised = [e * normlist2 for e in list2]
      return sum([(x * y)**0.5 for x, y in zip(list1_renormalised, list2_renormalised)])
   elif type_correlation ==  'Bhattacharyya coeff.':
      # following the definition in Greenbury et al. (2016)
      if max(list1) == 0 or max(list2) == 0:
         return np.nan
      norm_to_prob1, norm_to_prob2 = float(sum(list1)), float(sum(list2))
      list1_renormalised = [e/norm_to_prob1 for e in list1]
      list2_renormalised = [e/norm_to_prob2 for e in list2]
      return sum([(x * y)**0.5 for x, y in zip(list1_renormalised, list2_renormalised)])
   else:
      raise RuntimeError('type_correlation not defined: '+type_correlation)

def SpearmanR_nonzero_values(x_list, y_list):
   return spearmanr([x for i, x in enumerate(x_list) if x > 0 and y_list[i] > 0], [y for i, y in enumerate(y_list) if y > 0 and x_list[i] > 0])

def logPearsonR_nonzero_values(x_list, y_list):
   return pearsonr(np.log([x for i, x in enumerate(x_list) if x > 0 and y_list[i] > 0]), np.log([y for i, y in enumerate(y_list) if y > 0 and x_list[i] > 0]))

