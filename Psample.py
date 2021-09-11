import numpy as np
import functions.Psampling_shapes as inv
from functools import partial
import sys
from functions.rna_structural_functions import generate_all_allowed_dotbracket
import pandas as pd
from multiprocessing import Pool
from os.path import isfile
import parameters as param
from collections import Counter




  

Lmax_fullenumeration = 35
n_cpus = int(sys.argv[1])
L = int(sys.argv[2])
###################################################################################################
assert L <= Lmax_fullenumeration #prevent memory overflow in structure enumeration
if not isfile('./GPmapdata/Psample_seq_struct_data'+param.L_vsstring_parametersPsample[L]+'.csv'):
   ###############################################################################################
   print 'get shape sample' 
   ###############################################################################################
   filename_shape_list = './GPmapdata/allshapesL'+str(L)+'.txt'
   filename_fullshape_list = './GPmapdata/shapes_for_all_structuresL'+str(L)+'.txt'
   filename_structure_list = './GPmapdata/allstructuresL'+str(L)+'.txt'
   if isfile(filename_shape_list) and isfile(filename_fullshape_list) and isfile(filename_structure_list):
      shape_sample = inv.read_structure_list(filename_shape_list)
      full_shapes_list = inv.read_structure_list(filename_fullshape_list)
      full_secondary_str_list = inv.read_structure_list(filename_structure_list)
   else:
      full_secondary_str_list = generate_all_allowed_dotbracket(L, allow_isolated_bps=param.allow_isolated_bps)
      full_shapes_list = [param.cg_function(s) for s in full_secondary_str_list]
      shape_sample = sorted(list(set(full_shapes_list)), key=len)
      inv.save_structure_list(shape_sample, filename_shape_list)
      inv.save_structure_list(full_shapes_list, filename_fullshape_list)
      inv.save_structure_list(full_secondary_str_list, filename_structure_list)
   shape_vs_structure_list = {shape: [structure for structureindex, structure in enumerate(full_secondary_str_list) if full_shapes_list[structureindex] == shape] for shape in shape_sample}  
   ###############################################################################################
   print 'start sequences RNAinverse' 
   ###############################################################################################
   sequence_sample_list_inverse_file = './GPmapdata/sequence_sample_list_inverse'+param.L_vsstring_parametersPsample[L].split('sample')[0]+'.txt'
   inv_function_name = partial(inv.shape_inverse_RNAinverse, shape_vs_structure_list=shape_vs_structure_list, GPfunction=param.GPfunction, number_RNAinverse=param.number_RNAinverse)
   if not isfile(sequence_sample_list_inverse_file): # and L <= Lmax_fullenumeration:
      shape_list_with_duplicate_inverse_runs = [s for i in range(param.no_sitescans) for s in shape_sample]
      pool = Pool(processes = n_cpus)
      sequence_sample_list_inverse = pool.map(inv_function_name, shape_list_with_duplicate_inverse_runs)
      pool.close()
      pool.join()
      for i, seq in enumerate(sequence_sample_list_inverse):
         if seq:
            assert param.GPfunction(seq) == shape_list_with_duplicate_inverse_runs[i]
      inv.save_structure_list([s for s in sequence_sample_list_inverse if s], sequence_sample_list_inverse_file)
   else:
      sequence_sample_list_inverse = inv.read_structure_list(sequence_sample_list_inverse_file)
   ###############################################################################################
   print 'get sequence sample: sitescanning' 
   ###############################################################################################
   sequence_sample_list_filename = './GPmapdata/site_scanning_sample'+param.L_vsstring_parametersPsample[L]+'.txt'
   sitescanning_function_name = partial(inv.rw_sitescanning, length_per_walk=param.number_seq_per_str*L//param.no_sitescans, 
                                   subsample_size=param.number_seq_per_str//param.no_sitescans, GPfunction=param.GPfunction, shrep_function = param.shrep_function , string_parameters=param.L_vsstring_parametersPsample[L])
   if not isfile(sequence_sample_list_filename):
      pool = Pool(processes = n_cpus)
      sequence_sample_list_of_list = pool.map(sitescanning_function_name, [s for s in sequence_sample_list_inverse if s])
      pool.close()
      pool.join()
      for i, seq_list in enumerate(sequence_sample_list_of_list):
         for seq in seq_list:
            assert param.GPfunction(seq) == param.GPfunction([s for s in sequence_sample_list_inverse if s][i]) 
      sequence_sample_list = [s for l in sequence_sample_list_of_list for s in l]
      ###########################################################
      print 'keep only shapes with required number of sequences' 
      ###########################################################
      pool = Pool(processes = n_cpus)
      shape_list_sample = pool.map(param.GPfunction, sequence_sample_list)
      pool.close()
      pool.join()      
      shape_vs_number_seq = Counter(shape_list_sample)
      assert param.number_seq_per_str == max(shape_vs_number_seq.values())
      sequence_sample_list = [seq for seq, shape in zip(sequence_sample_list, shape_list_sample) if shape_vs_number_seq[shape] == param.number_seq_per_str]
      inv.save_structure_list(sequence_sample_list, sequence_sample_list_filename)
      del sequence_sample_list_of_list
   else:
      sequence_sample_list = inv.read_structure_list(sequence_sample_list_filename)
   ###############################################################################################
   print 'get data' 
   ###############################################################################################
   function_for_parallelised_calc = partial(inv.get_shape_robustness_from_sequence, GPfunction=param.GPfunction)
   # pool_result = [function_for_parallelised_calc(s) for s in sequence_sample_list]
   pool = Pool(processes = n_cpus)
   pool_result = pool.map(function_for_parallelised_calc, sequence_sample_list)
   pool.close()
   pool.join()
   shape_list_sample, rho_listS, rho_listI, rho_listD, phi_listS, phi_listI, phi_listD = zip(*pool_result)

   df_rho = pd.DataFrame.from_dict({'shape': shape_list_sample, 
                                'sequence': sequence_sample_list,
                                'robustness substitutions': rho_listS,
                                'robustness deletions': rho_listD,
                                'robustness insertions': rho_listI,
                                'phi substitutions': phi_listS,
                                'phi deletions': phi_listD,
                                'phi insertions': phi_listI})
   print 'number sequences found for each shape:',  {shape: shape_list_sample.count(shape) for shape in shape_sample}
   df_rho.to_csv('./GPmapdata/Psample_seq_struct_data'+param.L_vsstring_parametersPsample[L]+'.csv')
   del shape_list_sample, rho_listS, rho_listI, rho_listD, phi_listS, phi_listI, phi_listD

###############################################################################################
print 'quick tests on collected data' 
###############################################################################################
from functions.mutational_neighbourhood_functions import *
df_psample = pd.read_csv('./GPmapdata/Psample_seq_struct_data'+param.L_vsstring_parametersPsample[L]+'.csv')
shapes_for_test = np.random.choice(df_psample['shape'].tolist(), 4)
for shape in shapes_for_test:
   no_seq_analysed = 0
   for i, s in enumerate(df_psample['shape'].tolist()): # get local phipq for each genotype in neutral network
      if s == shape and no_seq_analysed < 15:
         no_seq_analysed += 1
         seq = df_psample['sequence'].tolist()[i]
         list_neighboursS = df_psample['phi substitutions'].astype('str').tolist()[i].split('o')
         list_neighboursD = df_psample['phi deletions'].astype('str').tolist()[i].split('o')
         list_neighboursI = df_psample['phi insertions'].astype('str').tolist()[i].split('o')
         list_neighboursS_test = [param.GPfunction(seq2) for seq2 in substitution_neighbours_letters(seq, L)]
         list_neighboursD_test = [param.GPfunction(seq2) for seq2 in neighbours_by_deletion_letters(seq, L)]
         list_neighboursI_test = [param.GPfunction(seq2) for seq2 in neighbours_by_insertion_letters(seq, L)]         
         for shape2 in list_neighboursS + list_neighboursS_test:
            if shape2 != shape and ('[' in shape2 or '_' in shape2 or '|' in shape2):
               assert list_neighboursS.count(shape2) == list_neighboursS_test.count(shape2)
         for shape2 in list_neighboursD + list_neighboursD_test:
            if shape2 != shape and ('[' in shape2 or '_' in shape2 or '|' in shape2):
               assert list_neighboursD.count(shape2) == list_neighboursD_test.count(shape2)
         for shape2 in list_neighboursI + list_neighboursI_test:
            if shape2 != shape and ('[' in shape2 or '_' in shape2 or '|' in shape2):
               assert list_neighboursI.count(shape2) == list_neighboursI_test.count(shape2)
         assert abs(df_psample['robustness substitutions'].tolist()[i] - list_neighboursS_test.count(shape)/float(len(list_neighboursS_test))) < 1/float(4 * L)
         assert abs(df_psample['robustness insertions'].tolist()[i] - list_neighboursI_test.count(shape)/float(len(list_neighboursI_test))) < 1/float(4 * L)
         assert abs(df_psample['robustness deletions'].tolist()[i] - list_neighboursD_test.count(shape)/float(len(list_neighboursD_test))) < 1/float(4 * L)
         print 'finished test for ', seq, shape
shape_count_sample = Counter(df_psample['shape'].tolist())
assert param.number_seq_per_str == min(shape_count_sample.values())
assert param.number_seq_per_str == max(shape_count_sample.values())
print 'number of shapes in sample:', len(set(df_psample['shape'].tolist()))