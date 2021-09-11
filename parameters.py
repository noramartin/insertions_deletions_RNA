import numpy as np
import functions.Vienna_cgshape_probs as shapes
from functools import partial
from functions.rna_structural_functions import dotbracket_to_coarsegrained_for_level
###############################################################################################
# 'general' 
###############################################################################################
L_list = list(range(10, 60))
shape_level, range_kbT = 2, 15
K = 4 #alphabet size (4 for RNA), following notation by Schaper and Louis 
dangling_ends_option = 2 # corresponds to overdangle grammar
folding_criterion = 1.1
allow_isolated_bps = False
param_general = '_' + str(int(10*folding_criterion)) + '_' + str(dangling_ends_option) + 'SL' + str(shape_level) + 'kbT' + str(range_kbT)
###############################################################################################
# 'G-sampling parameters for large G-sample (frequencies)' 
###############################################################################################
sample_sizeGsample = 10**7 #10**5
L_vsstring_parametersGsample = {L: 'L'+str(L)+'_gsample'+str(int(np.log10(sample_sizeGsample))) + param_general for L in L_list}
minimum_number_found = 10 
###############################################################################################
# 'P-sampling parameters' 
###############################################################################################
number_seq_per_str = 500
number_RNAinverse, no_sitescans = 5000, 10
L_vsstring_parametersPsample = {L: 'L' + str(L)  + param_general + '_'+str(no_sitescans)+'sample' + str(number_seq_per_str) for L in L_list}
assert number_seq_per_str % no_sitescans == 0 # lengths of independent random walks must be equal
###############################################################################################
# functions 
###############################################################################################
GPfunction = partial(shapes.find_most_freq_shape, shape_level=shape_level, range_kbT=range_kbT, folding_criterion=folding_criterion, allow_isolated_bps=allow_isolated_bps, dangling_ends_option=dangling_ends_option)
shrep_function = partial(shapes.get_shrep, shape_level=shape_level, range_kbT=range_kbT, allow_isolated_bps=allow_isolated_bps, dangling_ends_option=dangling_ends_option)
cg_function = partial(dotbracket_to_coarsegrained_for_level, shape_level=shape_level)
shape_and_prob_function = partial(shapes.get_shapes_prob_subopt, shape_level=shape_level, range_kbT=range_kbT, allow_isolated_bps=allow_isolated_bps, dangling_ends_option=dangling_ends_option)
###############################################################################################
# formatting 
###############################################################################################
type_m_vs_color = {'substitutions': 'lime', 'insertions': 'b', 'deletions': 'r'}
