import RNA
from os.path import isfile
from rna_structural_functions import get_basepair_indices_from_dotbracket
from mutational_neighbourhood_functions import sequence_robustness_and_phipq, substitution_neighbours_given_site_letters
import random # use random and not np.random because random uses a shared random state in multiprocessing

def save_structure_list(structure_sample, filename):
   """save list of structures to file, one structure per line"""
   with open(filename, 'w') as textfile:
      for structure in structure_sample:
         textfile.write(structure+'\n')

def read_structure_list(filename):
   """read list of structures from file"""
   with open(filename, 'r') as textfile:
      structure_sample = [str(line.strip()) for line in textfile.readlines()]
   return structure_sample

def get_shape_robustness_from_sequence(seq, GPfunction):
   """run robustness calculations for all types of mutations"""
   shape = GPfunction(seq)
   assert shape != '|'
   rho_s, phi_s = sequence_robustness_and_phipq(seq, type_m='substitutions', GPfunction=GPfunction)
   rho_d, phi_d = sequence_robustness_and_phipq(seq, type_m='deletions', GPfunction=GPfunction)
   rho_i, phi_i = sequence_robustness_and_phipq(seq, type_m='insertions', GPfunction=GPfunction)
   return shape, rho_s, rho_i, rho_d, phi_s, phi_i, phi_d

   
###############################################################################################
# RNAinverse
###############################################################################################
def run_RNAinverse(dotbracket_structure):
   """run inverse_fold for a random start sequence"""
   random_startseq = ''.join([random.choice(['A', 'U', 'C', 'G']) for c in range(len(dotbracket_structure))])
   (sequence_str, distance_from_target_str) = RNA.inverse_fold(random_startseq, dotbracket_structure) 
   return sequence_str, distance_from_target_str


def shape_inverse_RNAinverse(target_shape, shape_vs_structure_list, GPfunction, number_RNAinverse=10):
   """find a sequence folding into shape target_shape;
   need a list of full structures belonging to that shape - 
   one structure will be chosen at random and used for ViennaRNA's inverse fold;
   if the result of inverse fold does not meet the shape criteria, the process is repeated up to number_RNAinverse times;
   None is returned if no sequence found"""
   no_iterations = 0
   while no_iterations < number_RNAinverse:
      start_seq, dist_target = run_RNAinverse(random.choice(shape_vs_structure_list[target_shape]))
      if GPfunction(start_seq) == target_shape: #only accept sequences which are not perfect if no perfect sequences found repeatedly
         return start_seq
      no_iterations += 1
   return None

###############################################################################################
# site-scanning random walk
###############################################################################################
allowed_basepairs = ['AU', 'UA', 'GC', 'CG', 'GU', 'UG']

def bp_swaps_g_given_site(g_char, pos, paired_pos): 
   """perform all base pair swaps at the given position and its paired position"""
   g_list = []
   for new_bp in allowed_basepairs:
      if new_bp[0] != g_char[pos] or new_bp[1] != g_char[paired_pos]: # not identical bp 
         g_new = [x for x in g_char]
         g_new[pos] = new_bp[0]
         g_new[paired_pos] = new_bp[1]
         g_list.append(''.join(g_new))
   return g_list

def rw_sitescanning(startseq, length_per_walk, subsample_size, GPfunction, shrep_function, string_parameters=None):
   """perform a site-scanning random walk of length length_per_walk starting from startseq (tuple of ints) and subsample a sample of size subsample_size;
   site-scanning method adapted from Weiss and Ahnert (2020), Royal Society Interface - as described in Electronic Supplementary Material - IIA;
   here, this is applied with the RNAshapes GPfunction and base pair swaps are allowed at paired positions in the shrep"""
   if string_parameters:
      filename = './GPmapdata/site_scanning_sample'+string_parameters+'_'+startseq+'.txt' 
      if isfile(filename):
         print 'read site-scanning from file'
         return read_structure_list(filename)
   structure = GPfunction(startseq)
   seq_list_rw, current_seq = [], ''.join([g for g in startseq])
   L, K, site_to_scan = len(startseq), 4, 0
   seq_list_rw.append(current_seq)
   shrep = shrep_function(current_seq, structure)
   bp_dict = get_basepair_indices_from_dotbracket(shrep)
   print 'start site-scanning for '+structure+ '\n'
   while len(seq_list_rw) < length_per_walk:         
      if shrep[site_to_scan] != '.':
         neighbours_given_site = bp_swaps_g_given_site(current_seq, site_to_scan, bp_dict[site_to_scan])
      else:
         neighbours_given_site = substitution_neighbours_given_site_letters(current_seq, L, site_to_scan)
      for index in random.sample(xrange(len(neighbours_given_site)), k=len(neighbours_given_site)):
         g = ''.join([c for c in neighbours_given_site[index]])
         shape_g = GPfunction(g)
         if structure == shape_g:
            seq_list_rw.append(g)
            shrep = shrep_function(current_seq, shape_g)
            bp_dict = get_basepair_indices_from_dotbracket(shrep)
            current_seq = g
            break
      site_to_scan = (site_to_scan+1)%L  
      if site_to_scan == 0 and len(seq_list_rw) == 1: # no neutral neighbours
         print 'site-scanning unsuccessful for ', structure
         return []
      #print 'shape', structure, 'new seq', g
   assert len(seq_list_rw) == length_per_walk
   print 'finish site-scanning for ', structure
   seq_list_sitescanning = [seq_list_rw[random.randint(0, length_per_walk-1)] for i in range(subsample_size)]
   if string_parameters:
      save_structure_list(seq_list_sitescanning, filename)
   return seq_list_sitescanning

############################################################################################################
## test 
############################################################################################################
if __name__ == "__main__":
   print 'test reading/writing'
   l1 = ['abc', 'b', 'c', 'd']
   save_structure_list(l1, 'test.txt')
   l2 = read_structure_list('test.txt')
   assert len(l2) == len(l1)
   for i in range(len(l2)):
      assert l2[i] == l1[i]
   print 'base pair swaps'
   seq = 'AGAUUUUCU'
   structure = '(((...)))'
   bp_list = ['ACAUUUUGU', 'AAAUUUUUU', 'AUAUUUUAU', 'AUAUUUUGU', 'AGAUUUUUU']
   bp_list_test = bp_swaps_g_given_site(seq, 1, 7)
   assert len(bp_list_test) == len(bp_list) == len(set(bp_list_test)) == len(set(bp_list_test + bp_list))



