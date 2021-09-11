# the robustness and phi_{pq} definitions follow Wagner (2008) and Greenbury et al. (2016), as described in the main text of the paper
####################################################################################
# mutational neighbours of a genotype in letter notation
####################################################################################
base_list = ['A', 'U', 'C', 'G']

def neighbours_by_deletion_letters(g, L):  
   return [''.join([g[i] for i in range(L) if i!=position]) for position in range(L)]

def neighbours_by_insertion_letters(g, L):  
   return [insert_letters(position, new_base, g) for position in range(L+1) for new_base in base_list]

def insert_letters(position, new_base, g):
   return ''.join([b for i, b in enumerate(g) if i < position] + [new_base,] + [b for i, b in enumerate(g) if i >= position])

def substitution_neighbours_letters(g, L):  
   return [''.join([base if seqindex!=pos else new_base for seqindex, base in enumerate(g)]) for pos in range(L) for new_base in base_list if not g[pos] == new_base]

def substitution_neighbours_given_site_letters(g, L, pos):  
   return [''.join([base if seqindex!=pos else new_base for seqindex, base in enumerate(g)]) for new_base in base_list if not g[pos] == new_base]

####################################################################################
# phenotypic robustness
####################################################################################

def sequence_robustness_and_phipq(seq, type_m, GPfunction):
   """apply GPfunction to all sequences in the mutational neighbourhood and summarise
   neutral and non-neutral neighbours (genotype robustness and 
   concatenation of non-neutral shapes sererated by 'o' character"""
   assert type_m in ['substitutions', 'deletions', 'insertions']
   structure, L = GPfunction(seq), len(seq)
   if type_m == 'substitutions':
      g_list = substitution_neighbours_letters(seq, L)
   elif type_m == 'deletions':
      g_list = neighbours_by_deletion_letters(seq, L)
   elif type_m == 'insertions':
      g_list = neighbours_by_insertion_letters(seq, L)
   p_neighbours = [GPfunction(g) for g in g_list]
   rho = p_neighbours.count(structure)/float(len(g_list))
   phi = [s for s in p_neighbours if s != structure]
   return rho, 'o'+'o'.join(phi)


############################################################################################################
## test
############################################################################################################
if __name__ == "__main__":
   test_seq = 'ACGCUU'
   L = len(test_seq)
   neighboursS = substitution_neighbours_letters(test_seq, L)
   neighboursD = neighbours_by_deletion_letters(test_seq, L)
   neighboursI = neighbours_by_insertion_letters(test_seq, L)
   neighboursS_four = substitution_neighbours_given_site_letters(test_seq, L, 4)
   assert len(neighboursS) == 3 * L
   assert len(neighboursD) == L
   assert len(neighboursI) == 4 * (L + 1)
   assert len(neighboursS_four) == 3 
   assert 'ACGCAU' in neighboursS_four and 'ACGCCU' in neighboursS_four and 'ACGCGU' in neighboursS_four
   for seq2 in neighboursS:
      assert len(seq2) == L 
   for seq2 in neighboursI:
      assert len(seq2) == L + 1 
   for seq2 in neighboursD:
      assert len(seq2) == L - 1 
   assert len(set(neighboursS)) == len(neighboursS) 
   for i in range(L):
      for seq2 in substitution_neighbours_given_site_letters(test_seq, L, i):
         assert seq2 in neighboursS
   test_seq = 'ACUU'
   L = len(test_seq)
   neighboursD = neighbours_by_deletion_letters(test_seq, L)
   neighboursI = neighbours_by_insertion_letters(test_seq, L)
   assert neighboursD.count('ACU') == 2 and neighboursD.count('AUU') == 1 and neighboursD.count('CUU') == 1
   assert neighboursI.count('AACUU') == 2 and neighboursI.count('CACUU') == 1 and neighboursI.count('GACUU') == 1 and neighboursI.count('UACUU') == 1
   assert neighboursI.count('ACCUU') == 2 and neighboursI.count('AGCUU') == 1 and neighboursI.count('AUCUU') == 1
   assert neighboursI.count('ACAUU') == 1 and neighboursI.count('ACGUU') == 1 and neighboursI.count('ACUUU') == 3
   assert neighboursI.count('ACUAU') == 1 and neighboursI.count('ACUGU') == 1 and neighboursI.count('ACUCU') == 1
   assert neighboursI.count('ACUUA') == 1 and neighboursI.count('ACUUG') == 1 and neighboursI.count('ACUUC') == 1



      




