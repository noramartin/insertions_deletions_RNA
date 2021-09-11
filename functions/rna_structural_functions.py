import networkx as nx
import string

############################################################################################################
## manipulate structures
############################################################################################################

def generate_all_allowed_dotbracket(L, allow_isolated_bps=False):
   """recursively build up all possible structures of length L;
    - allow_isolated_bps (True/False) to control whether isolated base pairs are permitted"""
   assert L <= 35 # avoid overflow errors
   potential_structures = [s for s in return_all_allowed_structures_starting_with('', L, allow_isolated_bps=allow_isolated_bps) if is_likely_to_be_valid_structure(s, allow_isolated_bps=allow_isolated_bps)]
   return potential_structures

def return_all_allowed_structures_starting_with(structure, L, allow_isolated_bps=False):
   """recursively build up all possible structures of length L;
   already has some constraints of dot-bracket notations built in:
   only one closed bracket for each upstream open bracket, all brackets closed by the end, length of hairpin loops,
   optionally isolated base pairs"""
   if len(structure) == L:
      return [structure]
   assert len(structure) < L
   assert L-len(structure) <= 35
   structure_list = []
   if (len(structure) >= 2 and structure[-1] == '(' and structure[-2] != '(' and not allow_isolated_bps):
      structure_list += return_all_allowed_structures_starting_with(structure+'(', L, allow_isolated_bps) #need to add ( to avoid isolated bp
   elif (len(structure) == 1 and structure[-1] == '(' and not allow_isolated_bps):
      structure_list += return_all_allowed_structures_starting_with(structure+'(', L, allow_isolated_bps) #need to add ( to avoid isolated bp
   elif (len(structure) >= 2 and structure[-1] == ')' and structure[-2] != ')' and not allow_isolated_bps):
      structure_list += return_all_allowed_structures_starting_with(structure+')', L, allow_isolated_bps) #need to add ) to avoid isolated bp
   elif len(structure) <= 4: #at the beginning of a structure can only have opening brackets/loops
      structure_list += return_all_allowed_structures_starting_with(structure+'.', L, allow_isolated_bps) 
      structure_list += return_all_allowed_structures_starting_with(structure+'(', L, allow_isolated_bps) 
   elif structure.count('(') > structure.count(')')+L-len(structure):
      pass  #cannot close all base pairs anymore, return empty list
   elif structure.count('(') == structure.count(')')+L-len(structure) and '(' not in structure[-3:]: #need to close base pairs
      structure_list += return_all_allowed_structures_starting_with(structure+')', L, allow_isolated_bps)    
   elif structure.count('(') > structure.count(')') and allow_isolated_bps and '(' not in structure[-3:]: #upstream bps are open, so closing also allowed
      structure_list += return_all_allowed_structures_starting_with(structure+'.', L, allow_isolated_bps)
      if L-len(structure) > 3:
         structure_list += return_all_allowed_structures_starting_with(structure+'(', L, allow_isolated_bps)
      structure_list += return_all_allowed_structures_starting_with(structure+')', L, allow_isolated_bps)
   elif structure.count('(') > structure.count(')') and structure[-1] == ')' and not allow_isolated_bps: # everything allowed
      structure_list += return_all_allowed_structures_starting_with(structure+'.', L, allow_isolated_bps)
      if L-len(structure) > 3:
         structure_list += return_all_allowed_structures_starting_with(structure+'(', L, allow_isolated_bps)
      structure_list += return_all_allowed_structures_starting_with(structure+')', L, allow_isolated_bps)
   elif structure.count('(') > structure.count(')')+1 and not allow_isolated_bps and '(' not in structure[-3:]: #at least two upstream bps are open, so closing also allowed with non-isolated bps
      structure_list += return_all_allowed_structures_starting_with(structure+'.', L, allow_isolated_bps)
      if L-len(structure) > 3:
         structure_list += return_all_allowed_structures_starting_with(structure+'(', L, allow_isolated_bps)
      structure_list += return_all_allowed_structures_starting_with(structure+')', L, allow_isolated_bps)
   else: # can open new base pairs or introduce loop
      structure_list += return_all_allowed_structures_starting_with(structure+'.', L, allow_isolated_bps)
      if L-len(structure) > 3:
         structure_list += return_all_allowed_structures_starting_with(structure+'(', L, allow_isolated_bps)
   return structure_list

############################################################################################################
## coarse-grain structures
############################################################################################################

def cut_ends(dotbracketstring):
   """cut unpaired regions at either end of dotbacketstring"""
   while len(dotbracketstring) and dotbracketstring[0] in ['.', '_']:
      dotbracketstring = dotbracketstring[1:]
   while len(dotbracketstring) and dotbracketstring[-1] in ['.', '_']:
      dotbracketstring = dotbracketstring[:-1]
   if dotbracketstring:
      return dotbracketstring
   else:
      return '.' #for unfolded structure

def dotbracket_to_coarsegrained_for_level(dotbacketstring, shape_level):
   if shape_level == 1:
      return dotbracket_to_coarsegrained_lev1(dotbacketstring)
   elif shape_level == 2:
      return dotbracket_to_coarsegrained_lev2(dotbacketstring)
   elif shape_level == 5:
      return dotbracket_to_coarsegrained_lev5(dotbacketstring)
   else:
      raise ValueError('shape level not implemented:' +str(shape_level))


def dotbracket_to_coarsegrained_lev1(dotbracketstring):
   """transform a full dotbracket representation to a coarse-grained 
   (type 1 as defined in Janssen, Reeder, Giegerich (2008). BMC bioinformatics)"""
   fine_grained_to_coarse_grained_symbol = {'(': '[', ')': ']', '.': '_'}
   basepair_index_mapping = get_basepair_indices_from_dotbracket(dotbracketstring)
   coarse_grained_string_list = []
   for charindex, char in enumerate(dotbracketstring):
      if charindex == 0  or dotbracketstring[charindex-1] != dotbracketstring[charindex]:
         coarse_grained_string_list.append(fine_grained_to_coarse_grained_symbol[char])
      elif dotbracketstring[charindex-1] == dotbracketstring[charindex] and dotbracketstring[charindex] != '.': #two subsequent brackets of same type
         if not abs(basepair_index_mapping[charindex]-basepair_index_mapping[charindex-1])<1.5:
             coarse_grained_string_list.append(fine_grained_to_coarse_grained_symbol[char])
   return remove_hairpins_from_cg(coarse_grained_string_list)

def lev1_to_lev2(cg_string):
   """transform level-1 to level-2 coarse-grained shape,
    as defined in Janssen, Reeder, Giegerich (2008). BMC bioinformatics"""
   if len(cg_string) == 1:
      return cg_string
   tree_secondarystr = get_tree_bpmap_from_dotbracket(cg_string)[0]
   coarse_grained_string_list = []
   for charindex, char in enumerate(cg_string):
      if cg_string[charindex] != '_': 
         coarse_grained_string_list.append(char)
      elif not is_multiloop_from_tree(charindex, tree_secondarystr, cg_string) and not is_exterior_loop_cg(charindex, cg_string):
          coarse_grained_string_list.append(char) 
   return ''.join(coarse_grained_string_list)


def dotbracket_to_coarsegrained_lev2(dotbracketstring):
   """transform a full dotbracket representation to a coarse-grained 
   (type 2 as defined in Janssen, Reeder, Giegerich (2008). BMC bioinformatics)"""
   return lev1_to_lev2(dotbracket_to_coarsegrained_lev1(dotbracketstring))


def remove_hairpins_from_cg(cg_list):
   #cg_bp_map = get_basepair_indices_from_dotbracket(coarse_grained_string)
   #hairpin_indices = set([(i1+i2)//2 for i1, i2 in cg_bp_map.iteritems() if abs(i1-i2) == 2])
   hairpin_indices = [i-1 for i in range(2, len(cg_list)) if cg_list[i] == ']' and cg_list[i-1] == '_' and cg_list[i-2] == '[']
   return ''.join([c for i, c in enumerate(cg_list) if i not in hairpin_indices])

def dotbracket_to_coarsegrained_lev5(dotbracketstring):
   """transform a full dotbracket representation to a coarse-grained 
   (type 5 as defined in Janssen, Reeder, Giegerich (2008). BMC bioinformatics)"""
   square_bracket_to_round = {'[': '(', ']': ')'}
   simple_coarse_grained = dotbracket_to_coarsegrained_lev1(dotbracketstring)
   return lev1_or_2_to_lev5(simple_coarse_grained)

def lev1_or_2_to_lev5(simple_coarse_grained):
   """transform a full dotbracket representation to a coarse-grained 
   (type 5 as defined in Janssen, Reeder, Giegerich (2008). BMC bioinformatics)"""
   square_bracket_to_round = {'[': '(', ']': ')'}
   simple_coarse_grained_noloops = ''.join([square_bracket_to_round[c] for c in simple_coarse_grained if c != '_']) # remove all loops
   return dotbracket_to_coarsegrained_lev1(simple_coarse_grained_noloops) #shorten stacks

############################################################################################################
## extract base pairs
############################################################################################################

def get_basepair_indices_from_dotbracket(dotbracketstring):
   """extract a dictionary mapping each paired position with its partner:
   each base pair is represented twice: mapping from opening to closing bracket and vice versa"""
   base_pair_mapping = {}
   number_open_brackets = 0
   opening_level_vs_index = {}
   for charindex, char in enumerate(dotbracketstring):
      if char in ['(', '[']:
         number_open_brackets += 1
         opening_level_vs_index[number_open_brackets] = charindex
      elif char in [')', ']']:
         base_pair_mapping[charindex] = opening_level_vs_index[number_open_brackets]
         base_pair_mapping[opening_level_vs_index[number_open_brackets]] = charindex
         del opening_level_vs_index[number_open_brackets]
         number_open_brackets -= 1
      elif char in ['.', '_']:
         pass
      else:
         raise ValueError('invalid character in dot-bracket string')
      if number_open_brackets < 0:
         raise ValueError('invalid dot-bracket string')
   if number_open_brackets != 0:
      raise ValueError('invalid dot-bracket string')
   return base_pair_mapping

############################################################################################################
## functions supporting coarse-graining
############################################################################################################

def is_multiloop_from_tree(pos, tree_secondarystr, dotbracketstring):
   """is the unpaired site at pos a multiloop;
   which means an loop with one 'parent' stem and more than one 'child' stem
   this is tested by converting the secondary structure graph to a tree"""
   assert dotbracketstring[pos] in ['.', '_']
   parent_stem = [p for p in tree_secondarystr.predecessors(pos)][0]
   number_stems_below_parent_stem = len([x for x in tree_secondarystr.successors(parent_stem) if len([y for y in tree_secondarystr.successors(x)])>0])
   if number_stems_below_parent_stem > 1:
      return True
   else:
      return False

def is_exterior_loop_cg(pos, cg_string):
   """is the unpaired site at pos an exterior loop;
   which means either an unpaired stretch at either end of the chain
   or an unpaired region, which is not enclosed within the span of any base pair"""
   assert cg_string[pos] == '_'
   if cg_string[:pos].count('[') == cg_string[:pos].count(']'):
      return True 
   else:
      return False
      

def is_exterior_loop(pos, dotbracketstring):
   """is the unpaired site at pos an exterior loop;
   which means either an unpaired stretch at either end of the chain
   or an unpaired region, which is not enclosed within the span of any base pair"""
   assert dotbracketstring[pos] == '.'
   if '(' not in dotbracketstring[:pos] or ')' not in dotbracketstring[pos:]:
      return True 
   elif dotbracketstring[:pos].count('(') == dotbracketstring[:pos].count(')'):
      return True 
   else:
      return False


def get_tree_bpmap_from_dotbracket(dotbracketstring):
   """converting the secondary structure graph to a tree 
   (reference: https://www.tbi.univie.ac.at/RNA/ViennaRNA/doc/html/rna_structure_notations.html#sec_structure_representations_tree)
   each site is represented by a node whose integer reference is the site's sequence position;
   the tree is grounded in an additional node, '-1';
   each node is a child node of the next enclosing base pair"""
   tree, bp_mapping, current_node = nx.DiGraph(), {}, -1
   tree.add_node(current_node)
   for charindex, char in enumerate(dotbracketstring):
      if char in ['(', '[']:
         tree.add_edge(current_node, charindex)
         current_node = charindex
         if char == '[' and dotbracketstring[charindex+1] == ']':
            tree.add_edge(current_node, charindex + len(dotbracketstring) + 1) # add hairpin if not included
      elif char in [')', ']']:
         bp_mapping[charindex] = current_node
         bp_mapping[current_node] = charindex
         current_node = [n for n in tree.predecessors(current_node)][0]
      elif char in ['.', '_']:
         tree.add_edge(current_node, charindex)
      else:
         raise ValueError('invalid character in dot-bracket string')
   return tree, bp_mapping


############################################################################################################
## check viability of structure
############################################################################################################

def hairpin_loops_long_enough(structure):
   """check if any paired sites in the dot-bracket input structure
   are at least four sites apart"""
   bp_mapping = get_basepair_indices_from_dotbracket(structure)
   for baseindex1, baseindex2 in bp_mapping.iteritems():
      if abs(baseindex2-baseindex1) < 4:
         return False
   return True



def is_likely_to_be_valid_structure(structure, allow_isolated_bps=False):
   """tests if a structure in dotbracket format is likely to be a valid structure:
   basepairs closed, length of hairpin loops (>3), presence of basepairs and optionally isolated base pairs"""
   if not basepairs_closed(structure):
      return False
   if not hairpin_loops_long_enough(structure):
      return False
   if not structure.count(')') > 0:
      return False
   if not allow_isolated_bps and has_length_one_stack(structure):
      return False
   else:
      return True

def has_length_one_stack(dotbracketstring):
   """test if dotbracketstring has isolated base pairs"""
   for pos, char in enumerate(dotbracketstring):
      if char in [')', '('] and find_len_of_stack(pos, dotbracketstring) < 2:
         return 1
   return 0

def basepairs_closed(structure):
   """test if all brackets are closed correctly in a dot-bracket string"""
   try:
      bp_map = get_basepair_indices_from_dotbracket(structure)
      return True
   except (ValueError, KeyError):
      return False

def find_len_of_stack(pos, dotbracketstring):
   """ return the length of the stack at position pos in the structure given by the dot-bracket string
   bulges are defined as the end of the stack on both strands, i.e. if a base pair is at i, j, the base pair at i+1, j-2 would not belong to the same stack"""
   basepair_index_mapping = get_basepair_indices_from_dotbracket(dotbracketstring)
   assert pos in basepair_index_mapping
   node_of_basepair = min(pos, basepair_index_mapping[pos])
   ## make network of basepairs and connect adjacent basepairs in stacks - then stack size is size of component
   base_pair_neighbour_graph = nx.Graph()
   base_pair_neighbour_graph.add_nodes_from(set([min(b) for b in basepair_index_mapping.iteritems()])) # each base pair is represented by the pos of the opening bracket
   for b1, b2 in basepair_index_mapping.iteritems():
      for a1, a2 in basepair_index_mapping.iteritems():
         if b1<b2 and a1<a2: # both ordered from ( to )
            if b1 != a1: # distinct bas pairs
               if abs(b1-a1) == 1 and abs(b2-a2) == 1:
                  base_pair_neighbour_graph.add_edge(a1, b1)
   return len(nx.node_connected_component(base_pair_neighbour_graph, node_of_basepair))

############################################################################################################
## test
############################################################################################################
if __name__ == "__main__":
   teststructure1 = '...(((..((...)))))..((....)).'
   teststructure2 = '(((...))).'
   teststructure_ML1, teststructure_ML2, teststructure_ML3 = '.((((...))..((....))))..', '.((..((...))..((....))))..', '.((((...))..((....))..))..'
   ### test cutting ends
   assert cut_ends(teststructure1) == '(((..((...)))))..((....))'
   assert cut_ends(teststructure2) == '(((...)))'
   ### test stack of length one
   assert has_length_one_stack('.((.(...)))..') and has_length_one_stack('.((.(...).))..') and has_length_one_stack('.((.(...).))..')
   assert has_length_one_stack('.(...)...')
   assert not (has_length_one_stack(teststructure1) or has_length_one_stack(teststructure2))
   ### convert to coarsegrained structure - level 1
   assert dotbracket_to_coarsegrained_lev1(teststructure1) == '_[_[]]_[]_'
   assert dotbracket_to_coarsegrained_lev1(teststructure2) == '[]_'
   assert dotbracket_to_coarsegrained_lev1('.((.(...)))..') == '_[_[]]_' and dotbracket_to_coarsegrained_lev1('.((.(...).))..') == '_[_[]_]_' and dotbracket_to_coarsegrained_lev1('.(((...).))..') == '_[[]_]_'
   ### convert to coarsegrained structure -level 2
   assert dotbracket_to_coarsegrained_lev2(teststructure1) == '[_[]][]'
   assert dotbracket_to_coarsegrained_lev2(teststructure2) == '[]'
   assert dotbracket_to_coarsegrained_lev2('.((.(...)))..') == '[_[]]' and dotbracket_to_coarsegrained_lev2('.((.(...).))..') == '[_[]_]' and dotbracket_to_coarsegrained_lev2('.(((...).))..') == '[[]_]'
   assert dotbracket_to_coarsegrained_lev2(teststructure_ML1) == '[[][]]' and dotbracket_to_coarsegrained_lev2(teststructure_ML2) == '[[][]]' and dotbracket_to_coarsegrained_lev2(teststructure_ML3) == '[[][]]'
   ### convert to coarsegrained structure -level 5
   assert dotbracket_to_coarsegrained_lev5(teststructure1) == '[][]'
   assert dotbracket_to_coarsegrained_lev5(teststructure2) == '[]'
   assert dotbracket_to_coarsegrained_lev5('.((.(...)))..') == '[]' and dotbracket_to_coarsegrained_lev5('.((.(...).))..') == '[]' and dotbracket_to_coarsegrained_lev5('.(((...).))..') == '[]'
   assert dotbracket_to_coarsegrained_lev5('.((.((...)).((...))))..') == '[[][]]'
   assert dotbracket_to_coarsegrained_lev2(teststructure_ML1) == '[[][]]' and dotbracket_to_coarsegrained_lev2(teststructure_ML2) == '[[][]]' and dotbracket_to_coarsegrained_lev2(teststructure_ML3) == '[[][]]'
   ### test stack lengths
   assert find_len_of_stack(3, teststructure1) == find_len_of_stack(4, teststructure1) == find_len_of_stack(5, teststructure1) == find_len_of_stack(15, teststructure1) == find_len_of_stack(16, teststructure1) == 3
   assert find_len_of_stack(8, teststructure1) == find_len_of_stack(9, teststructure1) == find_len_of_stack(13, teststructure1) == find_len_of_stack(14, teststructure1) == 2
   assert find_len_of_stack(0, teststructure2) == find_len_of_stack(6, teststructure2) == find_len_of_stack(7, teststructure2) == find_len_of_stack(8, teststructure2) == 3
   ### test base pair extraction
   bp_mapping1 = get_basepair_indices_from_dotbracket(teststructure1)
   assert bp_mapping1[3] == 17 and bp_mapping1[4] == 16 and bp_mapping1[5] == 15 and bp_mapping1[8] == 14 and bp_mapping1[20] == 27
   assert bp_mapping1[17] == 3 and bp_mapping1[16] == 4 and bp_mapping1[15] == 5 and bp_mapping1[14] == 8 and bp_mapping1[27] == 20
   bp_mapping2 = get_basepair_indices_from_dotbracket(teststructure2)
   assert bp_mapping2[0] == 8 and bp_mapping2[1] == 7 and bp_mapping2[2] == 6 
   ### test length of hairpin loops
   assert hairpin_loops_long_enough(teststructure1) and hairpin_loops_long_enough(teststructure2)
   assert not hairpin_loops_long_enough('...(((..((..)))))..((....)).')
   ## test whether invalid structures are recognised correctly
   assert is_likely_to_be_valid_structure('...((..((...))))..((...))..', allow_isolated_bps=False)
   assert is_likely_to_be_valid_structure('...((..((...))))..(...)..', allow_isolated_bps=True) and not is_likely_to_be_valid_structure('...((..((...))))..(...)..', allow_isolated_bps=False)
   assert not is_likely_to_be_valid_structure('...((..((....)))))..((...)..', allow_isolated_bps=False)
   assert not is_likely_to_be_valid_structure('...((..((..))))..((...))..', allow_isolated_bps=False)
   assert not is_likely_to_be_valid_structure('...(((..((...))))..((...))..', allow_isolated_bps=False)
   teststructure3, teststructure4, teststructure5 = '.((((.((....)).))..((....))))..', '.((((((..((.....)))))((.....)))))..', '.((((...))..((....))..))..'
   for dotbracketstructure in [ teststructure1, teststructure2, teststructure3, teststructure4, teststructure5, teststructure_ML1, teststructure_ML2, teststructure_ML3]:
      assert lev1_to_lev2(dotbracket_to_coarsegrained_lev1(dotbracketstructure)) == dotbracket_to_coarsegrained_lev2(dotbracketstructure)
      assert lev1_or_2_to_lev5(dotbracket_to_coarsegrained_lev1(dotbracketstructure)) == lev1_or_2_to_lev5(dotbracket_to_coarsegrained_lev2(dotbracketstructure)) == dotbracket_to_coarsegrained_lev5(dotbracketstructure)
   ####### test structure generation
   allstructures_eleven = generate_all_allowed_dotbracket(L=11, allow_isolated_bps=False)
   assert len(set(allstructures_eleven)) == len(allstructures_eleven)
   test_allstructures_eleven = ['((((...))))', '(((...)))..',  '..(((...)))', '.(((...))).', '(((....))).',  '.(((....)))',  '(((.....)))',
                            '((...))....', '.((...))...', '..((...))..', '...((...)).', '....((...))',
                            '((....))...', '.((....))..', '..((....)).', '...((....))',
                            '((.....))..', '.((.....)).', '..((.....))', '((......)).', '.((......))', '((.......))']
   for structure in allstructures_eleven:
      assert structure in test_allstructures_eleven
   assert len(test_allstructures_eleven) == len(allstructures_eleven)
   print 'sucessfully finished tests'




