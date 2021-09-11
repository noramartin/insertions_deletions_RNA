import numpy as np
import random


base_list = ['A', 'C', 'U', 'G']


############################################################################################################
## plots
############################################################################################################

def decimal_to_scientific_notation(dec, n=1):
   """ transform floating point number dec to a string in scientific notation
   n is the number of digits for rounding (n=1 means that 233 will be 2.3 * 10**2)"""
   if dec == 0:
      return str(round(dec, n))
   try:
      b = int(np.floor(np.log10(abs(dec))))
      a = dec/float(10**b)
      if b == 0:
         return str(round(a, n))
      else:
         return str(round(a, n))+r'$\times 10^{{{0}}}$'.format(b)
   except OverflowError:
      return str(round(dec, n))
      
############################################################################################################
## sampling for shape frequency
############################################################################################################

def get_shape_from_random_sequence(repetition, L, GPfunction):
   """generate a random sequence of length L, 
      apply the GPfunction and return the sequence and GPfunction(sequence);
      repetition is an unused variable and introduced only so that pool.map can be used"""
   seq = ''.join([random.choice(['A', 'U', 'C', 'G']) for c in range(L)])
   return seq, GPfunction(seq)

############################################################################################################
## test 
############################################################################################################
if __name__ == "__main__":
   print 'test decimal_to_scientific_notation'
   assert decimal_to_scientific_notation(0.01, n=1) == '1.0' + r'$\times 10^{-2}$'
   assert decimal_to_scientific_notation(100, n=1) == '1.0' + r'$\times 10^{2}$'
   assert decimal_to_scientific_notation(123, n=1) == '1.2' + r'$\times 10^{2}$'
   assert decimal_to_scientific_notation(0.0, n=1) == '0.0'
   assert decimal_to_scientific_notation(1.2, n=1) == '1.2'
   print 'test get_shape_from_random_sequence - does it work for parallel'
   from multiprocessing import Pool
   from functools import partial
   from Vienna_cgshape_probs import find_most_freq_shape
   GPfunction = partial(find_most_freq_shape, shape_level=2, range_kbT=15, folding_criterion=1.1, allow_isolated_bps=False, dangling_ends_option=2)
   shape_funct_random_seq = partial(get_shape_from_random_sequence, L=30, GPfunction=GPfunction)
   pool = Pool(processes = 5)
   pool_result = pool.map(shape_funct_random_seq, np.arange(10))
   pool.close()
   pool.join()
   assert len(set(zip(*pool_result)[0])) == len(pool_result) # test that random seed applied correctly and no double sequences due to parallel sampling
   for seq, shape in pool_result:
      assert shape == GPfunction(seq)



