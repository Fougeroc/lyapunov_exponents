from lyapunov_exponent import *
from time import time
from sage.all import set_random_seed

#set_random_seed(1)

N = 10000


def main(R):
     #with open(str ing + 'dec_results_%s'%(N), 'w') as f:
     #     t0 = time()
     #     f.write("%s\n%s iterations\n\n"%(R,N))
     #     s = str(lyap_isotopic_decomposed(R,N))
     #     print s
     #     f.write(s)
     #     f.write("\n\ntime of process : %s\n"%(time()-t0))
     # f.closed
#     print R.intersection_forms_isotopic
     # for i in range(R.n_characters()):
     #      print R.isotopic_projection_matrix(i)
     # for i in range(R.n_characters()):
     #       print R.signatures_isotopic(i)          
     # print [R.signatures_isotopic(i) for i in range(R.n_characters())]
      print "Stratum :  ",
      print R.stratum()
      print R
      #print lyap(R, N)
      print R.lyapunov_exponents_H_plus(nb_iterations = 10000)
#main(cyclic_cover_iet(3, [5, 1, 3, 3]), 10000, "test")
#print "\n"

#Eierwohlmilchsau tous nuls
print "Eier.."
main(cyclic_cover_iet(4, [1, 1, 1, 1]))
print "\n"

print "Ornythorinque"
main(cyclic_cover_iet(6, [1, 1, 1, 3]))
print "\n"

intervals = [[Interval('a',1), Interval('a',-1), Interval('b',1)], [Interval('b',-1), Interval('c', -1), Interval('c',1)]]


sigma_a = Permutation('(1,2)(3,4)(5,6)')
sigma_b = Permutation('(1,3)(2,4)(5)(6)')
sigma_c = Permutation('(1,6)(2,5)(3,4)')

permutations = {'a': sigma_a, 'b': sigma_b, 'c': sigma_c}

R = IntExchange(intervals, permutations = permutations)

print "L a trois carreaux"
main(R)

# main(cyclic_cover_iet(5, [3, 1, 2, 2]), 10000, "test")
# print cyclic_cover_exact(5, [3, 1, 2, 2])
# print "\n"

# main(cyclic_cover_iet(5, [2, 2, 2, 2]), 10000, "test")
# print cyclic_cover_exact(5, [2, 2, 2, 2])
# print "\n"

#main(cyclic_cover_iet(5, [3, 1, 2, 2]), 1000, "test")
#print cyclic_cover_exact(5, [3, 1, 2, 2])
#print "\n"

# for i in range(2, 6):
#      main(pillow_case_iet(2*i+1), 1, "test")

# def main(N, a):
#      R = cyclic_cover_iet(N, a)
#      print [R.signatures_isotopic(i) for i in range(R.n_characters())]
#      print cyclic_cover_exact(N, a)
#      print "\n"

# main(7, [1, 5, 2, 4])









"""   
intervals = [[Interval('a',1), Interval('b',1), Interval('c',1), Interval('b',-1), Interval('d',1)], 
               [Interval('d',1), Interval('e',1), Interval('e',-1), Interval('c',1), Interval('a',1)]]

intervals = [[Interval('a',1), Interval('b',1), Interval('c',1), Interval('c',-1)], 
               [Interval('b',1), Interval('d',1), Interval('d',-1), Interval('a',1)]]


main(trivial_IntExchange(intervals).rand_lg(), N, "Q(2,-1,-1)" + "_base_")
#main(orientable_double_cover_IntExchange(intervals).rand_lg(), N, "Q(2,-1,-1)" + "_cover_")



#Quadratic[{-1, -1, -1, -1}]
intervals = [[Interval('a',1), Interval('a',-1), Interval('b',1)], [Interval('b',1), Interval('c', 1), Interval('c',-1)]]

#main(trivial_IntExchange(intervals).rand_lg(), N, "Q(-1,-1,-1,-1)" + "_base_")
main(orientable_double_cover_IntExchange(intervals).rand_lg(), N, "Q(-1,-1,-1,-1)" + "_cover_")

#Quadratic[{1,1,-1,-1}]
intervals = [[Interval('a',1), Interval('b',1), Interval('c',1), Interval('d',1), Interval('d',-1)], 
               [Interval('c',1), Interval('b',1), Interval('e',1), Interval('e',-1), Interval('a',1)]]

main(trivial_IntExchange(intervals).rand_lg(), N, "Q(1,1,-1,-1)" + "_base_")
main(orientable_double_cover_IntExchange(intervals).rand_lg(), N, "Q(1,1,-1,-1)" + "_cover_")

#Quadratic[{2, -1, -1}]
intervals = [[Interval('a',1), Interval('b',1), Interval('c',1), Interval('c',-1)], 
               [Interval('b',1), Interval('d',1), Interval('d',-1), Interval('a',1)]]

main(trivial_IntExchange(intervals).rand_lg(), N, "Q(2,-1,-1)" + "_base_")
main(orientable_double_cover_IntExchange(intervals).rand_lg(), N, "Q(2,-1,-1)" + "_cover_")

#Quadratic[{2,2}]
intervals = [[Interval('a',1), Interval('b',1), Interval('c',1), Interval('b',-1), Interval('d',1)], 
               [Interval('e',1), Interval('d',1), Interval('e',-1), Interval('c',1), Interval('a',1)]]

main(trivial_IntExchange(intervals).rand_lg(), N, "Q(2,2)" + "_base_")
main(orientable_double_cover_IntExchange(intervals).rand_lg(), N, "Q(2,2)" + "_cover_")


#Quadratic[{1,1,1,1}]
intervals = [[Interval('a',1), Interval('b',1), Interval('c',1), Interval('d',1), Interval('b',-1), Interval('e',1), Interval('f',1)], 
               [Interval('c',1), Interval('g',1), Interval('f',1), Interval('e',1), Interval('g',-1), Interval('d',1), Interval('a',1)]]

main(trivial_IntExchange(intervals).rand_lg(), N, "Q(1,1,1,1)" + "_base_")
main(orientable_double_cover_IntExchange(intervals).rand_lg(), N, "Q(1,1,1,1)" + "_cover_")

#Quadratic[{2,1,1}]
intervals = [[Interval('a',1), Interval('b',1), Interval('c',1), Interval('b',-1), Interval('d',1), Interval('e',1)], 
               [Interval('f',1), Interval('e',1), Interval('d',1), Interval('f',-1), Interval('c',1), Interval('a',1)]]

main(trivial_IntExchange(intervals).rand_lg(), N, "Q(2,1,1)" + "_base_")
main(orientable_double_cover_IntExchange(intervals).rand_lg(), N, "Q(2,1,1)" + "_cover_")

#Quadratic[{1,1,1,-1,-1,-1}]
intervals = [[Interval('a',1), Interval('b',1), Interval('c',1), Interval('d',1), Interval('d',-1), Interval('e',1), Interval('f',1), Interval('f',-1)], 
               [Interval('c',1), Interval('e',1), Interval('b',1), Interval('g',1), Interval('g',-1), Interval('a',1)]]

main(trivial_IntExchange(intervals).rand_lg(), N, "Q(1,1,1,-1,-1,-1)" + "_base_")

main(orientable_double_cover_IntExchange(intervals).rand_lg(), N, "Q(1,1,1,-1,-1,-1)" + "_cover_")

#Quadratic[{5,-1}]
intervals = [[Interval('a',1), Interval('b',1), Interval('c',1), Interval('b',-1), Interval('d',1)], 
               [Interval('d',1), Interval('e',1), Interval('e',-1), Interval('c',1), Interval('a',1)]]

main(trivial_IntExchange(intervals).rand_lg(), N, "Q(5,-1)" + "_base_")
main(orientable_double_cover_IntExchange(intervals).rand_lg(), N, "Q(5,-1)" + "_cover_")

#Quadratic[{4,1,-1}]
intervals = [[Interval('a',1), Interval('b',1), Interval('c',1), Interval('b',-1), Interval('d',1), Interval('e',1)], 
               [Interval('e',1), Interval('d',1), Interval('f',1), Interval('f',-1), Interval('c',1), Interval('a',1)]]

main(trivial_IntExchange(intervals).rand_lg(), N, "Q(4,1,-1)" + "_base_")
main(orientable_double_cover_IntExchange(intervals).rand_lg(), N, "Q(4,1,-1)" + "_cover_")
"""
