from lyapunov_exponent import *
from time import time
from sage.all import set_random_seed

#set_random_seed(1)

N = 10000

intervals = [[Interval('a',1), Interval('a',-1), Interval('b',1)], [Interval('b',1), Interval('c', 1), Interval('c',-1)]]


sigma_a = Permutation('(1,2)(3,4)(5,6)')
sigma_b = Permutation('(1,3)(2,4)(5)(6)')
sigma_c = Permutation('(1,6)(2,5)(3,4)')

permutations = {'a': sigma_a, 'b': sigma_b, 'c': sigma_c}

def main(R,N,string):
    with open(string + 'dec_results_%s'%(N), 'w') as f:
        t0 = time()
        f.write("%s\n%s iterations\n\n"%(R,N))
        f.write(str(lyap_isotopic_decomposed(R,N)))
        f.write("\n\ntime of process : %s\n"%(time()-t0))
    f.closed

def main(R, N, string):
    M = R.h_one_intersection_matrix()
    if M <> -M.transpose():
        print string
        print R._intervals
        print R.labels()
        print R.cycles
    T = R.intersection_isotopic()
    print 
    print string
    print R.signature_intersection
    print R.signatures_isotopic
#    for i in range(len(T)):
#        print T[i]
#        print R.signature(i)
#        print ","
    #print R._intervals
    #print R.labels()
    #v = R.canonical_VectPaths().copy_vector(R.h_one_to_generator[1])
    #w = R.canonical_VectPaths().copy_vector(R.h_one_to_generator[0])
    #print R.global_intersection(v, w)
    #print R.global_intersection(w,v)
    
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
