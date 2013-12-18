r"""
Python bindings for various computation of Lyapunov exponents.
"""

from libc.stdlib cimport malloc,free

cdef extern from "lyapunov_exponents.h":
    ctypedef struct quad_cyclic_cover:
        pass
    ctypedef struct generalized_permutation:
        pass
    ctypedef struct permutation:
        pass

    # initialisation/allocation/free
    generalized_permutation * new_generalized_permutation(int *perm, int *twin, int k, int n)
    quad_cyclic_cover * new_quad_cyclic_cover(generalized_permutation * gp, permutation ** sigma, size_t degree, size_t nb_vectors)
    void set_lengths(quad_cyclic_cover * qcc, long double	       *lengths)
    void set_random_lengths_quad_cyclic_cover(quad_cyclic_cover * qcc)
    void set_random_vectors(quad_cyclic_cover * qcc)
    void renormalize_length_quad_cyclic_cover(quad_cyclic_cover *qcc)
    permutation * new_permutation(size_t degree)
    void perm_copy_table(permutation *perm, size_t* tab)

    #int check_generalized_permutation(generalized_permutation *p)
    #int check_quad_cyclic_cover(quad_cyclic_cover * qcc)

    void free_generalized_permutation(generalized_permutation ** gp)
    void free_quad_cyclic_cover(quad_cyclic_cover ** qcc)
    void free_permutation(permutation **perm)

    # print
    #void print_generalized_permutation(generalized_permutation * p)
    void print_quad_cyclic_cover(quad_cyclic_cover * qcc)
    void print_vectors(quad_cyclic_cover * qcc)
    void print_permutation(permutation *perm)

    # algorithms
    #void renormalize_length_quad_cyclic_cover(quad_cyclic_cover * qcc)
    #void rauzy_induction_H_plus_quad_cyclic_cover(quad_cyclic_cover *qcc)

    int init_GS(size_t dim)
    void free_GS()
    #void orthogonalize_GS(quad_cyclic_cover * qcc, double * theta)

    void lyapunov_exponents_H_plus(quad_cyclic_cover *qcc, double *theta, size_t nb_induction)
    void top_lyapunov_exponents_H_plus(quad_cyclic_cover *qcc, double *theta, size_t nb_iterations)

def lyapunov_exponents_H_plus_cyclic_cover(
        gp, twin, k, n, sigma, degree,
        nb_vectors, nb_experiments, nb_iterations, lengths = None):
    r"""
    Compute the Lyapunov exponents of the H^+ part of the KZ-cocycle for
    covering locii.

    We assume that all the inputs are clean. If not, it may cause some SEGFAULT
    which would interrupt python!

    INPUT:

    - ``gp`` -- a generalized permutation given as a list of integers

    - ``twin`` -- the twin data of the gp

    - ``k`` -- the length of the top interval

    - ``n`` -- the length of gp

    - ``sigma`` -- covering data

    - ``nb_vectors`` -- the number of vectors to use

    - ``nb_experiments`` -- number of experimets

    - ``nb_iterations`` -- the number of iterations of the Rauzy-Zorich
      induction to perform

    - ``verbose`` -- if ``True`` print additional information concerning the
      mean and standard deviation
    """
    cdef int *p, *t   # permutation, twin, sigma
    cdef permutation **s
    cdef size_t *tab
    cdef generalized_permutation *gp_c
    cdef quad_cyclic_cover *qcc
    cdef double * theta
    


    # convert the data of into C values
    p = <int *> malloc(2 * n * sizeof(int))
    t = <int *> malloc(2 * n * sizeof(int))
    s = <permutation **> malloc(n * sizeof(permutation*))
    tab = <size_t *> malloc(degree * sizeof(size_t))

    for i from 0 <= i < n:
        p[i] = gp[i]
        t[i] = twin[i]
        for j in range(degree):
            tab[j] = sigma[j + degree * i]
        s[i] = new_permutation(degree)
        perm_copy_table(s[i], tab)
    for i from n <= i < 2*n:
        p[i] = gp[i]
        t[i] = twin[i]

    theta = <double *> malloc((nb_vectors+1) * sizeof(double))

    gp_c = new_generalized_permutation(p, t, k, n)
    qcc = <quad_cyclic_cover *> new_quad_cyclic_cover(gp_c,s,degree,nb_vectors)

    if lengths == None:
       set_random_lengths_quad_cyclic_cover(qcc)
    else:
        l = <long double *> malloc(n * sizeof(long double))	
        for i from 0 <= i < n:
            l[i] = <long double> lengths[i]
        set_lengths(qcc, l)
        free(l)

    free_generalized_permutation(&(gp_c))
    free(p)
    free(t)
    free(s)

    res = [[] for _ in xrange(nb_vectors+1)]
    if nb_vectors == 1:
        for i in xrange(nb_experiments):
            top_lyapunov_exponents_H_plus(qcc, theta, nb_iterations)
            for j in xrange(2):
                res[j].append(theta[j])

    else:
        init_GS(nb_vectors)
        for i in xrange(nb_experiments):
            lyapunov_exponents_H_plus(qcc, theta, nb_iterations)
            for j in xrange(nb_vectors+1):
                res[j].append(theta[j])

    free_quad_cyclic_cover(&qcc)
    free(theta)

    return res
