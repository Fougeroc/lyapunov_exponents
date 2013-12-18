from interval_exchange import *
from path_vector import *
from sage.all import ceil

def lyap(ie, it) :
    t = 0
    N = ie.genus()
    if N == 0:
        return []
    v =  [VectPaths(ie.labels(), ie.degrees()).random() for i in range(N)]
    theta = [0 for i in range(N)]
    if ie._lengths == None:
        ie.rand_lg()
    ie.normalise()
    for i in range(it) :
        #v =  [VectPaths(ie.labels(),ie.degrees()).id(n,a) for n in ie.degrees() for a in ie.labels()]
        (A, B, c, perm_one, perm_two) = ie.rauzy()
        
        for j in range(N):                                                    #first apply the change of B
            image = VectPaths(v[j]._labels, v[j]._degrees)
            for n, a in v[j]._iter:
                if a == B:                                                #return inversed and transposed matrix, thus  add (-1)
                    image._vect[n][a] = v[j].val(n,B) + c*v[j].val(perm_one(n),A)
                else:
                    image._vect[n][a] = v[j].val(n,a)
            v[j] = image
        for j in range(N):                                                    #then apply the permutation
            image = VectPaths(v[j]._labels, v[j]._degrees)
            for n, a in v[j]._iter:
                if a == A:
                    image._vect[n][a] = v[j].val(perm_two(n),A)
                else:
                    image._vect[n][a] = v[j].val(n,a)
            v[j] = image
        
        # for j in range(N):                                                    #first apply the change of B
        #     image = VectPaths(v[j]._labels, v[j]._degrees)
        #     for n, a in v[j]._iter:
        #         if a == B:                                                #return inversed and transposed matrix, thus  add (-1)
        #             image._vect[n][a] = v[j].val(n,B) + c*v[j].val(name_A(name_B.inverse()(n)),A)
        #         else:
        #             image._vect[n][a] = v[j].val(n,a)
        #     v[j] = image
        # for j in range(N):                                                    #then apply the permutation
        #     image = VectPaths(v[j]._labels, v[j]._degrees)
        #     for n, a in v[j]._iter:
        #         if a == A:
        #             if Ap_give_name:
        #                 image._vect[n][a] = v[j].val(name_A(ident_B.inverse()(n)),A)            #inverse and transpose of a permutation matrix are the same
        #             else:
        #                 image._vect[n][a] = v[j].val(name_A_twin(n),A)
        #         else:
        #             image._vect[n][a] = v[j].val(n,a)
        #     v[j] = image

        if ie.diff_sum() > 2**(-precision+20):
            t += -log(ie.normalise())
        if ie.min_length() < 2**(-6):
            nm = vect_paths_ortho_householder(v)
            t += -log(ie.normalise())
            for i in range(N):
                    theta[i] += log(nm[i])
    if t == 0:
        return [0 for i in xrange(N)]
    return [theta[i]/t for i in xrange(N)]


def lyap_isotopic_decomposed(ie, it) :
    t = 0
    isotopic_dimension = [ceil(ie.character_degree()[k]*ie.number_of_character_appeareance(k)/2) for k in range(ie.n_characters())]
    v_iter = [ (k, l) for k in range(ie.n_characters()) for l in range(isotopic_dimension[k])]
    v = [[ie.canonical_VectPaths().random() for l in range(isotopic_dimension[k])] for k in range(ie.n_characters())]
    theta = [[0 for l in range(isotopic_dimension[k])] for k in range(ie.n_characters())]
    def project():
        for k, l in v_iter:
            v[k][l] = ie.vector_isotopic_projection(k, v[k][l])
    for i in range(it) :
        project()
        (c,A,B, name_A, name_B, name_A_twin, Ap_give_name, ident_B) = ie.rauzy_rev()
        for k, l in v_iter:                                #first apply the change of B
            w = v[k][l]
            image = VectPaths(w._labels, w._degrees)
            for n, a in w._iter:
                if a == B:                                                #return inversed and transposed matrix, thus  add (-1)
                    image._vect[n][a] = CC(w.val(n,B) + c*w.val(name_A(name_B.inverse()(n)),A))
                else:
                    image._vect[n][a] = CC(w.val(n,a))
            v[k][l] = image
        for k, l in v_iter:                                #then apply the permutation
            w = v[k][l]
            image = VectPaths(w._labels, w._degrees)
            for n, a in w._iter:
                if a == A:
                    if Ap_give_name:
                        image._vect[n][a] = CC(w.val(name_A(ident_B.inverse()(n)),A))            #inverse and transpose of a permutation matrix are the same
                    else:
                        image._vect[n][a] = CC(w.val(name_A_twin(n),A))
                else:
                    image._vect[n][a] = CC(w.val(n,a))
            v[k][l] = image
        if ie.diff_sum() > diff_sum_seuil:
            t += -log(ie.normalise())
        if ie.min_length() < 2**(-10):
            t += -log(ie.normalise())
            for k in range(ie.n_characters()):
                nm = vect_paths_ortho_householder(v[k])
                for l in range(isotopic_dimension[k]):
                    theta[k][l] += log(nm[l])
    if t == 0:
        return [[0 for l in range(isotopic_dimension[k])] for k in range(ie.n_characters())]
    else:
        return [[theta[k][l]/t for l in range(isotopic_dimension[k])] for k in range(ie.n_characters())]
