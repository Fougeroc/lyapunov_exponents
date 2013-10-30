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
    ie.normalise()
    for i in range(it) :
        #v =  [VectPaths(ie.labels(),ie.degrees()).id(n,a) for n in ie.degrees() for a in ie.labels()]
        (c,A,B,name_A,name_B,name_A_twin,Ap_give_name,ident_B) = ie.rauzy_rev()
        for j in range(N):                                                    #first apply the change of B
            image = VectPaths(v[j]._labels, v[j]._degrees)
            for a in v[j]._labels:
                for n in v[j]._degrees:
                    if a == B:                                                #return inversed and transposed matrix, thus  add (-1)
                        image._vect[n][a] = v[j].val(n,B) - c*v[j].val(name_A(name_B.inverse()(n)),A)
                    else:
                        image._vect[n][a] = v[j].val(n,a)
            v[j] = image
        for j in range(N):                                                    #then apply the permutation
            image = VectPaths(v[j]._labels, v[j]._degrees)
            for a in v[j]._labels:
                for n in v[j]._degrees:
                    if a == A:
                        if Ap_give_name:
                            image._vect[n][a] = v[j].val(name_A(ident_B.inverse()(n)),A)            #inverse and transpose of a permutation matrix are the same
                        else:
                            image._vect[n][a] = v[j].val(name_A_twin(n),A)
                    else:
                        image._vect[n][a] = v[j].val(n,a)
            v[j] = image
        if ie.diff_sum() > 2**(-precision+20):
            t += -log(ie.normalise())
        if ie.min_length() < 2**(-6):
            nm = vect_paths_ortho_householder(v)
            t += -log(ie.normalise())
            for i in range(N):
                    theta[i] += log(nm[i])
    return [theta[i]/t for i in range(N)]


    
    


def lyap_isotopic_decomposed(ie, it) :
    t = 0
    isotopic_dimension = [ceil(ie.character_degree[j]*ie.number_of_character_appeareance[j]/2) for j in range(ie.n_characters)]
    v = [[ie.canonical_VectPaths().random() 
          for i in range(isotopic_dimension[j])] 
         for j in range(ie.n_characters)]
    theta = [[0 for i in range(isotopic_dimension[j])] for j in range(ie.n_characters)]
    v_iter = [(k,l) for k in range(ie.n_characters) for l in range(isotopic_dimension[k])]
    def project():
        for k, l in v_iter:
            v[k][l] = self.isotopic_projection(k, v[k][l])
    for i in range(it) :
        project()
        #v =  [VectPaths(ie.labels(),ie.degrees()).var(n,a) for n in ie.degrees() for a in ie.labels()]
        (c,A,B,name_A,name_B,name_A_twin,Ap_give_name,ident_B) = ie.rauzy_rev()
        for k in range(ie.n_characters):
            for l in range(isotopic_dimension[k]):                                #first apply the change of B
                w = v[k][l]
                image = VectPaths(w._labels, w._degrees)
                for a in w._labels:
                    for n in w._degrees:
                        if a == B:                                                #return inversed and transposed matrix, thus  add (-1)
                            image._vect[n][a] = w.val(n,B) - c*w.val(name_A(name_B.inverse()(n)),A)
                        else:
                            image._vect[n][a] = w.val(n,a)
                v[k][l] = image
        for k in range(ie.n_characters):
            for l in range(isotopic_dimension[k]):                                #then apply the permutation
                w = v[k][l]
                image = VectPaths(w._labels, w._degrees)
                for a in w._labels:
                    for n in w._degrees:
                        if a == A:
                            if Ap_give_name:
                                image._vect[n][a] = w.val(name_A(ident_B.inverse()(n)),A)            #inverse and transpose of a permutation matrix are the same
                            else:
                                image._vect[n][a] = w.val(name_A_twin(n),A)
                        else:
                            image._vect[n][a] = w.val(n,a)
                v[k][l] = image
        if ie.diff_sum() > diff_sum_seuil:
            t += -log(ie.normalise())
        if ie.min_length() < 2**(-10):
            t += -log(ie.normalise())
            for k in range(ie.n_characters):
                nm = vect_paths_ortho_householder(v[k])
                for i in range(isotopic_dimension[k]):
                    theta[k][i] += log(nm[i])
    return [[theta[k][i]/t for i in range(isotopic_dimension[k])] for k in range(ie.n_characters)]
