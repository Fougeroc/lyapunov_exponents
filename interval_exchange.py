#Here we give one function to get a random length vector for the IE
from sage.all import random, log, floor, Permutations, Permutation, copy, deepcopy, RealField, gap, SymmetricGroup, vector, matrix, ZZ, Integer, I, real_part
from path_vector import *

precision = 128
diff_sum_seuil = 2**(-precision+20)

R = RealField(precision)


class Interval(object) :
    def __init__(self, label, orientation):
        assert orientation == 1 or orientation == -1
        self.label = label
        self.orientation = orientation
    def __repr__(self):
        if self.orientation == 1:
            return " %s"%(self.label)
        else:
            return "-%s"%(self.label)
    def invert(self):
        if self.orientation == 1:
            self.orientation = -1
        else:
            self.orientation = 1


def index_cycle(i, cycles):
    """
    return index of the cycle in the array cycles where i appears
    """
    d = len(cycles)
    cursor = 0
    while cycles[cursor].count(t) == 0 and cursor <= d:
        cursor += 1
    if cursor == d:
        raise NameError('Index not in the cycle')
    return cursor

    
class IntExchange(object):
    r"""
    Create an interval exchange with possibility of considering a finite cover by giving
    the permutations associated to the natural generator of the  homotopy group of the suspension surface.

    We have a convention for labelling the saddle connections of the cover of the suspension surface :
    first we define a directed path between the two saddle connections, which will be a cycle in the ground surface,
    the direction will be choosen such that the intersection number with both saddle connections, directed according
    to the orientation is positive (egal to 1).
    The label bearing the number of the copy of the cover will be the one to which the path arrives.

    
    INPUT:

    - ``intervals`` -- a list of two lists of Interval class elements
    
    - ``lengths`` -- a dictionnary of positive floats for each label appearing in the intervals;
    the length of intervals of the corresponding label

    - ``permutations`` -- dictionnary of permutation for each label appearing in the intervals;
    The permutation associated to the label `a` gives for the saddle connection of label `a`
    bearing the number `i`, the number of the other copy on which there is the saddle connection
    to which it is identified.


    EXAMPLES::

        sage: intervals = [[Interval('a',1), Interval('a',-1), Interval('b',1)],
        ...   [Interval('b',1), Interval('c', 1), Interval('c',-1)]]

        sage: sigma_a = Permutation('(1,2)(3,4)(5,6)')
        sage: sigma_b = Permutation('(1,3)(2,4)(5)(6)')
        sage: sigma_c = Permutation('(1,6)(2,5)(3,4)')

        sage: permutations = {'a': sigma_a, 'b': sigma_b, 'c': sigma_c}

        sage: ie = IntExchange(intervals, None, permutations).rand_lg()

        sage: ie
        Intervals    : [ a, -a,  b]
                       [ b,  c, -c]
        Lengths      : {'a': 0.111439293741037, 'c': 0.111439293741037, 'b': 0.5143475134191677}
        Permutations : {'a': [2, 1, 4, 3, 6, 5], 'c': [6, 5, 4, 3, 2, 1], 'b': [3, 4, 1, 2, 5, 6]}

    """
    def __init__(self, intervals, lengths, permutations):
        self._intervals = deepcopy(intervals)
        self._lengths = deepcopy(lengths)
        self._permutations= deepcopy(permutations)
        self.degree = len(permutations[permutations.keys()[0]])
        self.number_of_intervals = [len(intervals[0]), len(intervals[1])]
        self._twins = [[self._twin(i, j) for j in range(self.number_of_intervals[i])] for i in range(2)]
        self.total_number_of_intervals_base = len(intervals[0]) + len(intervals[1])
        self.total_number_of_intervals_cover = (len(intervals[0]) + len(intervals[1]))*self.degree
        self.cycles = self.edge_cover_permutation().to_cycles()
        self.n_cycles = len(self.cycles)
        self.h_one_rel_to_generator, self.generator_to_h_one_rel = self._free_basis() 
        self.galois_group = self._galois_group()
        self.character_table, self.character_degree, self.galois_group_order, self.galois_group_permutation, self.n_characters  = self._characters()
#        self.h_one_projection_matrix = self._h_one_projection() #projection on the kernel inside h_one_rel
        self.h_one_to_h_one_rel, self.h_one_rel_to_h_one = self._h_one_to_h_one_rel(), self._h_one_rel_to_h_one()
        self.h_one_to_generator = self.h_one_to_h_one_rel*self.h_one_rel_to_generator
        self.generator_to_h_one = self.generator_to_h_one_rel*self.h_one_rel_to_h_one
        self.intersection_matrix = self.h_one_intersection_matrix()
        self.signature_intersection = self.signature(self.intersection_matrix)
        self.representation_matrix = self._representation_matrix()
        self.intersection_forms_isotopic = self.intersection_isotopic()
        self.signatures_isotopic = [self.signature(self.intersection_forms_isotopic[i]) for i in range(len(self.intersection_forms_isotopic))]
        self.number_of_character_appeareance = [self._number_of_character_appeareance(i) for i in range(self.n_characters)]
    def __repr__(self):
        return ("Intervals    : %s\n               %s\nLengths      : %s\nPermutations : %s\n"
                %(self._intervals[0], self._intervals[1], self._lengths,self._permutations))
    def _twin(self, i, j):
        r"""Give the position of the twin interval associated to the one in line i and place j
        """
        for n in range(self.number_of_intervals[i]):
            if n <> j and self._intervals[i][n].label == self._intervals[i][j].label:
                return(i,n)
        for n in range(self.number_of_intervals[(i+1)%2]):
            if self._intervals[(i+1)%2][n].label == self._intervals[i][j].label:
                return((i+1)%2,n)
    def degrees(self):
        r"""
        Give the list of numbers given to each copy of the suface in the cover.
        It goes from 1 to the degree of the cover.
        """
        return(range(1, self.degree + 1))
    def iter_lab(self):
        r"""
        Gives the list of the labels of intervals in order to iterate
        """
        return(self._permutations.iterkeys())
    def labels(self):
        r"""
        Gives the list of the labels of intervals
        """
        return(self._permutations.keys())
    def canonical_VectPaths(self):
        return(VectPaths(self.labels(), self.degrees()))
    def cover_generators(self):
        return([(d,a) for d in self.degrees() for a in self.labels()])
    def _free_basis(self):
        r"""
        Return matrix from free basis of H1(X, \Sigma) to R^(d*lab)
        and its pseudo inverse, from R^(d*lab) to H1(X, \Sigma)
        """
        bound_equation = {d: self.canonical_VectPaths() for d in self.degrees()}
        for d in self.degrees():
            for i in range(2):
                for j in range(self.number_of_intervals[i]):
                    bound_equation[d]._vect[self.name(i,j)(d)][self.label(i,j)] += ((-1)**i)*self._intervals[i][j].orientation
        s = []
        def bind(ind,deg,lab):
            line = deepcopy(bound_equation[ind])
            c = line.val(deg, lab)
            for i in self.degrees():
                coeff = bound_equation[i].val(deg,lab)
                for d,a in self.cover_generators():
                    bound_equation[i]._vect[d][a] -= line.val(d,a)*coeff/c
                    if projections[d][a].val(deg, lab) <> 0:
                        coeff_aux = projections[d][a].val(deg, lab)
                        for d_aux, a_aux in self.cover_generators():
                            projections[d][a]._vect[d_aux][a_aux] -= line.val(d_aux,a_aux)*coeff_aux/c
        def non_zero(deg, lab):
            for i in self.degrees():
                if bound_equation[i].val(deg,lab) <> 0:
                    return i
            return None
        projections = {d: {a: self.canonical_VectPaths().id(d,a) for a in self.labels()} for d in self.degrees()}
        for d,a in self.cover_generators():
            i = non_zero(d,a)
            if i == None:
                s.append((d,a))
            else:
                bind(i,d,a)
        def vect_paths_id((d,a)):
            return self.canonical_VectPaths().id(d,a).to_list()
        def eval(vect, (d,a)):
            return vect.val(d,a)
        return matrix(ZZ, [vect_paths_id(s[i]) for i in range(len(s))]), matrix(ZZ, [[eval(projections[d][a],s[i]) for i in range(len(s))] 
                                                                                     for d, a in self.cover_generators()])
    def _galois_group(self):
        r"""
        Return the galois group of the cover 
        """
        n = self.degree
        gap_list = ""
        for lab in self.iter_lab():
            perm = self._permutations[lab]
            perm_str = perm.cycle_string()
            gap_list += perm_str + ", "
        G = gap("Centralizer(SymmetricGroup(%s), "%n + "Group([" + gap_list + "]))") 
        return G
    def _characters(self):
        r"""
        We want to decompose homology space with the observation that it gives a representation of the galois group.
        
        RETURN
            - character: table of character, character[i][g] give the value of the i-th character on g
            g is given by a number
            - character_degree: list of degree of every character
            - g : Order of the group
            - perm : table s.t. perm[g] give the permutation associated to the group element g on the cover
            - n_cha : number of characters
        """
        G = self.galois_group
        G_order, T = gap.Order(G)._sage_(), gap.CharacterTable(G)
        irr_characters = gap.Irr(T)
        n_cha = len(irr_characters)
        character = [[irr_characters[i][j]._sage_() for j in range(1, G_order + 1)] for i in range(1, n_cha + 1)]
        character_degree = [gap.Degree(irr_characters[i])._sage_() for i in range(1, n_cha + 1)]
        gap_size_centralizers = gap.SizesCentralizers(T)
        gap_orders = gap.OrdersClassRepresentatives(T)
        def find_group_element(t):
            cursor = 1
            while gap.Order( gap.Centralizer(gap.SymmetricGroup(n), t) ) <> gap_size_centralizers[cursor] and gap.Order(t) <> gap_orders[cursor]:
                cursor += 1
            return cursor
        elements_group = gap.Elements(G)
        perm = [Permutation(str(elements_group[i]))*Permutations(self.degree).identity() for i in range(1, G_order + 1)] 
        #identity assures that the permutation has the right size
        return(character, character_degree, G_order, perm, n_cha)
    def _which_cycle(self, i_edge):
        c = self.cycles
        for k in range(len(c)):
            c_k = c[k]
            for s in range(len(c_k)):
                if c_k[s] == i_edge:
                    return k
        return None
    def _index_edge_cycle_orientated(self, d, a):
        r"""
        return index of the left edge, the right edge, depending on the orientation
        the departing and arriving edge
        """
        i, j = self.label_position(a)
        i_edge = self._index_edge(i, j, d) + 1
        if i == 1 and j == self.number_of_intervals[1] - 1:                                     #the edge with the biggest number
            i_next = self._index_edge(i, self.number_of_intervals[0], d) + 1
        else:
            i_next = i_edge + 1
        if self._intervals[i][j].orientation == 1:
            return self._which_cycle(i_edge), self._which_cycle(i_next)
        if self._intervals[i][j].orientation == -1:
            return self._which_cycle(i_next), self._which_cycle(i_edge)
    def _index_edge_cycle_non_orientated(self, d, a):
        i, j = self.label_position(a)
        i_left = self._index_edge(i, j, d) + 1
        if i == 1 and j == self.number_of_intervals[1] - 1:                                     #the edge with the biggest number
            i_right = self._index_edge(i, self.number_of_intervals[0], d) + 1
        else:
            i_right = i_left + 1
        return self._which_cycle(i_left), self._which_cycle(i_right)
    def _h_zero_sigma(self):
        r"""
        Return matrix of the border application
        """
        projected_vectors = self.canonical_VectPaths().copy(vector(ZZ, self.n_cycles))
        for d, a in self.cover_generators():
            i_depart_cycle, i_arrive_cycle = self._index_edge_cycle_orientated(d, a)
            projected_vectors.val(d,a)[i_arrive_cycle] += 1
            projected_vectors.val(d,a)[i_depart_cycle] -= 1
        return self.h_one_rel_to_generator*matrix(ZZ, [projected_vectors.val(d,a) for d, a in self.cover_generators()]) 
    def _h_one_projection(self):               #projection in h_one_rel on the kernel of the border application
        ker = self._h_zero_sigma().kernel() 
        d = ker.degree()
        r = ker.rank()
        def complete(basis):
            B = copy(basis)
            k = len(B) 
            c = 0
            def id(i):
                res = vector(ZZ, d)
                res[i] = 1
                return res
            while k < d:            
                if matrix(B + [id(c)]).rank() == k + 1:
                    B += [id(c)]
                    k += 1
                c += 1
            return matrix(B)
        B = complete(self._h_zero_sigma().kernel().basis())
        P = matrix(d)
        for i in range(r):
            P[i,i] = 1
        return B.inverse()*P*B
    def _h_one_rel_to_h_one(self):
        ker = self._h_zero_sigma().kernel() 
        d = ker.degree()
        r = ker.rank()
        def complete(basis):
            B = copy(basis)
            k = len(B) 
            c = 0
            def id(i):
                res = vector(ZZ, d)
                res[i] = 1
                return res
            while k < d:            
                if matrix(B + [id(c)]).rank() == k + 1:
                    B += [id(c)]
                    k += 1
                c += 1
            return matrix(B)
        B_c = complete(self._h_zero_sigma().kernel().basis())
        P = matrix(d,r)
        for i in range(r):
            P[i,i] = 1
        return B_c.inverse()*P
    def _h_one_to_h_one_rel(self):               #projection in h_one_rel on the kernel of the border application
        ker = self._h_zero_sigma().kernel()
        return matrix(ZZ, ker.dimension(), ker.degree(), ker.basis())
    def _representation_matrix(self):
        r"""
        Matrix from H1(X, Sigma) to H1(X, Sigma)
        """
        return [self.h_one_to_generator*matrix(ZZ, [[Integer(self.galois_group_permutation[t](n) == m and a == b) for n in self.degrees() for a in self.labels()] 
                            for m in self.degrees() for b in self.labels()])*self.generator_to_h_one
                for t in range(self.galois_group_order)]
    def _number_of_character_appeareance(self, ind_character):
        s = 0
        for t in range(self.galois_group_order):
            s += self.character_table[ind_character][t]*(self.representation_matrix[t].trace())
        return s/self.galois_group_order
    def dimension(self, i):
        s = 0
        for t in range(g):
            s += self.character[i][t]*self.perm[t].number_of_fixed_points()*self.degree
        return s/g            
    def _line_double(self,k):                                     #find a label in the line k which appears twice, return None if there is none
        line_k = self._intervals[k]
        n = len(line_k)
        for i in range(n):
            for j in range(i+1, n) :
                if line_k[i].label == line_k[j].label :
                    return j
        return(None)
    def _double(self):                                            #find a label which appears twice on the same line
        j_0 = self._line_double(0)
        if j_0 <> None:
            return (0, j_0)
        j_1 = self._line_double(1)
        if j_1 <> None:
            return (1, j_1)
        return (None, None)
    def _is_double(self, i, j):
        i_p, j_p = self._twins[i][j]
        return (i_p == i)
    def label(self, i, j):           # label of the i, j interval
        return(self._intervals[i][j].label)
    def label_position(self, lab):
        r"""
        Give position in polygone of the interval which will carry the name corresponding to the label and degree.
        """
        for j in range(self.number_of_intervals[0]):
            if self.label(0, j) == lab and self.give_name(0, j):
                return (0, j)
        for j in range(self.number_of_intervals[1]):
            if self.label(1, j) == lab  and self.give_name(1, j):
                return (1, j)
        return None
    def length(self, lab):
        return(self._lengths(lab))
    def permutation(self,i,j):
        r"""
        Returns the permutation associated to the label of the interval in line `i` and position `j`
        """
        return(self._permutations[self.label(i,j)])
    def min_length(self):
        r"""
        Gives the minimal length of all intervals
        """
        return(min(self._lengths.values()))
    def max_length(self):
        r"""
        Gives the maximal length of all intervals
        """
        return(max(self._lengths.values()))
    def nb_labels(self):
        return (len(self._lengths))
    def _sum(self,i,label):                                       #return(sum of lengths of intervals on a line, number of times label appears)
        lg = self._lengths
        line = self._intervals[i]
        n = len(line) 
        s, compt = 0, 0
        for k in range(n) :
            if line[k].label == label :
                compt += 1
            else :
                s += lg[line[k].label]
        return (s, compt)
    def rand_lg(self):
        r"""
        To get random length associated to your interval exchange. It modify directly your interval exchange
        and returns it. The total length of the two intervals is normalised to `1`.
        """
        lg = {self._intervals[i][j].label : R.random_element(0,1) for i in range(2) for j in range(self.number_of_intervals[i])}
        self._lengths = lg
        lab = lg.keys()
        i_double, j_double = self._double()
        if i_double == None:
            return(self)
        lab_double = self._intervals[i_double][j_double].label
        (s_0, l_0), (s_1, l_1) = self._sum(0,lab_double), self._sum(1,lab_double)
        l = (s_0 - s_1)/(l_1 - l_0)
        if l > 0:
            lg[lab_double] = l
        else:
            j_double = self._line_double((i_double + 1) % 2)
            if j_double == None:
                print self
                raise NameError('No double')
            lab_double = self._intervals[(i_double + 1) % 2][j_double].label
            (s_0, l_0), (s_1, l_1) = self._sum(0,lab_double), self._sum(1,lab_double)
            l = (s_0 - s_1)/(l_1 - l_0)
            if l <= 0:
                print self
                raise NameError('Valeur neg')
            lg[lab_double] = l
        return(self)
    def diff_sum(self):
        l_0, l_1 = self._intervals[0], self._intervals[1]
        lg = self._lengths
        s0 = sum([lg[l_0[i].label] for i in range(len(l_0))])
        s1 = sum([lg[l_1[i].label] for i in range(len(l_1))])
        return(abs(s0-s1))
    def normalise(self):
        r"""
        As the length goes to zero with the Rauzy transform, you sometime have to renormalise your total length to 1 in
        order to have non zero length (there is some troubles with computer approximation).
        This normalise the length of both intervals to one.
        """
        l_0, l_1 = self._intervals[0], self._intervals[1]
        lg = self._lengths
        s0 = sum([lg[l_0[i].label] for i in range(len(l_0))])
        s1 = sum([lg[l_1[i].label] for i in range(len(l_1))])
        for a in lg.iterkeys():
            lg[a] = 2*lg[a]/(s0+s1)
        if abs(s0 - s1) > diff_sum_seuil:
            (i,j) = self._double()
            label = self._intervals[i][j].label
            (s_0, l_0) = self._sum(0,label)
            (s_1, l_1) = self._sum(1,label)
            if l_0 <> l_1 :
                val = (s_0 - s_1)/(l_1 - l_0)
                if val > 0 :
                    lg[label] = val
                else :
                    print (self, val)
                    raise NameError('Valeur neg')
        return (s0)
    def rauzy_type(self) :
        r"""
        To apply Rauzy transform, you need to know which one of the two last intervals in the lines is the shortest.
        This returns the couple (lign where the last interval is the shortest, lign where it is the longest)

        EXAMPLES::
            sage: intervals = [[Interval('a',1), Interval('a',-1), Interval('b',1)],
            ...   [Interval('b',1), Interval('c', 1), Interval('c',-1)]]

            sage: lengths = {'a': 2, 'b': 3, 'c': 2}

            sage: permutations = {lab: Permutations(1).identity() for lab in lengths.iterkeys()}
        
            sage: IntExchange(intervals, lengths, permutations).rauzy_type()
            (1, 0)

            sage: lengths['b'] = 1
            sage: IntExchange(intervals, lengths, permutations).rauzy_type()
            (0, 1)

        """
        line_0, line_1 = self._intervals[0], self._intervals[1]
        a, b = line_0[len(line_0)-1].label, line_1[len(line_1)-1].label
        if self._lengths[a] <= self._lengths[b] :
            return (0,1)
        else :
            return (1,0)
#    def ident(self, i, j):                                             #return interval identified and change of orientation
#        line = self._intervals[i]
#        label = line[j].label
#        for k in range(len(line)):
#            if k<>j and line[k].label == label:                        #the third parameter is 1 if the orientation is changed
#                return (i,k,-1)
#        other_line = self._intervals[(i+1) % 2]
#        for k in range(len(other_line)):
#            if other_line[k].label == label:                           #it is 0 otherwise
#                return ((i+1)%2,k,1)
    def give_name(self, i, j) :
        r"""
        Answers if the interval correspond to the given convention and take the number of the corresponding copy.
        See IntExchange for more information about the convention.
        """
        if (self._intervals[i][j].orientation == 1 and i == 0) or (self._intervals[i][j].orientation == -1 and i == 1):
            return True
        else :
            return False
    def move(self, i_init, j_init, i_dest, j_dest, new_orientation):
        r"""
        Move the interval from position `i_{init}, j_{init}` to `i_{dest}, j_{dest}` with new_orientation as orientation

        EXAMPLES::
            sage: intervals = [[Interval('a',1), Interval('a',-1), Interval('b',1)],
            ...   [Interval('b',1), Interval('c', 1), Interval('c',-1)]]

            sage: lengths = {'a': 2, 'b': 3, 'c': 2}

            sage: permutations = {lab: Permutations(1).identity() for lab in lengths.iterkeys()}
        
            sage: ie = IntExchange(intervals, lengths, permutations)
            sage: ie
            Intervals    : [ a, -a,  b]
                           [ b,  c, -c]
            Lengths      : {'a': 2, 'c': 2, 'b': 3}
            Permutations : {'a': [1], 'c': [1], 'b': [1]}
            <BLANKLINE>

            sage: ie.move(0,1,0,0,-1)
            sage: ie
            Intervals    : [-a,  a,  b]
                           [ b,  c, -c]
            Lengths      : {'a': 2, 'c': 2, 'b': 3}
            Permutations : {'a': [1], 'c': [1], 'b': [1]}
            <BLANKLINE>

        """
        if i_init == i_dest:
            if j_init >= j_dest:
                self._intervals[i_dest].insert(j_dest, Interval(self._intervals[i_init][j_init].label, new_orientation))
                del self._intervals[i_init][j_init + 1]
            else:
                self._intervals[i_dest].insert(j_dest + 1, Interval(self._intervals[i_init][j_init].label, new_orientation))
                del self._intervals[i_init][j_init]
        else:
            self._intervals[i_dest].insert(j_dest, Interval(self._intervals[i_init][j_init].label, new_orientation))
            del self._intervals[i_init][j_init]
            self.number_of_intervals[i_init] -= 1
            self.number_of_intervals[i_dest] += 1
        self._twins = [[self._twin(i, j) for j in range(self.number_of_intervals[i])] for i in range(2)]
    def rauzy(self):
        r"""
        Apply the rauzy transform to the interval exchange.
        
        OUTPUT::

        EXAMPLES::
            sage: intervals = [[Interval('a',1), Interval('a',-1), Interval('b',1)],
            ...   [Interval('b',1), Interval('c', 1), Interval('c',-1)]]

            sage: lengths = {'a': 2, 'b': 3, 'c': 2}

            sage: permutations = {lab: Permutations(1).identity() for lab in lengths.iterkeys()}
        
            sage: ie = IntExchange(intervals, lengths, permutations)
            sage: ie
            Intervals    : [ a, -a,  b]
                           [ b,  c, -c]
            Lengths      : {'a': 2, 'c': 2, 'b': 3}
            Permutations : {'a': [1], 'c': [1], 'b': [1]}
            <BLANKLINE>


            sage: ie.rauzy()
            (-1, 'c', 'b')

            sage: ie
            Intervals    : [ a, -a,  b]
                           [ b, -c,  c]
            Lengths      : {'a': 2, 'c': 2, 'b': 1}
            Permutations : {'a': [1], 'c': [1], 'b': [1]}
            <BLANKLINE>

            sage: ie.rauzy()
            (1, 'b', 'c')

            sage: ie
            Intervals    : [ a, -a]
                           [ b, -b, -c,  c]
            Lengths      : {'a': 2, 'c': 1, 'b': 1}
            Permutations : {'a': [1], 'c': [1], 'b': [1]}
            <BLANKLINE>
    
        """
        (i_0, i_1) = self.rauzy_type()
        d_0, d_1 = self.number_of_intervals[i_0] - 1, self.number_of_intervals[i_1] - 1
        A = self._intervals[i_0][d_0]                                   #A is the shortest of both final intervals
        B = self._intervals[i_1][d_1]                                   #B is the longest
        self._lengths[B.label] = self._lengths[B.label] - self._lengths[A.label]                       #CHANGE LENGTH
        (i_1_twin,d_1_twin) = self._twins[i_1][d_1]                     #CHANGE PERMUTATION #compute image of the long interval to see first return
        if i_1 <> i_1_twin:
            new_orientation =  A.orientation
            dj = 1
        else:
            new_orientation = -A.orientation
            dj = 0                                 #the way we insert depends on the change of orientation #dj = 0 if the orientation changed 1 otherwise
        self.move(i_0, d_0, i_1_twin, d_1_twin + dj, new_orientation)
        return (A.orientation*B.orientation, A.label, B.label)
    def name(self,i, j) :
        r"""
        Give permutation which associates number of the copy from which interval takes his name.
        """
        if self.give_name(i,j) :
            return Permutations(self.degree).identity()
        else : 
            return self.permutation(i,j).inverse()
    def ident_rev(self,i,j):
        r"""
        Return the permutation for the given saddle connection of identifiction according to the level of the cover. 
        """
        if self.give_name(i,j):
            return self.permutation(i,j)
        else :
            return self.permutation(i,j).inverse()
    def rauzy_rev(self):
        
        r"""
        Apply the rauzy transform to the interval exchange.

        EXAMPLES::
        
            sage: intervals = [[Interval('a',1), Interval('a',-1), Interval('b',1)],
            ...   [Interval('b',1), Interval('c', 1), Interval('c',-1)]]

            sage: lengths = {'a': 2, 'b': 3, 'c': 2}

            sage: permutations = {lab: Permutations(1).identity() for lab in lengths.iterkeys()}
        
            sage: ie = IntExchange(intervals, lengths, permutations)
            sage: ie
            Intervals    : [ a, -a,  b]
                           [ b,  c, -c]
            Lengths      : {'a': 2, 'c': 2, 'b': 3}
            Permutations : {'a': [1], 'c': [1], 'b': [1]}
            <BLANKLINE>

            sage: ie.rauzy_rev()
            (-1, 'c', 'b', [1], [1], [1], True, [1])

            sage: ie
            Intervals    : [ a, -a,  b]
                           [ b, -c,  c]
            Lengths      : {'a': 2, 'c': 2, 'b': 1}
            Permutations : {'a': [1], 'c': [1], 'b': [1]}
            <BLANKLINE>
            sage: ie.rauzy_rev()
            (1, 'b', 'c', [1], [1], [1], True, [1])

            sage: ie
            Intervals    : [ a, -a]
                           [ b, -b, -c,  c]
            Lengths      : {'a': 2, 'c': 1, 'b': 1}
            Permutations : {'a': [1], 'c': [1], 'b': [1]}
            <BLANKLINE>

        
        Non trivial cover case
        ::

            sage: sigma_a = Permutation('(1,2)(3,4)(5,6)')
            sage: sigma_b = Permutation('(1,3)(2,4)(5)(6)')
            sage: sigma_c = Permutation('(1,6)(2,5)(3,4)')

            sage: permutations = {'a': sigma_a, 'b': sigma_b, 'c': sigma_c}
            sage: ie = IntExchange(intervals, lengths, permutations)

            sage: ie.rauzy_rev()
            (-1, 'c', 'b', [1, 2, 3, 4, 5, 6], [1, 2, 3, 4, 5, 6], [6, 5, 4, 3, 2, 1], True, [3, 4, 1, 2, 5, 6])
        
            sage: ie
            Intervals    : [ a, -a,  b]
                           [ b, -c,  c]
            Lengths      : {'a': 2, 'c': 2, 'b': 1}
            Permutations : {'a': [2, 1, 4, 3, 6, 5], 'c': [4, 3, 6, 5, 2, 1], 'b': [3, 4, 1, 2, 5, 6]}

        """
        (i_0, i_1) = self.rauzy_type()
        d_0, d_1 = self.number_of_intervals[i_0] - 1, self.number_of_intervals[i_1] - 1
        A, B = self._intervals[i_0][d_0], self._intervals[i_1][d_1]        #A is the shortest of both final intervals B the longest
        self._lengths[B.label] = self._lengths[B.label] - self._lengths[A.label]                       #CHANGE LENGTH
        (i_1_twin, d_1_twin) = self._twins[i_1][d_1]                        #CHANGE PERMUTATION #compute image of the long interval to see first return
        (i_0_twin, d_0_twin) = self._twins[i_0][d_0]                        #compute image of the short interval 
        name_A, name_A_twin, name_B = self.name(i_0, d_0), self.name(i_0_twin,d_0_twin), self.name(i_1,d_1)
        ident_A, ident_A_twin, ident_B = self.ident_rev(i_0,d_0), self.ident_rev(i_0_twin,d_0_twin), self.ident_rev(i_1,d_1)
        if i_1 <> i_1_twin:
            new_orientation =  A.orientation
            dj = 1
        else:
            new_orientation = -A.orientation
            dj = 0                                                          #the way we insert depends on the change of orientation
        self.move(i_0, d_0, i_1_twin, d_1_twin + dj, new_orientation)
        if self.give_name(i_1_twin, d_1_twin+dj) :                          #CHANGE PERMUTATIONS, if the image of A is the one who gives the name
            self._permutations[A.label] = ident_B.inverse()*ident_A         #beware, the sense of composition is inversed !!
        else :
            self._permutations[A.label] = ident_A_twin*ident_B
        return (A.orientation*B.orientation, A.label, B.label, name_A, name_B, name_A_twin, self.give_name(i_1_twin,d_1_twin+dj), ident_B)
    def _type(self, i, j):
        i_twin, _ = self._twins[i][j]
        if i_twin == i:
            return 0
        else:
            return 1
    def edge_permutation(self):
        r"""
        We number edges of the upper part of the polygon, and return the permutation of the edge you
        get while winding around it. It is numbered from 0 to n, the n-1 first edges numbered i correspond to the left edge of the
        i-th interval, the n-th is the rightmost one.
        """
        perm = range(n+1)
        for j in range (n+1):
            if j < self.number_of_intervals[0] :
                i_twin, j_twin = self._twins[0][j]
            else:
                i_twin, j_twin = self._twins[1][len(self._intervals[1])-1]
            if i_twin == 0:
                perm[j] = j_twin + 1
            else:
                while i_twin <> 0 and j_twin <> 0:
                    i_twin, j_twin = self._twins[i_twin][j_twin-1]
                if j_twin == 0 and i_twin == 1:
                    perm[j] = 0
                else:
                    perm[j] = j_twin + 1
        return (Permutation([perm[i] + 1 for i in range(n+1)]))
    def twin_cover(self, i, j, d):
        if j < self.number_of_intervals[i]:
            (i_t, j_t), s = self._twins[i][j], self.ident_rev(i,j)(d)
        if j == self.number_of_intervals[i]:
            i_t, j_t = self._twins[i][j-1]
            j_t += 1
            s = self.ident_rev(i, j-1)(d)
        return(i_t, j_t, s)
    def _index_edge(self, i, j, d):
        """
        Give the number given to the edge at the left of the interval at i, j in level c
        """
        n = self.total_number_of_intervals_base
        if j == 0:
            return((d-1)*self.total_number_of_intervals_base)
        if j == self.number_of_intervals[i]:
            return((d-1)*self.total_number_of_intervals_base + self.number_of_intervals[0])
        return( (d-1)*self.total_number_of_intervals_base + i*self.number_of_intervals[0] + j)
    def _next_edge(self, i, j, d):
        r"""
        return next edge if you go counterclockwise around the point at the left of i,j,c interval
        """
        if i == 0:
            i_t, j_t, d_t = self.twin_cover(0, j, d)
            j_t += 1 - i_t
        if i == 1:
            i_t, j_t, d_t = self.twin_cover(1, j-1, d)
            j_t += 1 - i_t
        return self._index_edge(i_t, j_t, d_t)
    def _next_path(self, i, j, d, side):
        r"""
        return next intersected path if you wind counterclockwise
        """
        if i == 0 and side == "Right":
            if j < self.number_of_intervals[0] - 1:
                return 0, j + 1, d, "Left"
            if j == self.number_of_intervals[0] - 1:
                return 1, self.number_of_intervals[1] - 1, d, "Right"
        if i == 0 and side == "Left":
            i_t, j_t, d_t = self.twin_cover(i, j, d)
            if i_t == 0:
                return i_t, j_t, d_t, "Right"
            if i_t == 1:
                return i_t, j_t, d_t, "Left"
        if i == 1 and side == "Right":
            i_t, j_t, d_t = self.twin_cover(i, j, d)
            if i_t == 0:
                return i_t, j_t, d_t, "Right"
            if i_t == 1:
                return i_t, j_t, d_t, "Left"
        if i == 1 and side == "Left":
            if j == 0:
                return 0, 0, d, "Left"
            if j > 0:
                return 1, j - 1, d, "Right"
    def _next_path_generator(self, d, a, side):
        r"""
        Return next path with format degree and label, st the end point of orientated path is the edge
        """
        i, j = self.label_position(a)
        i_aux, j_aux, d_aux, side_aux = self._next_path(i, j, d, side)
        while not self.give_name(i_aux, j_aux):
            i_aux, j_aux, d_aux, side_aux = self._next_path(i_aux, j_aux, d_aux, side_aux)
        return d_aux, self.label(i_aux, j_aux), side_aux
    def edge_cover_permutation(self):
        perm = range(self.total_number_of_intervals_cover)
        for c in self.degrees():
            for j in range(self.number_of_intervals[0]):
                perm[self._index_edge(0,j,c)] = self._next_edge(0, j, c) + 1    #index must be greater than 1 for permutations
            for j in range(1, self.number_of_intervals[1] + 1):
                perm[self._index_edge(1,j,c)] = self._next_edge(1, j, c) + 1
        return (Permutation(perm))
    def _is_extremal(self, i):
        l = [len(self._intervals[0]), len(self._intervals[1])]
        n = l[0] + l[1]
        i_base = (i-1) % n
        return( i_base == 0 or i_base == l[0] )
    def genus(self):
        sum_d = 0
        for k in range(self.n_cycles):
            c = self.cycles[k]
            d = 0
            for i in range(len(c)):
                if not self._is_extremal(c[i]):
                    d += 1
            sum_d += d - 2
        return((sum_d + 4)/4)
    def couple_paths(self, v, i_edge_cycle_ref):
        arrive, depart = [], []
        for d, a in self.cover_generators():
            if v.val(d,a) <> 0:
                i_left_cycle, i_right_cycle = self._index_edge_cycle_non_orientated(d, a)
                i_depart_cycle, i_arrive_cycle = self._index_edge_cycle_orientated(d, a)
                if v.val(d,a) < 0:
                    i_depart_cycle, i_arrive_cycle = i_arrive_cycle, i_depart_cycle
                if i_depart_cycle == i_edge_cycle_ref:
                    if i_left_cycle == i_right_cycle:
                        depart += [(d, a, "Left")]*Integer(abs(v.val(d,a)))
                    else:
                        if i_left_cycle == i_edge_cycle_ref:
                            depart += [(d, a, "Left")]*Integer(abs(v.val(d,a)))
                        if i_right_cycle == i_edge_cycle_ref:
                            depart += [(d, a, "Right")]*Integer(abs(v.val(d,a)))
                if i_arrive_cycle == i_edge_cycle_ref:
                    if i_left_cycle == i_right_cycle:
                        arrive += [(d, a, "Right")]*Integer(abs(v.val(d,a)))
                    else:
                        if i_left_cycle == i_edge_cycle_ref:
                            arrive += [(d, a, "Left")]*Integer(abs(v.val(d,a)))
                        if i_right_cycle == i_edge_cycle_ref:
                            arrive += [(d, a, "Right")]*Integer(abs(v.val(d,a)))
        return arrive, depart
    def interval_to_name(self, i, j, d):
        lab = self._intervals[i][j].label
        orientation = self._intervals[i][j].orientation
        deg = self.name(i,j)(d)
        return(deg, lab)
    def label_orientation(self, a):
        i, j = self.label_position(a)
        return self._intervals[i][j].orientation
    def is_cycle(self, d, a):
        left, right = self._index_edge_cycle_non_orientated(d, a)
        return left == right
    def intersection(self, (d_arr, a_arr, side_arr), (d_dep, a_dep, side_dep), v):
        r"""
        return the number of intersection of the intervals arriving to an edge and departing from it, with the intervals
        given by v.

        EXAMPLES::
            sage: intervals = [[Interval('a',1), Interval('a',-1), Interval('b',1)],
            ...   [Interval('b',1), Interval('c', 1), Interval('c',-1)]]
            sage: ie = orientable_double_cover_IntExchange(intervals)

            sage: v = ie.canonical_VectPaths()
            sage: v._vect[1]['a'] = 1
            sage: v._vect[2]['a'] = -1
            sage: ie.intersection((1, 'b'), (2, 'b'), v)
        """
        d_aux, a_aux, side_aux = self._next_path_generator(d_arr, a_arr, side_arr)
        c = 0
        s = 0
        #print (d_arr, a_arr, side_arr), (d_dep, a_dep, side_dep)
        while (d_aux, a_aux, side_aux) <> (d_dep, a_dep, side_dep):
            #print (d_aux, a_aux, side_aux)
#            if not self.is_cycle(d_aux, a_aux):
            if side_aux == "Right":
                c += v.val(d_aux, a_aux)*self.label_orientation(a_aux)
            if side_aux == "Left":
                c -= v.val(d_aux, a_aux)*self.label_orientation(a_aux)
            d_aux, a_aux, side_aux = self._next_path_generator(d_aux, a_aux, side_aux)
            s += 1
        #print "c: %s"%c
        return c
    def global_intersection(self, v, w):
        r"""
        give intersection of VectPath from the canonical generator family
        """
        res = 0
        for i in range(self.n_cycles):
            #print "cycle : %s"%i
            arrive, depart = self.couple_paths(v, i)
            #print arrive, depart, v, w
            if len(arrive) <> len(depart):
                print arrive, depart
                print v
                raise NameError('Not an element from kernel of border application')
            for j in range(len(arrive)):
                res += self.intersection(arrive[j], depart[j], w)
        return res
    def h_one_intersection_matrix(self):
        dim = self.h_one_to_generator.nrows()
        def inter(i, j):
            v_i = self.canonical_VectPaths().copy_vector(self.h_one_to_generator[i])
            v_j = self.canonical_VectPaths().copy_vector(self.h_one_to_generator[j])
            return self.global_intersection(v_i, v_j)
        return matrix(dim, inter)
    def isotopic_projection(self, i_character, v):
        r"""
        Return a new VectPaths corresponding to the projection of self to the isotropic space of the \chi_i character
        """
        res = self.canonical_VectPaths()
        for d, a in self.cover_generators():
            somme = 0
            for t in range(self.galois_group_order):
                somme += self.character_table[i_character][t]*v.val(self.galois_group_permutation[t].inverse()(d), a)
            res._vect[d][a] = self.character_degree[i_character]/self.galois_group_order*somme
        return(res)
    def isotopic_projection_matrix(self, i_character):
        return self.h_one_to_generator*matrix([self.isotopic_projection(i_character, self.canonical_VectPaths().id(d, a)).to_list() 
                                               for d, a in self.cover_generators()])*self.generator_to_h_one
    def intersection_isotopic(self):
        return [self.isotopic_projection_matrix(k)*self.intersection_matrix*self.isotopic_projection_matrix(k).transpose() 
                                      for k in range(self.n_characters)]
    def signature(self, M):
        p, q = 0, 0
        eigenvalues = M.eigenvalues()
        for i in range(len(eigenvalues)):
            ev = eigenvalues[i]*I/2
            if ev > 0:
                p += 1
            if ev < 0:
                q += 1
        return(p,q)


def trivial_IntExchange(intervals):
    return IntExchange(intervals, None, {intervals[i][j].label: Permutations(1).identity() for i in range(2) for j in range(len(intervals[i])) })

def orientable_double_cover_IntExchange(intervals):
    r"""
    Return the interval exchange of the orientated cover
    EXEMPLES ::
        sage: permutations = {'a': sigma_a, 'b': sigma_b, 'c': sigma_c}
        sage: orientable_double_cover_IntExchange(intervals)
        Intervals    : [ a, -a,  b]
                       [ b,  c, -c]
        Lengths      : None
        Permutations : {'a': [2, 1], 'c': [2, 1], 'b': [1, 2]}
        <BLANKLINE>

    """
    def twin(intervals, i, j):
        for n in range(len(intervals[i])):
            if n <> j and intervals[i][n].label == intervals[i][j].label:
                return(i,n)
        for n in range(len(intervals[(i+1)%2])):
            if intervals[(i+1)%2][n].label == intervals[i][j].label:
                return((i+1)%2,n)
    def orientable_double_cover_permutation(i, j):
        (i_twin, j_twin) = twin(intervals, i, j)
        if i_twin <> i :
            return Permutations(2).identity()
        else :
            return Permutation('(1,2)')
    return IntExchange(intervals, None, {intervals[i][j].label: orientable_double_cover_permutation(i, j) for i in range(2) for j in range(len(intervals[i])) })
