F.<x>=GF(2)[]
F_2.<mu> = GF(2)[]
F_4.<mu> =GF(2^2, name='mu', modulus=mu^2+mu+1)
Omega1.<x> = F_4[]

def d(f):
    return (f(x + mu) - f(x))*mu^(-1)
            
def partial_1(f):
    def helper(coeff):
        if coeff in [0, 1]:
            return coeff
        elif coeff == mu + 1:
            return 1
        elif coeff == mu:
            return 0
    l = list(f)
    return sum([helper(l[i]) * x^i for i in range(len(l))])
        
def partial_mu(f):
    def helper(coeff):
        if coeff in [mu, mu+1]:
            return 1
        elif coeff in [0, 1]:
            return 0
    l = list(f)
    return sum([helper(l[i]) * x^i for i in range(len(l))])

del_1 = lambda f: partial_1(d(f))
del_mu = lambda f: partial_mu(d(f))
import itertools
from multiprocessing import Pool
from functools import partial

def d0(f):
    return (del_1(f), del_mu(f))
    
def kerd1_filter(f4):
    f,g = f4
    if del_1(g) == del_mu(f):
        return f4
    else:
        return None

def get_H1(num_proc, deg, sub_om1 = None, sub_imd0=None, kerd1=None):
    pool = Pool(processes=num_proc)
    
    def get_Omega1(deg):
        polys = list(F.polynomials(max_degree=deg))
        return list(itertools.product(polys, polys))
        
    def get_Imd0(deg):
        polys = list(F.polynomials(max_degree=deg))
        ret = pool.map(d0, polys)
        return set(ret)
    
    if sub_om1 is None:
        print  "Getting a ", deg-2," degree subspace of Omega^1..."
        sub_om1 = get_Omega1(deg-2)
    print "Got Omega^1"
    if sub_imd0 is None:
        print "Getting a ", deg," degree subspace of Im(d^0)..."
        sub_imd0 = list(get_Imd0(deg))
    print "Got Im(d^0)"
    def get_kerd1():

        kerd1 = pool.map(kerd1_filter, sub_om1)
        kerd1 = filter(lambda f: f!=None, kerd1)
        return set(kerd1)
        #kerd1 = []
        #for f,g in sub_om1:
        #    if del_1(g) == del_mu(f):
        #        if (f,g) in sub_imd0:
        #            str = " Imd0"
        #        else:
        #            str = ""
        #        kerd1.append((f,g))
        #        print "(", f, ")*1 + (", g, ")*mu", str    
        #return kerd1
    
    if kerd1 is None:
        print "Getting a subset of Ker(d^1)"
        kerd1 = get_kerd1()
    print "Got ker(d^1)"

    classes = [sub_imd0]
    remainders = list(set(kerd1) - set(sub_imd0))
    while len(remainders) > 0:
        f1,fmu = remainders[0]
        print "looking for things equiv to: ", (f1, fmu)
        f_equivs = pool.map(partial(get_equiv_class, f1, fmu, sub_imd0), remainders)
        f_equivs = list(set(filter(lambda f:f!=None, f_equivs)))
        classes.append(f_equivs)
        if len(f_equivs) == 1:
            print "can't find any more ", len(remainders), " left..."
            return classes
        remainders = list(set(remainders)-set(f_equivs))
        print "found ", len(f_equivs), " in new class.  now looking through next ", len(remainders)
    return classes
    
def get_equiv_class(f1, fmu, sub_imd0, g):
    g1,gmu = g
    if (f1-g1, fmu-gmu) in sub_imd0:
        return g
    else:
        return None

h1_0 = get_H1(4,4)
