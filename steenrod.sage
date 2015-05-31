import itertools

A = SteenrodAlgebra(2,'adem')
A0 = sage.algebras.steenrod.steenrod_algebra.AA(0)
A1 = sage.algebras.steenrod.steenrod_algebra.AA(1)
A2 = sage.algebras.steenrod.steenrod_algebra.AA(2)
A3 = sage.algebras.steenrod.steenrod_algebra.AA(3)
A4 = sage.algebras.steenrod.steenrod_algebra.AA(4)
F2 = FiniteField(2)
Aw = SteenrodAlgebra(2,'woodz')
Ap = SteenrodAlgebra(2,'pst')
Apz = SteenrodAlgebra(2,'pst_z')
p01 = Apz.monomial(((0,1),))
p02 = Apz.monomial(((0,2),))
p03 = Apz.monomial(((0,3),))
p04 = Apz.monomial(((0,4),))
p11 = Apz.monomial(((1,1),))
p12 = Apz.monomial(((1,2),))
p13 = Apz.monomial(((1,3),))
p14 = Apz.monomial(((1,4),))
p21 = Apz.monomial(((2,1),))
p22 = Apz.monomial(((2,2),))
p23 = Apz.monomial(((2,3),))
p24 = Apz.monomial(((2,4),))
p31 = Apz.monomial(((3,1),))
p32 = Apz.monomial(((3,2),))
p33 = Apz.monomial(((3,3),))
p34 = Apz.monomial(((3,4),))
p41 = Apz.monomial(((4,1),))
p42 = Apz.monomial(((4,2),))
p43 = Apz.monomial(((4,3),))
p44 = Apz.monomial(((4,4),))

def ad(s):
    return s.change_basis('adem')

def mil(s):
    return s.change_basis('milnor')

def wall(s):
    return s.change_basis('wall_long')

def woody(s):
    return s.change_basis('wood_y')

def woodz(s):
    return s.change_basis('wood_z')

def pst(s):
    return s.change_basis('pst')

def powerset(iterable,start=0):
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return itertools.chain.from_iterable(itertools.combinations(s, r) for r in range(start, len(s)+1))

def alist_old(n,k=-1):
    basis=map(ad,list(sage.algebras.steenrod.steenrod_algebra.AA(n).basis()))
    if k<0:
        return basis
    else:
        s=sage.algebras.steenrod.steenrod_algebra.AA(k).top_class()
        basis2=[[] for i in range(sage.algebras.steenrod.steenrod_algebra.AA(n).top_class().degree())]

        # remove monomials that are zero in A(n)//A(k).  Note this doesn't currently remove polynomials that are zero, only if they're one term.
        keys=set()
        for t in basis:
            if t*s != 0:
                t2=sum(i for i in t.monomials() if i*s != 0)
                if t2 != 0 and t2*s not in keys:
                    basis2[t2.degree()].append(t2)
                    keys.add(t2*s)

        #reduce each degree of the basis - should make this just row-reduce
        for deg in basis2:
            deg.sort(key=len)
            #first just try to make elements simpler
            for i in range(len(deg)):
                for j in range(len(deg)):
                    if j!=i and len(deg[i]+deg[j])<len(deg[j]):
                        deg[j]+=deg[i]
            #then delete linearly dependent elements
            rem=set([0])
            while rem:
                for i in rem:
                    deg[:]=filter(lambda x: x!=i, deg)
                rem=set()
                for subset in powerset(deg,2):
                    if sum(subset)*s==0:
                        rem.add(subset[-1])
        return basis2

def alist(n,k=-1,module='right'):
    An=sage.algebras.steenrod.steenrod_algebra.AA(n)
    if k<0:
        return map(ad,list(An.basis()))
    else:
        s=sage.algebras.steenrod.steenrod_algebra.AA(k).top_class()
        basis2=[An.basis(i) for i in range(An.top_class().degree()+1)]
        if module=='right':
            def f(x): return x*s
        elif module=='left':
            def f(x): return s*x
        else:
            def f(x): return x
        return [reduceBasis(b,f) for b in basis2]

def head(i):
    return i.__iter__().next()

def generalizedInverse(m):
    i,p,q = m.smith_form()
    return q*i.T*p

def reduceBasis(basis,f,A=A):
    '''aka the coimage'''
    try:
        basis = map(A,basis)
    except:
        print basis
        raise
    monomials = set(map(head,itertools.chain.from_iterable(basis)))
    monomialsDict = {}
    monomialsTimesS = set()
    for m in monomials:
        r=set(map(head,f(A.monomial(m))))
        if r:
            monomialsDict[m]=r
            for i in r:
                monomialsTimesS.add(i)
    monomials = list(monomialsDict)
    monomialsTimesS = list(monomialsTimesS)
    # print monomials
    # print monomialsTimesS
    Basis = matrix(F2,[ [ 1 if m in map(head,b) else 0 for m in monomials] for b in basis])
    # print Basis
    Change = matrix(F2,[ [ 1 if n in monomialsDict[m] else 0 for n in monomialsTimesS] for m in monomials])
    # print Change
    NewBasis = ((Basis*Change).echelon_form()*generalizedInverse(Change)).echelon_form()
    # NewBasis = (Basis*Change).row_space().basis_matrix() #generates uglier representatives
    # print NewBasis
    newBasis = [ sum([A.monomial(monomials[i]) for i in range(len(monomials)) if b[i]==1]) for b in NewBasis.rows()]
    newBasis = filter(lambda x: x!=0, newBasis)
    # print newBasis
    return newBasis

def kernelBasis(basis,f,A=A):
    try:
        basis = map(A,basis)
    except:
        print basis
        raise
    monomials = set(map(head,itertools.chain.from_iterable(basis)))
    monomialsDict = {}
    monomialsTimesS = set()
    for m in monomials:
        r=set(map(head,f(A.monomial(m))))
        monomialsDict[m]=r
        for i in r:
            monomialsTimesS.add(i)
    monomials = list(monomialsDict)
    monomialsTimesS = list(monomialsTimesS)
    Basis = matrix(F2,[ [ 1 if m in map(head,b) else 0 for m in monomials] for b in basis])
    Change = matrix(F2,[ [ 1 if n in monomialsDict[m] else 0 for n in monomialsTimesS] for m in monomials])
    Kern=(Basis*Change).kernel().basis_matrix()
    kernBasis = [ sum([basis[i] for i in range(len(basis)) if b[i]==1]) for b in Kern.rows()]
    return kernBasis

def shorten(t,l=1):
    if t in ZZ:
        return t
    else:
        return sum([m for m in t.monomials() if len(list(m)[0][0])<=l])

def sqpr(*l,**kw):
    if 'A' not in kw:
        B=A
    else:
        B=kw['A']
    return prod((B.monomial((x,)) for x in l),B(1))
