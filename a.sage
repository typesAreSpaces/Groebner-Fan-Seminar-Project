import sage.all
import copy
import numpy as np
from scipy.optimize import linprog

debug = False

R.<x,y,z> = PolynomialRing(QQ)

# Given two monomials x^vec(a), x^vec(b)
# returns vec(a) - vec(b) as a tuple
def subtractExponents(monomial1, monomial2):
    exponents1 = monomial1.exponents(false)[0]
    exponents2 = monomial2.exponents(false)[0]
    return tuple(np.subtract(exponents1, exponents2))

# 'G' is an array of polynomials representing
# the Groebner basis
# 'Psi' is a map fro polynomials to initial
# polynomials
# 'h' is a polynomial
def minimalPolynomialCheck(G, Psi, h):
    if debug:
        print "G: ", G
        print "Psi: ", Psi
        print "h: ", h
    for t in h.monomials():
        for g in G:
            q, r = t.quo_rem(Psi[g])
            if r == 0:
                if debug:
                    print "t: ", t
                    print "Psi[g]: ", Psi[g], " g: ", g
                hCoeffT = h.monomial_coefficient(t)
                gCoeffPsiG = g.monomial_coefficient(Psi[g])
                return True, hCoeffT/gCoeffPsiG*t//Psi[g]*g
    return False, None

# 'inequalities' is an array of vectors
def isNotEmptyTO(inequalities):
    epsilon = 0.0001
    numInequalities = len(inequalities)
    numVariables = len(inequalities[0])
    # The upper bounds are epsilon values
    # so we can express strict inequalities
    zeros = [-epsilon for x in range(numInequalities)]
    ones = [1 for x in range(numVariables)]
    # Solutions cannot contain zeros
    result = linprog(ones, A_ub=inequalities, b_ub=zeros, bounds=(epsilon, None))
    if debug:
        print "Inequalities: ", inequalities
        print "Result x: ", result.x
        print "Status: ", result.status
    return result.success

# 'inputBasis' is an array of polynomials representing
# the Groebner basis
def groebnerFan(inputBasis):
    # Initialization
    L = ([], [], [], {}, [])
    Lnew = [L]
    for polynomial in inputBasis:
        Lold = Lnew
        Lnew = []
        for (G, M, E, Psi, B) in Lold:
            for leadingMonomial in polynomial.monomials():
                Gnew = G[:]
                Gnew.append(polynomial)
                Mnew = M[:]
                Mnew.append(leadingMonomial)
                Enew = E[:]
                for nonLeadingMonomial in polynomial.monomials():
                    if(nonLeadingMonomial != leadingMonomial):
                        # We substract the Leading Monomial to the
                        # Non Leading Monomial because the LP solver
                        # has <= as default inequalities
                        Enew.append(subtractExponents(nonLeadingMonomial, leadingMonomial))
                Psinew = copy.deepcopy(Psi)
                Psinew[polynomial] = leadingMonomial
                Bnew = B[:]
                for g in G:
                    Bnew.append((g, polynomial))
                if isNotEmptyTO(Enew):
                    L = (Gnew, Mnew, Enew, Psinew, Bnew)
                    Lnew.append(L)
    # Computation of the Groebner Bases
    Lwork = Lnew
    Lpartial = []
    while (Lwork != []):
        if debug:
            print "New Loop: -----------------------------------"
        G, M, E, Psi, B = Lwork.pop()
        f, g = B.pop()
        T = lcm(Psi[f], Psi[g])
        gCoeffPsiG = g.monomial_coefficient(Psi[g])
        fCoeffPsiF = f.monomial_coefficient(Psi[f])
        h = gCoeffPsiG*T*f//Psi[f] - fCoeffPsiF*T*g//Psi[g]
        check, subtract = minimalPolynomialCheck(G, Psi, h)
        if debug:
            print "Inequalities"
            print E
        while check:
            if debug:
                print "h before: ", h
            h = h - subtract
            if debug:
                print "h after: ", h
            check, subtract = minimalPolynomialCheck(G, Psi, h)
        if h == 0:
            if (B == []):
                Lpartial.append((G, M, E, Psi))
            else:
                Lwork.append((G, M , E, Psi, B))
        else:
            for leadingMonomial in h.monomials():
                Gnew = G[:]
                Gnew.append(h)
                Mnew = M[:]
                M.append(leadingMonomial)
                Enew = E[:]
                for nonLeadingMonomial in h.monomials():
                    if(nonLeadingMonomial != leadingMonomial):
                        # We substract the Leading Monomial to the
                        # Non Leading Monomial because the LP solver
                        # has <= as default inequalities
                        Enew.append(subtractExponents(nonLeadingMonomial, leadingMonomial))
                Psinew = copy.deepcopy(Psi)
                Psinew[h] = leadingMonomial
                Bnew = B[:]
                for g in G:
                    Bnew.append((g, h))
                if isNotEmptyTO(Enew):
                    Lwork.append((Gnew, Mnew, Enew, Psinew, Bnew))
    return Lpartial


result = groebnerFan([x^2 - y, y^2 - x*z - y*z])
for x in result:
    print x
print len(result)
