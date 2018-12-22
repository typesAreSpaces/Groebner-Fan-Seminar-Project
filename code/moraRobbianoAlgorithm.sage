import sage.all
import copy
import numpy as np
from scipy.optimize import linprog

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
    for t in h.monomials():
        for g in G:
            q, r = t.quo_rem(Psi[g])
            if r == 0 :
                hCoeffT = h.monomial_coefficient(t)
                gCoeffPsiG = g.monomial_coefficient(Psi[g])
                return True, hCoeffT/gCoeffPsiG*t//Psi[g]*g
    return False, None

def reducibilityCheck(G, M, Psi):
    for g in G:
        for m in M:
            if(m != Psi[g]):
                q, r = Psi[g].quo_rem(m)
                if r == 0 :
                    return g, True
    return None, False

# 'inequalities' is an array of vectors
def isNotEmptyTO(inequalities):
    epsilon = 0.0001
    numInequalities = len(inequalities)
    numVariables = len(inequalities[0])
    # The upper bounds are epsilon values
    # so we can express strict inequalities
    zeros = [-epsilon for x in range(numInequalities)]
    ones = [1  for x in range(numVariables)]
    # Solutions cannot contain zeros
    result = linprog(ones, A_ub=inequalities,
                     b_ub=zeros, bounds=(epsilon, None))
    return result.success

# 'G' is an array of polynomials
# (a Groebner basis)
# 'M' is an array of (initial) monomials
# 'g' is a polynomial
def redOneStep(G, M, g):
    index = 0
    for leadingMonomial in M:
        for gMonomial in g.monomials():
            q, r = gMonomial.quo_rem(leadingMonomial)
            if r == 0 :
                numerator = g.monomial_coefficient(gMonomial)
                denominator = G[index].monomial_coefficient(leadingMonomial)
                t = gMonomial//leadingMonomial
                return (g - numerator/denominator*t*G[index]), True
        index+=1
    return g, False

# 'G' is an array of polynomials
# (a Groebner basis)
# 'M' is an array of (initial) monomials
# 'g' is a polynomial
def red(G, M, g):
    polynomial = g
    polynomial, check = redOneStep(G, M, polynomial)
    while check:
        polynomial, check = redOneStep(G, M, polynomial)
    return polynomial

# M1, M2 are polynomial arrays
def containmentIdeals(M1, M2):
    flag = True
    for x in M1:
        flag = False
        for y in M2:
            q, r = x.quo_rem(y)
            if r == 0 :
                flag = True
                break
        if (not flag):
            return False
    return True

# M1, M2 are polynomial arrays
def equalIdeals(M1, M2):
    return containmentIdeals(M1, M2) and containmentIdeals(M2, M1)

# 'polynomials' is an array with polynomials
# 'Mon' is an array of arrays of polynomials
def membershipIdealArrayTest(polynomials, Mon):
    if (Mon == []):
        return False
    else:
        for ideal in Mon:
            if equalIdeals(polynomials, ideal):
                return True
        return False
            
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
                        Enew.append(subtractExponents(nonLeadingMonomial,
                                                      leadingMonomial))
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
        G, M, E, Psi, B = Lwork.pop()
        f, g = B.pop()
        T = lcm(Psi[f], Psi[g])
        gCoeffPsiG = g.monomial_coefficient(Psi[g])
        fCoeffPsiF = f.monomial_coefficient(Psi[f])
        h = gCoeffPsiG*T*f//Psi[f] - fCoeffPsiF*T*g//Psi[g]
        check, subtract = minimalPolynomialCheck(G, Psi, h)
        while check:
            h = h - subtract
            check, subtract = minimalPolynomialCheck(G, Psi, h)
        if h == 0 :
            if (B == []):
                Lpartial.append((G, M, E, Psi))
            else:
                Lwork.append((G, M , E, Psi, B))
        else:
            for leadingMonomial in h.monomials():
                Gnew = G[:]
                Gnew.append(h)
                Mnew = M[:]
                Mnew.append(leadingMonomial)
                Enew = E[:]
                for nonLeadingMonomial in h.monomials():
                    if(nonLeadingMonomial != leadingMonomial):
                        # We substract the Leading Monomial to the
                        # Non Leading Monomial because the LP solver
                        # has <= as default inequalities
                        Enew.append(subtractExponents(nonLeadingMonomial,
                                                      leadingMonomial))
                Psinew = copy.deepcopy(Psi)
                Psinew[h] = leadingMonomial
                Bnew = B[:]
                for g in G:
                    Bnew.append((g, h))
                if isNotEmptyTO(Enew):
                    Lwork.append((Gnew, Mnew, Enew, Psinew, Bnew))
    # Computation of the Reduced Groebner Bases
    # and of the Groebner Region
    Loutput = []
    Mon = []
    while (Lpartial != []):
        G, M, E, Psi = Lpartial.pop()
        if (not membershipIdealArrayTest(M, Mon)):
            polynomial, check = reducibilityCheck(G, M, Psi)
            while check:
                G.remove(polynomial)
                M.remove(Psi[polynomial])
                del Psi[polynomial]
                polynomial, check = reducibilityCheck(G, M, Psi)
            for g in G:
                G.remove(g)
                M.remove(Psi[g])
                gnew = red(G, M, g)
                coeffGNew = gnew.monomial_coefficient(Psi[g])
                gnew = 1/coeffGNew*gnew
                G.append(gnew)
                M.append(Psi[g])
                tempPsiG = Psi[g]
                del Psi[g]
                Psi[gnew] = tempPsiG
            E = []
            for g in G:
                for monomial in g.monomials():
                    if(monomial != Psi[g]):
                        E.append(subtractExponents(Psi[g],
                                                   monomial))
            Loutput.append((G, M, E, Psi))
            Mon.append(M)
    return Loutput


if __name__ == "__main__":
    R.<x,y,z> = PolynomialRing(QQ)
    examples = [[x^2  - y, y^2  - x*z - y*z],
                [y*z + x, x*y + z, x^2 -z^2],
                [x*y - x, x^2 + x*z, y^2*z + x]]
    
    for example in examples:
        result = groebnerFan(example)
        print "Groebner Fan of: ", example
        for (G, M, E, Psi) in result:
            print "G: ", G, " M: ", M
        print "Size: ", len(result)

