# Copyright (c) 2017 - 2023, Philip Herd
# This file is distributed under the BSD 2-Clause License

from cmath import sqrt
__doc__ = """
A collection of general purpose analytic formulas
for polynomials of degree 0 through 4

Polynomials with degrees higher than four do not have analytic solutions
"""

def root_0(a):
    """returns the roots for a constant equation
    a = 0, polynomial of degree 0"""
    return 0

def root_1(a, b):
    """returns the roots for a linear equation
    ax + b = 0,  polynomial of degree 1"""
    return -b / a

def root_2(a, b, c):
    """returns the roots for a quadratic equation
    ax^2 + bx + c = 0,  polynomial of degree 2"""
    p1 = sqrt(b * b - 4. * a * c)
    p2 = -2. * a

    x1 = (b - p1) / p2
    x2 = (b + p1) / p2

    return x1, x2

def root_3(a, b, c, d):
    """returns the roots for a cubic equation
    ax^3 + bx^2 + cx + d = 0,  polynomial of degree 3"""
    abc = a * b * c
    bbb = b * b * b
    aad = a * a * d

    dd = (18. * abc * d - 4. * bbb * d
          + b * b * c * c - 4. * a * c * c * c
          - 27. * aad * d)
    d0 = b * b - 3. * a * c

    # second and third cubic unity roots (first is just 1)
    cu2 = -0.5 + 0.86602540378443864676j
    cu3 = -0.5 - 0.86602540378443864676j

    if not dd and not d0:  # all real roots
        x1 = x2 = x3 = -b / (3. * a)
    elif not dd and d0:  # double root, simple root
        x1 = x2 = ((9. * a * d - b * c) / (2. * d0))
        x3 = (4. * abc - 9. * aad - bbb) / (a * d0)
    else:
        d1 = 2. * bbb - 9. * abc
        d1 = d1 + 27. * aad

        if not d0: cin = d1 + 0j # inner terms cancel
        else: cin = (d1 - sqrt(-27.0 * a * a * dd)) / 2.

        cc = cin ** (1. / 3.)
        p = (-1. / (3. * a))

        x1 = p * (b + cc + d0 / cc)
        x2 = p * (b + cu2 * cc + d0 / (cu2 * cc))
        x3 = p * (b + cu3 * cc + d0 / (cu3 * cc))

    return x1, x2, x3

def root_4(a, b, c, d, e):
    """returns the roots for a quartic equation
    ax^4 + bx^3 + cx^2 + dx + e = 0,  polynomial of degree 4"""
    aa, bb, cc, dd = b / a, c / a, d / a, e / a
    a2, b2 = aa * aa, bb * bb

    bq = (- (2. * b2 * bb) + 9. * aa * bb * cc
          - 27. * (cc * cc + a2 * dd)
          + 72. * bb * dd)
    c1 = (b2 - 3. * aa * cc + 12. * dd)
    cu2 = -0.5 + 0.86602540378443864676j

    p1 = sqrt(bq * bq - (4. * c1 * c1 * c1))
    v = (bq - p1) / -2.
    if not v: v = (bq + p1) / -2.  # choose non zero quad root

    u = a2 / 4. - (2. * bb) / 3.
    if not v: uu = u  # both quad roots zero, uu simplifies to u
    else:
        v3 = (v ** (1. / 3.)) * cu2
        uu = u + (1. / 3.) * (v3 + c1 / v3)

    p1 = - aa / 4.
    if not uu:  # degenerate, quadruple root
        x1 = x2 = x3 = x4 = p1
    else:
        p2 = 3. * a2 - 8. * bb - 4. * uu
        p3 = -(a2 * aa) + 4. * aa * bb - 8. * cc

        usq = sqrt(uu)
        usq2 = usq / 2.
        u4 = uu / 4.

        blkp = .25 * sqrt(p2 + p3 /  usq)
        blkm = .25 * sqrt(p2 + p3 / -usq)

        x1 = p1 + usq2 + blkp
        x2 = p1 - usq2 + blkm
        x3 = p1 + usq2 - blkp
        x4 = p1 - usq2 - blkm
    return x1, x2, x3, x4

def cons(a):
    "constant formula, returns root_0(a)"
    return root_0(a)

def lin(a, b):
    "linear formula, returns root_1(a, b)"
    return root_1(a, b)

def quad(a, b, c):
    "quadratic forumula, returns root_2(a, b, c)"
    return root_2(a, b, c)

def cubic(a, b, c, d):
    "cubic forumula, returns root_3(a, b, c, d)"
    return root_3(a, b, c, d)

def quartic(a, b, c, d, e):
    "quartic formula, returns root_4(a, b, c, d, e)"
    return root_4(a, b, c, d, e)
