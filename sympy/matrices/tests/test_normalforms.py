from sympy import Symbol, Poly
from sympy.polys.solvers import RawMatrix
from sympy.polys.polymatrix import PolyMatrix
from sympy.matrices.normalforms import invariant_factors, smith_normal_form
from sympy.polys.domains import ZZ, QQ

def test_smith_normal():
    m = RawMatrix([[12, 6, 4,8],[3,9,6,12],[2,16,14,28],[20,10,10,20]])
    setattr(m, 'ring', ZZ)
    smf = RawMatrix([[1, 0, 0, 0], [0, 10, 0, 0], [0, 0, -30, 0], [0, 0, 0, 0]])
    assert smith_normal_form(m) == smf

    x = Symbol('x')
    m = PolyMatrix([[Poly(x-1), Poly(1, x),Poly(-1,x)],
                [0, Poly(x), Poly(-1,x)],
                [Poly(0,x),Poly(-1,x),Poly(x)]])
    setattr(m, 'ring', QQ[x])
    invs = (Poly(1, x), Poly(x - 1), Poly(x**2 - 1))
    assert invariant_factors(m) == invs
    snf, s, t = smith_normal_form(m, full=True)
    assert snf == s*m*t

    m = RawMatrix([[2, 4]])
    setattr(m, 'ring', ZZ)
    smf = RawMatrix([[2, 0]])
    assert smith_normal_form(m) == smf

    m = RawMatrix([[QQ(2), QQ(1)/3], [QQ(4),QQ(2)/3]])
    setattr(m, "ring", QQ)
    snf, s, t = smith_normal_form(m, full=True)
    assert snf == s*m*t

