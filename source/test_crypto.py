import random

from crypto import *

def test_shamir():
    E = Curve("P-256")
    t = 3
    n = 5
    q = E.q

    S = random.randint(1, q - 1)
    f = Polynomial.random(S, t, q)

    shares = [(i, f(i)) for i in range(1, n + 1)]

    assert S == interpolate(shares[:t], q)
