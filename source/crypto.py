import math
import random

from Crypto.PublicKey import ECC
from typing import List, Tuple

class Curve():
    def __init__(self, curve: str):
        self._curve = ECC._curves[curve]

    @property
    def q(self):
        return int(self._curve.order)

    @property
    def G(self):
        return self._curve.G

class Polynomial():
    @staticmethod
    def random(S: int, t: int, q: int):
        return Polynomial([S] + [random.randint(1, q - 1) for i in range(t - 1)], q)

    def __init__(self, a: List[int], q: int):
        self.a = a
        self.q = q

    def __call__(self, x: int):
        a = self.a
        q = self.q
        t = len(a)
        return sum(a[i] * pow(x, i, q) for i in range(t)) % q

def interpolate(shares: List[Tuple[int, int]], q: int):
    t = len(shares)

    def x(i: int):
        return shares[i][0]

    def y(i: int):
        return shares[i][1]

    def ℓ(i: int):
        return math.prod(x(j) * pow(x(j) - x(i), -1, q) for j in range(t) if j != i)

    return sum(y(i) * ℓ(i) for i in range(t)) % q
