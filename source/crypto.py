import nacl.hash
import nacl.secret

from Crypto.PublicKey import ECC
from math import prod
from random import randint
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
        return Polynomial([S] + [randint(1, q) for i in range(t - 1)], q)

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
        return prod(x(j) * pow(x(j) - x(i), -1, q) for j in range(t) if j != i)

    return sum(y(i) * ℓ(i) for i in range(t)) % q

def derive_key(S: ECC.EccPoint):
    ikm = S.x.to_bytes(S.size_in_bytes(), byteorder='big')
    return nacl.hash.blake2b(ikm, digest_size=nacl.secret.SecretBox.KEY_SIZE, encoder=nacl.encoding.RawEncoder)

def encrypt(M: bytes, K: bytes):
    return nacl.secret.SecretBox(K).encrypt(M)

def decrypt(C: bytes, K: bytes):
    return nacl.secret.SecretBox(K).decrypt(C)
