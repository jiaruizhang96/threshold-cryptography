import nacl.hash
import nacl.secret

from Crypto.PublicKey import ECC
from math import prod
from random import randint
from typing import Dict, List

class Curve():
    def __init__(self, curve: str):
        self._curve = ECC._curves[curve]

    @property
    def order(self):
        return int(self._curve.order)

    @property
    def generator(self):
        return self._curve.G

class Polynomial():
    @staticmethod
    def Shamir(secret: int, threshold: int, order: int):
        S = secret
        t = threshold
        q = order
        return Polynomial([S] + [randint(1, q) for _ in range(t - 1)], q)

    def __init__(self, coefficients: List[int], order: int):
        self._coefficients = coefficients
        self._order = order

    def __call__(self, x: int):
        a = self._coefficients
        q = self._order
        t = len(a)
        return sum(a[i] * pow(x, i, q) for i in range(t)) % q

def interpolate_int(shares: Dict[int, int], order: int):
    q = order
    ℓ = lambda xᵢ : prod(xⱼ * pow(xⱼ - xᵢ, -1, q) for xⱼ in shares if xⱼ != xᵢ)

    return sum(yᵢ * ℓ(xᵢ) for xᵢ, yᵢ in shares.items()) % q

def interpolate_ecc(shares: Dict[int, ECC.EccPoint], order: int, start: ECC.EccPoint):
    q = order
    ℓ = lambda xᵢ : prod(xⱼ * pow(xⱼ - xᵢ, -1, q) for xⱼ in shares if xⱼ != xᵢ)

    return sum((yᵢ * ℓ(xᵢ) for xᵢ, yᵢ in shares.items()), start)

def derive_key(S: ECC.EccPoint):
    ikm = S.x.to_bytes(S.size_in_bytes(), byteorder='big')
    return nacl.hash.blake2b(ikm, digest_size=nacl.secret.SecretBox.KEY_SIZE, encoder=nacl.encoding.RawEncoder)

def encrypt(M: bytes, K: bytes):
    return nacl.secret.SecretBox(K).encrypt(M)

def decrypt(C: bytes, K: bytes):
    return nacl.secret.SecretBox(K).decrypt(C)
