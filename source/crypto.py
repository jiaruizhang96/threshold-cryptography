import base64
import json
import math
import nacl.hash
import nacl.secret
import random

from Crypto.PublicKey import ECC
from typing import Dict, List

class Point():
    def __init__(self, point: ECC.EccPoint):
        self._point = point

    @property
    def x(self):
        return int(self._point.x)

    @property
    def y(self):
        return int(self._point.y)

    @property
    def curve(self):
        return self._point.curve

    def __add__(self, point):
        return Point(self._point + point._point)

    def __mul__(self, scalar: int):
        return Point(self._point * scalar)

    def __rmul__(self, scalar: int):
        return Point(self._point * scalar)

    def to_dict(self):
        return {"x": self.x, "y": self.y, "curve": self.curve}

    def to_json(self):
        return json.dumps(self.to_dict())

    def to_bytes(self):
        return self.to_json().encode("utf-8")

    @classmethod
    def from_dict(cls, data):
        return cls(ECC.EccPoint(data["x"], data["y"], data["curve"]))

    @classmethod
    def from_json(cls, data):
        return cls.from_dict(json.loads(data))

class Curve():
    def __init__(self, curve: str):
        self._curve = ECC._curves[curve]

    @property
    def order(self):
        return int(self._curve.order)

    @property
    def generator(self):
        return Point(self._curve.G)

    @property
    def identity(self):
        return Point(self._curve.G.point_at_infinity())

class EncryptedValue():
    def __init__(self, public_key: Point, ciphertext: bytes):
        self.public_key = public_key
        self.ciphertext = ciphertext

    def to_dict(self):
        R = self.public_key
        C = self.ciphertext
        return {"publicKey": R.to_dict(), "ciphertext": base64.urlsafe_b64encode(C).decode("utf-8")}

    def to_json(self):
        return json.dumps(self.to_dict())

    def to_bytes(self):
        return self.to_json().encode("utf-8")

    @classmethod
    def from_dict(cls, data: Dict[str, str]):
        return cls(Point.from_dict(data["publicKey"]), base64.urlsafe_b64decode(data["ciphertext"].encode("utf-8")))

    @classmethod
    def from_json(cls, data):
        return cls.from_dict(json.loads(data))

class Polynomial():
    @staticmethod
    def shamir(secret: int, threshold: int, order: int):
        S = secret
        t = threshold
        q = order
        return Polynomial([S] + [random.randint(1, q) for _ in range(t - 1)], q)

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
    ℓ = lambda xᵢ : math.prod(xⱼ * pow(xⱼ - xᵢ, -1, q) for xⱼ in shares if xⱼ != xᵢ)

    return sum(yᵢ * ℓ(xᵢ) for xᵢ, yᵢ in shares.items()) % q

def interpolate(shares: Dict[int, Point], order: int, identity: Point):
    q = order
    I = identity
    ℓ = lambda xᵢ : math.prod(xⱼ * pow(xⱼ - xᵢ, -1, q) for xⱼ in shares if xⱼ != xᵢ)

    return sum((yᵢ * ℓ(xᵢ) for xᵢ, yᵢ in shares.items()), I)

def symmetric_derive_key(S: Point):
    return nacl.hash.blake2b(S.to_bytes(), digest_size=nacl.secret.SecretBox.KEY_SIZE, encoder=nacl.encoding.RawEncoder)

def symmetric_encrypt(M: str, K: bytes):
    return bytes(nacl.secret.SecretBox(K).encrypt(M.encode()))

def symmetric_decrypt(C: bytes, K: bytes):
    return nacl.secret.SecretBox(K).decrypt(C).decode()
