from crypto import *
from random import randint

def test_shamir():
    # Curve parameters.
    E = Curve("P-256")
    q = E.q

    # Threshold encryption parameters.
    t = 3
    n = 5

    # Shamir's secret sharing.
    S = randint(1, q)
    f = Polynomial.random(S, t, q)

    shares = [(i, f(i)) for i in range(1, n + 1)]

    # Check result is as expected.
    assert S == interpolate(shares[:t], q)

def test_ecies():
    # Curve parameters.
    E = Curve('P-256')
    q = E.q
    G = E.G

    # Message from Alice.
    M = b"Hello world"

    # Bob's key pair.
    k = randint(1, q)
    Q = k * G

    # Alice's ephemeral key pair.
    d = randint(1, q)
    R = d * G

    # Alice encrypts her message.
    S = d * Q
    K = derive_key(S)
    C = encrypt(M, K)

    # Bob decrypts Alice's message.
    S = k * R
    K = derive_key(S)
    M = decrypt(C, K)

    # Check result is as expected.
    assert M.decode('utf-8') == "Hello world"
