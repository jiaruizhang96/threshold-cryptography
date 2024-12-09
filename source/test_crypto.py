from crypto import *
from random import randint, sample

def test_shamir():
    # Curve parameters.
    E = Curve('P-256')
    q = E.order

    # Threshold encryption parameters.
    t = 3
    n = 5

    # Shamir's secret sharing.
    S = randint(1, q)
    f = Polynomial.Shamir(S, t, q)

    shares = {i: f(i) for i in range(1, n + 1)}
    subset = sample(range(1, n + 1), t)
    L_zero = interpolate_int({i: shares[i] for i in subset}, q)

    # Check result is as expected.
    assert S == L_zero

def test_ecies():
    # Curve parameters.
    E = Curve('P-256')
    q = E.order
    G = E.generator

    # Bob's key pair.
    k = randint(1, q)
    Q = k * G

    # Message from Alice.
    M = "Hello world"

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
    assert M == "Hello world"

def test_shamir_with_ecies():
    # Curve parameters.
    E = Curve('P-256')
    q = E.order
    G = E.generator

    # Threshold encryption parameters.
    t = 3
    n = 5

    # Bob's key pair, with the secret key split into shares.
    k = randint(1, q)
    Q = k * G

    f = Polynomial.Shamir(k, t, q)
    k = {i: f(i) for i in range(1, n + 1)}

    # Message from Alice.
    M = "Hello world"

    # Alice's ephemeral key pair.
    d = randint(1, q)
    R = d * G

    # Alice encrypts her message.
    S = d * Q
    K = derive_key(S)
    C = encrypt(M, K)

    # Bob decrypts Alice's message by interpolating decryption shares.
    subset = sample(range(1, n + 1), t)
    shares = {i: k[i] * R for i in subset}

    S = interpolate_ecc(shares, q, G.point_at_infinity())
    K = derive_key(S)
    M = decrypt(C, K)

    # Check result is as expected.
    assert M == "Hello world"
