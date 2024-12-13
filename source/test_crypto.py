from crypto import *

def test_shamir():
    # Elliptic-curve cryptography parameters.
    E = Curve('P-256')
    q = E.order

    # Threshold cryptography parameters.
    t = 3
    n = 5

    # Shamir's secret sharing.
    S = random.randint(1, q)
    f = Polynomial.shamir(S, t, q)

    subset = random.sample(range(1, n + 1), t)
    shares = {xᵢ: f(xᵢ) for xᵢ in subset}

    L0 = interpolate_int(shares, q)

    # Check result is as expected.
    assert L0 == S

def test_ecies():
    # Elliptic-curve cryptography parameters.
    E = Curve('P-256')
    q = E.order
    G = E.generator

    # Bob's key pair.
    k = random.randint(1, q)
    Q = k * G

    # Message from Alice.
    M = "Hello world"

    # Alice's ephemeral key pair.
    d = random.randint(1, q)
    R = d * G

    # Alice encrypts her message.
    S = d * Q
    K = symmetric_derive_key(S)
    C = symmetric_encrypt(M, K)

    # Bob decrypts Alice's message.
    S = k * R
    K = symmetric_derive_key(S)
    M = symmetric_decrypt(C, K)

    # Check result is as expected.
    assert M == "Hello world"

def test_shamir_with_ecies():
    # Elliptic-curve cryptography parameters.
    E = Curve('P-256')
    q = E.order
    G = E.generator
    I = E.identity

    # Threshold cryptography parameters.
    t = 3
    n = 5

    # Bob's key pair, with the secret key split into shares.
    k = random.randint(1, q)
    Q = k * G
    f = Polynomial.shamir(k, t, q)
    k = {xᵢ: f(xᵢ) for xᵢ in range(1, n + 1)}

    # Message from Alice.
    M = "Hello world"

    # Alice's ephemeral key pair.
    d = random.randint(1, q)
    R = d * G

    # Alice encrypts her message.
    S = d * Q
    K = symmetric_derive_key(S)
    C = symmetric_encrypt(M, K)

    # Bob decrypts Alice's message by interpolating decryption shares.
    subset = random.sample(range(1, n + 1), t)
    shares = {xᵢ: k[xᵢ] * R for xᵢ in subset}

    S = interpolate(shares, q, I)
    K = symmetric_derive_key(S)
    M = symmetric_decrypt(C, K)

    # Check result is as expected.
    assert M == "Hello world"
