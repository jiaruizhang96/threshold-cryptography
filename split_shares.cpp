#include <iostream>
#include <vector>
#include <random>
#include <ctime>
#include <algorithm> 

// Assuming a large prime number, for example, a 32-bit prime
const long long PRIME = 1000000007;

// Structure to store secret shares as pairs of (x, y)
struct SecretPair {
    int x;
    int y;
    SecretPair(int x, int y) : x(x), y(y) {}
};

// generate random number between min and max 
int genRandom(int min, int max) {
    static std::mt19937 rng((unsigned) std::time(nullptr)); // static to use seeding only once
    std::uniform_int_distribution<std::mt19937::result_type> dist(min, max);
    return dist(rng);
}

// generate polynomial coefficients
std::vector<int> genCoefficients(int k, int secret) {
    // A vector coefficients of size k: store the polynomial coefficients (int)
    std::vector<int> coefficients(k);
    // The constant term is the secret
    coefficients[0] = secret; 
    // Function = secret + a_1 * x + a_2 * x^2 + ...+ a_(k-1) * x^(k-1)
    // a_1, a_2 ... a_(k-1)are all random numbers 
    for (int i = 1; i < k; ++i) {
        coefficients[i] = genRandom(1, PRIME - 1);
    }
    return coefficients;
}

// gen secret pairs (shares)
std::vector<SecretPair> genSecretPairs(int n, const std::vector<int>& coefficients) {
    std::vector<SecretPair> shares;
    shares.reserve(n); // preallocates memory 
    int k = coefficients.size();

    for (int i = 1; i <= n; ++i) {
        // coefficients[0] -> constant, term = 0 
        long long y = coefficients[0];
        // get higher terms values , term = 1 ... (k-1)
        for (int j = 1; j < k; ++j) {
            long long term = 1;
            // term = i^j % PRIME
            // calculates the power of x (i^j) for the j-th degree term.
            for (int exp = 1; exp <= j; ++exp) {
                term = term * i % PRIME;
            }
            // (the coefficient for the x^j term) 
            // term = x^j (which we calculated above)
            y = (y + coefficients[j] * term) % PRIME;
        }
        // create (x,y) SecretPair in shares
        shares.emplace_back(i, static_cast<int>(y)); 
    }
    return shares;
}


// compute the modular inverse using Extended Euclidean Algorithm
// returns the modular inverse of 'a' under modulus 'm'
// it finds 'x' such that (a * x) % m == 1
// eg. modInverse(3, 11) = 4, (3*4)%11 = 1
long long modInverse(long long a, long long m) {
    long long m0 = m, t, q;
    long long x0 = 0, x1 = 1;
    // Base case: If m = 1, the modular inverse doesn't exist
    if (m == 1)
        return 0;
    // Loop continues until 'a' is reduced to 1 
    // modular inverse only exists when gcd(a, m) = 1)
    while (a > 1) {
        // q is quotient = q/m
        q = a / m;
        t = m;

        // euclidean update: m = remainder of 'a % m'
        m = a % m;
        // a = previous value of 'm' (stored in 't')
        a = t;
        t = x0;
        x0 = x1 - q * x0;
        //x1 is the new coefficient for a
        x1 = t;
    }

    // Make x1 positive
    if (x1 < 0)
        x1 += m0;

    return x1;
}

// reconstruct the secret using Lagrange interpolation
long long reconstructSecret(const std::vector<SecretPair>& shares) {
    long long secret = 0;
    int k = shares.size();
    for (int i = 0; i < k; i++) {
        long long li = 1;
        for (int j = 0; j < k; j++) {
            if (i != j) {
                // xi and xj are x-coordinates of the i-th and j-th shares
                long long xi = shares[i].x;
                long long xj = shares[j].x;
                // The numerator of Lagrange basis term: (0 - x_j) = -xj
                long long num = -xj % PRIME;
                // the denominator of Lagrange basis term: (x_i - x_j)
                long long den = (xi - xj) % PRIME;

                if (den < 0) {
                    den += PRIME; // correct negative modulus
                }
                // li *= (num * inverse(den)) % PRIME
                // 'modInverse' computes the modular inverse of den wrt PRIME
                li *= num * modInverse(den, PRIME) % PRIME;
                // make sure it's positive 
                li = (li % PRIME + PRIME) % PRIME; 
            }
        }
        long long yi = shares[i].y;
        secret = (secret + yi * li) % PRIME;
    }
    return secret;
}

// Function to select only 'k' shares from the total 'n' shares and recover the secret
long long thresholdRecover(int k, const std::vector<SecretPair>& shares) {
    // Create a copy of the original shares so we can shuffle them
    std::vector<SecretPair> selectedShares = shares;

    // shuffle the shares
    std::shuffle(selectedShares.begin(), selectedShares.end(), std::mt19937(std::random_device()()));

    // Select the first 'k' shares
    std::vector<SecretPair> thresholdShares(selectedShares.begin(), selectedShares.begin() + k);

    // Recover the secret using only the threshold number of shares
    return reconstructSecret(thresholdShares);
}

int main() {
    int n = 5; // Number of shares
    int k = 3; // Threshold: min number of shares needed to reconstruct the secret
    int secret = 123456789; 

    auto coefficients = genCoefficients(k, secret);
    auto shares = genSecretPairs(n, coefficients);
    // Recover the secret using only 'k' shares
    long long recovered_secret = thresholdRecover(k, shares);

    std::cout << "Original Secret: " << secret << std::endl;
    std::cout << "Recovered Secret using " << k << " shares: " << recovered_secret << std::endl;

    return 0;
}
