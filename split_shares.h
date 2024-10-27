// split_shares.h

#ifndef SPLIT_SHARES_H
#define SPLIT_SHARES_H

#include <vector>

// Structure to store secret shares as pairs of (x, y)
struct SecretPair {
    int x;
    int y;
    SecretPair(int x, int y) : x(x), y(y) {}
};

// Function declarations
int genRandom(int min, int max);
std::vector<int> genCoefficients(int k, int secret);
std::vector<SecretPair> genSecretPairs(int n, const std::vector<int>& coefficients);
long long modInverse(long long a, long long m);
long long reconstructSecret(const std::vector<SecretPair>& shares);
long long thresholdRecover(int k, const std::vector<SecretPair>& shares);

#endif // SPLIT_SHARES_H
