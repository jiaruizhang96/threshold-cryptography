#include <iostream>
#include <vector>
#include <string>
#include <utility>
#include "etcd.h"          
#include "split_shares.h"  

using namespace std;
using Share = SecretPair; 
/*
SecretKeyValueStore class:
A key-value store that uses Shamir's Secret Sharing and etcd. 
 */
class SecretKeyValueStore {
private:
    EtcdHandler& etcd; // reference to an EtcdHandler instance
    int n; // number of shares
    int k; // threshold to reconstruct the secret

public:
    // constructor 
    SecretKeyValueStore(int numShares, int threshold, EtcdHandler& etcdHandler) 
        : n(numShares), k(threshold), etcd(etcdHandler) {}

    // Put: store the secret using etcd
    void put(const string &key, int secret) {
        // 1. generate k coefficients and shares using Shamir's Secret Sharing
        auto coefficients = genCoefficients(k, secret);
        auto shares = genSecretPairs(n, coefficients);
        
        // 2. store each share in etcd with a unique identifier for the key
        for (int i = 0; i < shares.size(); i++) {
            // create a unique identifier share_key
            string share_key = key + "_share_" + to_string(shares[i].x); 
            // convert the secret share to string
            string share_value = to_string(shares[i].y);
            // put the key in etcd
            etcd.put(share_key, share_value); 
        }

        cout << "Secret stored securely under key: " << key;
    }

    // Get: retrieve the secret by reconstructing it from the shares
    long long get(const string &key) {
        // 1. retrieve at least k shares from etcd
        vector<Share> shares_to_recover;
        for (int i = 1; i <= k; i++) {
            // unique key for each share
            string share_key = key + "_share_" + to_string(i); 
            // get the value from etcd
            string share_value_str = etcd.get(share_key);

            if (!share_value_str.empty()) {
                // convert the value to int
                int y_value = stoi(share_value_str);
                // pair the value with key 
                shares_to_recover.emplace_back(i, y_value); 
            }
        }

        // check if we have enough shares
        if (shares_to_recover.size() < k) {
            cerr << "Error: Not enough shares to reconstruct the secret for key: " << key;
            return -1;
        }

        // 2. recover the secret using k shares
        long long recovered_secret = thresholdRecover(k, shares_to_recover); // Reconstruct secret

        cout << "Recovered Secret for key " << key << ": " << recovered_secret;
        return recovered_secret;
    }
};


int main() {
    // 1. start a single instance of EtcdHandler
    EtcdHandler etcd;
    etcd.startEtcdServer(); // Start the etcd server

    // 2. init some vars
    int n = 5; // Number of shares
    int k = 3; // Threshold to reconstruct the secret
    SecretKeyValueStore store(n, k, etcd);

    // 3. put 
    int secret = 123456789;
    string key = "skylerzhang";
    store.put(key, secret);

    // 4. get 
    long long recovered_secret = store.get(key);

    // Verify the recovered secret
    if (recovered_secret == secret) {
        cout << "Secret recovered successfully!\n";
    } else {
        cout << "Error: Secret recovery failed.\n";
    }

    return 0;
}