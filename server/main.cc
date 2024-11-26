#include <httplib.h>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <nlohmann/json.hpp>
#include "split_shares.cc"  

using namespace httplib;
using Share = SecretPair; 
std::vector<std::pair<std::string, int>> etcd_servers;


// This function gives a random etcd server.
std::pair<std::string, int> get_random_etcd_server() {
    return {"localhost", 8080};  // Always route requests through NGINX
}

// Helper function for PUT requests
bool put_helper(const std::string& key, const std::string& value) {
    auto [host, port] = get_random_etcd_server();
    Client etcd_server(host.c_str(), port);

    std::string put_path = "/v2/keys/" + key + "?value=" + value;
    std::cout << "PUT Path: " << put_path << "\n";

    auto etcd_res = etcd_server.Put(put_path.c_str(), "", "application/x-www-form-urlencoded");

    if (etcd_res && (etcd_res->status == 201 || etcd_res->status == 200)) {
        return true;
    }
    if (etcd_res) {
        std::cout << "Failed to PUT key: " << key << " with status: " << std::to_string(etcd_res->status) << "\n";
    } else {
        std::cout << "Failed to PUT key: " << key << " - No response" << "\n";
    }

    return false;
}

// Helper function for GET requests
std::string get_helper(const std::string& key) {
    auto [host, port] = get_random_etcd_server();
    Client etcd_server(host.c_str(), port);

    std::string get_path = "/v2/keys/" + key;
    std::cout << "GET Path: " << get_path << "\n";
    auto etcd_res = etcd_server.Get(get_path.c_str());

    if (etcd_res && etcd_res->status == 200) {
        std::cout << "Successfully GET key: " << key << " with status: " << etcd_res->status << "\n";
        auto response_json = nlohmann::json::parse(etcd_res->body);
        // Extract the value inside node
        std::string value = response_json["node"]["value"];
        std::cout << "Extracted value: " << value << "\n";
        return value;

    }
    if (etcd_res) {
        std::cout << etcd_res->body << "\n";
        std::cout << "Failed to GET key: " << key << " with status: " << std::to_string(etcd_res->status) << "\n";
    } else {
        std::cout << "Failed to GET key: " << key << " - No response" << "\n";
    }
    return "";
}

int main() {
    // Reading environment variables
    char* env_n = std::getenv("N");
    char* env_k = std::getenv("K");
    const int n = std::atoi(env_n);
    const int k = std::atoi(env_k);
    
    Server server;

    server.Post("/put", [&, n, k](const Request& req, Response& res) {
        std::string key = req.get_param_value("key");
        std::string value = req.get_param_value("value");

        // Convert the input string to an integer secret
        int secret = std::stoi(value);
        bool all_puts_successful = true;
        // 1. Generate k coefficients and shares using Shamir's Secret Sharing
        auto coefficients = genCoefficients(k, secret);
        auto shares = genSecretPairs(n, coefficients);
        // 2. Store each share in etcd with a unique identifier for the key
        for (int i = 0; i < shares.size(); i++) {
            // Create a unique identifier share_key
            std::string share_key = key + "_share_" + std::to_string(shares[i].x); 
            // Convert the secret share to string
            std::string share_value = std::to_string(shares[i].y);
            // Put the key in etcd
            if (!put_helper(share_key, share_value)) {
                res.status = 500;
                res.set_content("Failed to PUT key: " + share_key + " in etcd", "text/plain");
                all_puts_successful = false;        
                break;
            } 
        }
        // if all puts are success, status=200
        if (all_puts_successful) {
            res.status = 200;
            res.set_content("Successfully PUT key: " + key + " with value: " + value, "text/plain");
        }
    });

    server.Get("/get", [&, k](const Request& req, Response& res) {
        std::string key = req.get_param_value("key");
        std::cout << " current k value: " << k << "\n";
        // 1. Retrieve at least k shares from etcd
        std::vector<Share> shares_to_recover;
        for (int i = 1; i <= k; i++) {
            // Unique key for each share
            std::string share_key = key + "_share_" + std::to_string(i); 
            std::string share_value_str = get_helper(share_key);

            // check if we get the value 
            if (!share_value_str.empty()) {
                int y_value = std::stoi(share_value_str);
                shares_to_recover.emplace_back(i, y_value);
            } else {
                std::cout << "Share value string is empty for i = " << i << "\n";
            }
        }
        // Check if we have enough shares
        if (shares_to_recover.size() < k) {
            res.status = 500;
            res.set_content("Error: Not enough shares to reconstruct the secret for key in etcd", "text/plain");    
        }
        // 2. Recover the secret using k shares
        long long recovered_secret = thresholdRecover(k, shares_to_recover); // Reconstruct secret
        //std::cout << "line 156 main.cc server directroy recovered secret:" << recovered_secret << "\n";
            
        res.status = 200;
        res.set_content("Recovered Secret: " + std::to_string(recovered_secret), "text/plain");
    });

    std::cout << "Server is running on http://0.0.0.0:8081" << "\n";
    server.listen("0.0.0.0", 8081);

    return 0;
}
