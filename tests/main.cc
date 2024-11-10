#define CATCH_CONFIG_MAIN
#define CATCH_CONFIG_ENABLE_BENCHMARKING
#include <catch2/catch.hpp>
#include <httplib.h>
#include <string>
#include "../server/split_shares.cc"

using namespace httplib;

TEST_CASE("Benchmark PUT and GET key-value pair", "[put_get_benchmark]") {
    Client client("localhost", 8080);

    // key-value pair for testing
    std::string key = "benchmark-test-key";
    std::string value = "123456789";

    SECTION("Benchmark PUT request") {
        BENCHMARK("PUT request") {
            // PUT request to store the key-value pair
            std::string put_path = "/put?key=" + key + "&value=" + value;
            auto res = client.Post(put_path.c_str());

            REQUIRE(res != nullptr);
            REQUIRE((res->status == 200 || res->status == 201));
            REQUIRE(res->body == "Successfully PUT key: " + key + " with value: " + value);
        };
    }

    // wait 0.1s: PUT is complete before measuring GET
    std::this_thread::sleep_for(std::chrono::milliseconds(100));

    SECTION("Benchmark GET request") {
        BENCHMARK("GET request") {
            // GET request to retrieve the stored value
            std::string get_path = "/get?key=" + key;
            auto res = client.Get(get_path.c_str());
            REQUIRE(res != nullptr);
            REQUIRE(res->status == 200); // Check for 200 OK status
            std::cout << "line 40 main.cc tests directory recovered secret:" <<res->body << "\n";
            REQUIRE(res->body.find(value) != std::string::npos); // Check if response contains the value
        };
    }
}

TEST_CASE("Verify key splitting and reconstruction correctness", "[key_split_reconstruct]") {
    // Get the number of shares and threshold from environment variables
    const int n = 9;  // number of shares
    const int k = 5;  // threshold to reconstruct
    int secret = 123456789;  // Original secret to split and reconstruct

    SECTION("Key splitting and reconstruction correctness test") {
        // Generate coefficients and shares using Shamir's Secret Sharing
        auto coefficients = genCoefficients(k, secret);
        auto shares = genSecretPairs(n, coefficients);

        // Ensure the shares have been generated correctly
        REQUIRE(shares.size() == n);

        // Use the first k shares for reconstruction
        std::vector<SecretPair> shares_to_recover(shares.begin(), shares.begin() + k);

        // Reconstruct the secret
        long long reconstructed_secret = thresholdRecover(k, shares_to_recover);

        // Verify that the reconstructed secret matches the original
        REQUIRE(reconstructed_secret == secret);
    }
}

/*
TEST_CASE("Benchmark key splitting and reconstruction", "[key_split_reconstruct]") {
    //const int n = 3;  // number of shares
    //const int k = 2;  // threshold to reconstruct
    const int n = std::atoi(getenv("N"));  // number of shares
    const int k = std::atoi(getenv("K"));  // threshold to reconstruct
    int secret = 123456789;  
    SECTION("Benchmark key splitting") {
        BENCHMARK("Key splitting") {
            auto coefficients = genCoefficients(k, secret);
            auto shares = genSecretPairs(n, coefficients);
            return shares;  
        };
    }

    SECTION("Benchmark key reconstruction") {
        // generate the shares to use for reconstruction
        auto coefficients = genCoefficients(k, secret);
        auto shares = genSecretPairs(n, coefficients);

        // use the first k shares
        std::vector<SecretPair> shares_to_recover(shares.begin(), shares.begin() + k);

        BENCHMARK("Key reconstruction") {
            long long recovered_secret = thresholdRecover(k, shares_to_recover);
            return recovered_secret;  
        };
    }
}*/