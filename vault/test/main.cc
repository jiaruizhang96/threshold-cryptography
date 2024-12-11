#define CATCH_CONFIG_MAIN
#define CATCH_CONFIG_ENABLE_BENCHMARKING
#include <catch2/catch.hpp>
#include <httplib.h>
#include <string>
#include <iostream>
#include <thread>
#include <chrono>

class VaultClient {
private:
    std::string host;
    int port;
    std::string vault_token;

public:
    VaultClient(const std::string& host, int port, const std::string& token)
        : host(host), port(port), vault_token(token) {}

    bool put(const std::string& path, const std::string& key, const std::string& value) {
        httplib::Client client(host, port);
        std::string data = "{ \"data\": { \"" + key + "\": \"" + value + "\" } }";
        httplib::Headers headers = {{"X-Vault-Token", vault_token}};
        auto res = client.Post(path.c_str(), headers, data, "application/json");

        if (res && (res->status == 200 || res->status == 201)) {
            //std::cout << "PUT successful: " << res->body << std::endl;
            return true;
        }
        std::cerr << "PUT failed: " << (res ? res->status : 0) << std::endl;
        return false;
    }

    std::string get(const std::string& path) {
        httplib::Client client(host, port);
        httplib::Headers headers = {{"X-Vault-Token", vault_token}};
        auto res = client.Get(path.c_str(), headers);

        if (res && res->status == 200) {
            //std::cout << "GET successful: " << res->body << std::endl;
            return res->body;
        }
        std::cerr << "GET failed: " << (res ? res->status : 0) << std::endl;
        return "";
    }
};

TEST_CASE("Benchmark PUT and GET key-value pair", "[put_get_benchmark]") {
    std::string vault_token = "hvs.gPjyEYWIIKKkAOpbXgepv3Dk"; 
    VaultClient client("127.0.0.1", 8200, vault_token);

    std::string key = "benchmark-test-key";
    std::string value = "123456789";
    std::string put_path = "/v1/secret/data/" + key;
    std::string get_path = "/v1/secret/data/" + key;

    SECTION("Benchmark PUT request") {
        BENCHMARK("PUT request") {
            REQUIRE(client.put(put_path, "password", value));
        };
    }

    std::this_thread::sleep_for(std::chrono::milliseconds(100));

    SECTION("Benchmark GET request") {
        BENCHMARK("GET request") {
            std::string response = client.get(get_path);
            REQUIRE(!response.empty());
            REQUIRE(response.find(value) != std::string::npos);
        };
    }
}
