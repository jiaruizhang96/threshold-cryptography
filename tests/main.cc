#define CATCH_CONFIG_MAIN
#define CATCH_CONFIG_ENABLE_BENCHMARKING
#include <catch2/catch.hpp>
#include <httplib.h>
#include <string>

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
            REQUIRE(res->body.find(value) != std::string::npos); // Check if response contains the value
        };
    }
}
