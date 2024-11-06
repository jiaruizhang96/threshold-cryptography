#include <catch2/catch_test_macros.hpp>
#include <httplib.h>
#include <string>

using namespace httplib;

TEST_CASE("PUT and GET key-value pair", "[put_get]") {
    Client client("localhost", 8080);

    // Define key-value pair for testing
    std::string key = "test-key";
    std::string value = "123456789";

    SECTION("PUT value") {
        // PUT request to store the key-value pair
        std::string put_path = "/put?key=" + key + "&value=" + value;
        auto res = client.Post(put_path.c_str());

        REQUIRE(res != nullptr);
        // 200 OK status, 201 Created status 
        REQUIRE((res->status == 200 || res->status == 201));
        REQUIRE(res->body == "Successfully PUT key: " + key + " with value: " + value);
    }

    SECTION("GET value") {
        // GET request to retrieve the stored value
        std::string get_path = "/get?key=" + key;
        auto res = client.Get(get_path.c_str());

        REQUIRE(res != nullptr);
        REQUIRE(res->status == 200); // 200 OK status
        REQUIRE(res->body.find(value) != std::string::npos); // Check if response contains the value
    }
}
