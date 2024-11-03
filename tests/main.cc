#include <catch2/catch_test_macros.hpp>
#include <httplib.h>

using namespace httplib;

TEST_CASE("Server says hello", "[hi]") {
  Client client("localhost", 8080);

  auto res = client.Get("/hi");

  REQUIRE(res->status == 200);
  REQUIRE(res->body == "Hello World!");
}
