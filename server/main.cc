#include <httplib.h>

using namespace httplib;

int main() {
  Server server;

  server.Get("/hi", [](Request const&, Response& res) {
    res.set_content("Hello World!", "text/plain");
  });

  server.listen("0.0.0.0", 8080);
  return 0;
}
