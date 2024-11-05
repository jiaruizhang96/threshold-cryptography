#include <httplib.h>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>

using namespace httplib;

// List of etcd server addresses (host and port).
std::vector<std::pair<std::string, int>> etcd_servers = {
    {"localhost", 2379},
    {"localhost", 3379},
    {"localhost", 4379}
};

// This function will give you a etcd server randomly.
std::pair<std::string, int> get_random_etcd_server() {
    int index = rand() % etcd_servers.size();
    return etcd_servers[index];
}

int main() {
    srand(static_cast<unsigned>(time(0))); // Seed for random selection
    Server server;

    // PUT request
    server.Post("/put", [](const Request& req, Response& res) {
        std::string key = req.get_param_value("key");
        std::string value = req.get_param_value("value");

        // Select a random etcd server
        auto [host, port] = get_random_etcd_server();
        Client etcd_server(host.c_str(), port);

        std::cout << "Attempting to PUT key: " << key << " with value: " << value << " to etcd server at " << host << ":" << port << "\n";

        // PUT path for etcd (syntax from curl)
        std::string put_path = "/v2/keys/" + key + "?value=" + value;
        auto etcd_res = etcd_server.Put(put_path.c_str(), "", "application/x-www-form-urlencoded");

        if (etcd_res) {
            std::cout << "Received response from etcd" << "\n";
            std::cout << "Status: " << etcd_res->status << "\n";
            std::cout << "Body: " << etcd_res->body << "\n";
            // from etcd, the status could be 201 or 200 
            if ((etcd_res->status == 201) || (etcd_res->status == 200)) {
                res.status = etcd_res->status;
                res.set_content("Successfully PUT key: " + key + " with value: " + value, "text/plain");
            } else {
                res.status = etcd_res->status; // Reflect etcd status if not 200
                res.set_content("Failed to PUT key in etcd. Etcd response: " + etcd_res->body, "text/plain");
            }
        } else {
            std::cerr << "etcd_res is null. PUT request failed." << "\n";
            res.status = 500;
            res.set_content("Failed to PUT key in etcd. No response", "text/plain");
        }
    });

    // GET request
    server.Get("/get", [](const Request& req, Response& res) {
        std::string key = req.get_param_value("key");

        // Select a random etcd server
        auto [host, port] = get_random_etcd_server();
        Client etcd_server(host.c_str(), port);

        std::cout << "Attempting to GET key: " << key << " from etcd server at " << host << ":" << port << "\n";

        // Construct the GET path for etcd
        std::string get_path = "/v2/keys/" + key;
        auto etcd_res = etcd_server.Get(get_path.c_str());

        if (etcd_res) {
            std::cout << "Received response from etcd" << "\n";
            std::cout << "Status: " << etcd_res->status << "\n";
            std::cout << "Body: " << etcd_res->body << "\n";

            if (etcd_res->status == 200) {
                res.set_content("GET response: " + etcd_res->body, "application/json");
            } else {
                res.status = etcd_res->status; 
                res.set_content("Failed to GET key from etcd. Etcd response: " + etcd_res->body, "text/plain");
            }
        } else {
            std::cerr << "etcd_res is null. GET request failed." << "\n";
            res.status = 404;
            res.set_content("Key not found in etcd. No response", "text/plain");
        }
    });

    // Start the server on port 8080
    std::cout << "Server is running on http://0.0.0.0:8080" << "\n";
    server.listen("0.0.0.0", 8080);

    return 0;
}
