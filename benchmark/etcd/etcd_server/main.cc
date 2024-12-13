#include <httplib.h>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <nlohmann/json.hpp>

using namespace httplib;
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
    auto etcd_res = etcd_server.Get(get_path.c_str());

    if (etcd_res && etcd_res->status == 200) {
        auto response_json = nlohmann::json::parse(etcd_res->body);
        // Extract the value inside node
        std::string value = response_json["node"]["value"];
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

        if (!put_helper(key, value)) {
                res.status = 500;
                res.set_content("Failed to PUT key: " + key + " in etcd", "text/plain");
        } 
        else{
            res.status = 200;
            res.set_content("Successfully PUT key: " + key + " with value: " + value, "text/plain");
        }
        
    });

    server.Get("/get", [&, k](const Request& req, Response& res) {
        std::string key = req.get_param_value("key");
        std::string value = get_helper(key);
        res.status = 200;
        res.set_content("Recovered Secret: " + value, "text/plain");
    });

    std::cout << "Server is running on http://0.0.0.0:8081" << "\n";
    server.listen("0.0.0.0", 8081);

    return 0;
}
