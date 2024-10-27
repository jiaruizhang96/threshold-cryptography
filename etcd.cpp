#include <iostream>
#include <cstdlib>  // For system()
#include <cstdio>   // For popen(), fgets()
#include <memory>   // For std::unique_ptr
#include "etcd.h"

using namespace std;

void EtcdHandler::startEtcdServer() {
    // Start etcd in terminal 
    string command = "./etcd &"; 
    int result = system(command.c_str());
    if (result == 0) {
        cout << "etcd server started successfully.\n";
    } else {
        cerr << "Failed to start etcd server.\n";
    }
}

void EtcdHandler::put(const string& key, const string& value) {
    string command = "./etcdctl --endpoints=localhost:2379 put " + key + " " + value ;
    int result = system(command.c_str());
    if (result == 0) {
        cout << "Successfully wrote key-value pair to etcd: " << key << " = " << value << endl;
    } else {
        cerr << "Failed to write key-value pair to etcd.\n";
    }
}

string EtcdHandler::get(const string& key) {
    string command = "./etcdctl --endpoints=localhost:2379 get " + key;
    string result;

    // Use popen to execute the command and read the output
    unique_ptr<FILE, decltype(&pclose)> pipe(popen(command.c_str(), "r"), pclose);
    if (!pipe) {
        cerr << "Failed to run command.\n";
        return "";
    }

    char buffer[128];
    while (fgets(buffer, sizeof(buffer), pipe.get()) != nullptr) {
        result += buffer;
    }

    //cout << "Raw result from `etcdctl get " << key;
    //cout << result;

    // split the output to find the value part
    size_t firstNewline = result.find("\n");
    if (firstNewline != string::npos) {
        // the 1st new line is the key 
        // 2nd new line is the value 
        size_t secondNewline = result.find("\n", firstNewline + 1);
        if (secondNewline != string::npos) {
            // Extract everything after the first newline
            result = result.substr(firstNewline + 1, secondNewline - firstNewline - 1);
        } else {
            // If there's only one newline, take everything after it
            result = result.substr(firstNewline + 1);
        }
    }

    // trim trailing or leading whitespace
    result.erase(0, result.find_first_not_of(" \n\r\t"));
    result.erase(result.find_last_not_of(" \n\r\t") + 1);

    // Print the processed result after parsing
    cout << "Processed result for key '" << key << endl;
    cout << result << endl;

    if (result.empty()) {
        cerr << "Failed to retrieve value for key: " << key;
    } else {
        cout << "Retrieved value for key '" << key << endl;
        cout << result << endl;
    }

    return result;
}


/*

int main() {
    EtcdHandler etcd;

    // Start the etcd server (optional, if not already running)
    etcd.startEtcdServer();

    // Run an infinite loop for user interaction
    while (true) {
        string input;
        cout << "> "; // Command prompt symbol

        // Read the entire line of input from the user
        getline(cin, input);

        // Exit the loop if the user types 'exit'
        if (input == "exit") {
            cout << "Exiting...\n";
            break;
        }

        // Process the user input and execute the etcd commands
        processUserInput(etcd, input);
    }

    return 0;
}
*/
