// etcd.h

#ifndef ETCD_H
#define ETCD_H

#include <string>

class EtcdHandler {
public:
    // Start the etcd server
    void startEtcdServer();

    // Put a key-value pair in etcd
    void put(const std::string& key, const std::string& value);

    // Get a value by key from etcd
    std::string get(const std::string& key);
};

#endif // ETCD_H
