#!/bin/bash

# Function to generate docker-compose configuration
generate_compose() {
  local n=$1  # Number of etcd servers
  echo "services:"
  for ((i=0; i<n; i++)); do
    let "client_port=2379 + 1000 * i"
    let "peer_port=2380 + 1000 * i"
    if [ $i -eq 0 ]; then
      echo "  etcd$i: &etcd"
    else
      echo "  etcd$i:"
    fi
    echo "    image: quay.io/coreos/etcd:v3.5.16"
    echo "    ports:"
    echo "      - \"$client_port:2379\""
    echo "      - \"$peer_port:2380\""
    echo "    entrypoint:"
    echo "      - etcd"
    echo "      - --enable-v2"
    echo "      - --listen-client-urls=http://0.0.0.0:2379"
    echo "      - --listen-peer-urls=http://0.0.0.0:2380"
    echo "      - --initial-cluster=$(initial_cluster $n)"
    echo "      - --initial-cluster-state=new"
    echo "      - --initial-cluster-token=mys3cr3ttok3n"
    echo "    command:"
    echo "      - --name=etcd$i"
    echo "      - --advertise-client-urls=http://etcd$i:2379"
    echo "      - --initial-advertise-peer-urls=http://etcd$i:2380"
    echo "    volumes:"
    echo "      - etcd$i:/etcd_data"
  done

  echo "volumes:"
  for ((i=0; i<n; i++)); do
    echo "  etcd$i:"
  done
}

# Generate the initial cluster string required by etcd
initial_cluster() {
  local n=$1
  local cluster=""
  for ((i=0; i<n; i++)); do
    cluster+="etcd$i=http://etcd$i:$((2380)),"
  done
  echo "${cluster%,}"
}


# Test scenarios
declare -a tests=(7)

# Main loop through each test scenario
for n in "${tests[@]}"; do
  # Dynamically generate etcd server string
  ETCD_SERVERS=""
  for ((i=0; i<n; i++)); do
    let "port=2379 + 1000 * i"
    ETCD_SERVERS+="localhost:$port,"
  done
  ETCD_SERVERS=${ETCD_SERVERS%?}  
  echo "Running test with $n servers"
  generate_compose $n > docker-compose.yml
  docker-compose up -d
  sleep 30  # Wait for the servers to stabilize

  # Set environment variables
  export ETCD_SERVERS
  export N=$n
  export K=$(((n + 1) / 2))  # Set K as (n+1)/2

  # Build and run server
  mkdir build
  cd build
  cmake ..
  make all
  stdbuf -oL -eL ./server/server >> ../server_output.log 2>&1 & # Redirect stdout and stderr to a log file
  SERVER_PID=$!
  # Wait for the server to initialize
  sleep 30  # Adjust the sleep time
  # Run tests
  ./tests/tests
  cd ..
  # Cleanup
  rm -r build
  rm docker-compose.yml
  kill $SERVER_PID
  docker-compose down --remove-orphans
  docker stop $(docker ps -q)  
  docker rm $(docker ps -aq) 
done
