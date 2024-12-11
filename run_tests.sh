#!/bin/bash

# Function to generate docker-compose configuration
generate_compose() {
  local n=$1  # Number of etcd servers
  echo "version: '3.8'"
  echo "services:"
  for ((i=0; i<n; i++)); do
    let "client_port=2379 + 1000 * i"
    let "peer_port=2380 + 1000 * i"
    echo "  etcd$i:"
    echo "    image: quay.io/coreos/etcd:v3.5.16"
    echo "    container_name: etcd$i"
    echo "    ports:"
    echo "      - \"$client_port:2379\""
    echo "      - \"$peer_port:2380\""
    echo "    command: >"
    echo "      etcd"
    echo "      --name etcd$i"
    echo "      --data-dir /etcd-data"
    echo "      --listen-client-urls http://0.0.0.0:2379"
    echo "      --advertise-client-urls http://etcd$i:2379"
    echo "      --listen-peer-urls http://0.0.0.0:2380"
    echo "      --initial-advertise-peer-urls http://etcd$i:2380"
    echo "      --initial-cluster $(initial_cluster $n)"
    echo "      --initial-cluster-state new"
    echo "      --initial-cluster-token mys3cr3ttok3n"
    echo "      --enable-v2"
    echo "    volumes:"
    echo "      - etcd$i:/etcd_data"
  done

  echo "  nginx:"
  echo "    image: nginx:latest"
  echo "    container_name: nginx"
  echo "    ports:"
  echo "      - \"8080:8080\""
  echo "    volumes:"
  echo "      - ./nginx.conf:/etc/nginx/nginx.conf:ro"

  echo "volumes:"
  for ((i=0; i<n; i++)); do
    echo "  etcd$i:"
  done
}

# Function to generate the initial cluster string required by etcd
initial_cluster() {
  local n=$1
  local cluster=""
  for ((i=0; i<n; i++)); do
    cluster+="etcd$i=http://etcd$i:2380,"
  done
  echo "${cluster%,}"
}

# Function to generate NGINX configuration
generate_nginx_config() {
  local n=$1  # Number of etcd servers
  echo "events{}"
  echo "http {
    upstream etcd_cluster {
        least_conn;"  # Use least connected method
  for ((i=0; i<n; i++)); do
    let "port=2379 + i * 1000"  # Adjust for mapped ports
    echo "        server etcd$i:$port max_fails=3 fail_timeout=30s;"
  done
  echo "    }

    server {
        listen 8080;

        location / {
            proxy_pass http://etcd_cluster;
            proxy_set_header Host \$host;
            proxy_set_header X-Real-IP \$remote_addr;
            proxy_set_header X-Forwarded-For \$proxy_add_x_forwarded_for;
            proxy_http_version 1.1;
            proxy_set_header Connection '';
        }
    }
  }"
}

# Test scenarios
declare -a tests=(5)

# Main loop through each test scenario
for n in "${tests[@]}"; do
  # Dynamically generate etcd server string
  echo "Running test with $n servers"
  
  # Generate the docker-compose.yml file
  generate_compose $n > docker-compose.yml

  # Generate the nginx.conf file
  generate_nginx_config $n > nginx.conf

  # Start services
  docker-compose up -d
  sleep 30  # Wait for the services to stabilize

  # Set environment variables
  export N=$n
  export K=$(((n + 1) / 2))  # Set K as (n+1)/2

  # Build and run server
  mkdir -p build
  cd build
  cmake ..
  make all
    
  # Start the server in the background
  stdbuf -oL -eL ./server/server 2>&1 | tee server_output.log &
  SERVER_PID=$!
  sleep 30  # Wait for the server to initialize

  # Run tests
  ./tests/tests
  cd ..

  # Cleanup
  kill $SERVER_PID
  docker-compose down --remove-orphans
  docker stop $(docker ps -q)
  docker rm $(docker ps -aq)
done

# NGINX and the etcd instances are running inside the same Docker network 
# the internal ports (2379) are the ones to use in the upstream block.

# Inside the Docker network, the etcd0, etcd1, and etcd2 containers 
# use the default internal port 2379 for client communication.

# The external ports (2379, 3379, 4379) are host-mapped, not visible inside the Docker network
# NGINX must communicate with etcd containers on their internal 
# Docker network addresses and ports (etcd0:2379).
# if you call curl via command line, you can use external ports
# curl -X PUT "http://localhost:3379/v2/keys/hello" -d value=world
# curl "http://localhost:3379/v2/keys/hello"