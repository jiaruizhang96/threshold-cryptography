#!/bin/bash

# Set up certificate directory
CERTS_DIR=./certs
mkdir -p $CERTS_DIR

# Generate CA certificate
echo "Generating CA certificate..."
openssl genrsa -out $CERTS_DIR/ca.key 4096
openssl req -x509 -new -nodes -key $CERTS_DIR/ca.key -sha256 -days 3650 -out $CERTS_DIR/ca.crt -subj "/CN=Cluster CA"

# Generate certificates for each peer
for peer in peer1 peer2 peer3 peer4 peer5; do
    echo "Generating certificate for $peer..."
    openssl genrsa -out $CERTS_DIR/$peer.key 2048
    openssl req -new -key $CERTS_DIR/$peer.key -out $CERTS_DIR/$peer.csr -subj "/CN=$peer"
    openssl x509 -req -in $CERTS_DIR/$peer.csr -CA $CERTS_DIR/ca.crt -CAkey $CERTS_DIR/ca.key -CAcreateserial -out $CERTS_DIR/$peer.crt -days 365 -sha256
done

# Organize certificates into separate directories
echo "Organizing certificates..."
mkdir -p /peer1-certs /peer2-certs /peer3-certs /peer4-certs /peer5-certs

cp $CERTS_DIR/peer1.crt /peer1-certs/peer.crt
cp $CERTS_DIR/peer1.key /peer1-certs/peer.key
cp $CERTS_DIR/ca.crt /peer1-certs/ca.crt

cp $CERTS_DIR/peer2.crt /peer2-certs/peer.crt
cp $CERTS_DIR/peer2.key /peer2-certs/peer.key
cp $CERTS_DIR/ca.crt /peer2-certs/ca.crt

cp $CERTS_DIR/peer3.crt /peer3-certs/peer.crt
cp $CERTS_DIR/peer3.key /peer3-certs/peer.key
cp $CERTS_DIR/ca.crt /peer3-certs/ca.crt

echo "Certificate generation and organization complete."

