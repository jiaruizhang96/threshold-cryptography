FROM ubuntu:24.04

RUN apt-get update \
 && export DEBIAN_FRONTEND=noninteractive \
 && apt-get -y upgrade \
 && apt-get -y install ca-certificates clang clang-format cmake g++ gdb git lldb make openssh-client
