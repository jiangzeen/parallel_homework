cmake_minimum_required(VERSION 3.20)
project(parallel_mpi)

set(CMAKE_CXX_STANDARD 14)
link_libraries(
        /usr/local/lib64/libssl.so
        /usr/local/lib64/libcrypto.so
        /lib/x86_64-linux-gnu/libgmp.so
        libpthread.so.0
)
include_directories(
        links
        /usr/local/include/openssl
)
add_executable(parallel_mpi main.cpp)
