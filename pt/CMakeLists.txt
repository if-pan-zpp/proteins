cmake_minimum_required(VERSION 3.13)
project(prototype)
set(CMAKE_CXX_STANDARD 17)

find_package(Protobuf REQUIRED)
find_package(Eigen3 3.3 REQUIRED)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")

add_subdirectory(src)
