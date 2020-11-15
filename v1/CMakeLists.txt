cmake_minimum_required(VERSION 3.13)
project(v1 LANGUAGES CXX)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)

add_subdirectory(core)
add_subdirectory(eigen)
