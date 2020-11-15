cmake_minimum_required(VERSION 3.13)
project(proteins)

include(ExternalProject)
find_package(Git REQUIRED)
find_package(Eigen3 3.3)

add_subdirectory(v1)
