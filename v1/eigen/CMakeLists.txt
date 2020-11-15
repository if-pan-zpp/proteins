add_executable(v1.eigen)

target_sources(v1.eigen
    PRIVATE
    eigen/main.cpp)

target_link_libraries(v1.eigen
    PRIVATE
    Eigen3::Eigen)