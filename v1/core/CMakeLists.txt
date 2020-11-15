add_library(core STATIC)

target_sources(core
    PRIVATE core/residues.cpp
    PUBLIC core/pdb/pdb.h core/units.h core/residues.h)

target_include_directories(core
    PUBLIC core)

target_link_libraries(core
    PRIVATE Eigen3::Eigen)