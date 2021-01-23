add_executable(pt)

FILE(GLOB Sources *.cpp pdb/ir/*.cpp util/*.cpp potentials/*.cpp)
FILE(GLOB Headers *.hpp pdb/ir/*.hpp util/*.cpp conf/*.hpp def/*.hpp potentials/*.hpp)
set(ProtoFiles pdb/ir/ir.proto)
set(Resources 1ubq.pdb inputfile1)

foreach(Proto IN LISTS ProtoFiles)
    protobuf_generate_cpp(ProtoSources ProtoHeaders ${Proto})
    list(APPEND Sources ${ProtoSources})
    list(APPEND Headers ${ProtoHeaders})
endforeach()

target_sources(pt
    PRIVATE ${Sources}
    PUBLIC ${Headers})

target_link_libraries(pt
    PUBLIC Eigen3::Eigen protobuf::libprotobuf)

target_include_directories(pt
    PUBLIC ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_LIST_DIR})

foreach(Resource IN LISTS Resources)
    configure_file(${Resource} ${Resource} COPYONLY)
endforeach()
