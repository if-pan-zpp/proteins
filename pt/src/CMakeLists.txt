add_executable(pt)

set(Sources main.cpp pdb/ir/parser.cpp pdb/ir/auto_line_parser.cpp
    util/str.cpp conf/legacy/parser.cpp util/field.cpp)
set(Headers pdb/ir/parser.h pdb/ir/line_parser.h pdb/ir/auto_line_parser.h
    util/str.h conf/legacy/parser.h util/field.h def/types.h def/units.h conf/constants.h)
set(ProtoFiles pdb/ir/ir.proto conf/legacy/legacy.proto)
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
