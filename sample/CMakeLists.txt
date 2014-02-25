
add_library(loader ../src/plyfile.cpp loader.cpp)

add_executable(main main.cpp)
target_link_libraries(main ccd loader)
