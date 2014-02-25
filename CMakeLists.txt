cmake_minimum_required(VERSION 2.8.8)

project(CCD)

# FindOpenMP
find_package(OpenMP)
if(OPENMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

include_directories(inc)

add_subdirectory(src)

option(BUILD_SAMPLE_APP "Build Sample Application." OFF)
if(BUILD_SAMPLE_APP)
  add_subdirectory(sample)
endif(BUILD_SAMPLE_APP)
