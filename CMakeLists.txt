cmake_minimum_required(VERSION 3.16)
project(scip_branching_demo LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(scip_branching_demo src/branching_rule_demo.cpp)

# The official SCIP Optimization Suite .deb installs headers into
# /usr/include/scip and the versioned shared library
# /usr/lib/libscip.so.10.1.0 (there is no unversioned libscip.so symlink).
find_path(SCIP_INCLUDE_DIR NAMES scip/scip.h)
if(NOT SCIP_INCLUDE_DIR)
  message(FATAL_ERROR "SCIP headers (scip/scip.h) not found")
endif()

file(GLOB SCIP_LIB_CANDIDATES
     "/usr/lib/libscip.so*"
     "/usr/lib/*/libscip.so*")
if(NOT SCIP_LIB_CANDIDATES)
  message(FATAL_ERROR "libscip shared library not found")
endif()
list(GET SCIP_LIB_CANDIDATES 0 SCIP_LIBRARY)

message(STATUS "SCIP headers : ${SCIP_INCLUDE_DIR}")
message(STATUS "SCIP library : ${SCIP_LIBRARY}")
target_include_directories(scip_branching_demo PRIVATE "${SCIP_INCLUDE_DIR}")
target_link_libraries(scip_branching_demo PRIVATE "${SCIP_LIBRARY}")

include(CTest)
if(BUILD_TESTING)
  add_test(NAME scip_branching_demo COMMAND scip_branching_demo)
endif()
