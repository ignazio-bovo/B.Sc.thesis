CMAKE_MINIMUM_REQUIRED(VERSION 3.2)
PROJECT(significant_subgraph_mining VERSION 1.0 LANGUAGES CXX)
# INCLUDE_DIRECTORIES(include)
# INCLUDE_DIRECTORIES(${INCLUDE_DIRECTORIES} src)

# Set up directories
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${PROJECT_SOURCE_DIR}/src)
set (MAIN_FILE ${PROJECT_SOURCE_DIR}/main.cpp)

link_directories(${PROJECT_SOURCE_DIR}/include ${PROJECT_SOURCE_DIR}/src)

file(GLOB SOURCE_FILES ${MAIN_FILE} ${PROJECT_SOURCE_DIR}/include/*.h  ${PROJECT_SOURCE_DIR}/src/*.cpp ${PROJECT_SOURCE_DIR}/src/*.hpp)

# target_include_directories (testable PUBLIC "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>")
link_directories (${LINK_DIRECTORIES} "${INCLUDE_DIRECTORIES}")

add_executable (testable ${SOURCE_FILES})
# set (PROJECT_SOURCE_DECLARATION_DIRECTORY ${PROJECT_SOURCE_DIR}/include)
# set (PROJECT_SOURCE_DEFINITION_DIRECTORY ${PROJECT_SOURCE_DIR}/src)
set (MAIN_FILE ${PROJECT_SOURCE_DIR}/main.cpp)

# enviroinment
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/build)

# compiler-related
set (CMAKE_CXX_FLAGS "-Wall -O3 -g -std=c++1y")
# target_compile_features(lboost_program_options)

# Load settings for the external project "Boost".
# The components "system" and "filesystem" from "Boost" version 1.55.0 or
# greater are required.
find_package (Boost 1.64.0)
if (Boost_FOUND)
  INCLUDE_DIRECTORIES (${INCLUDE_DIRECTORIES} ${Boost_INCLUDE_DIR})
  LINK_DIRECTORIES (${LINK_DIRECTORIES} ${Boost_LIBRARY_DIR})
endif ()

# "Link a target to given libraries."
target_link_libraries (testable ${Boost_LIBRARIES})


