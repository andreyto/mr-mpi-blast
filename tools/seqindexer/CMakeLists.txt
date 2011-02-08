cmake_minimum_required(VERSION 2.6)
project(SEQINDEXER)

#set(CMAKE_C_FLAGS "-g -Wall")
add_executable(seqindexer seqindexer.c)

SET(CMAKE_BUILD_TYPE distribtion)
SET(CMAKE_CXX_FLAGS_DISTRIBUTION "-O3")
SET(CMAKE_C_FLAGS_DISTRIBUTION "-O3")

#SET(CMAKE_BUILD_TYPE Debug)
#SET(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -Wall")
#SET(CMAKE_C_FLAGS_DEBUG "-g -O0 -Wall")
