cmake_minimum_required(VERSION 3.16)
project(21s-pa03-christophera.exe)

set(CMAKE_CXX_STANDARD 14)
set(Boost_INCLUDE_DIR boost_1_66_0)

find_package(Boost)
include_directories(${Boost_INCLUDE_DIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall -Wextra")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

add_executable(21s-pa03-christophera src/main.cpp src/graphing.cpp src/graphing.h)