cmake_minimum_required(VERSION 3.14)
project(21s-pa03-christophera)
set (EXE_NAME 21s-pa03-christophera)

set (input01 "data/input.txt")
foreach(file IN LISTS input01)
    configure_file(${file} ${file} COPYONLY)
endforeach()

set(Boost_INCLUDE_DIR boost_1_66_0)
find_package(Boost)
include_directories(${Boost_INCLUDE_DIR})
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall -Wextra")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

add_executable(${EXE_NAME} src/main.cpp src/graphing.cpp src/graphing.h)