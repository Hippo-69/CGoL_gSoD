#autogenerated file
cmake_minimum_required(VERSION 3.5)
project(cpads)
set(CMAKE_CXX_STANDARD 17)

add_compile_options( -O3 -g )
add_link_options( -lpthread )

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-sign-compare -pthread -std=c++17")

include(CheckLanguage)
check_language(CUDA)
set(CMAKE_CUDA_HOST_COMPILER "${CMAKE_CXX_COMPILER}")

if(CMAKE_CUDA_COMPILER)
enable_language(CUDA)
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --ptxas-options=-v -lineinfo -std=c++17")
endif()

include_directories("include")
add_subdirectory(test)
