cmake_minimum_required(VERSION 3.14)

project(multifile LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++20")

add_executable(multifile     
  lfsr.hpp
  hamming.hpp  
  gf.hpp
  gf.cpp
  rsexh.hpp
  rsexh.cpp
  main.cpp
)
