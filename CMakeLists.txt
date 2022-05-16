cmake_minimum_required(VERSION 3.21)
project(math_kursovaya)

set(CMAKE_CXX_STANDARD 14)

add_executable(math_kursovaya main.cpp quanc8.h quanc8.cpp rkf45.h rkf45.cpp zeroin.h zeroin.cpp)
