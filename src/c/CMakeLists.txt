# src/c/CMakeLists.txt
cmake_minimum_required(VERSION 3.10)
project(tri_affine C)
set(CMAKE_C_STANDARD 11)
add_executable(tri_test tri_tri_affine.c tri_test.c)
