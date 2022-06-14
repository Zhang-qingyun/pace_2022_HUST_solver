cmake_minimum_required(VERSION 3.10.2)
project(PACE)

set(CMAKE_CXX_STANDARD 14)

set(CMAKE_EXE_LINKER_FLAGS "-static")

add_executable(PACE22_FVSP_HUST_NOSCP main.cpp SA.cpp SA.h LocalSearch.cpp LocalSearch.h maxVerVC.h minVerVC.h Utility.cpp Utility.h datastuct.h global.h config.cpp config.h)
