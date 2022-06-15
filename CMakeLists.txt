cmake_minimum_required(VERSION 3.10.2)
project(pace_2022_HUST_solver)

set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O1")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} - O3")

set(CMAKE_CXX_STANDARD 17)

link_directories(${PROJECT_SOURCE_DIR}/lib)
add_executable(PACE22_FVSP_HUST_NOSCP src/main.cpp src/SA.cpp src/SA.h src/LocalSearch.cpp src/LocalSearch.h src/maxVerVC.h src/minVerVC.h src/Utility.cpp src/Utility.h src/datastuct.h src/global.h src/config.cpp src/config.h)
target_link_libraries(${PROJECT_NAME} VerVCSolver.a)