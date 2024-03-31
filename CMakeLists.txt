cmake_minimum_required( VERSION 2.8 )
project( st_ekf_proj )

set( CMAKE_BUILD_TYPE "Release" )
#set( CMAKE_CXX_FLAGS "-O3" )
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
if(COMPILER_SUPPORTS_CXX11)
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11" )
endif()
# 添加Eigen头文件
find_package(Eigen3)
include_directories(${EIGEN3_INCLUDE_DIR})
# include_directories( "/usr/include/eigen3" )

add_executable(estimator estimator.cpp rotation.cpp)