cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(BSplineInSE3)
find_package(PCL REQUIRED)
find_package(OpenCV REQUIRED)

find_package(Sophus REQUIRED)
find_package(Eigen3 REQUIRED)

#include_directories(${Sophus_INCLUDE_DIRS})
#target_include_directories(BSplineInSE3 PUBLIC "D:/dependency/Sophus/include/sophus")
#message("sophus_include_dirs: ******************"${Sophus_INCLUDE_DIRS})
#include_directories(${OpenCV_INCLUDE_DIRS})
include_directories( ${PCL_INCLUDE_DIRS})
link_directories( ${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})
include_directories( ${EIGEN3_INCLUDE_DIR})



#set(Sophus_LIBRARIES "/usr/local/include/sophus/se3.hpp")

#set(CMAKE_CXX_FLAGS "-std=c++11 -O2" )

add_executable(BSplineInSE3 main.cpp)
target_link_libraries(BSplineInSE3 
        ${OpenCV_LIBS}
        ${PCL_LIBRARIES}
        Sophus::Sophus
        fmt::fmt

)
#link_directories(${OpenCV_LIBRARY_DIRS})

