cmake_minimum_required(VERSION 3.1)

project(mshw_ocgrid)

# https://eigen.tuxfamily.org/dox/TopicCMakeGuide.html
find_package (Eigen3 REQUIRED NO_MODULE)

include_directories(include ${EIGEN3_INCLUDE_DIRS} ../mshw_quatmath/include)

add_executable(mshw_ocgrid src/mshw_ocgridlib_test.cpp)
