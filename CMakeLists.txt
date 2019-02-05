cmake_minimum_required (VERSION 2.8.11)
project (HEPLike)

###########################
#
#    Packages search
#
###########################
# looking for yaml-cpp
find_package(Yaml-cpp)
# including yaml-cpp headers
include_directories(${CMAKE_SOURCE_DIR}/yaml-cpp/include)
# looking for boost
find_package(Boost 1.41)
find_package(GSL REQUIRED) 

include_directories (include)


# creating own libraries:
add_library(HEPStatistics src/HEPStats.cc)
add_library(HEPData src/HEPData.cc)
add_library(HEPBRLimit src/HEPBRLimit.cc) #



#############################################
# the executables:
add_executable(HEPLike main/main.cc)
target_link_libraries(HEPLike HEPBRLimit HEPData HEPStatistics)


#############################################
# Linking external libraries:
target_link_libraries(${PROJECT_NAME} yaml-cpp GSL::gsl)

