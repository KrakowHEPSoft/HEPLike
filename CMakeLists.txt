cmake_minimum_required (VERSION 2.8.11)
project (HEPLike)
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS}/build)


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
# GLS
find_package(GSL REQUIRED) 
# ROOT:
find_package(ROOT REQUIRED COMPONENTS MathCore RIO Hist Tree Net)


include(${ROOT_USE_FILE})
include_directories(${CMAKE_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})
add_definitions(${ROOT_CXX_FLAGS})


include_directories (include)


# creating own libraries:
add_library(HEPStatistics src/HEPStats.cc)
add_library(HEPData src/HEPData.cc)
add_library(HEPBR_Limit src/HEPBR_Limit.cc) #
add_library(HEPBR_Gaussian src/HEPBR_Gaussian.cc)
add_library(HEPBR_BifurGaussian src/HEPBR_BifurGaussian.cc)
add_library(HEPBR_ProfLikelihood src/HEPBR_ProfLikelihood.cc)

#############################################
# the executables:
add_executable(HEPLike main/main.cc)
target_link_libraries(HEPLike HEPBR_Limit HEPData HEPStatistics HEPBR_Gaussian HEPBR_BifurGaussian HEPBR_ProfLikelihood)



#############################################
# Linking external libraries:
target_link_libraries(${PROJECT_NAME} yaml-cpp GSL::gsl ${ROOT_LIBRARIES})

