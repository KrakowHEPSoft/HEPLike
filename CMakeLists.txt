cmake_minimum_required (VERSION 2.8.11)
project (HL_Like)
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS}/build)


###########################
#
#    Packages search
#
###########################
# looking for yaml-cpp
find_package(Yaml-cpp)
# including yaml-cpp headers
#include_directories(${CMAKE_SOURCE_DIR}/yaml-cpp/include)
# looking for boost
#include_directories(/opt/pkgs/boost_1_69/include)
#link_directories(/opt/pkgs/boost_1_69/lib)
#find_package(Boost 1.69 REQUIRED)
find_package(Boost REQUIRED)   
#include_directories(${Boost_INCLUDE_DIR})
# GLS
find_package(GSL REQUIRED) 
# ROOT:
find_package(ROOT REQUIRED COMPONENTS MathCore RIO Hist Tree Net)


include(${ROOT_USE_FILE})
include_directories(${CMAKE_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})
add_definitions(${ROOT_CXX_FLAGS})


include_directories (include)


# creating own libraries:
add_library(HL_Statistics src/HL_Stats.cc)
add_library(HL_Data src/HL_Data.cc)
add_library(HL_Limit src/HL_Limit.cc) #
add_library(HL_Gaussian src/HL_Gaussian.cc)
add_library(HL_BifurGaussian src/HL_BifurGaussian.cc)
add_library(HL_ProfLikelihood src/HL_ProfLikelihood.cc)
add_library(HL_nDimGaussian src/HL_nDimGaussian.cc)
add_library(HL_nDimBifurGaussian src/HL_nDimBifurGaussian.cc) 
add_library(HL_nDimLikelihood src/HL_nDimLikelihood.cc)

#############################################
# the executables:
add_executable(Br_example main/Br_example.cc)
add_executable(BrBifurGaussian_example main/BrBifurGaussian_example.cc)
add_executable(Limit_example main/Limit_example.cc)


target_link_libraries(Br_example HL_Limit HL_Data HL_Statistics HL_Gaussian HL_BifurGaussian HL_ProfLikelihood HL_nDimGaussian HL_nDimBifurGaussian HL_nDimLikelihood yaml-cpp GSL::gsl ${ROOT_LIBRARIES} ${Boost_LIBRARIES} )                         

target_link_libraries(BrBifurGaussian_example HL_Limit HL_Data HL_Statistics HL_Gaussian HL_BifurGaussian HL_ProfLikelihood HL_nDimGaussian HL_nDimBifurGaussian HL_nDimLikelihood yaml-cpp GSL::gsl ${ROOT_LIBRARIES} ${Boost_LIBRARIES} )

target_link_libraries(Limit_example HL_Limit HL_Data HL_Statistics HL_Gaussian HL_BifurGaussian HL_ProfLikelihood HL_nDimGaussian HL_nDimBifurGaussian HL_nDimLikelihood yaml-cpp GSL::gsl ${ROOT_LIBRARIES} ${Boost_LIBRARIES} )




#target_link_libraries(${PROJECT_NAME} HL_Limit HL_Data HL_Statistics HL_Gaussian HL_BifurGaussian HL_ProfLikelihood HL_nDimGaussian HL_nDimBifurGaussian HL_nDimLikelihood )


#############################################
# Linking external libraries:

#target_link_libraries(${PROJECT_NAME} yaml-cpp GSL::gsl ${ROOT_LIBRARIES} ${Boost_LIBRARIES})

