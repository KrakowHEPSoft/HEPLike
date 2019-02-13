cmake_minimum_required (VERSION 2.8.11)
project (HEPLike)
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS}/build)

###########################
# 
#    Compile seetings
#
###########################
# Add -fPIC for 64 bit systems
if(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
endif()
if(COMPILER_SUPPORTS_CXX11)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
endif()



###########################
#
#    Packages search
#
###########################
# looking for yaml-cpp
find_package(Yaml-cpp)
# looking for boost
find_package(Boost REQUIRED)   
# GSL
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
add_executable(Ndim_Gaussian main/Ndim_Gaussian.cc)
add_executable(ProfLikelihood_example main/ProfLikelihood_example.cc)
add_executable(Ndim_BifurGaussian_example main/Ndim_BifurGaussian_example.cc)
add_executable(Ndim_Likelihood_example main/Ndim_Likelihood_example.cc)



target_link_libraries(Br_example HL_Limit HL_Data HL_Statistics HL_Gaussian HL_BifurGaussian HL_ProfLikelihood HL_nDimGaussian HL_nDimBifurGaussian HL_nDimLikelihood yaml-cpp GSL::gsl ${ROOT_LIBRARIES} ${Boost_LIBRARIES} )                         

target_link_libraries(BrBifurGaussian_example HL_Limit HL_Data HL_Statistics HL_Gaussian HL_BifurGaussian HL_ProfLikelihood HL_nDimGaussian HL_nDimBifurGaussian HL_nDimLikelihood yaml-cpp GSL::gsl ${ROOT_LIBRARIES} ${Boost_LIBRARIES} )

target_link_libraries(Limit_example HL_Limit HL_Data HL_Statistics HL_Gaussian HL_BifurGaussian HL_ProfLikelihood HL_nDimGaussian HL_nDimBifurGaussian HL_nDimLikelihood yaml-cpp GSL::gsl ${ROOT_LIBRARIES} ${Boost_LIBRARIES} )

target_link_libraries(Ndim_Gaussian HL_Limit HL_Data HL_Statistics HL_Gaussian HL_BifurGaussian HL_ProfLikelihood HL_nDimGaussian HL_nDimBifurGaussian HL_nDimLikelihood yaml-cpp GSL::gsl ${ROOT_LIBRARIES} ${Boost_LIBRARIES} )     

target_link_libraries(ProfLikelihood_example HL_Limit HL_Data HL_Statistics HL_Gaussian HL_BifurGaussian HL_ProfLikelihood HL_nDimGaussian HL_nDimBifurGaussian HL_nDimLikelihood yaml-cpp GSL::gsl ${ROOT_LIBRARIES} ${Boost_LIBRARIES} )     

target_link_libraries(Ndim_BifurGaussian_example HL_Limit HL_Data HL_Statistics HL_Gaussian HL_BifurGaussian HL_ProfLikelihood HL_nDimGaussian HL_nDimBifurGaussian HL_nDimLikelihood yaml-cpp GSL::gsl ${ROOT_LIBRARIES} ${Boost_LIBRARIES} )  

target_link_libraries(Ndim_Likelihood_example  HL_Limit HL_Data HL_Statistics HL_Gaussian HL_BifurGaussian HL_ProfLikelihood HL_nDimGaussian HL_nDimBifurGaussian HL_nDimLikelihood yaml-cpp GSL::gsl ${ROOT_LIBRARIES} ${Boost_LIBRARIES} )  


#target_link_libraries(${PROJECT_NAME} HL_Limit HL_Data HL_Statistics HL_Gaussian HL_BifurGaussian HL_ProfLikelihood HL_nDimGaussian HL_nDimBifurGaussian HL_nDimLikelihood )


#############################################
# Linking external libraries:

#target_link_libraries(${PROJECT_NAME} yaml-cpp GSL::gsl ${ROOT_LIBRARIES} ${Boost_LIBRARIES})

