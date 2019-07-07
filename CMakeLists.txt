cmake_minimum_required (VERSION 2.8.11)
project (HEPLike)
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS}/build)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} $ENV{ROOTSYS}/etc/cmake)

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



add_library(HEPLike src/HL_Stats.cc src/HL_Data.cc src/HL_Limit.cc src/HL_Gaussian.cc src/HL_BifurGaussian.cc src/HL_ProfLikelihood.cc src/HL_nDimGaussian.cc src/HL_nDimBifurGaussian.cc src/HL_nDimLikelihood.cc src/HL_ExpPoints.cc)


add_library(HEPLikeSO SHARED src/HL_Stats.cc src/HL_Data.cc src/HL_Limit.cc src/HL_Gaussian.cc src/HL_BifurGaussian.cc src/HL_ProfLikelihood.cc src/HL_nDimGaussian.cc src/HL_nDimBifurGaussian.cc src/HL_nDimLikelihood.cc src/HL_ExpPoints.cc ${ROOT_LIBRARIES} ${Boost_LIBRARIES})


#############################################
# the executables:
#add_executable(Br_example main/Br_example.cc)
#add_executable(BrBifurGaussian_example main/BrBifurGaussian_example.cc)
#add_executable(Limit_example main/Limit_example.cc)
#add_executable(Ndim_Gaussian main/Ndim_Gaussian.cc)
#add_executable(ProfLikelihood_example main/ProfLikelihood_example.cc)
#add_executable(Ndim_BifurGaussian_example main/Ndim_BifurGaussian_example.cc)
#add_executable(Ndim_Likelihood_example main/Ndim_Likelihood_example.cc)
#add_executable(Data_Fit_example  main/Data_Fit_example.cc)
#add_executable(Test_YAML  main/Test_YAML.cc)


#target_link_libraries(Br_example HEPLike yaml-cpp GSL::gsl ${ROOT_LIBRARIES} ${Boost_LIBRARIES} )   
#target_link_libraries(BrBifurGaussian_example HEPLike yaml-cpp GSL::gsl ${ROOT_LIBRARIES} ${Boost_LIBRARIES} )
#target_link_libraries(Limit_example HEPLike yaml-cpp GSL::gsl ${ROOT_LIBRARIES} ${Boost_LIBRARIES} )
#target_link_libraries(Ndim_Gaussian HEPLike yaml-cpp GSL::gsl ${ROOT_LIBRARIES} ${Boost_LIBRARIES} )
#target_link_libraries(ProfLikelihood_example HEPLike yaml-cpp GSL::gsl ${ROOT_LIBRARIES} ${Boost_LIBRARIES} )
#target_link_libraries(Ndim_BifurGaussian_example HEPLike yaml-cpp GSL::gsl ${ROOT_LIBRARIES} ${Boost_LIBRARIES} )
#target_link_libraries(Ndim_Likelihood_example  HEPLike yaml-cpp GSL::gsl ${ROOT_LIBRARIES} ${Boost_LIBRARIES} )  
#target_link_libraries(Data_Fit_example HEPLike yaml-cpp GSL::gsl ${ROOT_LIBRARIES} ${Boost_LIBRARIES} )
#target_link_libraries(Test_YAML  HEPLike yaml-cpp GSL::gsl ${ROOT_LIBRARIES} ${Boost_LIBRARIES} ) 


target_link_libraries(HEPLikeSO PRIVATE GSL::gsl ${ROOT_LIBRARIES} ${Boost_LIBRARIES} )

