cmake_minimum_required (VERSION 2.8.11)
project (HEPLike)
enable_testing()


list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS}/build)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} $ENV{ROOTSYS}/etc/cmake)
set(CMAKE_MACOSX_RPATH TRUE) # use @rpath keyword for the macOS rpath

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
include_directories(${CMAKE_SOURCE_DIR}/include ${ROOT_INCLUDE_DIRS})
add_definitions(${ROOT_CXX_FLAGS})


###########################
#
#    Libraries
#
###########################
file(GLOB HL_SOURCES ${CMAKE_SOURCE_DIR}/src/*.cc)
add_library(HEPLike_static ${HL_SOURCES})
add_library(HEPLike_shared SHARED ${HL_SOURCES})
set_target_properties(HEPLike_static PROPERTIES OUTPUT_NAME HEPLike)
set_target_properties(HEPLike_shared PROPERTIES OUTPUT_NAME HEPLike)
target_link_libraries(HEPLike_shared PRIVATE ${ROOT_LIBRARIES} ${Boost_LIBRARIES} ${GSL_LIBRARIES} )


###########################
#
#    Executables
#
###########################

if(EXECUTABLE)
  
  execute_process (
    COMMAND bash -c "ln -s ../data_toy/ data"
    )

  add_executable(Br_example main/Br_example.cc)
  add_executable(BrBifurGaussian_example main/BrBifurGaussian_example.cc)
  add_executable(Limit_example main/Limit_example.cc)
  add_executable(Ndim_Gaussian main/Ndim_Gaussian.cc)
  add_executable(ProfLikelihood_example main/ProfLikelihood_example.cc)
  add_executable(Ndim_BifurGaussian_example main/Ndim_BifurGaussian_example.cc)
  add_executable(Ndim_Likelihood_example main/Ndim_Likelihood_example.cc)
  add_executable(Data_Fit_example  main/Data_Fit_example.cc)
  add_executable(Test_YAML  main/Test_YAML.cc)

  target_link_libraries(Br_example HEPLike_shared yaml-cpp ${ROOT_LIBRARIES} ${Boost_LIBRARIES} )
  target_link_libraries(BrBifurGaussian_example HEPLike_static yaml-cpp  ${ROOT_LIBRARIES} ${Boost_LIBRARIES} ${GSL_LIBRARIES} )
  target_link_libraries(Limit_example HEPLike_static yaml-cpp  ${ROOT_LIBRARIES} ${Boost_LIBRARIES} ${GSL_LIBRARIES})
  target_link_libraries(Ndim_Gaussian HEPLike_static yaml-cpp  ${ROOT_LIBRARIES} ${Boost_LIBRARIES}  ${GSL_LIBRARIES})  
  target_link_libraries(ProfLikelihood_example HEPLike_static yaml-cpp   ${ROOT_LIBRARIES} ${Boost_LIBRARIES}  ${GSL_LIBRARIES})  
  target_link_libraries(Ndim_BifurGaussian_example HEPLike_static yaml-cpp  ${ROOT_LIBRARIES} ${Boost_LIBRARIES}  ${GSL_LIBRARIES})  
  target_link_libraries(Ndim_Likelihood_example  HEPLike_static yaml-cpp  ${ROOT_LIBRARIES} ${Boost_LIBRARIES}  ${GSL_LIBRARIES})  
  target_link_libraries(Data_Fit_example HEPLike_static yaml-cpp  ${ROOT_LIBRARIES} ${Boost_LIBRARIES}  ${GSL_LIBRARIES})  
  target_link_libraries(Test_YAML  HEPLike_static yaml-cpp  ${ROOT_LIBRARIES} ${Boost_LIBRARIES}  ${GSL_LIBRARIES})  

  add_test(HL_Test_YAML Test_YAML ../data_toy/examples/test_3dimassym.yaml)
  add_test(HL_Limit Limit_example xterm)
  add_test(HL_Br_example Br_example)
  add_test(HL_BrBifurGaussian_example  BrBifurGaussian_example)
  add_test(HL_Ndim_Gaussian Ndim_Gaussian)
  add_test(HL_ProfLikelihood_example ProfLikelihood_example ${CMAKE_CURRENT_BINARY_DIR})
  add_test(HL_Ndim_BifurGaussian_example Ndim_BifurGaussian_example)
  

  
endif()
