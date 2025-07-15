cmake_minimum_required(VERSION 3.13)

set(CMAKE_CXX_STANDARD 14)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -fno-math-errno -fno-trapping-math")

project(quasar)

FIND_PACKAGE(Boost)
IF (Boost_FOUND)
    INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
    ADD_DEFINITIONS( "-DHAS_BOOST" )
ENDIF()

add_executable(quasar
  ${CMAKE_SOURCE_DIR}/src/Quasar.cpp
  ${CMAKE_SOURCE_DIR}/src/Data.cpp
  ${CMAKE_SOURCE_DIR}/src/Geno.cpp
  ${CMAKE_SOURCE_DIR}/src/Utils.cpp
  ${CMAKE_SOURCE_DIR}/src/Phi.cpp 
  ${CMAKE_SOURCE_DIR}/src/Residualise.cpp 
  ${CMAKE_SOURCE_DIR}/src/ModelFit.cpp 
  ${CMAKE_SOURCE_DIR}/src/ScoreTest.cpp 
  ${CMAKE_SOURCE_DIR}/src/QTLMappingUtils.cpp
)
target_include_directories(quasar PRIVATE ${CMAKE_SOURCE_DIR}/src)

set(EXTERN_LIBS_PATH "${CMAKE_SOURCE_DIR}/external_libs")

target_include_directories(quasar PRIVATE ${EXTERN_LIBS_PATH}/cxxopts/include/)
target_include_directories(quasar PRIVATE ${EXTERN_LIBS_PATH}/eigen-3.4.0/)
target_include_directories(quasar PRIVATE ${EXTERN_LIBS_PATH}/brent/)
