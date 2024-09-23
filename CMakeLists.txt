cmake_minimum_required(VERSION 3.13)

set(CMAKE_CXX_STANDARD 11)

project(quasar)

add_executable(quasar
  ${CMAKE_SOURCE_DIR}/src/Quasar.cpp
  ${CMAKE_SOURCE_DIR}/src/ProgramData.cpp
  ${CMAKE_SOURCE_DIR}/src/Data.cpp
  ${CMAKE_SOURCE_DIR}/src/Geno.cpp
  ${CMAKE_SOURCE_DIR}/src/Utils.cpp
)
target_include_directories(quasar PRIVATE ${CMAKE_SOURCE_DIR}/src)

set(EXTERN_LIBS_PATH "${CMAKE_SOURCE_DIR}/external_libs")

target_include_directories(quasar PRIVATE ${EXTERN_LIBS_PATH}/cxxopts/include/)
target_include_directories(quasar PRIVATE ${EXTERN_LIBS_PATH}/eigen-3.4.0/)

