cmake_minimum_required(VERSION 2.8.11)

enable_language(Fortran)
project(progs)

# Do not allow in-source builds:
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
  message(FATAL_ERROR "In-source builds are not allowed.")
endif()

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)

set(CMAKE_MODULE_PATH ${progs_SOURCE_DIR}/cmake)

set(CMAKE_INSTALL_PREFIX "/home/jvcosel/bin/progs")

add_subdirectory(source)


