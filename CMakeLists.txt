######################################################################
## CMakeLists.txt --- gsPreCICE
## This file is part of the G+Smo library.
##
######################################################################
if(POLICY CMP0076)# CMake 3.13
cmake_policy(SET CMP0076 NEW)
endif()


## gsPreCICE module
project(gsPreCICE)

# Apply G+Smo config
include(gsConfig)

# Collect file names
aux_header_directory(${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_NAME}_HEADERS)
#aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_NAME}_SOURCES)
#aux_tmpl_header_directory(${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_NAME}_HPPFILES)

# Apply same configuration as G+Smo
include(gsConfig)

find_package(precice REQUIRED CONFIG)

add_library(${PROJECT_NAME} INTERFACE)
# target_sources(${PROJECT_NAME} INTERFACE $<BUILD_INTERFACE:${${PROJECT_NAME}_HEADERS}>)
# target_link_libraries(${PROJECT_NAME} PRIVATE precice::precice)
set_property(TARGET ${PROJECT_NAME} APPEND PROPERTY INTERFACE_LINK_LIBRARIES precice::precice)

target_include_directories(${PROJECT_NAME} INTERFACE
 $<INSTALL_INTERFACE:gismo/gsPreCICE>
)

add_dependencies(${PROJECT_NAME} precice)

set(gismo_LINKER ${gismo_LINKER} ${PROJECT_NAME}
   CACHE INTERNAL "${PROJECT_NAME} extra linker objects")

install(DIRECTORY ${PROJECT_SOURCE_DIR}
        DESTINATION include/gismo
        FILES_MATCHING PATTERN "gsPreCICE/*.h")


# add filedata folder
add_definitions(-DELAST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/filedata/")


# add unittests
aux_gs_cpp_directory(${PROJECT_SOURCE_DIR}/unittests unittests_SRCS)
set(gismo_UNITTESTS ${gismo_UNITTESTS} ${unittests_SRCS}
  CACHE INTERNAL "gismo list of unittests")

# Look for CMakeLists.txt in examples
# add example files
if(GISMO_BUILD_EXAMPLES)
    add_custom_target(${PROJECT_NAME}-examples)
    add_subdirectory(examples)
else()
    add_subdirectory(examples EXCLUDE_FROM_ALL)
endif(GISMO_BUILD_EXAMPLES)

