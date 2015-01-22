cmake_minimum_required(VERSION 2.8)
PROJECT(actlearn)

#if you don't want the full compiler output, remove the following line
SET(CMAKE_VERBOSE_MAKEFILE OFF)

FIND_PACKAGE(ITK REQUIRED)
IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ENDIF(ITK_FOUND)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()


include_directories(${PROJECT_SOURCE_DIR})
set(Boost_DEBUG 0)
set(Boost_USE_STATIC_LIBS on)
set(Boost_USE_MULTITHREADED ON) 
set(Boost_USE_STATIC_RUNTIME OFF)
find_package(Boost 1.53.0 REQUIRED COMPONENTS program_options filesystem system)

if(Boost_FOUND)
  MESSAGE(STATUS "Project source dir: " ${PROJECT_SOURCE_DIR})
  include_directories(${Boost_INCLUDE_DIRS})
  add_executable(activeCutSeg
    grabcutseg.cxx
    loadfiles.cxx
    gmm.cxx
    query.cxx
    graphcuts.cxx
    graph.cpp
    maxflow.cpp
    utility.cxx
    )

  target_link_libraries(activeCutSeg ${ITK_LIBRARIES} ${Boost_LIBRARIES})

endif()

