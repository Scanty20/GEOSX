#
# Specify all headers
#
set(mainInterface_headers
    ProblemManager.hpp
    initialization.hpp
    GeosxState.hpp
   )


#
# Specify all sources
#

set(mainInterface_sources
    ProblemManager.cpp
    initialization.cpp
    GeosxState.cpp
   )

set( dependencyList events fileIO optionparser RAJA conduit)


if( ENABLE_MATHPRESSO )
   set( dependencyList ${dependencyList} mathpresso )
endif()

if ( ENABLE_OPENMP )
  set( dependencyList ${dependencyList} openmp )
endif()

if ( ENABLE_CUDA )
  set( dependencyList ${dependencyList} cuda )
endif()

if ( ENABLE_CALIPER )
  set( dependencyList ${dependencyList} caliper adiak )
endif()

if ( ENABLE_MKL )
  set( dependencyList ${dependencyList} mkl )
endif()

blt_add_library( NAME                  mainInterface
                 SOURCES               ${mainInterface_sources}
                 HEADERS               ${mainInterface_headers}
                 DEPENDS_ON            ${dependencyList}
                 OBJECT                ${GEOSX_BUILD_OBJ_LIBS}
               )

target_include_directories( mainInterface PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents)

geosx_add_code_checks( PREFIX mainInterface )
