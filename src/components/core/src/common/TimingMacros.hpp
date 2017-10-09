#ifndef SRC_COMPONENTS_CORE_SRC_COMMON_TIMING_MACROS_HPP_
#define SRC_COMPONENTS_CORE_SRC_COMMON_TIMING_MACROS_HPP_

#ifndef USE_CALIPER
#define USE_CALIPER 0
#endif

#if USE_CALIPER
#include <caliper/cali.h>

#define GEOS_MARK_FUNCTION CALI_CXX_MARK_FUNCTION

#define DO_STRINGIFY(arg) #arg
#define GEOS_CXX_MARK_LOOP_BEGIN(loop, loopName) CALI_CXX_MARK_LOOP_BEGIN(loop,DO_STRINGIFY(loopName))
#define GEOS_CXX_MARK_LOOP_END(loop) CALI_CXX_MARK_LOOP_END(loop)

#define GEOS_MARK_BEGIN(name) CALI_MARK_BEGIN(DO_STRINGIFY(name))
#define GEOS_MARK_END(name) CALI_MARK_END(DO_STRINGIFY(name))

#else

#define GEOS_CXX_MARK_LOOP_BEGIN(loop, loopName) 
#define GEOS_CXX_MARK_LOOP_END(loop)
#define GEOS_MARK_BEGIN(name) 
#define GEOS_MARK_END(name)
#endif





#endif


