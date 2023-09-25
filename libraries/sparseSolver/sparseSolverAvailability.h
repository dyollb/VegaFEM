#ifndef _SPARSESOLVER_AVAILABILITY_H_
#define _SPARSESOLVER_AVAILABILITY_H_

#include "vega-config.h"

#ifdef VEGA_USE_PARDISO
    #define PARDISO_SOLVER_IS_AVAILABLE
#endif

#ifdef VEGA_USE_SPOOLES
    #define SPOOLES_SOLVER_IS_AVAILABLE
#endif

#endif

