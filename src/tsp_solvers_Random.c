/*
 * \brief   Random solution.
 * \authors Francesco Cazzaro, Marco Cieno
 */
#include <errno.h>
#include <limits.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <cplex.h>

#include "tsp_solvers.h"
#include "logging.h"
#include "tsp.h"


void
Random_model ( instance *problem )
{
    struct timespec start, end;
    clock_gettime( CLOCK_MONOTONIC, &start );

    for ( size_t i = 0; i < problem->nnodes; ++i )
    {
        problem->solution[i][0] = i;
        problem->solution[i][1] = i + 1;
    }

    problem->solution[problem->nnodes - 1][1] = 0;

    clock_gettime( CLOCK_MONOTONIC, &end );

    problem->elapsedtime  = ( end.tv_sec - start.tv_sec ) + ( end.tv_nsec - start.tv_nsec ) / 1000000000.;
    problem->visitednodes = 0;
    problem->solcost      = compute_solution_cost( problem );
}
