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
#include <sys/timeb.h>

#include <cplex.h>

#include "tsp_solvers.h"
#include "logging.h"
#include "tsp.h"


void
Random_model ( instance *problem )
{
    struct timeb start, end;
    ftime( &start );

    for ( size_t i = 0; i < problem->nnodes; ++i )
    {
        problem->solution[i][0] = i;
        problem->solution[i][1] = i + 1;
    }

    problem->solution[problem->nnodes - 1][1] = 0;

    ftime( &end );

    problem->elapsedtime  = ( 1000. * ( end.time - start.time ) + end.millitm - start.millitm ) / 1000.;
    problem->visitednodes = 0;
    problem->solcost      = compute_solution_cost( problem );
}
