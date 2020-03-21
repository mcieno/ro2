/*
 * \brief   Implementation of tsp_solvers.h functions.
 * \authors Francesco Cazzaro, Marco Cieno
 */
#include <errno.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <ilcplex/cplex.h>

#include "logging.h"
#include "tsp.h"


void
dummy_solution ( instance *problem )
{
    for ( unsigned long i = 0; i < problem->nnodes; ++i ) {
        problem->solution[i] = i;
    }
}


void
dummy_cplex_solution ( instance *problem )
{
	int error;

	CPXENVptr env = CPXopenCPLEX( &error );
	CPXLPptr lp = CPXcreateprob( env, &error, problem->name ? problem->name : "TSP" );

    CPXfreeprob( env, &lp );
    CPXcloseCPLEX( &env );

}
