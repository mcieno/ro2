/*
 * \brief   Implementation of tsp_solvers.h functions.
 * \authors Francesco Cazzaro, Marco Cieno
 */
#include <errno.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <cplex.h>

#include "tsp_solvers.h"
#include "logging.h"
#include "tsp.h"


#define MAX_CNAME_LENGTH 128

/*!
 * \brief Retrieve the solution after CPXmimopt was run.
 *
 *
 * \param env
 *      cplex env parameter
 *
 * \param lp
 *      cplex lp parameter
 *
 * \param problem
 *     Pointer to the instance structure
 */
void
infer_cplex_solution ( CPXENVptr env, CPXLPptr lp, instance *problem )
{
	double sol[CPXgetnumcols( env, lp )];
	CPXsolution( env, lp, NULL, NULL, sol, NULL, NULL, NULL );

	unsigned long p = 0;

	for ( unsigned long i = 0; i < problem->nnodes; ++i ) {
		for ( unsigned long j = i + 1; j < problem->nnodes; ++j ) {
			if ( sol[xpos( i, j, problem )] == 1 )
			{
				problem->solution[p][0] = i;
				problem->solution[p][1] = j;
				++p;
			}
		}
	}
}


void
dummy_solution ( instance *problem )
{
    for ( unsigned long i = 0; i < problem->nnodes; ++i ) {
        problem->solution[i][0] = i;
		problem->solution[i][1] = i + 1;
    }

	problem->solution[problem->nnodes - 1][1] = 0;
}


int
xpos ( unsigned long i, unsigned long j, instance *problem )
{
	if ( i == j ) {
		perror( " i == j in xpos" );
		exit( EXIT_FAILURE );
	}

	if ( i > j ) return xpos( j, i, problem );

	return i * problem->nnodes + j - ( ( i + 1 ) * ( i + 2 ) / 2UL );
}


void
dummy_build_cplex_model ( instance *problem, CPXENVptr env, CPXLPptr lp )
{
	char binary = 'B';
	double lb = 0.0;
	double ub = 1.0;

	char *cname = calloc( MAX_CNAME_LENGTH, sizeof(  *cname ) );

	// add binary var.s x(i,j) for i < j
	for ( unsigned long i = 0; i < problem->nnodes; ++i )
	{
		for ( unsigned long j = i + 1; j < problem->nnodes; ++j )
		{
			snprintf( cname, MAX_CNAME_LENGTH, "x(%lu,%lu)", i + 1, j + 1 );
			double obj = _euclidean_distance(
				problem->xcoord[i],
				problem->ycoord[i],
				problem->xcoord[j],
				problem->ycoord[j]
			);

			if ( CPXnewcols( env, lp, 1, &obj, &lb, &ub, &binary, &cname ) ) {
				perror( "Wrong CPXnewcols on x var.s" );
				exit( EXIT_FAILURE );
			}

    		if ( CPXgetnumcols( env, lp ) - 1 != xpos( i, j, problem ) ) {
				perror( "Wrong position for x var.s" );
				exit( EXIT_FAILURE );
			}
		}
	}

	// add the degree constraints
	double rhs = 2.0;
	char sense = 'E';

	for ( unsigned long h = 0; h < problem->nnodes; ++h )
	{
		unsigned long lastrow = CPXgetnumrows( env, lp );

		snprintf( cname, MAX_CNAME_LENGTH, "degree(%lu)", h + 1 );
		if ( CPXnewrows( env, lp, 1, &rhs, &sense, NULL, &cname ) ) {
			perror( "Wrong CPXnewrows [degree]" );
			exit( EXIT_FAILURE );
		}

		for ( unsigned long i = 0; i < problem->nnodes; ++i )
		{
			if ( i == h ) continue;
			if ( CPXchgcoef( env, lp, lastrow, xpos( i, h, problem ), 1.0 ) ) {
				perror( "Wrong CPXchgcoef [degree]" );
				exit( EXIT_FAILURE );
			}
		}
	}

	free( cname );
}

void
dummy_cplex_solution ( instance *problem )
{
	int error;

	CPXENVptr env = CPXopenCPLEX( &error );
	CPXLPptr lp = CPXcreateprob( env, &error, problem->name ? problem->name : "TSP" );

    dummy_build_cplex_model( problem, env, lp );

	//CPXwriteprob (env, lp, "bin/myprob.lp", NULL);

    if ( CPXmipopt( env, lp ) ) {
		fprintf( stderr, "CPXmimopt true\n" );
	}

	infer_cplex_solution( env, lp, problem );

    CPXfreeprob( env, &lp );
    CPXcloseCPLEX( &env );
}
