/*
 * \brief   Basic MTZ model.
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
#include "tspconf.h"


/*!
 * \brief Get the position of variable x(i,j) in MTZ model. 0-indexed.
 *
 *
 * \param i
 *      i in x(i,j)
 *
 * \param j
 *      j in x(i,j)
 *
 * \param problem
 *     Pointer to the instance structure.
 */
size_t
_mtz_xpos ( size_t i, size_t j, const instance *problem )
{
    if ( i == j ) {
        errno = EFAULT;
        perror( CFATAL "_mtz_xpos: i == j" );
        exit( EXIT_FAILURE );
    }

    return i * ( problem->nnodes - 1UL ) + j - ( j > i ? 1UL : 0UL );
}


/*!
 * \brief Get the position of variable u(i) in MTZ model. 0-indexed.
 *
 *
 * \param i
 *      i in u(i)
 *
 * \param problem
 *     Pointer to the instance structure.
 */
size_t
_mtz_upos ( size_t i, const instance *problem )
{
    return problem->nnodes * ( problem->nnodes - 1UL ) + i - 1UL;
}


/*!
 * \brief Add MTZ constraints to the problem.
 *
 *
 * \param problem
 *     Pointer to the instance structure.
 *
 * \param env
 *     CPLEX environment.
 *
 * \param lp
 *     CPLEX problem.
 */
void
_add_constraints_mtz ( const instance *problem, CPXENVptr env, CPXLPptr lp )
{
    char ctype = CPX_BINARY;
    double lb = 0.0;
    double ub = 1.0;

    double obj;
    double rhs;
    char sense;
    size_t lastrow;

    char *cname = calloc( CPX_STR_PARAM_MAX, sizeof( *cname ) );

    // add binary vars x(i,j) for all i, j
    for ( size_t i = 0; i < problem->nnodes; ++i )
    {
        for ( size_t j = 0; j < problem->nnodes; ++j )
        {
            if( i == j ) continue;

            snprintf( cname, CPX_STR_PARAM_MAX, "x(%zu,%zu)", i + 1, j + 1 );
            obj = _euclidean_distance(
                problem->xcoord[i],
                problem->ycoord[i],
                problem->xcoord[j],
                problem->ycoord[j]
            );

            if ( CPXnewcols( env, lp, 1, &obj, &lb, &ub, &ctype, &cname ) ) {
                fprintf( stderr, CFATAL "_add_constraints_mtz: CPXnewcols [%s]\n", cname );
                exit( EXIT_FAILURE );
            }

            if ( CPXgetnumcols( env, lp ) - 1 != _mtz_xpos( i, j, problem ) ) {
                fprintf( stderr, CFATAL "_add_constraints_mtz: CPXgetnumcols [%s: x(%zu, %zu)]\n",
                    cname, i + 1, j + 1 );
                exit( EXIT_FAILURE );
            }
        }
    }

    // add degree constraints
    rhs = 1.0;
    sense = 'E';

    // Sum[ x(i,h) ]_{i | i != h} = 1  for all h
    for ( size_t h = 0; h < problem->nnodes; ++h )
    {
        snprintf( cname, CPX_STR_PARAM_MAX, "degree(%zu)", h + 1 );
        if ( CPXnewrows( env, lp, 1, &rhs, &sense, NULL, &cname ) ) {
            fprintf( stderr, CFATAL "_add_constraints_mtz: CPXnewrows [%s]\n", cname );
            exit( EXIT_FAILURE );
        }

        lastrow = CPXgetnumrows( env, lp ) - 1;

        for ( size_t i = 0; i < problem->nnodes; ++i )
        {
            if ( i == h ) continue;

            if ( CPXchgcoef( env, lp, lastrow, _mtz_xpos( i, h, problem ), 1.0 ) ) {
                fprintf( stderr, CFATAL "_add_constraints_mtz: CPXchgcoef [%s: x(%zu, %zu)]\n",
                    cname, i + 1, h + 1 );
                exit( EXIT_FAILURE );
            }
        }
    }

    // Sum[ x(h,i) ]_{i | i != h} = 1  for all h
    for ( size_t h = 0; h < problem->nnodes; ++h )
    {
        snprintf( cname, CPX_STR_PARAM_MAX, "degree(%zu)", h + 1 );
        if ( CPXnewrows( env, lp, 1, &rhs, &sense, NULL, &cname ) ) {
            fprintf( stderr, CFATAL "_add_constraints_mtz: CPXnewrows [%s]\n", cname );
            exit( EXIT_FAILURE );
        }

        lastrow = CPXgetnumrows( env, lp ) - 1;

        for ( size_t i = 0; i < problem->nnodes; ++i )
        {
            if ( i == h ) continue;

            if ( CPXchgcoef( env, lp, lastrow, _mtz_xpos( h, i, problem ), 1.0 ) ) {
                fprintf( stderr, CFATAL "_add_constraints_mtz: CPXchgcoef [%s: x(%zu, %zu)]\n",
                    cname, h + 1, i + 1 );
                exit( EXIT_FAILURE );
            }
        }
    }

    // add sequences vars u(i)
    ctype = CPX_CONTINUOUS;
    lb = 2.0;
    ub = problem->nnodes;
    obj = 0.;

    for ( size_t i = 1; i < problem->nnodes; ++i )
    {
        snprintf( cname, CPX_STR_PARAM_MAX, "u(%zu)", i + 1 );

        if( CPXnewcols( env, lp, 1, &obj, &lb, &ub, &ctype, &cname ) ) {
            fprintf( stderr, CFATAL "_add_constraints_mtz: CPXnewcols [%s]\n", cname );
            exit( EXIT_FAILURE );
        }

        if ( CPXgetnumcols( env, lp ) - 1 != _mtz_upos( i, problem ) ) {
            fprintf( stderr, CFATAL "_add_constraints_mtz: CPXgetnumcols [%s]\n", cname );
            exit( EXIT_FAILURE );
        }
    }

    // add sequences constraints
    // u(i) - u(j) + n x(i,j) <= n - 1  for all i, j | i != j, i != 0, j != 0
    rhs = problem->nnodes - 1;
    sense = 'L';

    for ( size_t j = 1; j < problem->nnodes; ++j )
    {
        for ( size_t i = 1; i < problem->nnodes; ++i )
        {
            if ( i == j ) continue;

            snprintf( cname, CPX_STR_PARAM_MAX, "sequence(%zu)", j + 1 );

            if ( CPXnewrows( env, lp, 1, &rhs, &sense, NULL, &cname ) ) {
                fprintf( stderr, CFATAL "_add_constraints_mtz: CPXnewrows [%s]\n", cname );
                exit( EXIT_FAILURE );
            }

            lastrow = CPXgetnumrows( env, lp ) - 1;

            if ( CPXchgcoef( env, lp, lastrow, _mtz_upos( j, problem ), -1.0 ) ) {
                fprintf( stderr, CFATAL "_add_constraints_mtz: CPXchgcoef [%s: u(%zu)]\n", cname, j + 1 );
                exit( EXIT_FAILURE );
            }

            if ( CPXchgcoef( env, lp, lastrow, _mtz_upos( i, problem ), 1.0 ) ) {
                fprintf( stderr, CFATAL "_add_constraints_mtz: CPXchgcoef [%s: u(%zu)]\n", cname, i + 1 );
                exit( EXIT_FAILURE );
            }

            if ( CPXchgcoef( env, lp, lastrow, _mtz_xpos( i, j, problem ), problem->nnodes ) ) {
                fprintf( stderr, CFATAL "_add_constraints_mtz: CPXchgcoef [%s: x(%zu,%zu)]\n",
                    cname, i + 1, j + 1);
                exit( EXIT_FAILURE );
            }
        }
    }

    free( cname );
}


void
mtz_model ( instance *problem )
{

    int error;

    CPXENVptr env = CPXopenCPLEX( &error );
    CPXLPptr lp = CPXcreateprob( env, &error, problem->name ? problem->name : "TSP" );

    /* CPLEX PARAMETERS */
    tspconf_apply( env );

    /* BUILD MODEL */
    _add_constraints_mtz( problem, env, lp );

    struct timeb start, end;
    ftime( &start );

    if ( CPXmipopt( env, lp ) ) {
        fprintf( stderr, CFATAL "mtz_model: CPXmimopt error\n" );
        exit( EXIT_FAILURE );
    }

    ftime( &end );

    double *xopt = malloc( CPXgetnumcols( env, lp ) * sizeof( *xopt ) );
    CPXsolution( env, lp, NULL, NULL, xopt, NULL, NULL, NULL );

    _xopt2solution( xopt, problem, &_mtz_xpos );

    free( xopt );

    problem->elapsedtime  = ( 1000. * ( end.time - start.time ) + end.millitm - start.millitm ) / 1000.;
    problem->visitednodes = CPXgetnodecnt( env, lp );
    problem->solcost      = compute_solution_cost( problem );

    CPXfreeprob( env, &lp );
    CPXcloseCPLEX( &env );
}
