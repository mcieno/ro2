/*
 * \brief   MTZ model with lazy constraints.
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
#include "tspconf.h"


/*!
 * \brief Get the position of variable x(i,j) in CPLEX internal state.
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
_LazyMTZ_xpos ( size_t i, size_t j, const instance *problem )
{
    if ( i == j ) {
        log_fatal( "i == j" );
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
_LazyMTZ_upos ( size_t i, const instance *problem )
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
_add_constraints_LazyMTZ ( const instance *problem, CPXENVptr env, CPXLPptr lp )
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
                log_fatal( "CPXnewcols [%s]", cname );
                exit( EXIT_FAILURE );
            }

            if ( CPXgetnumcols( env, lp ) - 1 != _LazyMTZ_xpos( i, j, problem ) ) {
                log_fatal( "CPXgetnumcols [%s: x(%zu, %zu)]",
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
            log_fatal( "CPXnewrows [%s]", cname );
            exit( EXIT_FAILURE );
        }

        lastrow = CPXgetnumrows( env, lp ) - 1;

        for ( size_t i = 0; i < problem->nnodes; ++i )
        {
            if ( i == h ) continue;

            if ( CPXchgcoef( env, lp, lastrow, _LazyMTZ_xpos( i, h, problem ), 1.0 ) ) {
                log_fatal( "CPXchgcoef [%s: x(%zu, %zu)]",
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
            log_fatal( "CPXnewrows [%s]", cname );
            exit( EXIT_FAILURE );
        }

        lastrow = CPXgetnumrows( env, lp ) - 1;

        for ( size_t i = 0; i < problem->nnodes; ++i )
        {
            if ( i == h ) continue;

            if ( CPXchgcoef( env, lp, lastrow, _LazyMTZ_xpos( h, i, problem ), 1.0 ) ) {
                log_fatal( "CPXchgcoef [%s: x(%zu, %zu)]",
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
            log_fatal( "CPXnewcols [%s]", cname );
            exit( EXIT_FAILURE );
        }

        if ( CPXgetnumcols( env, lp ) - 1 != _LazyMTZ_upos( i, problem ) ) {
            log_fatal( "CPXgetnumcols [%s]", cname );
            exit( EXIT_FAILURE );
        }
    }


    // Lazy constraints
    int rmatbeg = 0;
    int *rmatind = malloc( CPXgetnumcols( env, lp ) * sizeof( *rmatind ) );
    double *rmatval = malloc( CPXgetnumcols( env, lp ) * sizeof( *rmatval ) );

    if ( rmatind == NULL || rmatval == NULL ) {
        log_fatal( "Out of memory." );
        exit( EXIT_FAILURE );
    }

    // u(i) - u(j) + n x(i,j) <= n - 1  for all i, j | i != j, i != 0, j != 0
    rhs = problem->nnodes - 1;
    sense = 'L';

    for ( size_t j = 1; j < problem->nnodes; ++j )
    {
        for ( size_t i = 1; i < problem->nnodes; ++i )
        {
            if ( i == j ) continue;

            snprintf( cname, CPX_STR_PARAM_MAX, "sequence(%zu)", j + 1 );

            rmatind[0] = (int) _LazyMTZ_upos( i, problem );
            rmatval[0] = +1.0;

            rmatind[1] = (int) _LazyMTZ_upos( j, problem );
            rmatval[1] = -1.0;

            rmatind[2] = (int) _LazyMTZ_xpos( i, j, problem );
            rmatval[2] = (double) problem->nnodes;

            if ( CPXaddlazyconstraints( env, lp, 1, 3, &rhs, &sense, &rmatbeg, rmatind, rmatval, &cname ) ) {
                log_fatal( "CPXaddlazyconstraints [%s]", cname );
                exit( EXIT_FAILURE );
            }
        }
    }

    free( rmatind );
    free( rmatval );
    free( cname );
}


void
LazyMTZ_model ( instance *problem )
{

    int error;

    CPXENVptr env = CPXopenCPLEX( &error );
    CPXLPptr lp = CPXcreateprob( env, &error, problem->name ? problem->name : "TSP" );

    /* BUILD MODEL */
    log_info( "Adding constraints to the model." );
    _add_constraints_LazyMTZ( problem, env, lp );

    /* CPLEX PARAMETERS */
    tspconf_apply( env );

    struct timespec start, end;
    clock_gettime( CLOCK_MONOTONIC, &start );

    log_info( "Starting solver." );
    if ( CPXmipopt( env, lp ) ) {
        log_fatal( "CPXmipopt error." );
        exit( EXIT_FAILURE );
    }

    clock_gettime( CLOCK_MONOTONIC, &end );

    log_info( "Retrieving final solution." );
    double *xopt = malloc( CPXgetnumcols( env, lp ) * sizeof( *xopt ) );

    if ( xopt == NULL ) {
        log_fatal( "Out of memory." );
        exit( EXIT_FAILURE );
    }

    CPXsolution( env, lp, NULL, NULL, xopt, NULL, NULL, NULL );

    _xopt2solution( xopt, problem, &_LazyMTZ_xpos );

    free( xopt );

    problem->elapsedtime  = ( end.tv_sec - start.tv_sec ) + ( end.tv_nsec - start.tv_nsec ) / 1000000000.;
    problem->visitednodes = CPXgetnodecnt( env, lp );
    problem->solcost      = compute_solution_cost( problem );

    CPXfreeprob( env, &lp );
    CPXcloseCPLEX( &env );
}
