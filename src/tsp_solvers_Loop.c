/*
 * \brief   Branch and Cut model with SEC added at every restart.
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
#include "concorde.h"


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
_Loop_xpos ( size_t i, size_t j, const instance *problem )
{
    if ( i == j ) {
        log_fatal( "i == j" );
        exit( EXIT_FAILURE );
    }

    if ( i > j ) return _Loop_xpos( j, i, problem );

    return i * problem->nnodes + j - ( ( i + 1 ) * ( i + 2 ) / 2UL );
}


/*!
 * \brief Add Degree constraints to the problem.
 *
 *
 * \param problem
 *     Pointer to the instance structure.
 *
 * \param env
 *     CPLEX environment
 *
 * \param lp
 *     CPLEX problem.
 */
void
_add_constraints_Loop ( const instance *problem, CPXENVptr env, CPXLPptr lp )
{
    char ctype;
    double lb, ub, obj, rhs;
    char sense;

    char *cname = calloc( CPX_STR_PARAM_MAX, sizeof( *cname ) );

    // add binary vars x(i,j) for 0 <= i < j <= N
    ctype = CPX_BINARY;
    lb = 0.0;
    ub = 1.0;

    for ( size_t i = 0; i < problem->nnodes; ++i )
    {
        for ( size_t j = i + 1; j < problem->nnodes; ++j )
        {
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

            if ( CPXgetnumcols( env, lp ) - 1 != _Loop_xpos( i, j, problem ) ) {
                log_fatal( "x(%zu, %zu)]",
                    cname, i + 1, j + 1 );
                exit( EXIT_FAILURE );
            }
        }
    }

    // Degree constraints
    rhs = 2.0;
    sense = 'E';

    for ( size_t h = 0; h < problem->nnodes; ++h )
    {
        snprintf( cname, CPX_STR_PARAM_MAX, "degree(%zu)", h + 1 );
        if ( CPXnewrows( env, lp, 1, &rhs, &sense, NULL, &cname ) ) {
            log_fatal( "CPXnewrows [%s]", cname );
            exit( EXIT_FAILURE );
        }

        size_t lastrow = CPXgetnumrows( env, lp ) - 1;

        for ( size_t i = 0; i < problem->nnodes; ++i )
        {
            if ( i == h ) continue;
            if ( CPXchgcoef( env, lp, lastrow, _Loop_xpos( i, h, problem ), 1.0 ) ) {
                log_fatal( "x(%zu, %zu)]",
                    cname, i + 1, h + 1 );
                exit( EXIT_FAILURE );
            }
        }
    }


    free( cname );
}


void
_add_subtour_constraints ( const instance *problem,
                           CPXENVptr      env,
                           CPXLPptr       lp,
                           size_t         *next,
                           size_t         *comps,
                           size_t         ncomps )
{
    if ( ncomps == 1 ) {
        return;
    }

    int *cnodes;
    int *rmatind;
    double *rmatval;
    char *cname;
    void *memchunk = malloc(                     problem->nnodes * sizeof( *cnodes  )
                             + problem->nnodes * problem->nnodes * sizeof( *rmatind )
                             + problem->nnodes * problem->nnodes * sizeof( *rmatval )
                             +                 CPX_STR_PARAM_MAX * sizeof( *cname   ) );

    if ( memchunk == NULL ) {
        log_fatal( "Out of memory." );
        exit( EXIT_FAILURE );
    }

    cnodes  =           ( memchunk );
    rmatind =           ( cnodes + problem->nnodes );
    rmatval = (double*) ( rmatind + problem->nnodes * problem->nnodes );
    cname   = (char*)   ( rmatval + problem->nnodes * problem->nnodes );

    char sense = 'L';
    int rmatbeg[] = {0};

    /* Add constraint for k-th component */
    double rhs;
    for ( size_t k = 0; k < ncomps; ++k )
    {
        snprintf( cname, CPX_STR_PARAM_MAX, "SEC(%zu/%zu)", k + 1, ncomps);

        /* Find k-th component initial node */
        size_t j_start = 0;
        for ( j_start = 0; j_start < problem->nnodes; ++j_start ) {
            if ( comps[j_start] == k + 1 ) {
                break;
            }
        }

        int compsize = 1;
        cnodes[0] = j_start;
        for (size_t j = next[j_start]; j != j_start; cnodes[compsize++] = j, j = next[j])
            ;

        rhs = compsize - 1.0;

        /* Build rmatind/rmatval */
        int nzcnt = 0;
        for (size_t i = 0; i < compsize; ++i) {
            for (size_t j = i + 1; j < compsize; ++j) {
                rmatind[nzcnt] = _Loop_xpos( cnodes[i], cnodes[j], problem );
                rmatval[nzcnt] = 1.0;
                ++nzcnt;
            }
        }

        if ( CPXaddrows( env, lp, 0, 1, nzcnt, &rhs, &sense,
                         rmatbeg, rmatind, rmatval, NULL, &cname ) ) {
            log_fatal( "CPXaddrows [SEC(%zu/%zu)]", k + 1, ncomps);
            exit( EXIT_FAILURE );
        }
    }

    free( memchunk );
}


void
Loop_model ( instance *problem )
{
    int error;



    CPXENVptr env = CPXopenCPLEX( &error );
    CPXLPptr lp = CPXcreateprob( env, &error, problem->name ? problem->name : "TSP" );

    /* BUILD MODEL */
    log_info( "Adding constraints to the model." );
    _add_constraints_Loop(problem, env, lp);

    /* CPLEX PARAMETERS */
    tspconf_apply( env );

    size_t ncomps = 0;
    double *xopt  = malloc( CPXgetnumcols( env, lp ) * sizeof( *xopt ) );
    size_t *next =  calloc( problem->nnodes, sizeof( *next ) );
    size_t *comps = calloc( problem->nnodes, sizeof( *comps ) );

    if ( xopt == NULL || next == NULL || comps == NULL ) {
        log_fatal( "Out of memory." );
        exit( EXIT_FAILURE );
    }

    int visitednodes = 0;
    struct timespec start, end;
    clock_gettime( CLOCK_MONOTONIC, &start );

    log_info( "Starting solver main loop." );

    for ( size_t iter = 0; ncomps != 1; ++iter )
    {
        if ( CPXmipopt( env, lp ) ) {
            log_fatal( "CPXmipopt error." );
            exit( EXIT_FAILURE );
        }

        clock_gettime( CLOCK_MONOTONIC, &end );

        visitednodes += CPXgetnodecnt( env, lp ) + 1;
        CPXsolution( env, lp, NULL, NULL, xopt, NULL, NULL, NULL );
        _xopt2subtours( problem, xopt, next, comps, &ncomps, _Loop_xpos );

        log_debug( "Iteration %zu",                                     iter );
        log_debug( "    - Components: %zu",                           ncomps );
        log_debug( "    - Elapsed:    %lfs",
            ( end.tv_sec - start.tv_sec ) + ( end.tv_nsec - start.tv_nsec ) / 1000000000. );


        _add_subtour_constraints( problem, env, lp, next, comps, ncomps );
    }

    clock_gettime( CLOCK_MONOTONIC, &end );

    log_info( "Retrieving final solution." );
    _xopt2solution( xopt, problem, &_Loop_xpos );

    free( xopt );

    problem->elapsedtime  = ( end.tv_sec - start.tv_sec ) + ( end.tv_nsec - start.tv_nsec ) / 1000000000.;
    problem->visitednodes = visitednodes;
    problem->solcost      = compute_solution_cost( problem );

    CPXfreeprob( env, &lp );
    CPXcloseCPLEX( &env );
}
