/*
 * \brief   FLOW1 model with lazy constraints.
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
_LazyFlow1_xpos ( size_t i, size_t j, const instance *problem )
{
    if ( i == j ) {
        errno = EFAULT;
        perror( CFATAL "_Flow1_xpos: i == j" );
        exit( EXIT_FAILURE );
    }

    if ( i > j ) return _LazyFlow1_xpos( j, i, problem );

    return i * problem->nnodes + j - ( ( i + 1 ) * ( i + 2 ) / 2UL );
}


/*!
 * \brief Get the position of variable y(i,j) in Flow1 model. 0-indexed.
 *
 *
 * \param i
 *      i in y(i,j)
 *
 * \param j
 *      j in y(i,j)
 *
 * \param problem
 *     Pointer to the instance structure.
 */
size_t
_LazyFlow1_ypos ( size_t i, size_t j, const instance *problem )
{
    if ( i == j ) {
        errno = EFAULT;
        perror( CFATAL "_Flow1_ypos: i == j" );
        exit( EXIT_FAILURE );
    }

    size_t N = problem->nnodes;
    size_t padding = N * N - ( N * ( N + 1 ) ) / 2;
    return padding + i * ( N - 1 ) + j - ( j > i ? 1 : 0 );
}


/*!
 * \brief Add Flow1 constraints to the problem.
 *
 *
 * This function adds constraints of the Single Commodity Flow model,
 * that is both degree constraints and commodity and coupling constraints.
 *
 * The last two types of constraints are achieved by adding the following *continuous* variables:
 *
 * ```
 * y(i,j) = flow in arc (i,j)  for all i != j
 * ```
 *
 * The new constraince are then
 *
 * ```
 * (0)  Degree constraints
 * (1)  y(i,j) <= (n - 1) x(i,j)  for all i != j
 * (2)  Sum{ y(1,j) }_{ j | j != 1 } = n - 1
 * (3)  Sum{ y(i,h) }_{ i | i != h } - Sum{ y(h,j) }_{ j | j != h } = 1  for all h != 1
 * ```
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
_add_constraints_LazyFlow1( const instance *problem, CPXENVptr env, CPXLPptr lp )
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
                fprintf( stderr, CFATAL "_add_constraints_LazyFlow1: CPXnewcols [%s]\n", cname );
                exit( EXIT_FAILURE );
            }

            if ( CPXgetnumcols( env, lp ) - 1 != _LazyFlow1_xpos( i, j, problem ) ) {
                fprintf( stderr, CFATAL "_add_constraints_LazyFlow1: CPXgetnumcols [%s: x(%zu, %zu)]\n",
                    cname, i + 1, j + 1 );
                exit( EXIT_FAILURE );
            }
        }
    }

    // add constraints (0)
    // Degree constraints
    rhs = 2.0;
    sense = 'E';

    for ( size_t h = 0; h < problem->nnodes; ++h )
    {
        snprintf( cname, CPX_STR_PARAM_MAX, "degree(%zu)", h + 1 );
        if ( CPXnewrows( env, lp, 1, &rhs, &sense, NULL, &cname ) ) {
            fprintf( stderr, CFATAL "_add_constraints_LazyFlow1: CPXnewrows [%s]\n", cname );
            exit( EXIT_FAILURE );
        }

        size_t lastrow = CPXgetnumrows( env, lp ) - 1;

        for ( size_t i = 0; i < problem->nnodes; ++i )
        {
            if ( i == h ) continue;
            if ( CPXchgcoef( env, lp, lastrow, _LazyFlow1_xpos( i, h, problem ), 1.0 ) ) {
                fprintf( stderr, CFATAL "_add_constraints_LazyFlow1: CPXnewrows [%s: x(%zu, %zu)]\n",
                    cname, i + 1, h + 1 );
                exit( EXIT_FAILURE );
            }
        }
    }


    // add continuous vars y(i,j) for all i != j
    ctype = CPX_CONTINUOUS;
    lb = 0.0;
    ub = problem->nnodes;
    obj = 0.0;

    for ( size_t i = 0; i < problem->nnodes; ++i )
    {
        for ( size_t j = 0; j < problem->nnodes; ++j )
        {
            if ( i == j ) continue;

            snprintf( cname, CPX_STR_PARAM_MAX, "y(%zu,%zu)", i + 1, j + 1 );

            if ( CPXnewcols( env, lp, 1, &obj, &lb, &ub, &ctype, &cname ) ) {
                fprintf( stderr, CFATAL "_add_constraints_LazyFlow1: CPXnewcols [%s]\n", cname );
                exit( EXIT_FAILURE );
            }

            if ( CPXgetnumcols( env, lp ) - 1 != _LazyFlow1_ypos( i, j, problem ) ) {
                fprintf( stderr, CFATAL "_add_constraints_LazyFlow1: CPXgetnumcols [%s]\n", cname );
                exit( EXIT_FAILURE );
            }
        }
    }

    // Lazy constraints
    int rmatbeg = 0;
    int *rmatind = malloc( CPXgetnumcols( env, lp ) * sizeof( *rmatind ) );
    double *rmatval = malloc( CPXgetnumcols( env, lp ) * sizeof( *rmatval ) );

    if ( rmatind == NULL || rmatval == NULL ) {
        fprintf( stderr, CFATAL "_add_constraints_LazyFlow1: Not enough memory\n" );
        exit( EXIT_FAILURE );
    }

    // add constraints (1)
    // y(i,j) - (n - 1) x(i,j) <= 0  for all i != j
    rhs = .0;
    sense = 'L';

    for ( size_t i = 0; i < problem->nnodes; ++i )
    {
        for ( size_t j = 0; j < problem->nnodes; ++j )
        {
            if ( i == j ) continue;

            rmatind[0] = (int) _LazyFlow1_ypos( i, j, problem );
            rmatval[0] = +1.0;

            rmatind[1] = (int) _LazyFlow1_xpos( i, j, problem );
            rmatval[1] = +1.0 - (double) problem->nnodes;

            snprintf( cname, CPX_STR_PARAM_MAX, "LazyFlow1_1(%zu,%zu)", i + 1, j + 1 );
            if ( CPXaddlazyconstraints( env, lp, 1, 2, &rhs, &sense, &rmatbeg, rmatind, rmatval, &cname ) ) {
                fprintf( stderr, CFATAL "_add_constraints_LazyFlow1: CPXaddlazyconstraints [%s]\n", cname );
                exit( EXIT_FAILURE );
            }
        }
    }

    // add constraint (2)
    // Sum{ y(1,j) }_{ j | j != 1 } = n - 1
    rhs = (double) ( problem->nnodes - 1 );
    sense = 'E';
    snprintf( cname, CPX_STR_PARAM_MAX, "LazyFlow1_2" );

    for ( size_t j = 1; j < problem->nnodes; ++j )
    {
        rmatind[j] = (int) _LazyFlow1_ypos( 0, j, problem );
        rmatval[j] = +1.0;
    }

    if ( CPXaddlazyconstraints( env, lp, 1, problem->nnodes - 1, &rhs, &sense, &rmatbeg, rmatind, rmatval, &cname ) ) {
        fprintf( stderr, CFATAL "_add_constraints_LazyFlow1: CPXaddlazyconstraints [%s]\n", cname );
        exit( EXIT_FAILURE );
    }

    // add constraints (3)
    // Sum{ y(i,h) }_{ i | i != h } - Sum{ y(h,j) }_{ j | j != h } = 1  for all h != 1
    rhs = 1.0;
    sense = 'E';

    for ( size_t h = 1; h < problem->nnodes; ++h )
    {
        snprintf( cname, CPX_STR_PARAM_MAX, "LazyFlow1_3(%zu)", h + 1 );

        size_t ii = 0;
        for ( size_t i = 0; i < problem->nnodes; ++i )
        {
            if ( i == h ) continue;

            rmatind[2 * ii] = (int) _LazyFlow1_ypos( i, h, problem );
            rmatval[2 * ii] = +1.0;
            rmatind[2 * ii + 1] = (int) _LazyFlow1_ypos( h, i, problem );
            rmatval[2 * ii + 1] = -1.0;

            ++ii;
        }

        if ( CPXaddlazyconstraints( env, lp, 1, 2 * problem->nnodes - 2, &rhs, &sense, &rmatbeg, rmatind, rmatval, &cname ) ) {
            fprintf( stderr, CFATAL "_add_constraints_LazyFlow1: CPXaddlazyconstraints [%s]\n", cname );
            exit( EXIT_FAILURE );
        }
    }

    free( rmatind );
    free( rmatval );
    free( cname );
}


void
LazyFlow1_model ( instance *problem )
{
    int error;

    CPXENVptr env = CPXopenCPLEX( &error );
    CPXLPptr lp = CPXcreateprob( env, &error, problem->name ? problem->name : "TSP" );

    /* CPLEX PARAMETERS */
    tspconf_apply( env );

    /* BUILD MODEL */
    _add_constraints_LazyFlow1( problem, env, lp );

    struct timeb start, end;
    ftime( &start );

    if ( CPXmipopt( env, lp ) ) {
        fprintf( stderr, CFATAL "Flow1_model: CPXmimopt error\n" );
        exit( EXIT_FAILURE );
    }

    ftime( &end );

    double *xopt = malloc( CPXgetnumcols( env, lp ) * sizeof( *xopt ) );
    CPXsolution( env, lp, NULL, NULL, xopt, NULL, NULL, NULL );

    _xopt2solution( xopt, problem, &_LazyFlow1_xpos );

    free( xopt );

    problem->elapsedtime  = ( 1000. * ( end.time - start.time ) + end.millitm - start.millitm ) / 1000.;
    problem->visitednodes = CPXgetnodecnt( env, lp );
    problem->solcost      = compute_solution_cost( problem );

    CPXfreeprob(env, &lp);
    CPXcloseCPLEX(&env);
}
