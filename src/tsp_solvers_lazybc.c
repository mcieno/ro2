/*
 * \brief   Branch and Cut model with lazy constraint callback.
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


typedef struct
{
    instance *problem;
    int ncols;
}
cbinfo_t;


/*!
 * \brief Get the position of variable x(i,j) in B&B model with lazy constraint callback.
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
_lazyBC_xpos ( size_t i, size_t j, const instance *problem )
{
    if ( i == j ) {
        errno = EFAULT;
        perror( CFATAL "_lazyBC_xpos: i == j" );
        exit( EXIT_FAILURE );
    }

    if ( i > j ) return _lazyBC_xpos( j, i, problem );

    return i * problem->nnodes + j - ( ( i + 1 ) * ( i + 2 ) / 2UL );
}


/*!
 * \brief Add degree constraints to the problem.
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
_add_constraints_lazyBC ( const instance *problem, CPXENVptr env, CPXLPptr lp )
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
                fprintf( stderr, CFATAL "_add_constraints_lazyBC: CPXnewcols [%s]\n", cname );
                exit( EXIT_FAILURE );
            }

            if ( CPXgetnumcols( env, lp ) - 1 != _lazyBC_xpos( i, j, problem ) ) {
                fprintf( stderr, CFATAL "_add_constraints_lazyBC: CPXgetnumcols [%s: x(%zu, %zu)]\n",
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
            fprintf( stderr, CFATAL "_add_constraints_lazyBC: CPXnewrows [%s]\n", cname );
            exit( EXIT_FAILURE );
        }

        size_t lastrow = CPXgetnumrows( env, lp ) - 1;

        for ( size_t i = 0; i < problem->nnodes; ++i )
        {
            if ( i == h ) continue;
            if ( CPXchgcoef( env, lp, lastrow, _lazyBC_xpos( i, h, problem ), 1.0 ) ) {
                fprintf( stderr, CFATAL "_add_constraints_lazyBC: CPXchgcoef [%s: x(%zu, %zu)]\n",
                    cname, i + 1, h + 1 );
                exit( EXIT_FAILURE );
            }
        }
    }


    free( cname );
}


void
_add_subtour_constraints_lazyBC ( const instance *problem,
                                  CPXCENVptr     env,
                                  size_t         *next,
                                  size_t         *comps,
                                  size_t         ncomps,
                                  void           *cbdata,
                                  int            wherefrom )
{
    if ( ncomps == 1 ) {
        return;
    }

    int *cnodes;
    int *rmatind;
    double *rmatval;
    void *memchunk = malloc(                     problem->nnodes * sizeof( *cnodes  )
                             + problem->nnodes * problem->nnodes * sizeof( *rmatind )
                             + problem->nnodes * problem->nnodes * sizeof( *rmatval ) );

    if ( memchunk == NULL ) {
        fprintf( stderr, CFATAL "_add_subtour_constraints_lazyBC: out of memory\n" );
        exit( EXIT_FAILURE );
    }

    cnodes  =           ( memchunk );
    rmatind =           ( cnodes + problem->nnodes );
    rmatval = (double*) ( rmatind + problem->nnodes * problem->nnodes );

    char sense = 'L';
    int purgeable = CPX_USECUT_PURGE;

    /* Add constraint for k-th component */
    double rhs;
    for ( size_t k = 0; k < ncomps; ++k )
    {
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
                rmatind[nzcnt] = _lazyBC_xpos( cnodes[i], cnodes[j], problem );
                rmatval[nzcnt] = 1.0;
                ++nzcnt;
            }
        }

        if ( CPXcutcallbackadd( env, cbdata, wherefrom, nzcnt, rhs, sense, rmatind, rmatval, purgeable ) ) {
            fprintf( stderr, CFATAL "_add_subtour_constraints_lazyBC: CPXcutcallbackadd [SEC(%zu/%zu)]\n",
                k + 1, ncomps );
            exit( EXIT_FAILURE );
        }
    }

    free( memchunk );
}


static int CPXPUBLIC
_lazyconstraintcallback_lazyBC ( CPXCENVptr env, void *cbdata, int wherefrom, void *cbhandle, int *useraction_p )
{
    int status = 0;

    *useraction_p = CPX_CALLBACK_DEFAULT;
    cbinfo_t *info = (cbinfo_t *) cbhandle;

    size_t ncomps = 0;
    double *x     = malloc( info->ncols * sizeof( *x ) );
    size_t *next  = calloc( info->problem->nnodes, sizeof( *next ) );
    size_t *comps = calloc( info->problem->nnodes, sizeof( *comps ) );

    if ( x     == NULL ||
         next  == NULL ||
         comps == NULL  ) {
        fprintf(stderr, CERROR "_lazyconstraintcallback_lazyBC: Out of memory.\n");
        goto TERMINATE;
    }

    status = CPXgetcallbacknodex( env, cbdata, wherefrom, x, 0, info->ncols - 1 );

    if ( status ) {
        fprintf( stderr, CERROR "_lazyconstraintcallback_lazyBC: CPXgetcallbacknodex.\n" );
        goto TERMINATE;
    }

    _xopt2subtours( info->problem, x, next, comps, &ncomps, _lazyBC_xpos );

    if ( loglevel >= LOG_INFO ) {
        fprintf( stderr, CINFO "_lazyconstraintcallback_lazyBC: got %zu components.\n", ncomps );
    }

    if ( ncomps > 1 ) {
        _add_subtour_constraints_lazyBC( info->problem, env, next, comps, ncomps, cbdata, wherefrom );
        *useraction_p = CPX_CALLBACK_SET;
    }

TERMINATE :

    if (x     != NULL)  free( x     );
    if (next  != NULL)  free( next  );
    if (comps != NULL)  free( comps );

    return status;
}


void
lazyBC_model ( instance *problem )
{
    int error;

    CPXENVptr env = CPXopenCPLEX( &error );
    CPXLPptr  lp  = CPXcreateprob( env, &error, problem->name ? problem->name : "TSP" );

    /* BUILD MODEL */
    _add_constraints_lazyBC( problem, env, lp );

    cbinfo_t info = { problem, CPXgetnumcols( env, lp ) };
    CPXsetlazyconstraintcallbackfunc( env, _lazyconstraintcallback_lazyBC, &info );

    /* CPLEX PARAMETERS */
    tspconf_apply( env );
    CPXsetintparam( env, CPXPARAM_MIP_Strategy_CallbackReducedLP, CPX_OFF );


    struct timeb start, end;
    ftime( &start );

    if ( CPXmipopt( env, lp ) ) {
        fprintf( stderr, CFATAL "lazyBC_model: CPXmimopt error\n" );
        exit( EXIT_FAILURE );
    }

    ftime( &end );

    double *xopt  = malloc( CPXgetnumcols( env, lp ) * sizeof( *xopt ) );

    CPXsolution( env, lp, NULL, NULL, xopt, NULL, NULL, NULL );
    _xopt2solution( xopt, problem, &_lazyBC_xpos );

    free( xopt );

    problem->elapsedtime  = ( 1000. * ( end.time - start.time ) + end.millitm - start.millitm ) / 1000.;
    problem->visitednodes = CPXgetnodecnt( env, lp );
    problem->solcost      = compute_solution_cost( problem );

    CPXfreeprob( env, &lp );
    CPXcloseCPLEX( &env );
}
