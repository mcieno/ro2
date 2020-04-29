/*
 * \brief   Like GenericConcorde but only cuts nodes close to the root.
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
#include "concorde.h"


typedef struct
{
    instance *problem;
    int ncols;
}
cbinfo_t;

typedef struct
{
    const CPXCALLBACKCONTEXTptr context;
    cbinfo_t *info;
}
ccinfo_t;


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
_GenericConcordeRand_xpos ( size_t i, size_t j, const instance *problem )
{
    if ( i == j ) {
        log_fatal( "i == j" );
        exit( EXIT_FAILURE );
    }

    if ( i > j ) return _GenericConcordeRand_xpos( j, i, problem );

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
_add_constraints_GenericConcordeRand ( const instance *problem, CPXENVptr env, CPXLPptr lp )
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

            if ( CPXgetnumcols( env, lp ) - 1 != _GenericConcordeRand_xpos( i, j, problem ) ) {
                log_fatal( "CPXgetnumcols [%s: x(%zu, %zu)]",
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
            if ( CPXchgcoef( env, lp, lastrow, _GenericConcordeRand_xpos( i, h, problem ), 1.0 ) ) {
                log_fatal( "CPXchgcoef [%s: x(%zu, %zu)]",
                    cname, i + 1, h + 1 );
                exit( EXIT_FAILURE );
            }
        }
    }


    free( cname );
}


void
_add_subtour_constraints_GenericConcordeRand ( const instance       *problem,
                                               CPXCALLBACKCONTEXTptr context,
                                               size_t                *next,
                                               size_t                *comps,
                                               size_t                ncomps )
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
        log_fatal( "Out of memory." );
        CPXcallbackabort( context );
    }

    cnodes  =           ( memchunk );
    rmatind =           ( cnodes + problem->nnodes );
    rmatval = (double*) ( rmatind + problem->nnodes * problem->nnodes );

    char sense = 'L';
    int rmatbeg = 0;

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
                rmatind[nzcnt] = _GenericConcordeRand_xpos( cnodes[i], cnodes[j], problem );
                rmatval[nzcnt] = 1.0;
                ++nzcnt;
            }
        }

        if ( CPXcallbackrejectcandidate( context, 1, nzcnt, &rhs, &sense, &rmatbeg, rmatind, rmatval ) ) {
            log_fatal( "CPXcallbackaddusercuts [SEC(%zu/%zu)]",
                k + 1, ncomps );
            CPXcallbackabort( context );
        }
    }

    free( memchunk );
}


static int CPXPUBLIC
_candidatecutcallback_GenericConcordeRand ( CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void *userhandle )
{
    int status = 0;

    cbinfo_t *info = (cbinfo_t *) userhandle;

    size_t ncomps = 0;
    double *x     = malloc( info->ncols * sizeof( *x ) );
    size_t *next  = calloc( info->problem->nnodes, sizeof( *next ) );
    size_t *comps = calloc( info->problem->nnodes, sizeof( *comps ) );

    if ( x     == NULL ||
         next  == NULL ||
         comps == NULL  ) {
        log_fatal( "Out of memory.");
        goto TERMINATE;
    }

    int ispoint;
    status = CPXcallbackcandidateispoint( context, &ispoint );

    if ( status || !ispoint ) {
        /* Not a feasible solution */
        goto TERMINATE;
    }

    status = CPXcallbackgetcandidatepoint(context, x, 0, info->ncols - 1, NULL);

    if ( status ) {
        log_fatal( "CPXcallbackgetcandidatepoint returned %d.", status );
        goto TERMINATE;
    }

    _xopt2subtours( info->problem, x, next, comps, &ncomps, _GenericConcordeRand_xpos );

    log_debug( "Found %zu components.", ncomps );

    if ( ncomps > 1 ) {
        _add_subtour_constraints_GenericConcordeRand( info->problem, context, next, comps, ncomps );
    }

TERMINATE :

    if (x     != NULL)  free( x     );
    if (next  != NULL)  free( next  );
    if (comps != NULL)  free( comps );

    return status;
}


int
_concorde_callback_GenericConcordeRand ( double val, int cutcount, int *cut, void *userhandle )
{
    ccinfo_t *ccinfo = (ccinfo_t *) userhandle;

    log_debug( "Cut contains %d nodes.", cutcount );

    char sense      = 'G';
    double rhs      = 2;
    int purgeable   = CPX_USECUT_PURGE;
    int local       = 0;
    int rmatbeg     = 0;

    int *rmatind;
    double *rmatval;

    void *memchunk   = malloc( ccinfo->info->ncols * sizeof( *rmatind )
                               + ccinfo->info->ncols * sizeof( *rmatval ) );

    if ( memchunk == NULL ) {
        log_fatal( "Out of memory." );
        return 1;
    }

    rmatind = memchunk;
    rmatval = (double*) (rmatind + ccinfo->info->ncols );

    int nzcnt = 0;

    int i;
    for ( size_t l = 0; l < cutcount; ++l )
    {
        /* Node i is in S */
        i = cut[l];

        /* Find all nodes in V \ S and add the constraint */
        for ( size_t j = 0; j < ccinfo->info->problem->nnodes; ++j ) {
            if ( j == i ) {
                goto SKIPTHIS;
            }

            for ( size_t k = 0; k < cutcount; ++k ) {
                if ( cut[k] == j ) {
                    goto SKIPTHIS;
                }
            }

            /* Node j is in V \ S */
            rmatind[nzcnt] = _GenericConcordeRand_xpos(i, j, ccinfo->info->problem);
            rmatval[nzcnt] = 1.0;
            ++nzcnt;

        SKIPTHIS:
            continue;
        }
    }


    if ( CPXcallbackaddusercuts( ccinfo->context, 1, nzcnt, &rhs, &sense,
                                 &rmatbeg, rmatind, rmatval, &purgeable, &local ) )
    {
        log_fatal( "CPXcutcallbackadd");
        exit( EXIT_FAILURE );
    }

    free( memchunk );

    return 0;
}


static int CPXPUBLIC
_relaxationcutcallback_GenericConcordeRand ( CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void *userhandle )
{
    static unsigned int seed;

    CPXLONG nodedepth;
    int status = CPXcallbackgetinfolong( context, CPXCALLBACKINFO_NODEDEPTH, &nodedepth );

    if ( status ) {
        log_fatal( "CPXcallbackgetinfolong returned %d.", status);
        return status;
    }

    /* Cut node at depth i with probability 2^-i. Root is cut with probability 1. */
    int rnd = rand_r(&seed);
    if ( rnd > ( INT_MAX >> nodedepth ) ) return status;


    cbinfo_t *info = (cbinfo_t *) userhandle;
    ccinfo_t ccinfo = { context, info };

    int ncomp;
    int nedge = ( info->problem->nnodes * ( info->problem->nnodes - 1 ) ) / 2;

    int *elist;
    int *comps;
    int *compscount;
    double *x;

    void *memchunk = malloc(   nedge * 2             * sizeof( *elist      )
                             + info->problem->nnodes * sizeof( *comps      )
                             + info->problem->nnodes * sizeof( *compscount )
                             + info->ncols           * sizeof( *x          ) );

    if ( memchunk == NULL ) {
        log_fatal( "Out of memory.");
        goto TERMINATE;
    }

    elist      =           ( memchunk );
    comps      =           ( elist + nedge * 2 );
    compscount =           ( comps + info->problem->nnodes );
    x          = (double*) ( compscount + info->problem->nnodes );

    int loader = 0;

    for ( int i = 0; i < info->problem->nnodes; ++i ) {
        for ( int j = i + 1; j < info->problem->nnodes; ++j ) {
            elist[ loader++ ] = i;
            elist[ loader++ ] = j;
        }
    }

    status = CPXcallbackgetrelaxationpoint( context, x, 0, info->ncols - 1, NULL );

    if ( status ) {
        log_fatal( "CPXgetcallbacknodexreturned %d.", status );
        goto TERMINATE;
    }

    if ( CCcut_connect_components( info->problem->nnodes, nedge, elist, x, &ncomp, &compscount, &comps ) ) {
        log_fatal( "CCcut_connect_components" );
        status = 1;
        goto TERMINATE;
    }

    log_debug( "Relaxation graph is%s connected", ncomp == 1 ? "" : " not" );

    if ( ncomp == 1 ) {
        /* The solution is connected, search for violated cuts */

        if ( CCcut_violated_cuts( info->problem->nnodes, nedge, elist, x, 1.99,
                                  _concorde_callback_GenericConcordeRand, &ccinfo) )
        {
            log_fatal( "CCcut_violated_cuts" );
            status = 1;
            goto TERMINATE;
        }

    } else {
        /* The solution has subtours and we can add the corresponding SEC's */

        int *rmatind;
        double *rmatval;
        void *_memchunk = malloc(   info->problem->nnodes * info->problem->nnodes * sizeof( *rmatind )
                                 + info->problem->nnodes * info->problem->nnodes * sizeof( *rmatval ) );

        if ( _memchunk == NULL ) {
            log_fatal( "Out of memory." );
            CPXcallbackabort( context );
        }

        rmatind =           ( _memchunk );
        rmatval = (double*) ( rmatind + info->problem->nnodes * info->problem->nnodes );

        char sense = 'L';
        int purgeable = CPX_USECUT_PURGE;
        int rmatbeg = 0;
        int local = 0;

        int i = 0;
        int compend = i;
        double rhs;
        int nzcnt;

        for ( size_t k = 0; k < ncomp; ++k ) {
            rhs = compscount[k] - 1.0;
            compend += compscount[k];

            nzcnt = 0;

            for ( ; i < compend; ++i ) {
                for ( int j = i + 1; j < compend; ++j ) {
                    rmatind[nzcnt] = _GenericConcordeRand_xpos( comps[i], comps[j], info->problem );
                    rmatval[nzcnt] = 1.0;
                    ++nzcnt;
                }
            }

            if ( CPXcallbackaddusercuts( context, 1, nzcnt, &rhs, &sense,
                                         &rmatbeg, rmatind, rmatval, &purgeable, &local ) )
            {
                log_fatal( "CPXcutcallbackadd[SEC(%zu/%d)]", k + 1, ncomp );
                CPXcallbackabort( context );
            }

            i = compend;
        }

        free( _memchunk );
    }


TERMINATE:

    if ( memchunk != NULL )  free( memchunk );

    return status;
}


static int CPXPUBLIC
_callbackfunc_GenericConcordeRand ( CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void *userhandle )
{
    int status = 1;

    switch ( contextid )
    {
        case CPX_CALLBACKCONTEXT_RELAXATION:
            status = _relaxationcutcallback_GenericConcordeRand( context, contextid, userhandle );
            break;

        case CPX_CALLBACKCONTEXT_CANDIDATE:
            status = _candidatecutcallback_GenericConcordeRand( context, contextid, userhandle );
            break;

        default:
            CPXcallbackabort( context );
    }

    return status;
}


void
GenericConcordeRand_model ( instance *problem )
{
    int error;

    CPXENVptr env = CPXopenCPLEX( &error );
    CPXLPptr  lp  = CPXcreateprob( env, &error, problem->name ? problem->name : "TSP" );

    /* BUILD MODEL */
    log_info( "Adding constraints to the model." );
    _add_constraints_GenericConcordeRand( problem, env, lp );

    log_info( "Setting up callbacks." );
    cbinfo_t info = { problem, CPXgetnumcols( env, lp ) };
    CPXcallbacksetfunc( env, lp, CPX_CALLBACKCONTEXT_RELAXATION | CPX_CALLBACKCONTEXT_CANDIDATE,
                        _callbackfunc_GenericConcordeRand, &info );

    /* CPLEX PARAMETERS */
    tspconf_apply( env );
    CPXsetintparam( env, CPXPARAM_MIP_Strategy_CallbackReducedLP, CPX_OFF );


    struct timeb start, end;
    ftime( &start );

    log_info( "Starting solver." );
    if ( CPXmipopt( env, lp ) ) {
        log_fatal( "CPXmipopt error." );
        exit( EXIT_FAILURE );
    }

    ftime( &end );

    log_info( "Retrieving final solution." );
    double *xopt = malloc( CPXgetnumcols( env, lp ) * sizeof( *xopt ) );

    if ( xopt == NULL ) {
        log_fatal( "Out of memory." );
        exit( EXIT_FAILURE );
    }

    CPXsolution( env, lp, NULL, NULL, xopt, NULL, NULL, NULL );
    _xopt2solution( xopt, problem, &_GenericConcordeRand_xpos );

    free( xopt );

    problem->elapsedtime  = ( 1000. * ( end.time - start.time ) + end.millitm - start.millitm ) / 1000.;
    problem->visitednodes = CPXgetnodecnt( env, lp );
    problem->solcost      = compute_solution_cost( problem );

    CPXfreeprob( env, &lp );
    CPXcloseCPLEX( &env );
}
