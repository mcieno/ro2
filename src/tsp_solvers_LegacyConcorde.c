/*
 * \brief   Like Legacy, but also cuts the relaxation using Concorde routines.
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
    const CPXCENVptr env;
    void *cbdata;
    int wherefrom;
    int *useraction_p;
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
_LegacyConcorde_xpos ( size_t i, size_t j, const instance *problem )
{
    if ( i == j ) {
        errno = EFAULT;
        perror( CFATAL "_LegacyConcorde_xpos: i == j" );
        exit( EXIT_FAILURE );
    }

    if ( i > j ) return _LegacyConcorde_xpos( j, i, problem );

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
_add_constraints_LegacyConcorde ( const instance *problem, CPXENVptr env, CPXLPptr lp )
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
                fprintf( stderr, CFATAL "_add_constraints_LegacyConcorde: CPXnewcols [%s]\n", cname );
                exit( EXIT_FAILURE );
            }

            if ( CPXgetnumcols( env, lp ) - 1 != _LegacyConcorde_xpos( i, j, problem ) ) {
                fprintf( stderr, CFATAL "_add_constraints_LegacyConcorde: CPXgetnumcols [%s: x(%zu, %zu)]\n",
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
            fprintf( stderr, CFATAL "_add_constraints_LegacyConcorde: CPXnewrows [%s]\n", cname );
            exit( EXIT_FAILURE );
        }

        size_t lastrow = CPXgetnumrows( env, lp ) - 1;

        for ( size_t i = 0; i < problem->nnodes; ++i )
        {
            if ( i == h ) continue;
            if ( CPXchgcoef( env, lp, lastrow, _LegacyConcorde_xpos( i, h, problem ), 1.0 ) ) {
                fprintf( stderr, CFATAL "_add_constraints_LegacyConcorde: CPXchgcoef [%s: x(%zu, %zu)]\n",
                    cname, i + 1, h + 1 );
                exit( EXIT_FAILURE );
            }
        }
    }


    free( cname );
}


void
_add_subtour_constraints_LegacyConcorde ( const instance *problem,
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
        fprintf( stderr, CFATAL "_add_subtour_constraints_LegacyConcorde: out of memory\n" );
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
                rmatind[nzcnt] = _LegacyConcorde_xpos( cnodes[i], cnodes[j], problem );
                rmatval[nzcnt] = 1.0;
                ++nzcnt;
            }
        }

        if ( CPXcutcallbackadd( env, cbdata, wherefrom, nzcnt, rhs, sense, rmatind, rmatval, purgeable ) ) {
            fprintf( stderr, CFATAL "_add_subtour_constraints_LegacyConcorde: CPXcutcallbackadd [SEC(%zu/%zu)]\n",
                k + 1, ncomps );
            exit( EXIT_FAILURE );
        }
    }

    free( memchunk );
}


static int CPXPUBLIC
_lazyconstraintcallback_LegacyConcorde ( CPXCENVptr env,
                                         void       *cbdata,
                                         int        wherefrom,
                                         void       *cbhandle,
                                         int        *useraction_p )
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
        fprintf(stderr, CERROR "_lazyconstraintcallback_LegacyConcorde: Out of memory.\n");
        goto TERMINATE;
    }

    status = CPXgetcallbacknodex( env, cbdata, wherefrom, x, 0, info->ncols - 1 );

    if ( status ) {
        fprintf( stderr, CERROR "_lazyconstraintcallback_LegacyConcorde: CPXgetcallbacknodex.\n" );
        goto TERMINATE;
    }

    _xopt2subtours( info->problem, x, next, comps, &ncomps, _LegacyConcorde_xpos );

    if ( loglevel >= LOG_INFO ) {
        fprintf( stderr, CINFO "_lazyconstraintcallback_LegacyConcorde: got %zu components.\n", ncomps );
    }

    if ( ncomps > 1 ) {
        _add_subtour_constraints_LegacyConcorde( info->problem, env, next, comps, ncomps, cbdata, wherefrom );
        *useraction_p = CPX_CALLBACK_SET;
    }

TERMINATE :

    if (x     != NULL)  free( x     );
    if (next  != NULL)  free( next  );
    if (comps != NULL)  free( comps );

    return status;
}


int
_concorde_callback_LegacyConcorde( double val, int cutcount, int *cut, void *userhandle )
{
    ccinfo_t *ccinfo = (ccinfo_t *) userhandle;

    if ( loglevel >= LOG_DEBUG ) {
        fprintf( stderr, CDEBUG "_concorde_callback_LegacyConcorde: %d nodes in the cut\n", cutcount );
    }

    char sense    = 'G';
    int purgeable = CPX_USECUT_PURGE;
    double rhs    = 2.0;

    int *cutind;
    double *cutval;

    void *memchunk   = malloc( ccinfo->info->ncols * sizeof( *cutind )
                               + ccinfo->info->ncols * sizeof( *cutval ) );

    if ( memchunk == NULL ) {
        fprintf( stderr, CFATAL "_concorde_callback_GenericConcorde: out of memory\n" );
        return 1;
    }

    cutind = memchunk;
    cutval = (double*) (cutind + ccinfo->info->ncols );

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
            cutind[nzcnt] = _LegacyConcorde_xpos( i, j, ccinfo->info->problem );
            cutval[nzcnt] = 1.0;
            ++nzcnt;

        SKIPTHIS:
            continue;
        }
    }


    if ( CPXcutcallbackadd( ccinfo->env, ccinfo->cbdata, ccinfo->wherefrom,
                            nzcnt, rhs, sense, cutind, cutval, purgeable ) )
    {
        fprintf( stderr, CFATAL "_concorde_callback_LegacyConcorde: CPXcutcallbackadd \n");
        exit( EXIT_FAILURE );
    }

    free( memchunk );
    *ccinfo->useraction_p = CPX_CALLBACK_SET;

    return 0;
}


int
_usercutcallback_LegacyConcorde( CPXCENVptr env,
                                 void       *cbdata,
                                 int        wherefrom,
                                 void       *cbhandle,
                                 int        *useraction_p )
{
    *useraction_p = CPX_CALLBACK_DEFAULT;
    int status = 0;

    cbinfo_t *info = (cbinfo_t *) cbhandle;
    ccinfo_t ccinfo = { env, cbdata, wherefrom, useraction_p, info };

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
        fprintf(stderr, CERROR "_relaxationcutcallback_GenericConcorde: out of memory.\n");
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

    status = CPXgetcallbacknodex( env, cbdata, wherefrom, x, 0, info->ncols - 1 );

    if ( status ) {
        fprintf( stderr, CERROR "_usercutcallback_LegacyConcorde: CPXgetcallbacknodex.\n" );
        goto TERMINATE;
    }


    if ( CCcut_connect_components( info->problem->nnodes, nedge, elist, x, &ncomp, &compscount, &comps ) ) {
        fprintf( stderr, CERROR "_usercutcallback_LegacyConcorde: CCcut_connect_components.\n" );
        goto TERMINATE;
    }

    if ( loglevel >= LOG_DEBUG ) {
        fprintf( stderr, CDEBUG "_usercutcallback_LegacyConcorde: relaxation graph is%s connected\n",
            ncomp == 1 ? "" : " NOT" );
    }

    if ( ncomp == 1 ) {
        /* The solution is connected, search for violated cuts */

        if ( CCcut_violated_cuts( info->problem->nnodes, nedge, elist, x, 1.99,
                                  _concorde_callback_LegacyConcorde, &ccinfo ) )
        {
            fprintf( stderr, CERROR "_usercutcallback_LegacyConcorde: CCcut_violated_cuts.\n" );
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
            fprintf( stderr, CFATAL "_usercutcallback_LegacyConcorde: out of memory.\n" );
            exit( EXIT_FAILURE );
        }

        rmatind =           ( _memchunk );
        rmatval = (double*) ( rmatind + info->problem->nnodes * info->problem->nnodes );

        char sense = 'L';
        int purgeable = CPX_USECUT_PURGE;

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
                    rmatind[nzcnt] = _LegacyConcorde_xpos( comps[i], comps[j], info->problem );
                    rmatval[nzcnt] = 1.0;
                    ++nzcnt;
                }
            }

            if ( ( status = CPXcutcallbackadd( env, cbdata, wherefrom, nzcnt, rhs,
                                               sense, rmatind, rmatval, purgeable ) ) )
            {
                fprintf( stderr, CFATAL
                    "_usercutcallback_LegacyConcorde: CPXcutcallbackadd [SEC(%zu/%d)]\n", k + 1, ncomp );
                goto TERMINATE;
            }

            i = compend;
        }

        free( _memchunk );
    }


TERMINATE :

    if ( memchunk != NULL )  free( memchunk );

    return status;
}


void
LegacyConcorde_model ( instance *problem )
{
    int error;

    CPXENVptr env = CPXopenCPLEX( &error );
    CPXLPptr  lp  = CPXcreateprob( env, &error, problem->name ? problem->name : "TSP" );

    /* BUILD MODEL */
    _add_constraints_LegacyConcorde( problem, env, lp );

    cbinfo_t info = { problem, CPXgetnumcols( env, lp ) };
    CPXsetlazyconstraintcallbackfunc( env, _lazyconstraintcallback_LegacyConcorde, &info );
    CPXsetusercutcallbackfunc(env, _usercutcallback_LegacyConcorde, &info);

    /* CPLEX PARAMETERS */
    tspconf_apply( env );
    CPXsetintparam( env, CPXPARAM_MIP_Strategy_CallbackReducedLP, CPX_OFF );


    struct timeb start, end;
    ftime( &start );

    if ( CPXmipopt( env, lp ) ) {
        fprintf( stderr, CFATAL "LegacyConcorde_model: CPXmimopt error\n" );
        exit( EXIT_FAILURE );
    }

    ftime( &end );

    double *xopt  = malloc( CPXgetnumcols( env, lp ) * sizeof( *xopt ) );

    CPXsolution( env, lp, NULL, NULL, xopt, NULL, NULL, NULL );
    _xopt2solution( xopt, problem, &_LegacyConcorde_xpos );

    free( xopt );

    problem->elapsedtime  = ( 1000. * ( end.time - start.time ) + end.millitm - start.millitm ) / 1000.;
    problem->visitednodes = CPXgetnodecnt( env, lp );
    problem->solcost      = compute_solution_cost( problem );

    CPXfreeprob( env, &lp );
    CPXcloseCPLEX( &env );
}
