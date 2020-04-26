/*
 * \brief   Hard-Fix heuristic method
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
 * \brief Get the position of variable x(i,j) in B&B model.
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
_hardfix_xpos ( size_t i, size_t j, const instance *problem )
{
    if ( i == j ) {
        errno = EFAULT;
        perror( CFATAL "_hardfix_xpos: i == j" );
        exit( EXIT_FAILURE );
    }

    if ( i > j ) return _hardfix_xpos( j, i, problem );

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
_add_constraints_hardfix ( const instance *problem, CPXENVptr env, CPXLPptr lp )
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
                fprintf( stderr, CFATAL "_add_constraints_hardfix: CPXnewcols [%s]\n", cname );
                exit( EXIT_FAILURE );
            }

            if ( CPXgetnumcols( env, lp ) - 1 != _hardfix_xpos( i, j, problem ) ) {
                fprintf( stderr, CFATAL "_add_constraints_hardfix: CPXgetnumcols [%s: x(%zu, %zu)]\n",
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
            fprintf( stderr, CFATAL "_add_constraints_hardfix: CPXnewrows [%s]\n", cname );
            exit( EXIT_FAILURE );
        }

        size_t lastrow = CPXgetnumrows( env, lp ) - 1;

        for ( size_t i = 0; i < problem->nnodes; ++i )
        {
            if ( i == h ) continue;
            if ( CPXchgcoef( env, lp, lastrow, _hardfix_xpos( i, h, problem ), 1.0 ) ) {
                fprintf( stderr, CFATAL "_add_constraints_hardfix: CPXchgcoef [%s: x(%zu, %zu)]\n",
                    cname, i + 1, h + 1 );
                exit( EXIT_FAILURE );
            }
        }
    }


    free( cname );
}


void
_add_subtour_constraints_hardfix ( const instance *problem,
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

    if (memchunk == NULL) {
        fprintf( stderr, CFATAL "_add_subtour_constraints: out of memory\n" );
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
                rmatind[nzcnt] = _hardfix_xpos( cnodes[i], cnodes[j], problem );
                rmatval[nzcnt] = 1.0;
                ++nzcnt;
            }
        }

        if ( CPXaddrows( env, lp, 0, 1, nzcnt, &rhs, &sense,
                         rmatbeg, rmatind, rmatval, NULL, &cname ) ) {
            fprintf( stderr, CFATAL "_add_subtour_constraints: CPXaddrows [SEC(%zu/%zu)]\n", k + 1, ncomps);
            exit( EXIT_FAILURE );
        }
    }

    free( memchunk );
}

void
hardfix_solve(instance *problem, CPXENVptr env,  CPXLPptr lp, size_t ncomps, double *xopt, size_t *next, size_t *comps, int visitednodes )
{
    //First start with a dummy solution
    for(int i=0; i<problem->nnodes-1;i++){
        problem->solution[i][0] = i;
        problem->solution[i][1] = i+1;
    }
    problem->solution[problem->nnodes-1][0] = problem->nnodes-1;
    problem->solution[problem->nnodes-1][1] = 0;

    int fix_ind[problem->nnodes];
    char fix_lu[problem->nnodes];
    double fix_bd[problem->nnodes];

    for(int k =0; k<100; k++) //intanto facciamolo cento volte
    {

        //fissa circa 90% dei lati
        int counter_fixed =0;
        for(int i =0; i<problem->nnodes;i++){
            int rnd_values = rand() % 10 + 1;
            if(rnd_values<=9){
                fix_ind[counter_fixed] = _hardfix_xpos(problem->solution[i][0], problem->solution[i][1], problem);
                fix_lu[counter_fixed] = 'L';
                fix_bd[counter_fixed] = 1.0;
                counter_fixed++;
            }
        }
        if(CPXchgbds(env, lp, counter_fixed, fix_ind, fix_lu, fix_bd)){
             fprintf( stderr, CFATAL "hardfix_solve: CPXchgbds");
                exit( EXIT_FAILURE );
        }

        //Solve
        for (size_t iter = 0; ncomps != 1; ++iter)
        {
            if ( CPXmipopt( env, lp ) ) {
                fprintf( stderr, CFATAL "loopBC_model: CPXmimopt error\n" );
                exit( EXIT_FAILURE );
            }


            visitednodes += CPXgetnodecnt( env, lp ) + 1;
            CPXsolution( env, lp, NULL, NULL, xopt, NULL, NULL, NULL );
            _xopt2subtours( problem, xopt, next, comps, &ncomps, _hardfix_xpos );

            _add_subtour_constraints_hardfix( problem, env, lp, next, comps, ncomps );
        }
        ncomps=0;

        _xopt2solution( xopt, problem, &_hardfix_xpos );
        //Sfissare lati
        for(int i=0; i<counter_fixed; i++){
            fix_bd[i]=0.0;
        }
        if(CPXchgbds(env, lp, counter_fixed, fix_ind, fix_lu, fix_bd)){
             fprintf( stderr, CFATAL "hardfix_solve: CPXchgbds");
                exit( EXIT_FAILURE );
        }
        double cost = compute_solution_cost( problem );
        fprintf(stderr, "solcost: %f\n", cost);

    }


}


void
hardfix_model ( instance *problem )
{
    int error;

    CPXENVptr env = CPXopenCPLEX( &error );
    CPXLPptr lp = CPXcreateprob( env, &error, problem->name ? problem->name : "TSP" );

    /* CPLEX PARAMETERS */
    tspconf_apply( env );

    /* BUILD MODEL */
    _add_constraints_hardfix(problem, env, lp);

    size_t ncomps = 0;
    double *xopt  = malloc( CPXgetnumcols( env, lp ) * sizeof( *xopt ) );
    size_t *next =  calloc( problem->nnodes, sizeof( *next ) );
    size_t *comps = calloc( problem->nnodes, sizeof( *comps ) );

    int visitednodes = 0;
    struct timeb start, end;
    ftime( &start );

    hardfix_solve(problem, env, lp, ncomps, xopt, next, comps, visitednodes);

    ftime( &end );

    _xopt2solution( xopt, problem, &_hardfix_xpos );

    free( xopt );

    problem->elapsedtime  = ( 1000. * ( end.time - start.time ) + end.millitm - start.millitm ) / 1000.;
    problem->visitednodes = visitednodes;
    problem->solcost      = compute_solution_cost( problem );

    CPXfreeprob( env, &lp );
    CPXcloseCPLEX( &env );
}
