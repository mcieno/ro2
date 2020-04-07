/*
 * \brief   Basic Branch and Bound model.
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
#include "tspplot.h"


/*!
 * \brief Get the position of variable x(i,j) in dummy model.
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
_branch_and_bound_xpos ( size_t i, size_t j, const instance *problem )
{
    if ( i == j ) {
        errno = EFAULT;
        perror( CFATAL "_flow1_xpos: i == j" );
        exit( EXIT_FAILURE );
    }

    if ( i > j ) return _branch_and_bound_xpos( j, i, problem );

    return i * problem->nnodes + j - ( ( i + 1 ) * ( i + 2 ) / 2UL );
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
_add_constraints_branch_and_bound( const instance *problem, CPXENVptr env, CPXLPptr lp )
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
                fprintf( stderr, CFATAL "_add_constraints_flow1: CPXnewcols [%s]\n", cname );
                exit( EXIT_FAILURE );
            }

            if ( CPXgetnumcols( env, lp ) - 1 != _branch_and_bound_xpos( i, j, problem ) ) {
                fprintf( stderr, CFATAL "_add_constraints_flow1: CPXgetnumcols [%s: x(%zu, %zu)]\n",
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
            fprintf( stderr, CFATAL "_add_constraints_flow1: CPXnewrows [%s]\n", cname );
            exit( EXIT_FAILURE );
        }

        size_t lastrow = CPXgetnumrows( env, lp ) - 1;

        for ( size_t i = 0; i < problem->nnodes; ++i )
        {
            if ( i == h ) continue;
            if ( CPXchgcoef( env, lp, lastrow, _branch_and_bound_xpos( i, h, problem ), 1.0 ) ) {
                fprintf( stderr, CFATAL "_add_constraints_flow1: CPXnewrows [%s: x(%zu, %zu)]\n",
                    cname, i + 1, h + 1 );
                exit( EXIT_FAILURE );
            }
        }
    }


    free( cname );
}

void
add_soubtour_constraints(const instance *problem, CPXENVptr env, CPXLPptr lp, size_t *ncomps, size_t *comps)
{
    if(*ncomps==1){
        return;
    }


    size_t nodes[problem->nnodes]; //nodes in a subtour
    size_t num_nodes; //number of nodes in a subtour
    char *cname = calloc( CPX_STR_PARAM_MAX, sizeof( *cname ) );
    double rhs;
    char sense;
    
   

    for(int i=0; i<*ncomps; i++){
        num_nodes =0;
        for(size_t j =0; j<problem->nnodes; j++){ //fetch all nodes in a comp
            if(comps[j]==i+1){
                nodes[num_nodes]=j;
                num_nodes++;
            }
        }


        //add subtour constraint of the comp
        sense = 'L';
        rhs = num_nodes-1;

        snprintf( cname, CPX_STR_PARAM_MAX, "subtour constraint");
        if ( CPXnewrows( env, lp, 1, &rhs, &sense, NULL, &cname ) ) {
            fprintf( stderr, CFATAL "_add_soubtour_constraints in branch and bound: CPXnewrows [%s]\n", cname );
            exit( EXIT_FAILURE );
            }

        size_t lastrow = CPXgetnumrows( env, lp ) - 1;


        for ( size_t j = 0; j < num_nodes-1; ++j )
        {
            for(size_t k =j+1; k<num_nodes; ++k){
            
                if ( CPXchgcoef( env, lp, lastrow, _branch_and_bound_xpos( nodes[j], nodes[k], problem ), 1.0 ) ) {
                    fprintf( stderr, CFATAL "_add_constraints_subtours branch and bound: CPXnewrows [%s: x(%zu, %zu)]\n",
                        cname, j + 1, k + 1 );
                    exit( EXIT_FAILURE );
                }
            }
        }


    }


    free(cname);

}

double
branch_and_bound_model  ( instance *problem )
{
    int error;

    CPXENVptr env = CPXopenCPLEX( &error );
    CPXLPptr lp = CPXcreateprob( env, &error, problem->name ? problem->name : "TSP" );

    _add_constraints_branch_and_bound( problem, env, lp );

    size_t ncomps = 9999;
    double *xopt = malloc( CPXgetnumcols( env, lp ) * sizeof( *xopt ) );
    size_t *comps = calloc( problem->nnodes, sizeof( *comps ) );

    size_t *next =  calloc(problem->nnodes, sizeof(*next));


    struct timeb start, end;
    ftime( &start );

    while(ncomps!=1){
        if ( CPXmipopt( env, lp ) ) {
            fprintf( stderr, CFATAL "flow1_model: CPXmimopt error\n" );
            exit( EXIT_FAILURE );
        }
        CPXsolution( env, lp, NULL, NULL, xopt, NULL, NULL, NULL );
        _xopt2subtours(problem, xopt, next, comps, &ncomps, _branch_and_bound_xpos);

       
        add_soubtour_constraints(problem, env, lp, &ncomps, comps);

        //_xopt2solution( xopt, problem, &_branch_and_bound_xpos );
        //plot_solution(problem);
    
    }

    ftime( &end );

    _xopt2solution( xopt, problem, &_branch_and_bound_xpos );

    free( xopt );

    CPXfreeprob(env, &lp);
    CPXcloseCPLEX(&env);

    return ( 1000. * ( end.time - start.time ) + end.millitm - start.millitm ) / 1000.;
}