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
#define EPS 1e-5


/* UTILITIES */

/*!
 * \brief Get the position of variable x(i,j)
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
int
xpos ( unsigned long i, unsigned long j, const instance *problem )
{
    if ( i == j ) {
        errno = EDOM;
        perror( "i == j in xpos" );
        exit( EXIT_FAILURE );
    }

    if ( i > j ) return xpos( j, i, problem );

    return i * problem->nnodes + j - ( ( i + 1 ) * ( i + 2 ) / 2UL );
}


/*!
 * \brief Get the position of variable y(i,j) of Flow1 model. 0-indexed.
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
unsigned long
ypos_flow1 ( unsigned long i, unsigned long j, const instance *problem )
{
    if ( i == j ) {
        errno = EDOM;
        perror( "i == j in ypos_flow1" );
        exit( EXIT_FAILURE );
    }

    unsigned long padding = xpos( problem->nnodes - 1, problem->nnodes, problem );
    return padding + i * (problem->nnodes - 1) + j - ( ( j > i ) ? 1 : 0 );
}


/*!
 * \brief Retrieve the solution after CPXmimopt was run.
 *
 *
 * \param xopt
 *      CPLEX incumbent solution.
 *
 * \param problem
 *     Pointer to the instance structure
 */
void
_xopt2solution ( double *xopt, instance *problem )
{
    unsigned long p = 0;

    for ( unsigned long i = 0; i < problem->nnodes; ++i ) {
        for ( unsigned long j = i + 1; j < problem->nnodes; ++j ) {
            if ( xopt[xpos( i, j, problem )] > .5 )
            {
                problem->solution[p][0] = i;
                problem->solution[p][1] = j;
                ++p;
            }
        }
    }
}


/*!
 * \brief Given a CPLEX-generated solution, create a more convenient representation.
 *
 *
 * \param problem
 *     Pointer to the instance structure.
 *
 * \param xopt
 *     CPLEX incumbent solution.
 *     `xstar[xpos(i, j)] == 1` iff the edge was selected.
 *
 * \param next
 *     Array of adjacencies to be filled.
 *     `next[i] = j` means that there is an arc going from node `i` to node `j`.
 *
 * \param comps
 *     Array of components indices to be filled.
 *     `comps[i] = k` means that node `i` belongs to connected component `k`.
 *
 * \param ncomps
 *     Pointer to an integer where to store the number of connected components in the solution.
 *     If 1, the solution is a tour.
 */
void
_xopt2subtours ( const instance *problem,
                 const double *xopt,
                 unsigned long *next,
                 unsigned long *comps,
                 unsigned long *ncomps )
{
    /*
     * First we parse array `xstar` and build the following adjacency list:
     *
     * ```
     *    0     1             i
     * +-----+-----+-------+-----+-------+
     * | a_0 | a_1 |  ...  | a_i |  ...  |
     * +-----+-----+-------+-----+-------+
     * | b_0 | b_1 |  ...  | b_i |  ...  |
     * +-----+-----+-------+-----+-------+
     * ```
     *
     * Where the i-th element contains the two adjacient nodes to `i`.
     * That is, node `i` is linked to nodes `a_i` and `b_i`.
     * We know there always are two and only two nodes because of the constraints on the solution.
     */

    // Create adjacency list and initialize its values to nonsense.
    unsigned long adj[problem->nnodes][2];
    for ( unsigned long i = 0; i < problem->nnodes; ++i ) {
        adj[i][0] = adj[i][1] = ULONG_MAX;
		comps[i] = 0;
    }

    // Fill adjacency list
    for ( unsigned long i = 0; i < problem->nnodes; ++i ) {
        if ( adj[i][1] != ULONG_MAX ) continue;
        for ( unsigned long j = i + 1; j < problem->nnodes; ++j ) {
            if ( xopt[xpos(i, j, problem)] > .5 ) {
                // Fill the free spot of adj[i] and adj[j].
                *( ( adj[i][0] == ULONG_MAX ) ? &adj[i][0] : &adj[i][1] ) = j;
                *( ( adj[j][0] == ULONG_MAX ) ? &adj[j][0] : &adj[j][1] ) = i;

                if ( adj[i][1] != ULONG_MAX ) {
                    break;
                }
            }
        }
    }

    // Fill `next` and `comps`
	*ncomps = 0;
    for ( unsigned long start = 0; start < problem->nnodes; ++start ) {
        if ( comps[start] != 0 ) continue;

        ++*ncomps;

        unsigned long from = start;
        unsigned long to = adj[start][0];

        do {
            next[from] = to;
            comps[from] = *ncomps;
            // One edge of `adj[to]` is equal to `from`. We care about the other.
            to = ( adj[to][0] ^ adj[to][1] ^ from );
            from = next[from];
        }
        while ( from != start );
    }
}


/*!
 * \brief Build a dummy CPLEX model. Utility for dummy_cplex_solution().
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
_dummy_cplex_instance2model ( instance *problem, CPXENVptr env, CPXLPptr lp )
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


/* SOLVERS */

void
dummy_solution ( instance *problem )
{
    for ( unsigned long i = 0; i < problem->nnodes; ++i ) {
        problem->solution[i][0] = i;
        problem->solution[i][1] = i + 1;
    }

    problem->solution[problem->nnodes - 1][1] = 0;
}


void
dummy_cplex_solution ( instance *problem )
{
    int error;

    CPXENVptr env = CPXopenCPLEX( &error );
    CPXLPptr lp = CPXcreateprob( env, &error, problem->name ? problem->name : "TSP" );

    _dummy_cplex_instance2model( problem, env, lp );

    //CPXwriteprob (env, lp, "bin/myprob.lp", NULL);

    if ( CPXmipopt( env, lp ) ) {
        fprintf( stderr, "CPXmimopt error\n" );
    }

    double *xopt = malloc( CPXgetnumcols( env, lp ) * sizeof( *xopt ) );
    CPXsolution( env, lp, NULL, NULL, xopt, NULL, NULL, NULL );

    _xopt2solution( xopt, problem );

    free(xopt);

    CPXfreeprob( env, &lp );
    CPXcloseCPLEX( &env );
}


/*!
 * \brief Add Flow1 commodity and coupling constraints.
 *
 *
 * This function adds the commodity and coupling constraints to the model.
 * That is, it introduces the continuous variables
 *
 * ```
 * y(i,j) = flow in arc (i,j)  for all i != j
 * ```
 *
 * The following constraints are added to the model:
 *
 * ```
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
_add_constraints_flow1( const instance *problem, CPXENVptr env, CPXLPptr lp )
{
    char continuous = 'C';
    double lb = 0.0;
    double ub = 1.0;
    double obj = 0.;
    double rhs;
    char sense;
    unsigned long lastrow;

    char *cname = calloc( MAX_CNAME_LENGTH, sizeof( *cname ) );

    // add binary var.s x(i,j) for i < j
    for ( unsigned long i = 0; i < problem->nnodes; ++i )
    {
        for ( unsigned long j = 0; j < problem->nnodes; ++j )
        {
            if ( i == j ) continue;

            snprintf( cname, MAX_CNAME_LENGTH, "y(%lu,%lu)", i + 1, j + 1 );

            if ( CPXnewcols( env, lp, 1, &obj, &lb, &ub, &continuous, &cname ) ) {
                fprintf( stderr,  "_add_constraints_flow1: CPXnewcols [%s]\n", cname );
                exit( EXIT_FAILURE );
            }

            if ( CPXgetnumcols( env, lp ) - 1 != ypos_flow1( i, j, problem ) ) {
                fprintf( stderr,  "_add_constraints_flow1: CPXgetnumcols [%s]\n", cname );
                exit( EXIT_FAILURE );
            }
        }
    }

    // add constraints (1)
    // y(i,j) - (n - 1) x(i,j) <= 0  for all i != j
    rhs = .0;
    sense = 'L';

    for ( unsigned long i = 0; i < problem->nnodes; ++i )
    {
        for ( unsigned long j = 0; j < problem->nnodes; ++j )
        {
            if ( i == j ) continue;

            lastrow = CPXgetnumrows( env, lp );

            snprintf( cname, MAX_CNAME_LENGTH, "flow1_1(%lu,%lu)", i + 1, j + 1 );
            if ( CPXnewrows( env, lp, 1, &rhs, &sense, NULL, &cname ) ) {
                fprintf( stderr,  "_add_constraints_flow1: CPXnewrows [%s]\n", cname );
                exit( EXIT_FAILURE );
            }

            if ( CPXchgcoef( env, lp, lastrow, xpos( i, j, problem ), - (double) ( problem->nnodes - 1 ) ) ) {
                fprintf( stderr,  "_add_constraints_flow1: CPXchgcoef [%s: x]\n", cname);
                exit( EXIT_FAILURE );
            }

            if ( CPXchgcoef( env, lp, lastrow, ypos_flow1( i, j, problem ), 1.0 ) ) {
                fprintf( stderr,  "_add_constraints_flow1: CPXchgcoef [%s: y]\n", cname );
                exit( EXIT_FAILURE );
            }
        }
    }

    // add constraint (2)
    // Sum{ y(1,j) }_{ j | j != 1 } = n - 1
    rhs = (double) ( problem->nnodes - 1 );
    sense = 'E';

    lastrow = CPXgetnumrows( env, lp );

    snprintf( cname, MAX_CNAME_LENGTH, "flow1_2" );
    if ( CPXnewrows( env, lp, 1, &rhs, &sense, NULL, &cname ) ) {
        fprintf( stderr,  "_add_constraints_flow1: CPXnewrows [%s]\n", cname );
        exit( EXIT_FAILURE );
    }

    for ( unsigned long j = 0; j < problem->nnodes; ++j )
    {
        if ( 0 == j ) continue;

        if ( CPXchgcoef( env, lp, lastrow, ypos_flow1( 0, j, problem ), 1.0 ) ) {
            fprintf( stderr,  "_add_constraints_flow1: CPXchgcoef [%s: %lu]\n", cname, j );
            exit( EXIT_FAILURE );
        }
    }

    // add constraints (3)
    // Sum{ y(i,h) }_{ i | i != h } - Sum{ y(h,j) }_{ j | j != h } = 1  for all h != 1
    rhs = 1.0;
    sense = 'E';

    for ( unsigned long h = 1; h < problem->nnodes; ++h )
    {
        snprintf( cname, MAX_CNAME_LENGTH, "flow1_3(%lu)", h );
        if ( CPXnewrows( env, lp, 1, &rhs, &sense, NULL, &cname ) ) {
            fprintf( stderr,  "_add_constraints_flow1: CPXnewrows [%s]\n", cname );
            exit( EXIT_FAILURE );
        }

        lastrow = CPXgetnumrows( env, lp );

        for ( unsigned long i = 0; i < problem->nnodes; ++i )
        {
            if ( i == h ) continue;

            if ( CPXchgcoef( env, lp, lastrow, ypos_flow1( i, h, problem ), 1.0 ) ) {
                fprintf( stderr,  "_add_constraints_flow1: CPXchgcoef [%s: i=%lu]\n", cname, i );
                exit( EXIT_FAILURE );
            }

            if ( CPXchgcoef( env, lp, lastrow, ypos_flow1( h, i, problem ), -1.0 ) ) {
                fprintf( stderr,  "_add_constraints_flow1: CPXchgcoef [%s: j=%lu]\n", cname, i );
                exit( EXIT_FAILURE );
            }
        }
    }

    free( cname );
}


void
flow1_solution ( instance *problem )
{
    int error;

    CPXENVptr env = CPXopenCPLEX( &error );
    CPXLPptr lp = CPXcreateprob( env, &error, problem->name ? problem->name : "TSP" );

    /* Add conventional constraints */
    _dummy_cplex_instance2model( problem, env, lp );
    _add_constraints_flow1( problem, env, lp );

    if ( CPXmipopt( env, lp ) ) {
        fprintf( stderr, "CPXmimopt error\n" );
        exit( EXIT_FAILURE );
    }

    double *xopt = malloc( CPXgetnumcols( env, lp ) * sizeof( *xopt ) );
    CPXsolution( env, lp, NULL, NULL, xopt, NULL, NULL, NULL );

    _xopt2solution( xopt, problem );

    free(xopt);

    CPXfreeprob( env, &lp );
    CPXcloseCPLEX( &env );
}