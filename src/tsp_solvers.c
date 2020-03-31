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


#define EPS 1e-5


/* UTILITIES */

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
unsigned long
xpos_dummy ( unsigned long i, unsigned long j, const instance *problem )
{
    if ( i == j ) {
        errno = EDOM;
        perror( "i == j in xpos" );
        exit( EXIT_FAILURE );
    }

    if ( i > j ) return xpos_dummy( j, i, problem );

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
unsigned long
ypos_flow1 ( unsigned long i, unsigned long j, const instance *problem )
{
    if ( i == j ) {
        errno = EDOM;
        perror( "i == j in ypos_flow1" );
        exit( EXIT_FAILURE );
    }

    unsigned long padding = xpos_dummy( problem->nnodes - 1, problem->nnodes, problem );
    return padding + i * (problem->nnodes - 1) + j - ( ( j > i ) ? 1UL : 0UL );
}


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
unsigned long
xpos_mtz ( unsigned long i, unsigned long j, const instance *problem )
{
    return i * problem->nnodes + j;
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
int
upos_mtz ( unsigned long i, const instance *problem )
{
    return  problem->nnodes * problem->nnodes + i;
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
 *
 * \param pos
 *     Pointer to a function that given coordinates `i` and `j` returns the position in \p xopt fo `x(i,j)`.
 */
void
_xopt2solution ( const double  *xopt,
                 instance      *problem,
                 unsigned long (*pos)(unsigned long, unsigned long, const instance *) )
{
    unsigned long p = 0;

    for ( unsigned long i = 0; i < problem->nnodes; ++i ) {
        for ( unsigned long j = i + 1; j < problem->nnodes; ++j ) {
            // Try both (i,j) and (j,i) because some models may use oriented arcs
            if ( xopt[(*pos)( i, j, problem )] > .5 || xopt[(*pos)( j, i, problem )] > .5 )
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
 * Given the incumbent solution \p xopt, where `xopt[e] = 1 <==> edge e was selected`, populate provided
 *  arrays \p next and \p comps so that `next[i] = j <==> the tour goes from node i to node j` and
 * `comps[i] = k <==> node i is part of the k-th subtour`.
 *
 * The number of subtours is written to \p ncomps, hence \p xopt is a valid TSP solution iff `ncomps == 1`.
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
                 const double   *xopt,
                 unsigned long  *next,
                 unsigned long  *comps,
                 unsigned long  *ncomps )
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
            if ( xopt[xpos_dummy(i, j, problem)] > .5 ) {
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
 * \brief Add degree constraints to the problem.
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
_add_constraints_dummy ( const instance *problem, CPXENVptr env, CPXLPptr lp )
{
    char ctype = CPX_BINARY;
    double lb = 0.0;
    double ub = 1.0;

    char *cname = calloc( CPX_STR_PARAM_MAX, sizeof(  *cname ) );

    // add binary vars x(i,j) for 0 <= i < j <= N
    for ( unsigned long i = 0; i < problem->nnodes; ++i )
    {
        for ( unsigned long j = i + 1; j < problem->nnodes; ++j )
        {
            snprintf( cname, CPX_STR_PARAM_MAX, "x(%lu,%lu)", i + 1, j + 1 );
            double obj = _euclidean_distance(
                problem->xcoord[i],
                problem->ycoord[i],
                problem->xcoord[j],
                problem->ycoord[j]
            );

            if ( CPXnewcols( env, lp, 1, &obj, &lb, &ub, &ctype, &cname ) ) {
                fprintf( stderr,  "_add_constraints_dummy: CPXnewcols [%s]\n", cname );
                exit( EXIT_FAILURE );
            }

            if ( CPXgetnumcols( env, lp ) - 1 != xpos_dummy( i, j, problem ) ) {
                fprintf( stderr,  "_add_constraints_dummy: CPXgetnumcols [%s: x(%lu, %lu)]\n", cname, i + 1, j + 1 );
                exit( EXIT_FAILURE );
            }
        }
    }

    // add degree constraints
    double rhs = 2.0;
    char sense = 'E';

    for ( unsigned long h = 0; h < problem->nnodes; ++h )
    {
        snprintf( cname, CPX_STR_PARAM_MAX, "degree(%lu)", h + 1 );
        if ( CPXnewrows( env, lp, 1, &rhs, &sense, NULL, &cname ) ) {
            fprintf( stderr,  "_add_constraints_dummy: CPXnewrows [%s]\n", cname );
            exit( EXIT_FAILURE );
        }

        unsigned long lastrow = CPXgetnumrows( env, lp ) - 1;

        for ( unsigned long i = 0; i < problem->nnodes; ++i )
        {
            if ( i == h ) continue;
            if ( CPXchgcoef( env, lp, lastrow, xpos_dummy( i, h, problem ), 1.0 ) ) {
                fprintf( stderr,  "_add_constraints_dummy: CPXnewrows [%s: x(%lu, %lu)]\n", cname, i + 1, h + 1 );
                exit( EXIT_FAILURE );
            }
        }
    }

    free( cname );
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

    char *cname = calloc( CPX_STR_PARAM_MAX, sizeof( *cname ) );

    // add binary var x(i,j) for all i, j
    for ( unsigned long i = 0; i < problem->nnodes; ++i )
    {
        for ( unsigned long j = 0; j < problem->nnodes; ++j )
        {
            snprintf( cname, CPX_STR_PARAM_MAX, "x(%lu,%lu)", i + 1, j + 1 );
            double obj = _euclidean_distance(
                problem->xcoord[i],
                problem->ycoord[i],
                problem->xcoord[j],
                problem->ycoord[j]
            );

            if( i != j ) {
                if ( CPXnewcols( env, lp, 1, &obj, &lb, &ub, &ctype, &cname ) ) {
                    fprintf( stderr,  "_add_constraints_mtz: CPXnewcols [%s]\n", cname );
                    exit( EXIT_FAILURE );
                }
            } else {
                double temp_ub = 0.0;
                if ( CPXnewcols( env, lp, 1, &obj, &lb, &temp_ub, &ctype, &cname ) ) {
                    fprintf( stderr,  "_add_constraints_mtz: CPXnewcols [%s]\n", cname );
                    exit( EXIT_FAILURE );
                }
            }

            if ( CPXgetnumcols( env, lp ) - 1 != xpos_mtz( i, j, problem ) ) {
                fprintf( stderr,  "_add_constraints_mtz: CPXgetnumcols [%s: x(%lu, %lu)]\n", cname, i + 1, j + 1 );
                exit( EXIT_FAILURE );
            }
        }
    }

    // add degree constraints
    double rhs = 1.0;
    char sense = 'E';

    // Sum[ x(i,h) ]_{i | i != h} = 1  for all h
    for ( unsigned long h = 0; h < problem->nnodes; ++h )
    {
        unsigned long lastrow = CPXgetnumrows( env, lp );

        snprintf( cname, CPX_STR_PARAM_MAX, "degree(%lu)", h + 1 );
        if ( CPXnewrows( env, lp, 1, &rhs, &sense, NULL, &cname ) ) {
            fprintf( stderr,  "_add_constraints_mtz: CPXnewrows [%s]\n", cname );
            exit( EXIT_FAILURE );
        }

        for ( unsigned long i = 0; i < problem->nnodes; ++i )
        {
            if ( i == h ) continue;

            if ( CPXchgcoef( env, lp, lastrow, xpos_mtz( i, h, problem ), 1.0 ) ) {
                fprintf( stderr,  "_add_constraints_mtz: CPXchgcoef [%s: x(%lu, %lu)]\n", cname, i + 1, h + 1 );
                exit( EXIT_FAILURE );
            }
        }
    }

    // Sum[ x(h,i) ]_{i | i != h} = 1  for all h
    for ( unsigned long h = 0; h < problem->nnodes; ++h )
    {
        unsigned long lastrow = CPXgetnumrows( env, lp );

        snprintf( cname, CPX_STR_PARAM_MAX, "degree(%lu)", h + 1 );
        if ( CPXnewrows( env, lp, 1, &rhs, &sense, NULL, &cname ) ) {
            fprintf( stderr,  "_add_constraints_mtz: CPXnewrows [%s]\n", cname );
            exit( EXIT_FAILURE );
        }

        for ( unsigned long i = 0; i < problem->nnodes; ++i )
        {
            if ( i == h ) continue;

            if ( CPXchgcoef( env, lp, lastrow, xpos_mtz( h, i, problem ), 1.0 ) ) {
                fprintf( stderr,  "_add_constraints_mtz: CPXchgcoef [%s: x(%lu, %lu)]\n", cname, h + 1, i + 1 );
                exit( EXIT_FAILURE );
            }
        }
    }

    // add continuous vars u(i)
    ctype = CPX_CONTINUOUS;
    lb = 2.0;
    ub = problem->nnodes;

    for ( unsigned long i = 0; i < problem->nnodes; ++i )
    {

        snprintf( cname, CPX_STR_PARAM_MAX, "u(%lu)", i + 1 );
        double obj = 0.0;

        if( i == 0 ) {
            double temp_lb = 1.0;
            double temp_ub = 1.0;

            if( CPXnewcols( env, lp, 1, &obj, &temp_lb, &temp_ub, &ctype, &cname ) ) {
                fprintf( stderr,  "_add_constraints_mtz: CPXnewcols [%s]\n", cname );
                exit( EXIT_FAILURE );
            }

        } else {

            if( CPXnewcols( env, lp, 1, &obj, &lb, &ub, &ctype, &cname ) ) {
                fprintf( stderr,  "_add_constraints_mtz: CPXnewcols [%s]\n", cname );
                exit( EXIT_FAILURE );
            }

        }

        if ( CPXgetnumcols( env, lp ) - 1 != upos_mtz( i, problem ) ) {
            fprintf( stderr,  "_add_constraints_mtz: CPXgetnumcols [%s]\n", cname );
            exit( EXIT_FAILURE );
        }
    }

    // add sequencs constraints
    // u(i) - u(j) + n x(i,j) <= n - 1
    rhs = problem->nnodes - 1;
    sense = 'L';

    for ( unsigned long j = 1; j < problem->nnodes; ++j )
    {
        for ( unsigned long i = 1; i < problem->nnodes; ++i ){

            if ( i == j ) continue;

            snprintf( cname, CPX_STR_PARAM_MAX, "degree(%lu)", j + 1 );

            if ( CPXnewrows( env, lp, 1, &rhs, &sense, NULL, &cname ) ) {
                fprintf( stderr,  "_add_constraints_mtz: CPXnewrows [%s]\n", cname );
                exit( EXIT_FAILURE );
            }

            unsigned long lastrow = CPXgetnumrows( env, lp ) - 1;

            if ( CPXchgcoef( env, lp, lastrow, upos_mtz( j, problem ), -1.0 ) ) {
                fprintf( stderr,  "_add_constraints_mtz: CPXchgcoef [%s: u(%lu)]\n", cname, j + 1 );
                exit( EXIT_FAILURE );
            }

            if ( CPXchgcoef( env, lp, lastrow, upos_mtz( i, problem ), 1.0 ) ) {
                fprintf( stderr,  "_add_constraints_mtz: CPXchgcoef [%s: u(%lu)]\n", cname, i + 1 );
                exit( EXIT_FAILURE );
            }

            if ( CPXchgcoef( env, lp, lastrow, xpos_mtz( i, j, problem ), problem->nnodes ) ) {
                fprintf( stderr,  "_add_constraints_mtz: CPXchgcoef [%s: x(%lu,%lu)]\n", cname, i + 1, j + 1 );
                exit( EXIT_FAILURE );
            }
        }
    }

    free( cname );
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
    // Add degree constraints
    _add_constraints_dummy( problem, env, lp );

    char ctype = CPX_CONTINUOUS;
    double lb = 0.0;
    double ub = CPX_INFBOUND;
    double obj = 0.;

    double rhs;
    char sense;
    unsigned long lastrow;

    char *cname = calloc( CPX_STR_PARAM_MAX, sizeof( *cname ) );

    // add binary var.s x(i,j) for i < j
    for ( unsigned long i = 0; i < problem->nnodes; ++i )
    {
        for ( unsigned long j = 0; j < problem->nnodes; ++j )
        {
            if ( i == j ) continue;

            snprintf( cname, CPX_STR_PARAM_MAX, "y(%lu,%lu)", i + 1, j + 1 );

            if ( CPXnewcols( env, lp, 1, &obj, &lb, &ub, &ctype, &cname ) ) {
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

            snprintf( cname, CPX_STR_PARAM_MAX, "flow1_1(%lu,%lu)", i + 1, j + 1 );
            if ( CPXnewrows( env, lp, 1, &rhs, &sense, NULL, &cname ) ) {
                fprintf( stderr,  "_add_constraints_flow1: CPXnewrows [%s]\n", cname );
                exit( EXIT_FAILURE );
            }

            lastrow = CPXgetnumrows( env, lp ) - 1;

            if ( CPXchgcoef( env, lp, lastrow, xpos_dummy( i, j, problem ), - (double) ( problem->nnodes - 1 ) ) ) {
                fprintf( stderr,  "_add_constraints_flow1: CPXchgcoef [%s: x(%lu,%lu)]\n", cname, i + 1, j + 1);
                exit( EXIT_FAILURE );
            }

            if ( CPXchgcoef( env, lp, lastrow, ypos_flow1( i, j, problem ), 1.0 ) ) {
                fprintf( stderr,  "_add_constraints_flow1: CPXchgcoef [%s: y(%lu,%lu)]\n", cname, i + 1, j + 1 );
                exit( EXIT_FAILURE );
            }
        }
    }

    // add constraint (2)
    // Sum{ y(1,j) }_{ j | j != 1 } = n - 1
    rhs = (double) ( problem->nnodes - 1 );
    sense = 'E';

    lastrow = CPXgetnumrows( env, lp );

    snprintf( cname, CPX_STR_PARAM_MAX, "flow1_2" );
    if ( CPXnewrows( env, lp, 1, &rhs, &sense, NULL, &cname ) ) {
        fprintf( stderr,  "_add_constraints_flow1: CPXnewrows [%s]\n", cname );
        exit( EXIT_FAILURE );
    }

    for ( unsigned long j = 1; j < problem->nnodes; ++j )
    {
        if ( CPXchgcoef( env, lp, lastrow, ypos_flow1( 0, j, problem ), 1.0 ) ) {
            fprintf( stderr,  "_add_constraints_flow1: CPXchgcoef [%s: y(1,%lu)]\n", cname, j + 1 );
            exit( EXIT_FAILURE );
        }
    }

    // add constraints (3)
    // Sum{ y(i,h) }_{ i | i != h } - Sum{ y(h,j) }_{ j | j != h } = 1  for all h != 1
    rhs = 1.0;
    sense = 'E';

    for ( unsigned long h = 1; h < problem->nnodes; ++h )
    {
        snprintf( cname, CPX_STR_PARAM_MAX, "flow1_3(%lu)", h + 1 );
        if ( CPXnewrows( env, lp, 1, &rhs, &sense, NULL, &cname ) ) {
            fprintf( stderr,  "_add_constraints_flow1: CPXnewrows [%s]\n", cname );
            exit( EXIT_FAILURE );
        }

        lastrow = CPXgetnumrows( env, lp ) - 1;

        for ( unsigned long i = 0; i < problem->nnodes; ++i )
        {
            if ( i == h ) continue;

            if ( CPXchgcoef( env, lp, lastrow, ypos_flow1( i, h, problem ), 1.0 ) ) {
                fprintf( stderr,  "_add_constraints_flow1: CPXchgcoef [%s: y(%lu,%lu)]\n", cname, i + 1, h + 1);
                exit( EXIT_FAILURE );
            }

            if ( CPXchgcoef( env, lp, lastrow, ypos_flow1( h, i, problem ), -1.0 ) ) {
                fprintf( stderr,  "_add_constraints_flow1: CPXchgcoef [%s: y(%lu,%lu)]\n", cname, h + 1, i + 1 );
                exit( EXIT_FAILURE );
            }
        }
    }

    free( cname );
}


/* SOLVERS */


void
random_model ( instance *problem )
{
    for ( unsigned long i = 0; i < problem->nnodes; ++i )
    {
        problem->solution[i][0] = i;
        problem->solution[i][1] = i + 1;
    }

    problem->solution[problem->nnodes - 1][1] = 0;
}


void
dummy_model ( instance *problem )
{
    int error;

    CPXENVptr env = CPXopenCPLEX( &error );
    CPXLPptr lp = CPXcreateprob( env, &error, problem->name ? problem->name : "TSP" );

    _add_constraints_dummy( problem, env, lp );

    if ( CPXmipopt( env, lp ) ) {
        fprintf( stderr, "dummy_cplex_solution CPXmimopt error\n" );
    }

    double *xopt = malloc( CPXgetnumcols( env, lp ) * sizeof( *xopt ) );
    CPXsolution( env, lp, NULL, NULL, xopt, NULL, NULL, NULL );

    _xopt2solution( xopt, problem, &xpos_dummy );

    free( xopt );

    CPXfreeprob( env, &lp );
    CPXcloseCPLEX( &env );
}


void
mtz_model ( instance *problem )
{

    int error;

    CPXENVptr env = CPXopenCPLEX( &error );
    CPXLPptr lp = CPXcreateprob( env, &error, problem->name ? problem->name : "TSP" );

    _add_constraints_mtz( problem, env, lp );

    if ( CPXmipopt( env, lp ) ) {
        fprintf( stderr, "miller_tucker_cplex_solution CPXmimopt error\n" );
    }

    double *xopt = malloc( CPXgetnumcols( env, lp ) * sizeof( *xopt ) );
    CPXsolution( env, lp, NULL, NULL, xopt, NULL, NULL, NULL );

    _xopt2solution(xopt, problem, &xpos_mtz);

    free( xopt );

    CPXfreeprob( env, &lp );
    CPXcloseCPLEX( &env );
}


void
flow1_model ( instance *problem )
{
    int error;

    CPXENVptr env = CPXopenCPLEX( &error );
    CPXLPptr lp = CPXcreateprob( env, &error, problem->name ? problem->name : "TSP" );

    /* Add conventional constraints */
    _add_constraints_flow1( problem, env, lp );

    if ( CPXmipopt( env, lp ) ) {
        fprintf( stderr, "CPXmimopt error\n" );
        exit( EXIT_FAILURE );
    }

    double *xopt = malloc( CPXgetnumcols( env, lp ) * sizeof( *xopt ) );
    CPXsolution( env, lp, NULL, NULL, xopt, NULL, NULL, NULL );

    _xopt2solution( xopt, problem, &xpos_dummy );

    free( xopt );

    CPXfreeprob(env, &lp);
    CPXcloseCPLEX(&env);
}
