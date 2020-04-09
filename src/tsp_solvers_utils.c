/*
 * \brief   Utility functions for TSP solvers
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


#define EPS 1e-5
double timelimit = __DBL_MAX__;


void
_xopt2solution ( const double *xopt,
                 instance     *problem,
                 size_t       (*pos)(size_t, size_t, const instance *) )
{
    size_t p = 0;

    for ( size_t i = 0; i < problem->nnodes; ++i ) {
        for ( size_t j = i + 1; j < problem->nnodes; ++j ) {
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


void
_xopt2subtours ( const instance *problem,
                 const double   *xopt,
                 size_t         *next,
                 size_t         *comps,
                 size_t         *ncomps,
                 size_t         (*pos)(size_t, size_t, const instance *))
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
    size_t adj[problem->nnodes][2];
    for ( size_t i = 0; i < problem->nnodes; ++i ) {
        adj[i][0] = adj[i][1] = SIZE_MAX;
        comps[i] = 0;
    }

    // Fill adjacency list
    for ( size_t i = 0; i < problem->nnodes; ++i ) {
        if ( adj[i][1] != SIZE_MAX ) continue;
        for ( size_t j = i + 1; j < problem->nnodes; ++j ) {
            if ( xopt[pos(i, j, problem)] > .5 || xopt[pos(j, i, problem)] > .5 ) {
                // Fill the free spot of adj[i] and adj[j].
                *( ( adj[i][0] == SIZE_MAX ) ? &adj[i][0] : &adj[i][1] ) = j;
                *( ( adj[j][0] == SIZE_MAX ) ? &adj[j][0] : &adj[j][1] ) = i;

                if ( adj[i][1] != SIZE_MAX ) {
                    break;
                }
            }
        }
    }

    // Fill `next` and `comps`
    *ncomps = 0;
    for ( size_t start = 0; start < problem->nnodes; ++start ) {
        if ( comps[start] != 0 ) continue;

        ++*ncomps;

        size_t from = start;
        size_t to = adj[start][0];

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
