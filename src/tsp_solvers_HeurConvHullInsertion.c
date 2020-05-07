/*
 * \brief   Convex Hull Insertion heuristic method.
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

#include "tsp_solvers.h"
#include "logging.h"
#include "tsp.h"
#include "tspconf.h"
#include "tspplot.h"


static unsigned int __SEED;

typedef struct
{
    int available;
    int index;
    double x;
    double y;
}
node_t;

/**
 * Comparator for the node_t type.
 *
 * @param a void pointer to an `node_t` element.
 * @param b void pointer to an `node_t` element.
 * @returns negative if `x(a) < x(b)` or `x(a) == x(b) AND y(a) < y(b)`,
 *          positive if `x(a) > x(b)` or `x(a) == x(b) AND y(a) > y(b)`,
 *          0 if `x(a) == x(b) AND y(a) == y(b)`.
 */
int
_cmp_HeurConvHullInsertion( const void *a, const void *b ) {
    const node_t *na = (const node_t *) a;
    const node_t *nb = (const node_t *) b;

    return na->x < nb->x
               ? -1
               : na->x > nb->x
                     ? +1
                     : na->y > nb->y
                           ? -1
                           : na->y > nb->y
                                 ? +1
                                 : 0;
}


/*********************************************************************************************************************/
// Copyright 2001 softSurfer, 2012 Dan Sunday
// This code may be freely used and modified for any purpose
// providing that this copyright notice is included with it.
// SoftSurfer makes no warranty for this code, and cannot be held
// liable for any real or imagined damage resulting from its use.
// Users of this code must verify correctness for their application.


// isLeft(): tests if a point is Left|On|Right of an infinite line.
//    Input:  three points P0, P1, and P2
//    Return: >0 for P2 left of the line through P0 and P1
//            =0 for P2 on the line
//            <0 for P2 right of the line
//    See: Algorithm 1 on Area of Triangles
double
isLeft( size_t P0, size_t P1, size_t P2, instance *problem )
{
    return ( problem->xcoord[P1] - problem->xcoord[P0] ) * ( problem->ycoord[P2] - problem->ycoord[P0] )
            - ( problem->xcoord[P2] - problem->xcoord[P0] ) * ( problem->ycoord[P1] - problem->ycoord[P0] );
}

// chainHull_2D(): Andrew s monotone chain 2D convex hull algorithm
//     Input:  P[] = an array of 2D points
//                  presorted by increasing x and y-coordinates
//             n =  the number of points in P[]
//     Output: H[] = an array of the convex hull vertices (max is n)
//     Return: the number of points in H[]
size_t
chainHull_2D( size_t *P, size_t n, size_t *H, instance *problem )
{
    // indices for bottom and top of the stack H
    size_t    bot = 0, top = SIZE_MAX;
    size_t    i;

    // Get the indices of points with min x-coord and min|max y-coord
    size_t minmin = 0, minmax;
    double xmin = problem->xcoord[P[0]];
    for ( i = 1; i < n && problem->xcoord[P[i]] != xmin; ++i )
        ;
    minmax = i - 1;

    if ( minmax == n - 1 ) {
        // degenerate case: all x-coords == xmin
        H[++top] = P[minmin];
        // nontrivial segment
        if ( problem->ycoord[P[minmax]] != problem->ycoord[P[minmin]] )  H[++top] =  P[minmax];
        // add polygon endpoint
        H[++top] = P[minmin];
        return top + 1;
    }

    // Get the indices of points with max x-coord and min|max y-coord
    size_t maxmin, maxmax = n-1;
    double xmax = problem->xcoord[P[n - 1]];
    for ( i = n - 2; i >= 0 && problem->xcoord[P[i]] != xmax; --i )
        ;
    maxmin = i + 1;

    // Compute the lower hull on the stack H
    H[++top] = P[minmin];      // push  minmin point onto stack
    i = minmax;

    while ( ++i <= maxmin ) {
        // the lower line joins P[minmin]  with P[maxmin]

        // ignore P[i] above or on the lower line
        if ( isLeft( P[minmin], P[maxmin], P[i], problem )  >= 0 && i < maxmin )  continue;

        while (top > 0) {
            // there are at least 2 points on the stack

            // test if  P[i] is left of the line at the stack top
            if ( isLeft( H[top-1], H[top], P[i], problem ) > 0 )
                 break;         // P[i] is a new hull  vertex
            else
                 top--;         // pop top point off  stack
        }

        H[++top] = P[i];        // push P[i] onto stack
    }

    // Next, compute the upper hull on the stack H above  the bottom hull
    if ( maxmax != maxmin )     // if  distinct xmax points
        H[++top] = P[maxmax];   // push maxmax point onto stack
    bot = top;                  // the bottom point of the upper hull stack
    i = maxmin;

    --i;
    while ( i >= minmax && i <= maxmin )
    {

        // the upper line joins P[maxmax]  with P[minmax]

        if ( isLeft( P[maxmax], P[minmax], P[i], problem ) >= 0 && i > minmax ) {
            --i;
            continue;           // ignore P[i] below or on the upper line
        }

        // at least 2 points on the upper stack
        while ( top > bot ) {
            // test if  P[i] is left of the line at the stack top
            if ( isLeft(  H[top - 1], H[top], P[i], problem ) > 0 )  break;  // P[i] is a new hull vertex
            else                                                     --top;  // pop top point off stack
        }

        // push P[i] onto stack
        H[++top] = P[i];
        --i;
    }

    // push  joining endpoint onto stack
    if ( minmax != minmin )  H[++top] = P[minmin];

    return top + 1;
}
/*********************************************************************************************************************/


void
HeurConvHullInsertion_solve( instance *problem, double *bestcost, size_t *H, size_t k )
{

    //initialize edges and nodes array
    size_t edges[problem->nnodes][2];
    for ( size_t i = 0; i < k - 1; ++i ) {
        edges[i][0] = H[i];
        edges[i][1] = H[i + 1];
    }
    size_t edges_counter = k - 1;

    size_t nodes[problem->nnodes];
    size_t nodes_counter = problem->nnodes;
    for ( size_t i = 0; i < problem->nnodes; ++i ) {
        nodes[i] = i;
    }
    for ( size_t i = 0; i < k - 1; ++i ) {
        for ( size_t j = 0; j < nodes_counter; ++j ) {
            if ( nodes[j] == H[i] ) {
                --nodes_counter;
                nodes[j] = nodes[nodes_counter - 1];
                continue;
            }
        }
    }

    size_t x;

    /* Make insertions */
    size_t first_insertion = 0;
    for ( size_t iter = 0; iter < problem->nnodes - k + 1; ++iter ) {
        /* Select a random node */
        x           = rand_r(&__SEED) % nodes_counter;
        size_t node = nodes[x];
        nodes[x]    = nodes[nodes_counter - 2];
        --nodes_counter;

        /* Find the shortest insertion path */
        size_t best_edge = 0;
        double best_path_cost = _euclidean_distance(
                                    problem->xcoord[edges[0][0]],
                                    problem->ycoord[edges[0][0]],
                                    problem->xcoord[node],
                                    problem->ycoord[node]) +
                                _euclidean_distance(
                                    problem->xcoord[edges[0][1]],
                                    problem->ycoord[edges[0][1]],
                                    problem->xcoord[node],
                                    problem->ycoord[node]);

        for ( size_t i = 1; i < edges_counter; ++i ) {
            double cost = _euclidean_distance(
                              problem->xcoord[edges[i][0]],
                              problem->ycoord[edges[i][0]],
                              problem->xcoord[node], problem->ycoord[node]) +
                          _euclidean_distance(
                              problem->xcoord[edges[i][1]],
                              problem->ycoord[edges[i][1]],
                              problem->xcoord[node],
                              problem->ycoord[node]);

            if ( cost < best_path_cost ) {
                best_edge = i;
                best_path_cost = cost;
            }
        }

        /* Insert node */
        if ( first_insertion ) {
            edges[edges_counter][0]   = node;
            edges[edges_counter++][1] = edges[best_edge][0];
            edges[edges_counter][0]   = node;
            edges[edges_counter++][1] = edges[best_edge][1];
            first_insertion           = 0;
        } else {
            edges[edges_counter][0]   = node;
            edges[edges_counter++][1] = edges[best_edge][1];
            edges[best_edge][1]       = node;
        }

    }

    /* Compute solution cost */
    double currentcost = 0;
    for ( size_t i = 0; i < problem->nnodes; ++i ) {

        size_t node_1 = edges[i][0];
        size_t node_2 = edges[i][1];
        currentcost  += _euclidean_distance( problem->xcoord[node_1], problem->ycoord[node_1],
                                             problem->xcoord[node_2], problem->ycoord[node_2] );
    }

    if ( currentcost < *bestcost ) {
        log_info( "Heuristic solution improved (%.3e < %.3e).", currentcost, bestcost );

        for ( size_t i = 0; i < problem->nnodes; ++i ) {
            problem->solution[i][0] = edges[i][0];
            problem->solution[i][1] = edges[i][1];
        }

        *bestcost = currentcost;
    }


}


void
HeurConvHullInsertion_model ( instance *problem )
{
    struct timeb start, end;

    node_t *nodes = malloc( problem->nnodes * sizeof( *nodes ) );
    size_t *P     = malloc( problem->nnodes * sizeof( *P     ) );
    size_t *H     = malloc( problem->nnodes * sizeof( *H     ) );

    if ( nodes == NULL || P == NULL || H == NULL ) {
        log_fatal( "Out of memory." );
        exit( EXIT_FAILURE );
    }

    /* Sort the nodes and put their indices in P.  */
    ftime( &start );
    log_debug( "Sorting nodes by X-coordinates." );

    for ( size_t i = 0; i < problem->nnodes; ++i ) {
        nodes[i] = (node_t) {1, i, problem->xcoord[i], problem->ycoord[i] };
    }

    qsort( nodes, problem->nnodes, sizeof( *nodes ), _cmp_HeurConvHullInsertion );

    for ( size_t i = 0; i < problem->nnodes; ++i ) {
        P[i] = nodes[i].index;
    }

    ftime( &end );
    log_debug( "Done sorting in %.3lf seconds.",
               ( 1000. * ( end.time - start.time ) + end.millitm - start.millitm ) / 1000. );

    /* Compute the convex hull.  */
    size_t k = chainHull_2D( P, problem->nnodes, H, problem );

    double bestcost = __DBL_MAX__;
    __SEED = conf.seed;
    double elapsedtime = 0;

    for ( size_t j = 0; elapsedtime + 1e-3 < conf.heurtime; ++j ) {
        HeurConvHullInsertion_solve( problem, &bestcost, H, k );

        ftime( &end );
        elapsedtime = ( 1000. * ( end.time - start.time ) + end.millitm - start.millitm ) / 1000.;

        log_debug( "Found heuristic solution #%zu. Still %.3lf seconds remaining.",
                   j + 1, conf.heurtime - elapsedtime );
    }

    problem->elapsedtime  = elapsedtime;
    problem->solcost = compute_solution_cost( problem );

    free( nodes );
    free( P     );
    free( H     );
}
