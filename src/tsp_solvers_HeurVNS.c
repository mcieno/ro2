/*
 * \brief   VNS heuristic method.
 * \authors Francesco Cazzaro, Marco Cieno
 */

#include <errno.h>
#include <limits.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

#include "tsp_solvers.h"
#include "logging.h"
#include "tsp.h"
#include "tspconf.h"


static unsigned int __SEED;

typedef struct
{
    int available;
    int v;
    int u;
    double cost;
}
edge_t;


/**
 * Comparator for the edge_t type.
 *
 * @param a void pointer to an `edge_t` element.
 * @param b void pointer to an `edge_t` element.
 * @returns negative if `cost(a) < cost(b)`, positive if `cost(a) > cost(b)`,
 *          0 if `cost(a) == cost(b)`.
 */
int
_cmp_VNS( const void *a, const void *b ) {
    const edge_t *ea = (const edge_t *)a;
    const edge_t *eb = (const edge_t *)b;

    return ea->cost < eb->cost
               ? -1
               : ea->cost == eb->cost
                     ? 0
                     : +1;
}

/**
 * Invert symmetrically all elements in the array from index a to index b included.
 *
 * @param ordered_nodes size_t pointer to array to be inverted.
 * @param a index of node a.
 * @param b index of node b.
 * @param problem instance pointer to problem instance.
 */
void
invert_array(size_t *ordered_nodes, size_t a, size_t b, instance *problem)
{
    if ( a == b ) {
        return;
    }
    size_t iter;
    if ( a>b ) {
        iter = ((b+problem->nnodes)-a+1)/2;
    }else{ iter = (b-a+1)/2;}

    for ( int i = 0; i<iter; ++i ) {
        if ( a>b ) {
            b = b+problem->nnodes;
        }
        size_t temp = ordered_nodes[(a+i)%problem->nnodes];
        ordered_nodes[(a+i)%problem->nnodes] = ordered_nodes[(b-i)%problem->nnodes];
        ordered_nodes[(b-i)%problem->nnodes]= temp;
    }

}

void
swap_blocks(size_t *ordered_nodes, size_t a, size_t b, size_t c, size_t d, instance *problem){

    if ( a == c && b == d ) {
        return;
    }

    if ( b == c ) {
        return;
    }

    size_t *supp_array = calloc( problem->nnodes, sizeof( *supp_array ) );
    for ( int i = 0; i<problem->nnodes; ++i ) {
        supp_array[i] = ordered_nodes[i];
    }

    size_t iter = 0;
    if ( c>d ) {
        iter = (d+problem->nnodes)-c+1;
    }else{
        iter = d-c+1;
    }

    size_t index = a;
    for ( size_t i = 0; i<iter; ++i ) {
        ordered_nodes[index%problem->nnodes] = supp_array[(c+i)%problem->nnodes];
        index++;
    }

    if ( b>c ) {
        iter = (c+problem->nnodes)-b-1;
    }else{
        iter = c-b-1;
    }

    for ( size_t i = 0; i<iter; ++i ) {
        ordered_nodes[index%problem->nnodes] = supp_array[(b+i+1)%problem->nnodes];
        index++;
    }


    if ( a>b ) {
        iter = (b+problem->nnodes)-a+1;
    }else{
        iter = b-a+1;
    }


    for ( size_t i = 0; i<iter; ++i ) {
        ordered_nodes[index%problem->nnodes] = supp_array[(a+i)%problem->nnodes];
        index++;
    }



    free(supp_array);
}

/**
 * \brief Apply a random 5-OPT MOVE to diversificate the provided solution.
 *
 *
 * \param currentsol The solution to be diversificated. Will be diversificated in place.
 * \param problem Pointer to the problem instance object.
 *
 * \warning The value of \p currentsol must be a feasible solution.
 */
void _5opt_diversificate_VNS(size_t **currentsol, instance *problem)
{
    int num_opt = 5;

    size_t replaced_m_edges = 0;
    size_t *ind_m_edges = calloc(num_opt, sizeof(*ind_m_edges));
    size_t *nodes_m_edges = calloc(2 * num_opt, sizeof(*nodes_m_edges));
    size_t *phantom_edges = calloc(2 * num_opt, sizeof(*phantom_edges));
    if ( ind_m_edges == NULL )
    {
        log_fatal("Out of memory.");
        exit(EXIT_FAILURE);
    }

    for ( size_t i = 0; i < num_opt; ++i )
    {
        int found = 0;
        while ( found == 0 )
        {
            int unique = 1;
            size_t e = rand_r(&__SEED) % problem->nnodes;
            for ( size_t j = 0; j < i; ++j )
            {
                if ( ind_m_edges[j] == e )
                {
                    unique = 0;
                }
            }
            if ( unique == 1 )
            {
                found = 1;
                ind_m_edges[i] = e;
                nodes_m_edges[2 * i] = currentsol[e][0];
                nodes_m_edges[2 * i + 1] = currentsol[e][1];
                phantom_edges[2 * i] = currentsol[e][0];
                phantom_edges[2 * i + 1] = currentsol[e][1];
            }
        }
    }

    size_t *ordered_nodes = calloc(problem->nnodes, sizeof(*ordered_nodes));
    /* Build sequence of nodes of the solution to do so in the process, we
     * convert `currentsol` to `next/prev` representation */
    size_t *next = malloc(problem->nnodes * sizeof(*next));
    size_t *prev = malloc(problem->nnodes * sizeof(*prev));

    if ( next == NULL || prev == NULL )
    {
        log_fatal("Out of memory.");
        exit(EXIT_FAILURE);
    }

    size_t u = currentsol[0][0];
    size_t v = currentsol[0][1];
    for ( size_t k = 0; k < problem->nnodes; ++k )
    {
        next[u] = v;
        prev[v] = u;

        for ( size_t kk = 0; kk < problem->nnodes; ++kk )
        {
            if ( (currentsol[kk][0] == v && currentsol[kk][1] != u ) || (currentsol[kk][1] == v && currentsol[kk][0] != u))
            {
                /* kk-th edge in the solution is the next edge going from v to some other node. */
                u = v;
                v ^= currentsol[kk][0] ^ currentsol[kk][1];
                break;
            }
        }
    }

    size_t pointer = 0;
    for ( int i = 0; i < problem->nnodes; ++i )
    {
        ordered_nodes[i] = next[pointer];
        pointer = next[pointer];
    }

    size_t length_nodes_m_edges = 2 * num_opt;
    while ( length_nodes_m_edges > 2 )
    {

        /* Pick two random nodes from m-opt list */
        size_t ind_a = rand_r(&__SEED) % length_nodes_m_edges;
        size_t node_a = nodes_m_edges[ind_a];
        nodes_m_edges[ind_a] = nodes_m_edges[length_nodes_m_edges - 1];
        nodes_m_edges[length_nodes_m_edges - 1] = node_a;
        length_nodes_m_edges--;
        size_t ind_b = rand_r(&__SEED) % length_nodes_m_edges;
        size_t node_b = nodes_m_edges[ind_b];
        nodes_m_edges[ind_b] = nodes_m_edges[length_nodes_m_edges - 1];
        nodes_m_edges[length_nodes_m_edges - 1] = node_b;
        length_nodes_m_edges--;

        if ( node_a == node_b )
        {
            length_nodes_m_edges += 2;
            continue;
        }

        int new_edge_status = -1;
        for ( int i = 0; i < problem->nnodes; ++i )
        {
            if ( ordered_nodes[i] == node_a )
            {
                ind_a = i;
            }
            if ( ordered_nodes[i] == node_b )
            {
                ind_b = i;
            }
        }

        if ( ordered_nodes[(ind_a + 1 ) % problem->nnodes] == node_b && length_nodes_m_edges > 2)
        {
            length_nodes_m_edges += 2;

            continue;
        }
        if ( ordered_nodes[(ind_b + 1 ) % problem->nnodes] == node_a && length_nodes_m_edges > 2)
        {
            length_nodes_m_edges += 2;

            continue;
        }

        if ( length_nodes_m_edges == 2 && nodes_m_edges[0] == nodes_m_edges[1] )
        {

            length_nodes_m_edges += 2;
            continue;
        }

        if ( (ind_a + 1 ) % problem->nnodes == ind_b)
        {
            size_t ind_temp = ind_b;
            size_t node_temp = node_b;
            ind_b = ind_a;
            node_b = node_a;
            ind_a = ind_temp;
            node_a = node_temp;
        }

        size_t temp_index = ind_a + 1;
        while ( new_edge_status == -1 )
        {

            for ( int i = 0; i < length_nodes_m_edges; ++i )
            {
                if ( ordered_nodes[temp_index % problem->nnodes] == nodes_m_edges[i] )
                {
                    new_edge_status = 1;
                    break;
                }
            }
            if ( ordered_nodes[temp_index % problem->nnodes] == node_b )
            {
                /* If find node_b before other nodes of nodes_m_edges then we
                 * cannot connect, except in this same iteration the node_b
                 * was found. This means both edges of node_b have been
                 * unattached so we can connect.  */
                if ( new_edge_status != 1 )
                {
                    new_edge_status = 0;
                }

                break;
            }

            temp_index++;
        }

        if ( new_edge_status == 1 && !(length_nodes_m_edges == 2 && nodes_m_edges[0] == nodes_m_edges[1] ) )
        {

            new_edge_status = -1;
            temp_index = problem->nnodes + ind_a - 1;
            while ( new_edge_status == -1 )
            {

                for ( int i = 0; i < length_nodes_m_edges; ++i )
                {
                    if ( ordered_nodes[temp_index % problem->nnodes] == nodes_m_edges[i] )
                    {
                        new_edge_status = 1;
                        break;
                    }
                }
                if ( ordered_nodes[temp_index % problem->nnodes] == node_b )
                {
                    /* If find node_b before other nodes of nodes_m_edges then we
                     * cannot connect, except in this same iteration the node_b
                     * was found. This means both edges of node_b have been
                     * unattached so we can connect.  */
                    if ( new_edge_status != 1 )
                    {
                        new_edge_status = 0;
                    }

                    break;
                }

                temp_index--;
            }
        }

        if ( new_edge_status == 0 )
        {
            length_nodes_m_edges += 2;
        }
        else
        {

            if ( node_a < node_b )
            {
                currentsol[ind_m_edges[replaced_m_edges]][0] = node_a;
                currentsol[ind_m_edges[replaced_m_edges]][1] = node_b;
            }
            else
            {
                currentsol[ind_m_edges[replaced_m_edges]][0] = node_b;
                currentsol[ind_m_edges[replaced_m_edges]][1] = node_a;
            }

            if ( length_nodes_m_edges == 2 )
            {
                replaced_m_edges++;
                break;
            }

            /* reorder ordered_nodes appropriately first find index and edges
             * of phantom_edges involved */
            size_t a_y = 0, b_y = 0, edge_ind_a = 0, edge_ind_b = 0;
            for ( int i = replaced_m_edges; i < num_opt; ++i )
            {
                if ( phantom_edges[2 * i] == node_a )
                {
                    edge_ind_a = i;
                    a_y = phantom_edges[2 * i + 1];
                }
                else if ( phantom_edges[2 * i + 1] == node_a )
                {
                    edge_ind_a = i;
                    a_y = phantom_edges[2 * i];
                }

                if ( phantom_edges[2 * i] == node_b )
                {
                    edge_ind_b = i;
                    b_y = phantom_edges[2 * i + 1];
                }
                else if ( phantom_edges[2 * i + 1] == node_b )
                {
                    edge_ind_b = i;
                    b_y = phantom_edges[2 * i];
                }
            }

            for ( int i = 0; i < problem->nnodes; ++i )
            {
            }
            size_t candidate_a_y = ordered_nodes[(ind_a + 1) % problem->nnodes];
            size_t candidate_b_y = ordered_nodes[(ind_b + 1) % problem->nnodes];

            for ( size_t i = 0; i < problem->nnodes; ++i )
            {
                int skip = 0;
                for ( int j = replaced_m_edges; j < num_opt; ++j )
                {
                    if ( i == ind_m_edges[j] )
                    {

                        skip = 1;
                    }
                }

                if ( skip == 1 )
                {

                    continue;
                }

                if ( currentsol[i][0] == node_a && currentsol[i][1] == candidate_a_y )
                {
                    candidate_a_y = ordered_nodes[(ind_a - 1 + problem->nnodes) % problem->nnodes];
                }
                if ( currentsol[i][1] == node_a && currentsol[i][0] == candidate_a_y )
                {
                    candidate_a_y = ordered_nodes[(ind_a - 1 + problem->nnodes) % problem->nnodes];
                }
                if ( currentsol[i][0] == node_b && currentsol[i][1] == candidate_b_y )
                {
                    candidate_b_y = ordered_nodes[(ind_b - 1 + problem->nnodes) % problem->nnodes];
                }
                if ( currentsol[i][1] == node_b && currentsol[i][0] == candidate_b_y )
                {
                    candidate_b_y = ordered_nodes[(ind_b - 1 + problem->nnodes) % problem->nnodes];
                }
            }

            a_y = candidate_a_y;
            b_y = candidate_b_y;

            length_nodes_m_edges += 2;

            size_t blk_a_start = 0, blk_a_end = 0, blk_b_start = 0, blk_b_end = 0, blk_next_start = 0, blk_next_end = 0;
            int is_start_iteration = 0;
            int status = 0;
            int b_status = 0;
            int is_a_last = 0;
            size_t last = 0;
            int deattached_a_count = 0;

            /* Check if a is a node with two deattached edges, special case */
            for ( int j = 0; j < length_nodes_m_edges; ++j )
            {

                if ( ordered_nodes[ind_a % problem->nnodes] == nodes_m_edges[j] )
                {
                    deattached_a_count++;
                }
                if ( deattached_a_count == 2 )
                {
                    blk_a_start = ind_a;
                    blk_a_end = ind_a;
                    status = 1;
                    is_start_iteration = 1;
                }
            }
            for ( int i = ind_a + 1; i < ind_a + problem->nnodes; ++i )
            {
                for ( int j = 0; j < length_nodes_m_edges; ++j )
                {
                    if ( ordered_nodes[i % problem->nnodes] == nodes_m_edges[j] )
                    {

                        /* Determine if ind_a is start or end */
                        if ( status == 0 )
                        {
                            if ( ordered_nodes[i % problem->nnodes] != a_y )
                            {
                                blk_a_start = ind_a;
                                blk_a_end = i % problem->nnodes;
                                status = 10;
                                is_start_iteration = 0;
                            }
                            else
                            {
                                blk_a_end = ind_a;
                                is_a_last = 1;
                                blk_next_start = i % problem->nnodes;
                                status = 2;
                                is_start_iteration = 1;
                            }
                        }
                        else if ( status == 1 )
                        {
                            blk_next_start = i % problem->nnodes;
                            status = 2;
                        }
                        else if ( status == 2 )
                        {
                            blk_next_end = i % problem->nnodes;
                            status = 3;
                        }

                        if ( ordered_nodes[i % problem->nnodes] == node_b && b_status != 4 )
                        {
                            if ( is_start_iteration )
                            {
                                blk_b_start = i % problem->nnodes;

                                b_status = 4;
                                blk_b_end = last;
                            }
                            else
                            {

                                blk_b_end = i % problem->nnodes;
                                blk_b_start = last;
                                if ( ordered_nodes[last] == b_y )
                                {
                                    blk_b_start = blk_b_end;
                                    b_status = 4;
                                }
                            }
                        }
                        else if ( b_status == 4 )
                        {
                            if ( ordered_nodes[i % problem->nnodes] == b_y )
                            {
                                size_t tmp = blk_b_end;
                                blk_b_end = blk_b_start;
                                blk_b_start = tmp;
                                b_status = 3;
                            }
                            else
                            {
                                blk_b_end = i % problem->nnodes;
                                b_status = 3;
                            }
                        }

                        last = i % problem->nnodes;

                        if ( is_start_iteration )
                        {
                            is_start_iteration = 0;
                        }
                        else
                        {
                            is_start_iteration = 1;
                        }
                    }
                }
                if ( status == 10 )
                {
                    status = 1;
                }
            }
            if ( is_a_last )
            {
                blk_a_start = last;
            }

            length_nodes_m_edges -= 2;

            if ( !is_a_last )
            {
                invert_array(ordered_nodes, blk_a_start, blk_a_end, problem);
            }
            if ( ordered_nodes[blk_b_end] == node_b )
            {
                invert_array(ordered_nodes, blk_b_start, blk_b_end, problem);
            }

            swap_blocks(ordered_nodes, blk_next_start, blk_next_end, blk_b_start, blk_b_end, problem);

            /* Update phantom_edges */
            phantom_edges[2 * edge_ind_b] = a_y + edge_ind_a;
            phantom_edges[2 * edge_ind_b + 1] = b_y;
            phantom_edges[2 * edge_ind_a] = phantom_edges[2 * replaced_m_edges];
            phantom_edges[2 * edge_ind_a + 1] = phantom_edges[2 * replaced_m_edges + 1];

            replaced_m_edges++;
        }
    }

    /* Connect the last remaining edge */
    if ( nodes_m_edges[0] < nodes_m_edges[1] )
    {
        currentsol[ind_m_edges[replaced_m_edges]][0] = nodes_m_edges[0];
        currentsol[ind_m_edges[replaced_m_edges]][1] = nodes_m_edges[1];
    }
    else
    {
        currentsol[ind_m_edges[replaced_m_edges]][0] = nodes_m_edges[1];
        currentsol[ind_m_edges[replaced_m_edges]][1] = nodes_m_edges[0];
    }

    free(next);
    free(prev);

    free(ind_m_edges);
    free(ordered_nodes);
    free(nodes_m_edges);
}

/**
 * \brief Repeatedly apply 2-OPT MOVE for refininig the provided solution.
 *
 *
 * \param currentsol The solution to be improved. Will be modified in place.
 * \param problem Pointer to the problem instance object.
 *
 * \warning The value of \p currentsol must be a feasible solution.
 */
void
_2opt_refine_VNS( size_t **currentsol, instance *problem )
{
    /* Convert `currentsol` to `next/prev` representation */
    size_t *next = malloc( problem->nnodes * sizeof( *next ) );
    size_t *prev = malloc( problem->nnodes * sizeof( *prev ) );

    if ( next == NULL || prev == NULL ) {
        log_fatal( "Out of memory." );
        exit( EXIT_FAILURE );
    }

    size_t u = currentsol[0][0];
    size_t v = currentsol[0][1];
    for ( size_t k = 0; k < problem->nnodes; ++k ) {
        next[u] = v;
        prev[v] = u;

        for ( size_t kk = 0; kk < problem->nnodes; ++kk ) {
            if ( (currentsol[kk][0] == v && currentsol[kk][1] != u )
                 || (currentsol[kk][1] == v && currentsol[kk][0] != u) ) {
                /* kk-th edge in the solution is the next edge going from v to some other node. */
                u = v;
                v ^= currentsol[kk][0] ^ currentsol[kk][1];
                break;
            }
        }
    }


    /* Weights of the edges under scrutiny.  */
    double wuv, wu_v_, wuu_, wvv_;

    /* Nodes under scrutiny. */
    size_t u_, v_;

    /* Boolean flag exiting refinement loop in case local optimality was reached.  */
    int wasrefined;

    do {
        wasrefined = 0;

        for ( u = 0, v = next[u]; v != 0; u = v, v = next[u] )
        {
            wuv = _euclidean_distance( problem->xcoord[u], problem->ycoord[u], problem->xcoord[v], problem->ycoord[v] );

            for ( u_ = next[v], v_ = next[u_]; v_ != u; u_ = v_, v_ = next[u_] )
            {
                wu_v_ = _euclidean_distance( problem->xcoord[u_], problem->ycoord[u_], problem->xcoord[v_], problem->ycoord[v_] );

                /* Calculate cost of 2-OPT move, paying attention to not disconnecting the graph */
                wuu_ = _euclidean_distance( problem->xcoord[u], problem->ycoord[u], problem->xcoord[u_], problem->ycoord[u_] );
                wvv_ = _euclidean_distance( problem->xcoord[v], problem->ycoord[v], problem->xcoord[v_], problem->ycoord[v_] );

                if ( wuu_ + wvv_ < wuv + wu_v_ )
                {
                    log_trace( "2-OPT MOVE: (%zu, %zu) X (%zu, %zu)", u, v, u_, v_ );
                    /* Swap all the rout from v to u_ */
                    for ( size_t vv = v, t = next[next[vv]]; vv != u_; vv = prev[t], t = next[t] ) {
                        prev[vv] = prev[t];
                        next[prev[t]] = vv;
                    }

                    /* Swap edges if any improvement.  */
                    next[u] = u_;
                    next[v] = v_;
                    prev[u_] = u;
                    prev[v_] = v;

                    wasrefined = 1;
                    break;
                }
            }
        }

    } while ( wasrefined ) ;

    /* Update solution */
    size_t k;
    for ( k = 0, u = 0, v = next[u]; v != 0; ++k, u = v, v = next[u] ) {
        currentsol[k][0] = u;
        currentsol[k][1] = v;
    }
    currentsol[k][0] = u;
    currentsol[k][1] = v;

    free(next);
    free(prev);
}


void
HeurVNS_solve ( instance *problem )
{
    struct timespec start, end;

    size_t nedges = ( problem->nnodes * ( problem->nnodes + 1 ) ) / 2 - problem->nnodes;

    /* Build a list with all edges and sort them.  */
    edge_t *edges = malloc( nedges * sizeof( *edges ) );
    if ( edges == NULL ) {
        log_fatal( "Out of memory." );
        exit( EXIT_FAILURE );
    }

    size_t pos = 0;
    for ( size_t i = 0; i < problem->nnodes; ++i ) {
        for ( size_t j = i + 1; j < problem->nnodes; ++j ) {
            edges[pos++] = (edge_t) {
                1, i, j,
                _euclidean_distance(
                    problem->xcoord[i],
                    problem->ycoord[i],
                    problem->xcoord[j],
                    problem->ycoord[j]
                )
            };
        }
    }

    log_debug( "Sorting edges by cost." );
    clock_gettime( CLOCK_MONOTONIC, &start );

    qsort( edges, nedges, sizeof( *edges ), _cmp_VNS );

    clock_gettime( CLOCK_MONOTONIC, &end );
    log_debug( "Done sorting in %.3lf seconds.",
               ( end.tv_sec - start.tv_sec ) + ( end.tv_nsec - start.tv_nsec ) / 1000000000. );

    size_t from;

    /* `currentsol` will contain the solution obtained starting the nearest
     * neighbor heuristic from `startnode`.  */
    double bestcost = __DBL_MAX__;
    size_t **bestsol = problem->solution;
    size_t **currentsol = malloc( problem->nnodes * sizeof( *currentsol ) );
    if ( currentsol == NULL ) {
        log_fatal( "Out of memory." );
        exit( EXIT_FAILURE );
    }
    for ( size_t i = 0; i < problem->nnodes; ++i ) {
        currentsol[i] = malloc( 2 * sizeof( *currentsol[i] ) );
        if ( currentsol[i] == NULL ) {
            log_fatal( "Out of memory." );
            exit( EXIT_FAILURE );
        }
    }


    double elapsedtime = 0;

    clock_gettime( CLOCK_MONOTONIC, &start );

    /* Initialize our solution starting from a solution computed with the neirest neighborood method  */
    size_t startnode = rand_r(&__SEED) % problem->nnodes;
    from = startnode;
    for ( size_t k = 0; k < problem->nnodes - 1; ++k )
    {
        /* Search the almost-shortest edge where `from` occurs.  */
        for ( pos = 0; ; pos = (pos + 1 ) % nedges ) {
            /* Flip a coin and decide whether to stop here or keep going.
             * If we run out of edges prior to successfull coinflip,
             * simply restart the loop. Note that this method terminates
             * with probability 1.  */
            if ( edges[pos].available && ( edges[pos].v == from || edges[pos].u == from ) )
            {
                if ( rand_r( &__SEED ) >= RAND_MAX / 4 ) {
                    currentsol[k][0] = edges[pos].v;
                    currentsol[k][1] = edges[pos].u;
                    break;
                }
            }
        }

        log_trace("Greedy choice #%zu was %zu - %zu", k, currentsol[k][0], currentsol[k][1]);

        /* Avoid subtours removing all other edges where `from` occurs.  */
        for ( pos = 0; pos < nedges; ++pos ) {
            if ( edges[pos].available && ( edges[pos].v == from || edges[pos].u == from ) )
            {
                edges[pos].available = 0;
            }
        }

        /* Update `from` to start from the just added neighbor.  */
        from ^= currentsol[k][0] ^ currentsol[k][1];
    }

    /* Set last edge.  */
    currentsol[problem->nnodes - 1][0] = from;
    currentsol[problem->nnodes - 1][1] = startnode;

    for ( int i = 0; i<problem->nnodes; ++i ) {
        problem->solution[i][0] = currentsol[i][0];
        problem->solution[i][1] = currentsol[i][1];
    }

    /* Run iteratively 2-opt and a random 5-opt jump */
    for ( size_t iter = 0; elapsedtime + 1e-3 < conf.heurtime; iter++ )
    {
        log_debug( "Applying 5-OPT refinement" );
        _5opt_diversificate_VNS(currentsol, problem);

        log_debug( "Applying 2-OPT refinement" );
        _2opt_refine_VNS( currentsol, problem );

        clock_gettime( CLOCK_MONOTONIC, &end );

        elapsedtime = ( end.tv_sec - start.tv_sec ) + ( end.tv_nsec - start.tv_nsec ) / 1000000000.;

        /* Compute solution cost and update `bestsol` and `bestcost`
         * accordingly.  */
        problem->solution = currentsol;
        problem->solcost = compute_solution_cost( problem );

        if ( problem->solcost < bestcost ) {
            /* Swap `currentsol` and `bestsol`, so we can reuse the arrays
             * during the next itaration.  */
            log_info( "Heuristic solution improved at #%d (%.3e < %.3e).", startnode + 1, problem->solcost, bestcost );
            currentsol = bestsol;
            bestsol    = problem->solution;
            bestcost   = problem->solcost;
        }
    }

    problem->solution = bestsol;
    problem->solcost  = bestcost;

    /* Free `currentsol`, which may have been swapped in the mean time.  */
    for ( size_t i = 0; i < problem->nnodes; ++i ) free( currentsol[i] );
    free( currentsol );
}


void
HeurVNS_model ( instance *problem )
{
    struct timespec start, end;
    clock_gettime( CLOCK_MONOTONIC, &start );
    __SEED = conf.seed;

    log_debug( "Starting solver." );
    HeurVNS_solve( problem );

    clock_gettime( CLOCK_MONOTONIC, &end );

    problem->elapsedtime  = ( end.tv_sec - start.tv_sec ) + ( end.tv_nsec - start.tv_nsec ) / 1000000000.;
    problem->visitednodes = 0;
}