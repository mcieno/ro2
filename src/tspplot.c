/*
 * \brief   Implementation of tspplot.h functions.
 * \authors Francesco Cazzaro, Marco Cieno
 */
#include <errno.h>
#include <limits.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "logging.h"
#include "tsp.h"
#include "tspplot.h"


char *tspplot_tmpfile = "/tmp/tspplot.dat";


void
instance_to_plot_dat ( instance *problem, char *outputfile )
{
    FILE * fd = fopen( outputfile, "w" );

    for ( size_t i = 0; i < problem->nnodes; ++i ) {
        fprintf( fd, "%zu %lf %lf \n", i + 1, problem->xcoord[i], problem->ycoord[i] );
    }

    fclose( fd );
}


void
solution_to_plot_dat ( instance *problem, char *outputfile )
{
    FILE * fd = fopen( outputfile, "w" );


    for ( size_t k = 0; k < problem->nnodes; ++k )
    {
        long unsigned i = problem->solution[k][0];
        long unsigned j = problem->solution[k][1];

        // idx_from idx_to x_from y_from x_to y_to
        fprintf(
            fd, "%zu %zu %lf %lf %lf %lf \n", i + 1, j + 1,
            problem->xcoord[i],
            problem->ycoord[i],
            problem->xcoord[j],
            problem->ycoord[j]
        );
    }

    fclose( fd );
}


void
plot_instance ( instance *problem )
{
    instance_to_plot_dat( problem, tspplot_tmpfile );

    FILE * gnuplot = popen( "gnuplot -persistent", "w" );

    fprintf( gnuplot, "set title '%s' \n", problem->name );
    fprintf( gnuplot, "plot '%s' u 2:3 pt 2 notitle,"
                      " '' u 2:3:1 w labels offset 0.7,0.7 notitle \n", tspplot_tmpfile );

    fflush( gnuplot );
    fclose( gnuplot );

    remove( tspplot_tmpfile );
}


void
plot_solution ( instance *problem )
{
    solution_to_plot_dat( problem, tspplot_tmpfile );

    FILE * gnuplot = popen( "gnuplot -persistent", "w" );
    fprintf( gnuplot, "set style arrow 1 nohead \n" );

    fprintf( gnuplot, "set title '%s' \n", problem->name );
    fprintf( gnuplot, "plot '%s' u 3:4:($5-$3):($6-$4) w vectors arrowstyle 1 notitle, "
                      "'' u 3:4:1 w labels offset 0.7,0.7 notitle, "
                      "'' u 5:6:2 w labels offset 0.7,0.7 notitle \n", tspplot_tmpfile );

    fflush( gnuplot );
    fclose( gnuplot );

    remove( tspplot_tmpfile );
}
