/*
 * \brief   Implementation of tspplot.h functions.
 * \authors Francesco Cazzaro, Marco Cieno
 */
#include <errno.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>

#include "logging.h"
#include "tsp.h"
#include "tspplot.h"


void
instance_to_plot_dat ( instance *problem, char *outputfile )
{
    FILE * fd = fopen( outputfile, "w" );

    for ( unsigned long i = 0; i < problem->nnodes; ++i ) {
        fprintf( fd, "%lu %lf %lf \n", i + 1, problem->xcoord[i], problem->ycoord[i] );
    }

    fclose( fd );
}


void
solution_to_plot_dat ( instance *problem, char *outputfile )
{
    FILE * fd = fopen( outputfile, "w" );

    for ( unsigned long i = 0; i < problem->nnodes; ++i ) {
        fprintf(
            fd, "%lu %lf %lf \n", i + 1,
            problem->xcoord[problem->solution[i]],
            problem->ycoord[problem->solution[i]]
        );
    }

    // Repeat starting to close the tour
    fprintf(
        fd, "%lu %lf %lf \n", 0,
        problem->xcoord[problem->solution[0]],
        problem->ycoord[problem->solution[0]]
    );

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

    fprintf( gnuplot, "set title '%s' \n", problem->name );
    fprintf( gnuplot, "plot '%s' u 2:3 pt 2 w lines notitle,"
                      " '' using 2:3:1 w labels offset 0.7,0.7 notitle \n", tspplot_tmpfile );

    fflush( gnuplot );
    fclose( gnuplot );

    remove( tspplot_tmpfile );
}
