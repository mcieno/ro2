/*!
 * \file    tspplot.h
 * \brief   Type and function definitions for a TSP instance.
 * \authors Francesco Cazzaro, Marco Cieno
 */
#ifndef TSPPLOT_H
#define TSPPLOT_H

#include "tsp.h"


static char tspplot_tmpfile[] = "/tmp/tspplot.dat";  /*!< Temporary file to store information for GnuPlot. */


/**
 * \brief Write to file the information to plot an instance.
 *
 * \param problem
 *     Pointer to the instance structure.
 *
 * \param outputfile
 *     Name of the file where to save the plot data.
 */
void
instance_to_plot_dat ( instance *problem, char *outputfile );


/**
 * \brief Write to file the information to plot the solution of an instance.
 *
 * \param problem
 *     Pointer to the instance structure.
 *
 * \param outputfile
 *     Name of the file where to save the plot data.
 */
void
solution_to_plot_dat ( instance *problem, char *outputfile );


/**
 * \brief Plot a TSP instance.
 *
 * This function will generate a temporary file in the relative path defined by \p tspplot_tmpfile.
 * It will then spawn a gnuplot process and draw on it.
 *
 * \param problem
 *     Pointer to the instance structure.
 *
 * \note In order for this function to work properly, GnuPlot binary should be compiled with --enable-datastrings
 *       and available in the system `PATH`.
 */
void
plot_instance ( instance *problem );


/**
 * \brief Plot a solved TSP instance.
 *
 * This function will generate a temporary file in the relative path defined by \p tspplot_tmpfile.
 * It will then spawn a gnuplot process and draw on it.
 *
 * \param problem
 *     Pointer to the instance structure.
 *
 * \note In order for this function to work properly, GnuPlot binary should be compiled with --enable-datastrings
 *       and available in the system `PATH`.
 */
void
plot_solution ( instance *problem );


#endif
