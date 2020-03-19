/*!
 * \file    tsp.h
 * \brief   Type and function definitions for a TSP instance.
 * \authors Francesco Cazzaro, Marco Cieno
 */
#ifndef TSPPLOT_H
#define TSPPLOT_H

#include "tsp.h"


static char tspplotfile[] = "tspplot.dat";  /*!< Temporary file to store the GnuPlot. */


/**
 * \brief Plot a TSP instance.
 *
 * This function will generate a file in the relative path defined by \p tspplotfile.
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

#endif
