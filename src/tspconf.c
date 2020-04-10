/*
 * \brief   Implementation of tspconf.h functions.
 * \authors Francesco Cazzaro, Marco Cieno
 */
#include <cplex.h>

#include "logging.h"
#include "tspconf.h"


/* Global configuration */
tspconf_t conf;

void
tspconf_init ( char      *filename,
               instance  *problem,
               int       shouldplot,
               model_t   solving_method,
               size_t    threads,
               size_t    memory,
               size_t    nodelimit,
               double    timelimit,
               double    epgap )
{
    conf.filename       = filename;
    conf.problem        = problem;
    conf.shouldplot     = shouldplot;
    conf.solving_method = solving_method;
    conf.threads        = threads;
    conf.memory         = memory;
    conf.timelimit      = nodelimit;
    conf.nodelimit      = timelimit;
    conf.epgap          = epgap;
}


void
tspconf_apply ( CPXENVptr env )
{
    if (loglevel >= LOG_INFO) {
        fprintf( stderr, CINFO "Applying configuration:\n"                                      );
        fprintf( stderr, CINFO "    Time limit          : %zu hours %zu minutes %zu seconds\n",
                                                              ((size_t) conf.timelimit) / 3600,
                                                         ((size_t) conf.timelimit) % 3600 / 60,
                                                                 ((size_t) conf.timelimit) % 60 );
        fprintf( stderr, CINFO "    Maximum memory      : %zu MB\n",                conf.memory );
        conf.threads == 0 ?
        fprintf( stderr, CINFO "    Use multithread     : yes\n"                               ):
        fprintf( stderr, CINFO "    Use multithread     : %s (%zu)\n",
                                                 conf.threads > 1 ? "yes" : "no", conf.threads );
    }
    if ( conf.timelimit > 0. ) CPXsetdblparam( env, CPX_PARAM_TILIM,   conf.timelimit );
    if (     conf.epgap > 0. ) CPXsetdblparam( env, CPX_PARAM_EPGAP,   conf.epgap     );
    if ( conf.nodelimit > 0  ) CPXsetintparam( env, CPX_PARAM_NODELIM, conf.nodelimit );
    if ( conf.threads   > 0  ) CPXsetintparam( env, CPXPARAM_Threads,  conf.threads   );
}
