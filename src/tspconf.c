#include <cplex.h>

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
    conf.filename       =  filename;
    conf.problem        =  problem;
    conf.shouldplot     =  shouldplot;
    conf.solving_method =  solving_method;
    conf.threads        =  threads;
    conf.memory         =  memory;
    conf.timelimit      =  nodelimit;
    conf.nodelimit      =  timelimit;
    conf.epgap          = epgap;
}


void
tspconf_apply ( CPXENVptr env )
{
    if ( conf.timelimit > 0. ) CPXsetdblparam( env, CPX_PARAM_TILIM,   conf.timelimit );
    if (     conf.epgap > 0. ) CPXsetdblparam( env, CPX_PARAM_EPGAP,   conf.epgap     );
    if ( conf.nodelimit > 0  ) CPXsetintparam( env, CPX_PARAM_NODELIM, conf.nodelimit );
}
