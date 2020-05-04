/*
 * \brief   Implementation of tspconf.h functions.
 * \authors Francesco Cazzaro, Marco Cieno
 */
#include <cplex.h>

#include "logging.h"
#include "tspconf.h"
#include <sys/sysinfo.h>


/* Global configuration */
tspconf_t conf;

void
tspconf_init ()
{
    conf.filename       = NULL;
    conf.name           = NULL;
    conf.shouldplot     = 1;
    conf.solving_method = TSP_SOLVER_Generic;
    conf.threads        = 0;
    conf.memory         = 0;
    conf.timelimit      = 0.;
    conf.heurtime       = 600.;
    conf.nodelimit      = 0;
    conf.cutup          = 0.;
    conf.epgap          = 0.;
    conf.scrind         = CPX_OFF;
    conf.seed           = 0;

    log_set_level( loglevel );
}


void
tspconf_apply ( CPXENVptr env )
{
    log_info( "Applying configuration:"                                              );

    if ( conf.timelimit <= 0 )
    log_info( "    Time limit                 : No limit"                            );
    else
    log_info( "    Time limit                 : %zu hours %zu minutes %zu seconds",
                                                   ((size_t) conf.timelimit) / 3600,
                                              ((size_t) conf.timelimit) % 3600 / 60,
                                                      ((size_t) conf.timelimit) % 60 );

    if ( conf.nodelimit == 0 )
    log_info( "    MIP node limit             : No limit"                            );
    else
    log_info( "    MIP node limit             : %zu",                 conf.nodelimit );

    if ( conf.memory == 0 )
    log_info( "    Maximum working memory     : Automatic"                           );
    else
    log_info( "    Maximum working memory     : %zu MB",                 conf.memory );

    if ( conf.threads == 0 )
    log_info( "    Use multithread            : yes (ALL:%d)", get_nprocs()          );
    else
    log_info( "    Use multithread            : %s (%zu)",
                                       conf.threads > 1 ? "yes" : "no", conf.threads );

    if ( conf.cutup == 0 )
    log_info( "    Upper cutoff               : Don't cut"                           );
    else
    log_info( "    Upper cutoff               : %lf",                     conf.cutup );

    if ( conf.epgap <= 0 )
    log_info( "    Relative MIP gap tolerance : Automatic"                           );
    else
    log_info( "    Relative MIP gap tolerance : %lf",                     conf.epgap );

    if ( conf.scrind == CPX_OFF )
    log_info( "    Display messages on screen : No"                                  );
    else
    log_info( "    Display messages on screen : Yes"                                 );

    if ( conf.seed == 0 )
    log_info( "    Random seed                : Automatic"                           );
    else
    log_info( "    Random seed                : %i",                       conf.seed );

    if ( conf.timelimit > 0. ) CPXsetdblparam( env, CPXPARAM_TimeLimit,                  conf.timelimit );
    if ( conf.nodelimit > 0  ) CPXsetintparam( env, CPXPARAM_MIP_Limits_Nodes,           conf.nodelimit );
    if (    conf.memory > 0  ) CPXsetdblparam( env, CPXPARAM_WorkMem,                    conf.memory    );
    if (   conf.threads > 0  ) CPXsetintparam( env, CPXPARAM_Threads,                    conf.threads   );
    else                       CPXsetintparam( env, CPXPARAM_Threads,                    get_nprocs()   );
    if (     conf.cutup > 0. ) CPXsetdblparam( env, CPXPARAM_MIP_Tolerances_UpperCutoff, conf.cutup     );
    if (     conf.epgap > 0. ) CPXsetdblparam( env, CPXPARAM_MIP_Tolerances_MIPGap,      conf.epgap     );
    if (    conf.scrind > 0  ) CPXsetintparam( env, CPXPARAM_ScreenOutput,               conf.scrind    );
    if (      conf.seed > 0  ) CPXsetintparam( env, CPXPARAM_RandomSeed,                 conf.seed      );
}


void
tspconf_destroy()
{
    if ( conf.filename != NULL)  free( conf.filename );
    if ( conf.name     != NULL)  free( conf.name     );
    tspconf_init();
}
