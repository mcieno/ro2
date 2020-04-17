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
    conf.solving_method = TSP_SOLVER_LOOPBB;
    conf.threads        = 0;
    conf.memory         = 0;
    conf.timelimit      = 0.;
    conf.nodelimit      = 0;
    conf.cutup          = 0.;
    conf.epgap          = 0.;
    conf.scrind         = CPX_OFF;
    conf.seed           = 0;
}


void
tspconf_apply ( CPXENVptr env )
{
    if (loglevel >= LOG_INFO) {
        fprintf( stderr, CINFO "Applying configuration:\n"                                             );

        conf.timelimit <= 0 ?
        fprintf( stderr, CINFO "    Time limit                 : No limit\n"                           ):
        fprintf( stderr, CINFO "    Time limit                 : %zu hours %zu minutes %zu seconds\n",
                                                                     ((size_t) conf.timelimit) / 3600,
                                                                ((size_t) conf.timelimit) % 3600 / 60,
                                                                        ((size_t) conf.timelimit) % 60 );

        conf.nodelimit == 0 ?
        fprintf( stderr, CINFO "    MIP node limit             : No limit\n"                           ):
        fprintf( stderr, CINFO "    MIP node limit             : %zu\n",                conf.nodelimit );

        conf.memory == 0 ?
        fprintf( stderr, CINFO "    Maximum working memory     : Automatic\n"                          ):
        fprintf( stderr, CINFO "    Maximum working memory     : %zu MB\n",                conf.memory );

        conf.threads == 0 ?
        fprintf( stderr, CINFO "    Use multithread            : yes (ALL)\n"                          ):
        fprintf( stderr, CINFO "    Use multithread            : %s (%zu)\n",
                                                         conf.threads > 1 ? "yes" : "no", conf.threads );

        conf.cutup == 0 ?
        fprintf( stderr, CINFO "    Upper cutoff               : Don't cut\n"                          ):
        fprintf( stderr, CINFO "    Upper cutoff               : %lf\n",                    conf.cutup );

        conf.epgap <= 0 ?
        fprintf( stderr, CINFO "    Relative MIP gap tolerance : Automatic\n"                          ):
        fprintf( stderr, CINFO "    Relative MIP gap tolerance : %lf\n",                    conf.epgap );

        conf.scrind == CPX_OFF ?
        fprintf( stderr, CINFO "    Display messages on screen : No\n"                                  ):
        fprintf( stderr, CINFO "    Display messages on screen : Yes\n"                                 );

        conf.seed == 0 ?
        fprintf( stderr, CINFO "    Random seed                : Automatic\n"                          ):
        fprintf( stderr, CINFO "    Random seed                : %i\n",                      conf.seed );

    }

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
