/*!
 * \file    tsp.c
 * \brief   Parse a TSP problem into a convenient datastructure.
 * \authors Francesco Cazzaro, Marco Cieno
 */
#include <argp.h>
#include <errno.h>
#include <float.h>
#include <limits.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include "logging.h"
#include "tsp.h"


/* Function declarations */

/*!
 * \fn    parse_opt
 * \brief Option parser.
 *
 * This function will be the callback used by argp_parse() to check command line arguments consistency and populate
 * the instance structure.
 *
 *
 * \param key
 *     Value of the current option being parsed.
 *
 * \param arg
 *     Pointer to the value of such option if any, NULL otherwise.
 *
 * \param argp_state
 *     Pointer to the current parser state. Will hold a pointer to the instance structure in the `input` field.
 *
 * \returns Numeric value of the error occurred. 0 on success.
 */
static error_t
parse_opt ( int key, char *arg, struct argp_state *state );

/*!
 * \fn    _print_parsed_args
 * \brief Debugging function to inspect parsed command line arguments.
 *
 *
 * \param problem
 *     Pointer to the instance structure.
 */
void
_print_parsed_args ( instance *problem );

/*!
 * \fn    parse_tsp_file
 * \brief Read input TSP file into instance data structure.
 *
 *
 * \param problem
 *     Pointer to the problem being parsed.
 */
void
parse_tsp_file ( instance *problem );


/* Argp setup */

const char *argp_program_version     = "tsp_parser brought to you by Francesco Cazzaro and Marco Cieno";
const char *argp_program_bug_address = "{marco.cieno, francesco.cazzaro}@studenti.unipd.it";
static char doc[]                    = "Parse a TSP problem file into a convenient data structure.";
static char args_doc[]               = "TSP_FILE";
static struct argp_option options[]  =
{
    { "cutoff",    'c',     "VALUE",   OPTION_NO_USAGE,                       "Master cutoff value."              },
    { "threads",   'j',     "N",       OPTION_NO_USAGE | OPTION_ARG_OPTIONAL, "Use multithread. Default 4"        },
    { "memory",    'm',     "AVAIL",   OPTION_NO_USAGE,                       "Available memory (size in MB)."    },
    { "timelimit", 't',     "SECONDS", OPTION_NO_USAGE,                       "Maximum time the program may run." },

    { "verbose",   LOG_VBS, NULL,      OPTION_NO_USAGE,                       "Set program logging level."        },
    { "debug",     LOG_DBG, NULL,      OPTION_ALIAS,                          NULL                                },
    { "hidebug",   LOG_HID, NULL,      OPTION_ALIAS,                          NULL                                },

    { NULL },
};

static struct argp argp = { options, parse_opt, args_doc, doc, 0, 0, 0 };


int
main ( int argc, char *argv[] )
{


    instance problem;


    init_instance( &problem );

    argp_parse( &argp, argc, argv, 0, 0, &problem );
    _print_parsed_args( &problem );

    parse_tsp_file( &problem );

    // ...

    destroy_instance( &problem );
}


void
_print_parsed_args ( instance *problem )
{
    if ( problem->loglevel >= LOG_VBS ) {
        fprintf( stderr, "Arguments parsed:\n" );
        fprintf( stderr, "  * TSP file            : %s\n",                                   problem->filename );
        fprintf( stderr, "  * Time limit          : %llu hours %llu minutes %llu seconds\n",
            problem->timelimit / 3600ULL, problem->timelimit % 3600ULL / 60ULL, problem->timelimit % 60        );
        fprintf( stderr, "  * Master cutoff value : %e\n",                                   problem->cutoff   );
        fprintf( stderr, "  * Maximum memory      : %llu MB\n",                              problem->memory   );
        fprintf( stderr, "  * Use multithread     : %s (%lu)\n\n",
            problem->threads > 1 ? "yes" : "no", problem->threads                                              );
    }
}


static error_t
parse_opt ( int key, char *arg, struct argp_state *state )
{
    instance *problem = state->input;

    switch ( key )
    {
        case 'c':
            problem->cutoff = strtod( arg, NULL );
            if ( errno != 0 || problem->cutoff == 0. ) {
                argp_error(
                    state,
                    "Bad value for option -c --cutoff: %s", strerror( errno ? errno : EDOM )
                );
            }

            break;


        case 'j':
            if ( arg == NULL ) {
                /* Unless differently specified, `threads` will be set to 4. */
                problem->threads = 4UL;
            } else {
                problem->threads = strtoul( arg, NULL, 10 );
                if ( errno != 0 || problem->threads == 0UL ) {
                    argp_error(
                        state,
                        "Bad value for option -j --threads: %s", strerror( errno ? errno : EDOM )
                    );
                }
            }

            break;


        case 'm':
            problem->memory = strtoull( arg, NULL, 10 );
            if ( errno != 0 || problem->timelimit == 0ULL ) {
                argp_error(
                    state,
                    "Bad value for option -m --memory: %s", strerror( errno ? errno : EDOM )
                );
            }

            break;


        case 't':
            problem->timelimit = strtoull( arg, NULL, 10 );
            if ( errno != 0 || problem->timelimit == 0ULL ) {
                argp_error(
                    state,
                    "Bad value for option -t --timelimit: %s", strerror( errno ? errno : EDOM )
                );
                return ERANGE;
            }

            break;


        case LOG_VBS:
            problem->loglevel = LOG_VBS;

            break;


        case LOG_DBG:
            problem->loglevel = LOG_DBG;

            break;


        case LOG_HID:
            problem->loglevel = LOG_HID;

            break;


        case ARGP_KEY_ARG:
            problem->filename = calloc( strlen( arg ), sizeof( *arg ) );
            if ( problem->filename == NULL ) {
                argp_failure(
                    state,
                    1, errno,
                    "Could store filename: %s", arg
                );
            }

            problem->filename = strcpy( problem->filename, arg );

            break;


        case ARGP_KEY_END:
            if ( problem->filename == NULL ) {
                argp_error( state, "Missing TSP file name" );
            }

            break;


        default:
            return ARGP_ERR_UNKNOWN;
    }

    return 0;
}


void
parse_tsp_file ( instance *problem )
{
    FILE *fin = fopen( problem->filename, "r" );

    if ( fin == NULL ) {
        perror( "parse_tsp_file" );
        exit( EXIT_FAILURE );
    }

    problem->nnodes = 0;

    // TODO
    char line[180];
	char *par_name;   
	char *token1;
	char *token2;
	
	int active_section = 0; // =1 NODE_COORD_SECTION 
	

	while ( fgets(line, sizeof(line), fin) != NULL ) 
	{

		if ( strlen(line) <= 1 ) continue; 
	    par_name = strtok(line, " :");
        
		if ( strncmp(par_name, "NAME", 4) == 0 ) 
		{
			active_section = 0;
			continue;
		}

		if ( strncmp(par_name, "COMMENT", 7) == 0 ) 
		{
			active_section = 0;   
			continue;
		}   
		
		if ( strncmp(par_name, "TYPE", 4) == 0 ) 
		{
			token1 = strtok(NULL, " :");  
            //fprintf(stderr, "TRESD: %d", *token1);
			if ( strncmp(token1, "TSP",3) != 0 )  {perror( "format error: only TYPE == TSP \n"); exit( EXIT_FAILURE ); }
			active_section = 0;
			continue;
		}
		

		if ( strncmp(par_name, "DIMENSION", 9) == 0 ) 
		{
			if ( problem->nnodes > 0 ) {perror( "Repeated DIMENSION error in file\n" ); exit( EXIT_FAILURE ); }
			token1 = strtok(NULL, " :");
			problem->nnodes = strtoull(token1, NULL, 0);	 
			problem->xcoord = (double *) calloc(problem->nnodes, sizeof(double)); 	 
			problem->ycoord = (double *) calloc(problem->nnodes, sizeof(double));    
			active_section = 0;  
			continue;
		}


		if ( strncmp(par_name, "EDGE_WEIGHT_TYPE", 16) == 0 ) 
		{ 
			active_section = 0;
			continue;
		}            
		
		if ( strncmp(par_name, "NODE_COORD_SECTION", 18) == 0 ) 
		{
			if ( problem->nnodes <= 0 ) {perror( "DIMENSION section should appear before NODE_COORD_SECTION\n" ); exit( EXIT_FAILURE ); }
			active_section = 1;   
			continue;
		}
		
		if ( strncmp(par_name, "EOF", 3) == 0 ) 
		{
			active_section = 0;
			break;
		}
		
			
		if ( active_section == 1 ) // within NODE_COORD_SECTION
		{
			int i = atoi(par_name) - 1; 
			if ( i < 0 || i >= problem->nnodes ) {perror( "node out of range in NODE_COORD_SECTION\n" ); exit( EXIT_FAILURE ); }     
			token1 = strtok(NULL, " :,");
			token2 = strtok(NULL, " :,");
			problem->xcoord[i] = atof(token1);
			problem->ycoord[i] = atof(token2);
			continue;
		}    
		  
		
        {perror( "UNKNOWN FORMAT ARGUMENT IN FILE\n" ); exit( EXIT_FAILURE ); }     
		    
	}                


    fclose(fin);

    if(problem->loglevel >=LOG_VBS){
        fprintf( stderr, "TSP file parameters: \n");
        fprintf( stderr, "  * Number of nodes     : %llu\n\n", problem->nnodes );

    }

    if(problem->loglevel >= LOG_DBG){
        fprintf(stderr, "TSP nodes: \n");
        fprintf(stderr, "  * node 1 : %f, %f \n", problem->xcoord[0], problem->ycoord[0]);
        fprintf(stderr, "  * node 48: %f, %f \n\n", problem->xcoord[47], problem->ycoord[47]);
    }



}
