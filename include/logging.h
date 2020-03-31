/*!
 * \file    logging.h
 * \brief   Define logging levels.
 * \authors Francesco Cazzaro, Marco Cieno
 */
#ifndef LOGGING_H
#define LOGGING_H


#define LOG_OFF   0U    /*!< No logging. */
#define LOG_FATAL 100U  /*!< Fatals-only logging. */
#define LOG_ERROR 200U  /*!< Errors logging. */
#define LOG_WARN  300U  /*!< Verbose logging. */
#define LOG_INFO  400U  /*!< Informative logging. */
#define LOG_DEBUG 500U  /*!< Debug logging. */
#define LOG_TRACE 600U  /*!< Verbose debug logging. */

#define CSUCC  "\033[92m[+]\033[0m  "  /*!< Prefix for SUCCESS logs.  */
#define CERROR "\033[91m[-]\033[0m  "  /*!< Prefix for ERROR logs.  */
#define CFATAL "\033[31m[!]\033[0m  "  /*!< Prefix for FATAL logs.  */
#define CINFO  "\033[94m[*]\033[0m  "  /*!< Prefix for INFO logs.  */
#define CDEBUG "\033[96m[.]\033[0m  "  /*!< Prefix for INFO logs.  */
#define CTRACE "\033[95m[^]\033[0m  "  /*!< Prefix for INFO logs.  */

typedef unsigned loglevel_t;

extern loglevel_t loglevel;


#endif
