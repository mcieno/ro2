/*!
 * \file    debug.h
 * \brief   Definitions of logging levels.
 * \authors Francesco Cazzaro, Marco Cieno
 */
#ifndef LOGGING_H
#define LOGGING_H

#define LOG_ERR 600U  /*!< Errors-only logging level. */
#define LOG_VBS 700U  /*!< Verbose logging level. */
#define LOG_DBG 800U  /*!< Debug logging level. */
#define LOG_HID 900U  /*!< Highly verbose debug logging level. */

typedef unsigned loglevel_t;

#endif
