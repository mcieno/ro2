/**
 * Copyright (c) 2017 rxi
 *
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the MIT license. See `log.c` for details.
 */

#ifndef LOG_H
#define LOG_H
#define LOG_USE_COLOR

#include <stdio.h>
#include <stdarg.h>

#define LOG_VERSION "0.1.0"

typedef void (*log_LockFn)(void *udata, int lock);

typedef enum
{
    LOG_TRACE,
    LOG_DEBUG,
    LOG_INFO,
    LOG_WARN,
    LOG_ERROR,
    LOG_FATAL,
    LOG_OUT
}
loglevel_t;

extern loglevel_t loglevel;

#ifdef NO_LOG

#define log_trace(...) do{}while(0)
#define log_debug(...) do{}while(0)
#define log_info(...)  do{}while(0)
#define log_warn(...)  do{}while(0)
#define log_error(...) do{}while(0)
#define log_fatal(...) do{}while(0)
#define log_out(...) do{}while(0)

#else

#define log_trace(...) log_log(LOG_TRACE, __FILE__, __func__, __LINE__, __VA_ARGS__)
#define log_debug(...) log_log(LOG_DEBUG, __FILE__, __func__, __LINE__, __VA_ARGS__)
#define log_info(...) log_log(LOG_INFO, __FILE__, __func__, __LINE__, __VA_ARGS__)
#define log_warn(...) log_log(LOG_WARN, __FILE__, __func__, __LINE__, __VA_ARGS__)
#define log_error(...) log_log(LOG_ERROR, __FILE__, __func__, __LINE__, __VA_ARGS__)
#define log_fatal(...) log_log(LOG_FATAL, __FILE__, __func__, __LINE__, __VA_ARGS__)
#define log_out(...) log_log(LOG_OUT, __FILE__, __func__, __LINE__, __VA_ARGS__)

#endif

void log_set_udata(void *udata);
void log_set_lock(log_LockFn fn);
void log_set_fp(FILE *fp);
void log_set_level(int level);
void log_set_quiet(int enable);

void log_log(int level, const char *file, const char *func, int line, const char *fmt, ...);


#endif
