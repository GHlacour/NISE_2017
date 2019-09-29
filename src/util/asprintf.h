/**
 * This file implements asprintf for the MSVC compiler. This function is already part of the GNU
 * extensions of the library.
 */

#ifndef ASPRINTF_H
#define ASPRINTF_H

#if defined(__GNUC__) && ! defined(_GNU_SOURCE)
#define _GNU_SOURCE
#endif

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

#ifdef _MSC_VER
int vasprintf(char **strp, const char* format, va_list ap) {
    int len = _vscprintf(format, ap);
    if (len == -1) return -1;

    char* str = malloc((size_t)len + 1);
    if (!str) return -1;

    int retval = vsnprintf(str, len + 1, format, ap);
    if(retval == -1) {
        free(str);
        return -1;
    }

    *strp = str;
    return retval;
}

int asprintf(char **strp, const char* format, ...) {
    va_list ap;
    va_start(ap, format);

    int retval = vasprintf(strp, format, ap);

    va_end(ap);
    return retval;
}
#endif
#endif