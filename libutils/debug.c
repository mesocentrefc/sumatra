/*
 * debug.c
 *
 *  Created on: 4 sept. 2012
 *      Author: coissac
 */

#include <stdlib.h>
#include <inttypes.h>
#include "debug.h"

char* int2bin(int64_t i,size_t bits)
{
    static char str[65];
    uint64_t u;

    if (bits > 64)
    	return NULL;

    str[bits] = 0;

    // type punning because signed shift is implementation-defined
    u = *(unsigned *)&i;

    for(; bits--; u >>= 1)
        str[bits] = u & 1 ? '1' : '0';

    return str;
}



