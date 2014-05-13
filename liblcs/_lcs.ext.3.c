#include "_lcs.h"

#include <string.h>
#include <stdlib.h>
#include <limits.h>

#include <stdio.h>




#define VSIZE (16)
#define VTYPE vInt8
#define STYPE int8_t
#define CMENB byte
#define VMODE true
#define FASTLCSSCORE fastLCSScore8
#define INSERT_REG _MM_INSERT_EPI8
#define EXTRACT_REG _MM_EXTRACT_EPI8
#define EQUAL_REG  _MM_CMPEQ_EPI8
#define GREATER_REG  _MM_CMPGT_EPI8
#define SMALLER_REG  _MM_CMPLT_EPI8
#define ADD_REG    _MM_ADD_EPI8
#define SUB_REG    _MM_SUB_EPI8
#define AND_REG    _MM_AND_SI128
#define ANDNOT_REG    _MM_ANDNOT_SI128
#define OR_REG    _MM_OR_SI128
#define SET_CONST  _MM_SET1_EPI8
#define GET_MAX    _MM_MAX_EPI8
#define GET_MIN    _MM_MIN_EPI8
#define MIN_SCORE  INT8_MIN
#define MAX_SCORE  127

#include "_lcs_fast.h"
