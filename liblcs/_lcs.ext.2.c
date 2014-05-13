#include "_lcs.h"

#include <string.h>
#include <stdlib.h>
#include <limits.h>

#include <stdio.h>




#define VSIZE (8)
#define VTYPE vInt16
#define STYPE int16_t
#define CMENB shrt
#define VMODE false
#define FASTLCSSCORE fastLCSScore16
#define INSERT_REG _MM_INSERT_EPI16
#define EXTRACT_REG _MM_EXTRACT_EPI16
#define EQUAL_REG  _MM_CMPEQ_EPI16
#define GREATER_REG  _MM_CMPGT_EPI16
#define SMALLER_REG  _MM_CMPLT_EPI16
#define ADD_REG    _MM_ADD_EPI16
#define SUB_REG    _MM_SUB_EPI16
#define AND_REG    _MM_AND_SI128
#define ANDNOT_REG    _MM_ANDNOT_SI128
#define OR_REG    _MM_OR_SI128
#define SET_CONST  _MM_SET1_EPI16
#define GET_MAX    _MM_MAX_EPI16
#define GET_MIN    _MM_MIN_EPI16
#define MIN_SCORE  INT16_MIN
#define MAX_SCORE  32000

#include "_lcs_fast.h"
