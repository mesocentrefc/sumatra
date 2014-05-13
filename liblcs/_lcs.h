#include "../libsse/_sse.h"

#define bool char
#define false (1==0)
#define true  (1==1)

typedef struct {
	int16_t    size;

	union { int16_t *shrt;
	        int8_t  *byte;
	      } data;

	union { int16_t *shrt;
			int8_t  *byte;
		  } score;


} column_t, **column_pp, *column_p;

column_p allocateColumn(int length,column_t *column, bool mode8bits);

void freeColumn(column_p column);

int fastLCSScore16(const char* seq1, const char* seq2,column_pp ppcolumn,int32_t* lpath);
int fastLCSScore8(const char* seq1, const char* seq2,column_pp ppcolumn,int32_t* lpath);
int simpleLCS(const char* seq1, const char* seq2,column_pp ppcolumn,int32_t* lpath);

int fastLCSScore(const char* seq1, const char* seq2,column_pp column,int32_t* lpath);
