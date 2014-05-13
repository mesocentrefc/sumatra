/*
 * mtcompare_sumatra.h
 *
 *  Created on: 12 mars 2013
 *      Author: celinemercier
 */

#ifndef MTCOMPARE_SUMATRA_H_
#define MTCOMPARE_SUMATRA_H_

int mt_compare_sumatra(fastaSeqCount *db1, fastaSeqCount *db2, double threshold, BOOL normalize, int reference, BOOL lcsmode, BOOL extradata, int n);

#endif /* MTCOMPARE_SUMATRA_H_ */
