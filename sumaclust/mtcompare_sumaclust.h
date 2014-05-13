/*
 * mtcompare.h
 *
 *  Created on: 12 mars 2013
 *      Author: celinemercier
 */

#ifndef MTCOMPARE_H_
#define MTCOMPARE_H_

int mt_compare_sumaclust(fastaSeqPtr* db, int n, BOOL fast, double threshold, BOOL normalize, int reference, BOOL lcsmode, int threads_number);

#endif /* MTCOMPARE_H_ */
