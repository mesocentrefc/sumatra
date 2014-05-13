/*
 * sumatra.h
 *
 *  Created on: 8 aožt 2010
 *      Author: coissac
 */

#ifndef SUMATRA_H_
#define SUMATRA_H_

#include "../libfasta/sequence.h"

void printResults(fastaSeqPtr seq1, fastaSeqPtr seq2, double score, BOOL extradata, int64_t pairs, BOOL print);

#endif /* SUMATRA_H_ */
