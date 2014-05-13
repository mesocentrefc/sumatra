
/*
 * Print a SSE register for debug purpose
 */

#ifdef __SSE2__

static void  printreg(VTYPE r)
{
	STYPE a0,a1,a2,a3,a4,a5,a6,a7;
#if VMODE
	STYPE a8,a9,a10,a11,a12,a13,a14,a15;
#endif

	a0= EXTRACT_REG(r,0);
	a1= EXTRACT_REG(r,1);
	a2= EXTRACT_REG(r,2);
	a3= EXTRACT_REG(r,3);
	a4= EXTRACT_REG(r,4);
	a5= EXTRACT_REG(r,5);
	a6= EXTRACT_REG(r,6);
	a7= EXTRACT_REG(r,7);
#if VMODE
	a8= EXTRACT_REG(r,8);
	a9= EXTRACT_REG(r,9);
	a10= EXTRACT_REG(r,10);
	a11= EXTRACT_REG(r,11);
	a12= EXTRACT_REG(r,12);
	a13= EXTRACT_REG(r,13);
	a14= EXTRACT_REG(r,14);
	a15= EXTRACT_REG(r,15);
#endif

printf( "a00 :-> %7d  %7d  %7d  %7d "
		" %7d  %7d  %7d  %7d "
#if VMODE
		"%7d  %7d  %7d  %7d "
		" %7d  %7d  %7d  %7d "
#endif
		"\n"
		, a0,a1,a2,a3,a4,a5,a6,a7
#if VMODE
		, a8,a9,a10,a11,a12,a13,a14,a15
#endif
);
}

/*
 * set position p of a SSE register with the value v
 */

static inline VTYPE insert_reg(VTYPE r, STYPE v, int p)
{
	switch (p) {
	case 0: return INSERT_REG(r,v,0);
	case 1: return INSERT_REG(r,v,1);
	case 2: return INSERT_REG(r,v,2);
	case 3: return INSERT_REG(r,v,3);
	case 4: return INSERT_REG(r,v,4);
	case 5: return INSERT_REG(r,v,5);
	case 6: return INSERT_REG(r,v,6);
	case 7: return INSERT_REG(r,v,7);
#if VMODE
	case 8: return INSERT_REG(r,v,8);
	case 9: return INSERT_REG(r,v,9);
	case 10: return INSERT_REG(r,v,10);
	case 11: return INSERT_REG(r,v,11);
	case 12: return INSERT_REG(r,v,12);
	case 13: return INSERT_REG(r,v,13);
	case 14: return INSERT_REG(r,v,14);
	case 15: return INSERT_REG(r,v,15);
#endif
	}
	return _MM_SETZERO_SI128();
}

static inline STYPE extract_reg(VTYPE r, int p)
{
	switch (p) {
	case 0: return EXTRACT_REG(r,0);
	case 1: return EXTRACT_REG(r,1);
	case 2: return EXTRACT_REG(r,2);
	case 3: return EXTRACT_REG(r,3);
	case 4: return EXTRACT_REG(r,4);
	case 5: return EXTRACT_REG(r,5);
	case 6: return EXTRACT_REG(r,6);
	case 7: return EXTRACT_REG(r,7);
#if VMODE
	case 8: return EXTRACT_REG(r,8);
	case 9: return EXTRACT_REG(r,9);
	case 10: return EXTRACT_REG(r,10);
	case 11: return EXTRACT_REG(r,11);
	case 12: return EXTRACT_REG(r,12);
	case 13: return EXTRACT_REG(r,13);
	case 14: return EXTRACT_REG(r,14);
	case 15: return EXTRACT_REG(r,15);
#endif
	}
	return 0;
}

#define GET_H_SYMBOLE(s,p) ((p && p < lseq1) ? (s)[(p)-1]:255)
#define GET_V_SYMBOLE(s,p) ((p && p < lseq2) ? (s)[(p)-1]:0)

#define LSHIFT_SCORE(r)      { r = _MM_SLLI_SI128((r),sizeof(STYPE)); }
#define SET_H_SYMBOLE(r,p,s) { r = insert_reg((r),(STYPE)GET_H_SYMBOLE(seq1,(s)),(p)); }
#define PUSH_V_SYMBOLE(r,s)  { r = insert_reg(_MM_SLLI_SI128((r),sizeof(STYPE)),(STYPE)GET_V_SYMBOLE(seq2,(s)),0); }
#define EQUAL(f1,f2)         _MM_AND_SI128(EQUAL_REG((f1),(f2)),SET_CONST(1))

int FASTLCSSCORE(const char* seq1, const char* seq2,column_pp ppcolumn,int32_t* lpath)
{
	int lseq1,lseq2;		// length of the both sequences

	int itmp;				// tmp variables for swap
	const char* stmp;		//

	int nbands;				// Number of bands of width eight in the score matrix
	int lastband;			// width of the last band

							// Register for scanning the score matrix
	VTYPE  minus1;
	VTYPE  minus2;
	VTYPE  current;

	VTYPE  left;
	VTYPE  top;
	VTYPE  diag;


	VTYPE  sminus1;
	VTYPE  sminus2;
	VTYPE  scurrent;

	VTYPE  sleft;
	VTYPE  stop;
	VTYPE  sdiag;

	VTYPE  way;
	VTYPE  onevect;
	VTYPE  maxvect;

	VTYPE  fhseq;          	// The fragment of the horizontal sequence
							// to consider for aligment
	VTYPE  fvseq;			// The fragment of the horizontal sequence
							// to consider for aligment
	VTYPE  match;

	int band;
	int line;
	int limit;

	int lcs;

	int h;
	int i;

	column_t *column;


		// Made seq1 the longest sequences
	lseq1=strlen(seq1);
	lseq2=strlen(seq2);

	if (lseq1 < 10 || lseq2 < 10)
		return simpleLCS(seq1,seq2,ppcolumn,lpath);

	if (lseq1 < lseq2)
	{
		itmp=lseq1;
		lseq1=lseq2;
		lseq2=itmp;

		stmp=seq1;
		seq1=seq2;
		seq2=stmp;
	}

							// we add one to both lengths for taking into
							// account the extra line and column in the score
							// matrix

	lseq1++;
	lseq2++;

							// a band sized to the smallest sequence is allocated

	if (ppcolumn)
		column = *ppcolumn;
	else
		column=NULL;

	column = allocateColumn(lseq2,column,VMODE);

							// Check memory allocation
	if (column == NULL)
		return -1;

	for (i=0; i<lseq2;i++)
	{
		column->data.CMENB[i]=MIN_SCORE;
		column->score.CMENB[i]=-1;
	}

	nbands = lseq1 / VSIZE;					// You have VSIZE element in one SSE register
											// Alignment will be realized in nbands

	lastband = lseq1 - (nbands * VSIZE);	// plus one of width lastband except if
											// lastband==0

	if (lastband) nbands++;
	else lastband=VSIZE;

	lastband--;

//	printf("seq1 : %s  seq2 : %s\n",seq1,seq2);


	minus2 = SET_CONST(MIN_SCORE);
	minus1 = _MM_SETZERO_SI128();

	sminus1= _MM_SETZERO_SI128();
	sminus2= _MM_SETZERO_SI128();
	onevect= SET_CONST(1);
	maxvect= SET_CONST(MAX_SCORE);

	h=0;

	fhseq = _MM_SETZERO_SI128();
	fvseq = _MM_SETZERO_SI128();

					//
					// Beginning of the first band
					//

	for (line = 0; line < VSIZE; line++,h++) // avant VSIZE - 1
	{
//		printf("line= %4d   h= %4d\n",line,h);
		SET_H_SYMBOLE(fhseq,line,h)
		PUSH_V_SYMBOLE(fvseq,line)
		minus2 = insert_reg(minus2,0,h);
		minus1 = insert_reg(minus1,MIN_SCORE,line); // 0 avant
		match = EQUAL(fhseq,fvseq);

		if (lpath)
		{
			sminus2 = insert_reg(sminus2,line-1,line);  // Je ne suis pas certain de l'initialisation
			sminus1 = insert_reg(sminus1,0,line);
		}

//		printreg(fvseq);
//		printreg(fhseq);
//		printreg(match);
//		printf("================================\n");

		current = minus1;      // The best score is the upper one
							   // It cannot be the best as set to MIN_SCORE

		left = minus1;

//		printf("Vert = "); printreg(current);


		LSHIFT_SCORE(minus1)    // I shift minus1 so now I'll compare with the left position
		minus1=insert_reg(minus1,(column)->data.CMENB[line],0);

		top=minus1;

		if (lpath)
		{
			sleft=sminus1;  // I store the path length corresponding to the upper path
			LSHIFT_SCORE(sminus1)  // I shift to prepare the score coming from the left side
			sminus1=insert_reg(sminus1,(column)->score.CMENB[line],0);
			stop=sminus1;
			sdiag=sminus2;

		}

//		printf("Horz = "); printreg(minus1);

		current = GET_MAX(current,minus1); // Look for the best between upper and left

//		printf("BstHV= "); printreg(current);
//
//		printf("Diag = "); printreg(ADD_REG(minus2,match));

		diag=minus2;

		// minus2 = ;	// Minus2 contains the diagonal score, so I add the match reward
		                // Diag score are setup to 0 so this one will win on the first iteration
		current = GET_MAX(current,ADD_REG(minus2,match));

		if (lpath)
		{
//			printf("\n");
//			printf("current: ");
//			printreg(current);
//			printf("current: ");
//			printreg(SUB_REG(current,match));
//			printf("diag   : ");
//			printreg(diag);
//			printf("left   : ");
//			printreg(left);
//			printf("top    : ");
//			printreg(top);


			way     = EQUAL_REG(SUB_REG(current,match),diag);
			scurrent= OR_REG(AND_REG(way,sdiag),
					         ANDNOT_REG(way,maxvect));
//			printf("sdiag  : ");
//			printreg(scurrent);
			way     = EQUAL_REG(current,left);
			scurrent= GET_MIN(scurrent,OR_REG(AND_REG(way,sleft),
							  ANDNOT_REG(way,maxvect)));

//			printf("sleft  : ");
//			printreg(scurrent);
			way     = EQUAL_REG(current,top);
			scurrent= GET_MIN(scurrent,OR_REG(AND_REG(way,stop),
							  ANDNOT_REG(way,maxvect)));
//			printf("stop   : ");
//			printreg(scurrent);

			scurrent= ADD_REG(scurrent,onevect);

			sminus2=sminus1;
			sminus1=scurrent;
		}
//		printf("line %d :Best = ",line); printreg(current);
//
//		printf("================================\n");

		minus2=minus1;
		minus1=current;

//		printf("min2 = "); printreg(minus2);
//		printf("min1 = "); printreg(minus1);
//		printf("================================\n");

//		printf("\n");
//		printf("sdiag  : ");
//		printreg(sminus2);
//		printf("scur   : ");
//		printreg(scurrent);
//		printf("current: ");
//		printreg(current);
//		printf("%8s\n",seq1);
//		printf("%8s\n",seq2);
//		printf("================================\n");


	}  ///// <<<<<<<<------- Fin du debut de la premiere bande


//		printf("================================\n");

	(column)->data.CMENB[lseq2-VSIZE+line]=EXTRACT_REG(current,VSIZE-1);


	if (lpath)
		(column)->score.CMENB[lseq2-VSIZE+line]=EXTRACT_REG(scurrent,VSIZE-1);



	for (band=0; band < nbands; band++)
	{
//		SET_H_SYMBOLE(fhseq,line,h)
//		minus2 = insert_reg(minus2,0,line);
//		minus1 = insert_reg(minus1,MIN_SCORE,line); // 0 avant
//		h++;

		for (; line < lseq2; line++)
		{
//			printf("Je tourne avec line= %d \n",line);
			PUSH_V_SYMBOLE(fvseq,line)

			match = EQUAL(fhseq,fvseq);

//			printreg(fvseq);
//			printreg(fhseq);
//			printreg(match);
//			printf("================================\n");

			current = minus1;

			left = minus1;

			// Store the last current score in extra column
			(column)->data.CMENB[line-VSIZE]=EXTRACT_REG(current,VSIZE-1);
			LSHIFT_SCORE(minus1)
			minus1=insert_reg(minus1,(column)->data.CMENB[line],0);

			top = minus1;

//			printf("Vert = "); printreg(current);

			if (lpath)
			{
				sleft= sminus1;
				(column)->score.CMENB[line-VSIZE]=EXTRACT_REG(scurrent,VSIZE-1);
				LSHIFT_SCORE(sminus1)
				sminus1=insert_reg(sminus1,(column)->score.CMENB[line],0);
				stop=sminus1;
				sdiag=sminus2;
			}

//			printf("line = %d --> get = %d\n",line,(column)->data.CMENB[line]);

//			printf("Horz = "); printreg(minus1);

			current = GET_MAX(current,minus1);

			diag=minus2;

			current = GET_MAX(current,ADD_REG(minus2,match));

			if (lpath)
			{
//				printf("\n");
//				printf("current: ");
//				printreg(current);
//				printf("current: ");
//				printreg(SUB_REG(current,match));
//				printf("diag   : ");
//				printreg(diag);
//				printf("left   : ");
//				printreg(left);
//				printf("top    : ");
//				printreg(top);

				way     = EQUAL_REG(SUB_REG(current,match),diag);
				scurrent= OR_REG(AND_REG(way,sdiag),
						         ANDNOT_REG(way,maxvect));

//				printf("sdiag  : ");
//				printreg(scurrent);

				way     = EQUAL_REG(current,left);
				scurrent= GET_MIN(scurrent,OR_REG(AND_REG(way,sleft),
								  ANDNOT_REG(way,maxvect)));

//				printf("sleft  : ");
//				printreg(scurrent);

				way     = EQUAL_REG(current,top);
				scurrent= GET_MIN(scurrent,OR_REG(AND_REG(way,stop),
								  ANDNOT_REG(way,maxvect)));

//				printf("stop   : ");
//				printreg(scurrent);

				scurrent= ADD_REG(scurrent,onevect);

				sminus2=sminus1;
				sminus1=scurrent;
			}

			minus2=minus1;
			minus1=current;

//			printf("\n");
//			printf("sdiag  : ");
//			printreg(sminus2);
//			printf("scur   : ");
//			printreg(scurrent);
//			printf("current: ");
//			printreg(current);
//			printf("%8s\n",seq1);
//			printf("%8s\n",seq2);
		}
//		printf("================================\n");

								// end of the band and beginnig of the next one

		limit=(band==(nbands-1)) ? lastband:VSIZE;

		for (line = 0; line < limit; line++,h++)
		{
//			printf("Je fini avec line= %d \n",line);

			SET_H_SYMBOLE(fhseq,line,h)
			PUSH_V_SYMBOLE(fvseq,line)


			minus2 = insert_reg(minus2,MIN_SCORE,line);
			minus1 = insert_reg(minus1,MIN_SCORE,line);
			current = minus1;
			left=minus1;

			match = EQUAL(fhseq,fvseq);

			if (lpath)
			{
				sminus2 = insert_reg(sminus2,lseq2-VSIZE+line,line);
				sminus1 = insert_reg(sminus1,h,line);
				sleft= sminus1;
			}


//			printf("\n");
//			printf("fhseq = "); printreg(fhseq);
//			printf("fvseq = "); printreg(fvseq);
//			printf("----------------------------------------------------------------\n");
//			printf("match = "); printreg(match);


			(column)->data.CMENB[lseq2-VSIZE+line]=EXTRACT_REG(current,VSIZE-1);
			LSHIFT_SCORE(minus1)
			minus1=insert_reg(minus1,(column)->data.CMENB[line],0);
			top=minus1;

			current = GET_MAX(current,minus1);

			if (lpath)
			{
				(column)->score.CMENB[lseq2-VSIZE+line]=EXTRACT_REG(scurrent,VSIZE-1);
				LSHIFT_SCORE(sminus1)
				sminus1=insert_reg(sminus1,(column)->score.CMENB[line],0);
				stop=sminus1;
				sdiag=sminus2;

				way     = EQUAL_REG(current,minus1);

				scurrent= OR_REG(AND_REG(way,sminus1),
						         ANDNOT_REG(way,scurrent));
			}


			diag=minus2;

			current = GET_MAX(current,ADD_REG(minus2,match));

			if (lpath)
			{
				way     = EQUAL_REG(SUB_REG(current,match),diag);
				scurrent= OR_REG(AND_REG(way,sdiag),
						         ANDNOT_REG(way,maxvect));

				way     = EQUAL_REG(current,left);
				scurrent= GET_MIN(scurrent,OR_REG(AND_REG(way,sleft),
								  ANDNOT_REG(way,maxvect)));

				way     = EQUAL_REG(current,top);
				scurrent= GET_MIN(scurrent,OR_REG(AND_REG(way,stop),
								  ANDNOT_REG(way,maxvect)));

				scurrent= ADD_REG(scurrent,onevect);

				sminus2=sminus1;
				sminus1=scurrent;
			}

//			printf("currt = "); printreg(current);

			minus2=minus1;
			minus1=current;

//			printf("\n");
//			printf("sdiag  : ");
//			printreg(sminus2);
//			printf("scur   : ");
//			printreg(scurrent);
//			printf("current: ");
//			printreg(current);
//			printf("%8s\n",seq1);
//			printf("%8s\n",seq2);

//			printf("Je stocke line= %d la valeur %d\n",lseq2-VSIZE+line,(column)->data.CMENB[lseq2-VSIZE+line]);
		}

	}

//	printf("\n");
//	printf("line = %d, h= %d, lastband = %d\n",line,h,lastband);
//	printf("currt = "); printreg(current);
	lcs  = extract_reg(current,lastband);

	if(lpath)
		*lpath= extract_reg(scurrent,lastband);
//	printf("lastband = %d (%d) lcs = %d\n",lastband,lseq2,lcs);

	if (ppcolumn)
		*ppcolumn=column;
	else
		freeColumn(column);

	return lcs;
}

#else
int FASTLCSSCORE(const char* seq1, const char* seq2,column_pp ppcolumn,int32_t* lpath)
{
	return simpleLCS(seq1,seq2,ppcolumn,lpath);
}

#endif /* __SSE2__ */

