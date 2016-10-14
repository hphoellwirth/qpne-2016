/* qpform.c
 * 18 Mar 2016
 */

#include <stdio.h>
#include <stdlib.h>
	/* free()       */
#include "alloc.h"
	/* CALLOC(n,s), TALLOC(n,type), T2ALLOC(ptr,nrows,ncols,type)   */
	/* FREE2(ptr,nrows)                                             */
#include "col.h"
#include "rat.h"
	/* typedef Rat                  */
	/* ratadd, ratfromi, ratmult    */
#include "lemke.h"
#include "treedef.h"
#include "sfnf.h"
#include "seqform.h"

#include "qpform.h"

/* global variables for qpe sequence form   */

Payvec **sfpay;
int **sfconstr[PLAYERS];
int *sfconsvec[PLAYERS];
Rat *sfperturb[PLAYERS];

void allocqpf(void)
{
    static int oldnseqs1 = 0;
    static int oldconstrows[PLAYERS] = {0, 0, 0};
    static int oldconvecrows[PLAYERS] = {0, 0, 0};
    int pl, i, j;
    int nrows;
    
    /* payoff matrices, two players only here, init to pay 0        */
    FREE2(sfpay, oldnseqs1);
    oldnseqs1 = nseqs[1];
    T2ALLOC (sfpay, nseqs[1], nseqs[2], Payvec);
    for (i=0; i<nseqs[1]; i++)
	for (j=0; j<nseqs[2]; j++)
	    for (pl=1; pl < PLAYERS; pl++)
		sfpay[i][j][pl-1] = ratfromi(0);
		
    /* constraint matrices, any number of players           */
    /* sfconstr[0] stays unused                             */
    for (pl=1; pl < PLAYERS; pl++)
	{
	FREE2(sfconstr[pl], oldconstrows[pl]);
	oldconstrows[pl] = nrows = nisets[pl]+1;   /* extra row for seq 0  */
	T2ALLOC (sfconstr[pl], nrows, nseqs[pl], int);
	}
	
    /* constraint vectors, any number of players            */
    for (pl=1; pl < PLAYERS; pl++)
	{
	    free(sfconsvec[pl]);
	    sfconsvec[pl] = TALLOC(nisets[pl]+1, int);	
	}
			
    /* perturbing vectors, any number of players            */
    for (pl=1; pl < PLAYERS; pl++)
	{
	    free(sfperturb[pl]);
	    sfperturb[pl] = TALLOC(nseqs[pl], Rat);
	}	
}       /* end of allocsf()     */


void genqpf(int eps)
{
    int pl, i, j;
    Outcome z; 
    allocqpf();
    
    behavtorealprob(0);     /* get realization probabilities of leaves      */
    
    /* sf payoff matrices                   */
    for (z=outcomes; z < lastoutcome; z++)
	{
	Node u = z->whichnode;
	i = u->defseq[1] - firstmove[1];
	j = u->defseq[2] - firstmove[2];
	for (pl=1; pl < PLAYERS; pl++)
	    sfpay[i][j][pl-1] = ratadd(sfpay[i][j][pl-1],
		    ratmult(u->defseq[0]->realprob, z->pay[pl-1]) );
	}
	
    /* sf constraint matrices, sparse fill  */
    for (pl=1; pl < PLAYERS; pl++)
	{
	sfconstr[pl][0][0] = 1;     /* empty sequence                       */
	for (i=0; i < nisets[pl]; i++)
	    sfconstr[pl][i+1][(firstiset[pl]+i)->seqin - firstmove[pl]] = -1;
	for (j=1; j < nseqs[pl]; j++)
	    sfconstr[pl][(firstmove[pl]+j)->atiset - firstiset[pl]+1][j] = 1;
	}
	
	if (eps > 0)
	{
        /* sf constraint vertices  */ 
        for (pl=1; pl < PLAYERS; pl++)
	    {
	        sfconsvec[pl][0] = 1;       
	        for (i=1; i < nisets[pl]+1; i++)
	            sfconsvec[pl][i] = 0;	 
	    }		
	
        /* sf perturbing vectors  */
        for (pl=1; pl < PLAYERS; pl++)
	    {
	        sfperturb[pl][0] = ratfromi(0); 
            for (i=1; i < nseqs[pl]; i++)
	        { 
	            int seqlen = seqtoalen(firstmove[pl] + i, pl);
	            sfperturb[pl][i] = ratinv(ratfromi(pow(eps,seqlen)));     
	        }  
	    }
	}
}       /* end of  gensf()              */


void sfqpelcp(int eps)
{
    int i;

    genqpf(eps);
    setlcp( nseqs[1] + 2*(nisets[2]+1) + nseqs[2] + 2*(nisets[1]+1) );
    /* fill  M  */
    /* -A       */
    payratmatcpy(sfpay, 0, 1, 0, nseqs[1], nseqs[2], 
	lcpM, 0, nseqs[1]);
    /* E\T      */
    intratmatcpy(sfconstr[1], 0, 1, nisets[1]+1, nseqs[1], 
	lcpM, 0, nseqs[1] + nseqs[2]);	
    /* -E\T     */
    intratmatcpy(sfconstr[1], 1, 1, nisets[1]+1, nseqs[1], 
	lcpM, 0, nseqs[1] + nseqs[2] + nisets[1]+1);

    /* -B\T     */
    payratmatcpy(sfpay, 1, 1, 1, nseqs[1], nseqs[2], 
	lcpM, nseqs[1], 0);
    /* F\T      */
    intratmatcpy(sfconstr[2], 0, 1, nisets[2]+1, nseqs[2], 
	lcpM, nseqs[1], nseqs[1] + nseqs[2] + 2*(nisets[1]+1));	
    /* -F\T     */
    intratmatcpy(sfconstr[2], 1, 1, nisets[2]+1, nseqs[2], 
	lcpM, nseqs[1], nseqs[1] + nseqs[2] + 2*(nisets[1]+1) + nisets[2]+1);

    /* -E       */
    intratmatcpy(sfconstr[1], 1, 0, nisets[1]+1, nseqs[1], 
	lcpM, nseqs[1] + nseqs[2], 0 );	
    /* E        */
    intratmatcpy(sfconstr[1], 0, 0, nisets[1]+1, nseqs[1], 
	lcpM, nseqs[1] + nseqs[2] + nisets[1]+1, 0 );
	
    /* -F       */
    intratmatcpy(sfconstr[2], 1, 0, nisets[2]+1, nseqs[2], 
	lcpM, nseqs[1] + nseqs[2] + 2*(nisets[1]+1), nseqs[1]);	
    /* F        */
    intratmatcpy(sfconstr[2], 0, 0, nisets[2]+1, nseqs[2], 
	lcpM, nseqs[1] + nseqs[2] + 2*(nisets[1]+1) + nisets[2]+1, nseqs[1]);	
	
    /* define RHS q,  using special shape of SF constraints RHS e,f     */
    for (i = 0; i < lcpdim; i++)
	    rhsq[i] = ratfromi(0);
	
	int pos;	
    /* -Al        */ 
	Rat sfpayvec1[nseqs[2]];    
    for (i = 0; i < nseqs[1]; i++)
    {
        for (int j = 0; j < nseqs[2]; j++)
            sfpayvec1[j] = sfpay[i][j][0];          
        pos = i;    
        ratscalarprod(nseqs[2], sfpayvec1, sfperturb[2], &(rhsq[pos]));
        rhsq[pos] = ratneg(rhsq[pos]);
    }  
    /* -B\Tk       */
	Rat sfpayvec2[nseqs[1]];    
    for (i = 0; i < nseqs[2]; i++)
    {
        for (int j = 0; j < nseqs[1]; j++)
            sfpayvec2[j] = sfpay[j][i][1];          
        pos = nseqs[1] + i;    
        ratscalarprod(nseqs[1], sfpayvec2, sfperturb[1], &(rhsq[pos]));
        rhsq[pos] = ratneg(rhsq[pos]);
    }    
    /* e - Ek     */    
    for (i = 0; i < nisets[1]+1; i++)
    {
    pos = nseqs[1] + nseqs[2] + i;
    intratscalarprod(nseqs[1], sfconstr[1][i], sfperturb[1], &(rhsq[pos]));
    rhsq[pos] = ratadd(ratfromi(sfconsvec[1][i]), ratneg(rhsq[pos]));
	}
	/* -e + Ek    */
    for (i = 0; i < nisets[1]+1; i++)
    {
    pos = nseqs[1] + nseqs[2] + nisets[1]+1 + i;
    intratscalarprod(nseqs[1], sfconstr[1][i], sfperturb[1], &(rhsq[pos]));
    rhsq[pos] = ratadd(ratneg(ratfromi(sfconsvec[1][i])), rhsq[pos]);
	}	
	/* f - Fl     */
    for (i = 0; i < nisets[2]+1; i++)
    {
    pos = nseqs[1] + nseqs[2] + 2*(nisets[1]+1) + i;
    intratscalarprod(nseqs[2], sfconstr[2][i], sfperturb[2], &(rhsq[pos]));
    rhsq[pos] = ratadd(ratfromi(sfconsvec[2][i]), ratneg(rhsq[pos]));
	} 	
	/* -f + Fl    */
    for (i = 0; i < nisets[2]+1; i++)
    {
    pos = nseqs[1] + nseqs[2] + 2*(nisets[1]+1) + nisets[2]+1 + i;
    intratscalarprod(nseqs[2], sfconstr[2][i], sfperturb[2], &(rhsq[pos]));
    rhsq[pos] = ratadd(ratneg(ratfromi(sfconsvec[2][i])), rhsq[pos]);
	} 		   
}

void transperturbsol(int pl, Rat *rplan)
{
    char s[MAXSTRL];
    int i,j=0,k=0;
    Move c;
    Iset h;
    Rat bprob, predprob, redistr_sum, nzero_sum, redistr_share;
    Rat cpplan[nseqs[pl]];
   
    for (i = 0; i < nseqs[pl]; i++) 
    {
        cpplan[i] = rplan[i];       
    }

    for (h = firstiset[pl]; h < firstiset[pl+1]; h++)
    {
        /* identify actions with minimum (epsilon) probability */
        /* and prepare re-distribution of probability to       */
        /* actions with non-minimum probability                */    
        redistr_sum = ratfromi(0);
        nzero_sum   = ratfromi(0);
        for (c = h->move0, i=0; i < h->nmoves; c++, i++)
            {
            j++;
            bprob = rplan[ c - firstmove[pl] ];
            
            if ( bprob.num == 0)
                redistr_sum = ratadd(redistr_sum, sfperturb[pl][j]);
            else
                nzero_sum = ratadd(nzero_sum, ratadd(bprob,sfperturb[pl][j]));
            }       
          
        /* adjust non-minimum probabilities */   
        for (c = h->move0, i=0; i < h->nmoves; c++, i++)
            {
            k++;  
            bprob = rplan[ c - firstmove[pl] ];
                         
            if ( bprob.num != 0)
                {     
                /* add minimum probability */        
                bprob = ratadd(bprob, sfperturb[pl][k]);               
                cpplan[ c - firstmove[pl] ] = bprob;
                
                /* add share of redistributable probability */
                redistr_share = ratmult(redistr_sum, ratdiv(bprob, nzero_sum));
                bprob = ratadd(bprob, redistr_share);
                predprob = cpplan[ h->seqin - firstmove[pl]];                                                                     
                
                if (ratiseq(rplan[ h->seqin - firstmove[pl]], ratfromi(0))) 
                    rplan[ c - firstmove[pl] ] = ratfromi(0);
                else             
                    rplan[ c - firstmove[pl] ] = ratround(ratdiv(bprob, predprob));
                }
            }
    }
} 

int  properqpmixisets(int pl, Rat *bplan)
{
    int mix = 0;
    int i;
    Move c;
    Iset h;

    for (h = firstiset[pl]; h < firstiset[pl+1]; h++)
        for (c = h->move0, i=0; i < h->nmoves; c++, i++)
            if ( bplan[ c - firstmove[pl] ].num != 0  &&
                !ratiseq( bplan[ c - firstmove[pl] ],
                          ratfromi(1)) )
                {
                mix++ ;
                break ;
                }
    return mix;
}

void outqpbehavstrat(int pl, Rat *bplan, Bool bnewline)
{
    char s[MAXSTRL];
    int i;
    Move c;
    Iset h;
    Rat bprob;

    for (h = firstiset[pl]; h < firstiset[pl+1]; h++)
        for (c = h->move0, i=0; i < h->nmoves; c++, i++)
            {
            bprob = bplan[ c - firstmove[pl] ];
            if ( bprob.num != 0)
                {
                movetoa(c, pl, s);
                printf(" %s", s);
                if (!ratiseq(bprob, ratfromi(1) ) )
                    {
                    rattoa(bprob, s);
                    printf(":%s", s);
                    }
                }           
            }
    if (bnewline)
        printf("\n");
}



