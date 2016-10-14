/* gmp-wrap.c
 * 13 July 2000
 * wrapper for GMP functions similar to  mp.h
 */ 

#include <stdio.h>
#include <stdlib.h>
        /* atoi()       */
#include <limits.h>
        /* INT_MAX, INT_MIN     */

#include "gmp.h"
#include "gmpwrap.h"

void greduce(mpz_t a, mpz_t b)
{
    mpz_t tmp;
    mpz_init(tmp);
    mpz_gcd(tmp, a, b);
    mpz_divexact(a, a, tmp);
    mpz_divexact(b, b, tmp);
    mpz_clear(tmp);
}

int gmptoi(mpz_t a, int *result, int bcomplain) 
{
    char smp [MAXGMPCHARS];       /* string to print  mp  into */
    
    gmptoa(a, smp);
    *result = atoi(smp);
    if (*result == INT_MAX || *result == INT_MIN)
        {
        if (bcomplain)
            {
            printf("Warning: Long integer %s ", smp);
            printf("overflown, replaced by %d\n", *result);
            }
        return 1;
        }
    else
        return 0;
}

void gmptorat(mpz_t num, mpz_t den, int *result_num, int *result_den) 
{
    char nummp [MAXGMPCHARS], trnumsmp [MAXGMPCHARS]; 
    char denmp [MAXGMPCHARS], trdensmp [MAXGMPCHARS];
    int trsize = 8;     
    
    memset(nummp, '\0', sizeof(nummp));
    memset(denmp, '\0', sizeof(denmp));
    memset(trnumsmp, '\0', sizeof(trnumsmp));
    memset(trdensmp, '\0', sizeof(trdensmp));
    
    gmptoa(den, denmp);
    gmptoa(num, nummp);
    
    if (strlen(denmp) >= strlen(nummp))
    {
        if (strlen(denmp) <= trsize)
        {
            *result_den = atoi(denmp);
            *result_num = atoi(nummp);
        }
        else
        {
            strncpy(trdensmp, denmp, trsize);
            *result_den = atoi(trdensmp);
            if (strlen(denmp) - strlen(nummp) > trsize) 
                *result_num = 0;
            else
            {
                strncpy(trnumsmp, nummp, trsize+strlen(nummp)-strlen(denmp));
                *result_num = atoi(trnumsmp);
            }                   
        }        
    }
    else
    {
        if (strlen(nummp) <= trsize)
        {
            *result_den = atoi(denmp);
            *result_num = atoi(nummp);
        }
        else
        {
            strncpy(trnumsmp, nummp, trsize);
            *result_num = atoi(trnumsmp);
            if (strlen(nummp) - strlen(denmp) > trsize) 
                *result_den = 0;
            else
            {
                strncpy(trdensmp, denmp, trsize+strlen(denmp)-strlen(nummp));
                *result_den = atoi(trdensmp);
            }                   
        }        
    }    
}

