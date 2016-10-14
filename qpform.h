/* qpform.h
 * 18 Mar 2016
 */

/* #include before: "treedef.h"         */

/* global variables for sequence form   */

/* sf payoffs [row][col]                */
extern  Payvec **sfpay;

/* constraint matrix [player][row][col]
 * sacrifice one unused pointer for player 0
 */
extern  int **sfconstr[PLAYERS];

/* allocates sequence form payoff and constraint matrices
 * sets  all  sf  payoffs to 0
 * allocate  realplan[0..PLAYERS-1]
 * assumes  nseqs[], nisets[]  are set via  genseqin()
 * frees old data and keeps track of old dimensions
 */
void allocqpf(void);

/* allocate & generate sequence form:  payoff and constraint matrices
 * h->seqin must be defined via  genseqin()
 */
void genqpf(int eps);

/* LCP  (M,q) for the QPE sequence form
 * not the covering vector  d
 * generates  sf,  allocates  LCP,  fills LCP
 */
void sfqpelcp(int eps);

/* translate LCP solution of perturbed game 
 * to quasi-perfect equilibrium solution of the original game
 * and returns behaviour strategy for player pl
 */
void transperturbsol(int pl, Rat *rplan);

/* how many isets of player  pl  have nondeterministic moves
 * in the behaviour plan  bplan  (e.g. if 0: pure strategy)
 * of the quasi-perfect equilibrium
 */
int  propermixisetsqe(int pl, Rat *bplan);

/* gives the behaviour plan  bplan  of player  pl
 * to stdout, in one line
 * bnewline: terminate with \n  
 */
void outqpbehavstrat(int pl, Rat *bplan, Bool bnewline);



