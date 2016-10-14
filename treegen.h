/* treegen.h
 * 13 Apr 2000
 * generating examples of game trees
 */

#define PAYOFFSEEDPLUS  200 /* added to  seed  when creating random payoffs */

/* create a binary game tree with two-element info sets 
 * returns root = nodes + ROOT
 * seed for custom random generator
 */
void createbintree (int levels, int seed);

/* create tree with example in BvS/Elzen/Talman	 */
void tracingexample(void);

/* create tree with example Hurkens/Hauk */
void forwardexample(void);

/* create tree with Mertens simple example */
void mertensexample(void);

/* create tree with Mertens voting example */
void mertensvotingexample(void);

/* create tree with example van Damme */
void vandammeexample(void);

/* create tree for variant of matching pennies */
void matchingpennies(void);


