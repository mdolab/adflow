static char dpSid[]="$Id: dpStack.c 1663 2007-02-27 10:42:54Z llh $";

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#define ONE_BLOCK_SIZE 16384
#define CHUNK_SIZE 4096

/* The main stack is a double-chain of DoubleChainedBlock objects.
 * Each DoubleChainedBlock holds an array[ONE_BLOCK_SIZE] of char. */
typedef struct _doubleChainedBlock{
  struct _doubleChainedBlock *prev ;
  char                       *stackBottom ;
  char                       *stackTop ;
  char                       *stackLimit ;
  struct _doubleChainedBlock *next ;
} DoubleChainedBlock ;

/* Globals that define the current position in the stack: */
static DoubleChainedBlock *dpcurStack = NULL ;
static char               *dpcurStackBottom = NULL ;
static char               *dpcurStackTop    = NULL ;
static char               *dpcurStackLimit  = NULL ;
/* Globals that define the current LOOKing position in the stack: */
static DoubleChainedBlock *dplookStack = NULL ;
static char               *dplookStackBottom = NULL ;
static char               *dplookStackTop    = NULL ;
static char               *dplookStackLimit  = NULL ;

/* Used before PUSHing.
 * Resets the LOOKing position if it was active.
 * Checks that there is enough space left to hold "nbChars" chars.
 * Otherwise, allocates the necessary space. */
void dpcheck(unsigned int nbChars) {
  if (dplookStack) dplookStack = NULL ;
  if (dpcurStackTop+nbChars > dpcurStackLimit) {
    if ((dpcurStack == NULL) || (dpcurStack->next == NULL)) {
      DoubleChainedBlock *newStack ;
      char *contents = (char*)malloc(ONE_BLOCK_SIZE*sizeof(char)) ;
      newStack = (DoubleChainedBlock*)malloc(sizeof(DoubleChainedBlock)) ;
      newStack->prev = dpcurStack ;
      if (dpcurStack != NULL) {
	dpcurStack->stackTop = dpcurStackTop ;
	dpcurStack->next = newStack ;
      }
      newStack->next = NULL ;
      newStack->stackBottom = contents ;
      newStack->stackTop = contents ;
      newStack->stackLimit = contents+ONE_BLOCK_SIZE ;
      dpcurStack = newStack ;
    } else {
      dpcurStack->stackTop = dpcurStackTop ;
      dpcurStack = dpcurStack->next ;
    }
    dpcurStackBottom = dpcurStack->stackBottom ;
    dpcurStackTop    = dpcurStackBottom ;
    dpcurStackLimit  = dpcurStack->stackLimit ;
  }
}
/* Used before POPping.
 * Resets the LOOKing position if it was active.
 * Checks that the current block is not empty.
 * If it is, put the pointer back to the previous block. */
void dpcheckBack() {
  if (dplookStack) dplookStack = NULL ;
  if (dpcurStackTop == dpcurStackBottom) {
    dpcurStack->stackTop = dpcurStackBottom ;
    dpcurStack = dpcurStack->prev ;
    dpcurStackBottom = dpcurStack->stackBottom ;
    dpcurStackTop    = dpcurStack->stackTop ;
    dpcurStackLimit  = dpcurStack->stackLimit ;
  }
}
/* Used before LOOKing.
 * Activates the LOOKing position if it was reset.
 * Checks that the current LOOKing block is not empty.
 * If it is, put the LOOK pointer back to the previous block. */
void dpcheckLookBack() {
  if (dplookStack == NULL) {
    dplookStack = dpcurStack ;
    dplookStackBottom = dpcurStackBottom ;
    dplookStackTop = dpcurStackTop ;
    dplookStackLimit = dpcurStackLimit ;
  }
  if (dplookStackTop == dplookStackBottom) {
    dplookStack->stackTop = dplookStackBottom ;
    dplookStack = dplookStack->prev ;
    dplookStackBottom = dplookStack->stackBottom ;
    dplookStackTop    = dplookStack->stackTop ;
    dplookStackLimit  = dplookStack->stackLimit ;
  }
}

/* PUSHes "nbChars" consecutive chars,
 * from a location starting at address "x".
 * nbChars is assumed no larger than CHUNK_SIZE */
void dppushN(char *x, unsigned int nbChars) {
  dpcheck(nbChars) ;
  memcpy(dpcurStackTop,x,nbChars);
  dpcurStackTop+=nbChars ;
}
/* POPs "nbChars" consecutive chars,
 * to a location starting at address "x".
 * nbChars is assumed no larger than CHUNK_SIZE */
void dppopN(char *x, unsigned int nbChars) {
  dpcheckBack() ;
  if (dpcurStackTop-dpcurStackBottom < nbChars) {
    printf("Help! DP stack corrupted, height:%li !!",dpcurStackTop-dpcurStackBottom-nbChars) ;
    exit(0) ;
  }
  dpcurStackTop-=nbChars ;
  memcpy(x,dpcurStackTop,nbChars);
}
/* LOOKs "nbChars" consecutive chars,
 * to a location starting at address "x".
 * LOOKing is just like POPping, except that the main pointer
 * remains in place, so that the value is not POPped.
 * Further PUSHs or POPs will start from the same place as if
 * no LOOK had been made.
 * nbChars is assumed no larger than CHUNK_SIZE */
void dplookN(char *x, unsigned int nbChars) {
  dpcheckLookBack() ;
  dplookStackTop-=nbChars ;
  memcpy(x,dplookStackTop,nbChars);
}

/* PUSHes a large number "n" of consecutive chars,
 * from a location starting at address "x".
 * This "n"-sized array is cut into pieces no larger than CHUNK_SIZE */
void dppushNarray(char *x, unsigned int n) {
  unsigned int tailSize = n%CHUNK_SIZE ;
  char *xmax = x+n-tailSize ;
  char *xin = x ;
  while(xin<xmax) {
    dppushN(xin, CHUNK_SIZE) ;
    xin += CHUNK_SIZE ;
  }
  if (tailSize>0) dppushN(xin, tailSize) ;
}
/* POPs a large number "n" of consecutive chars,
 * to a location starting at address "x".
 * This "n"-sized array is cut into pieces no larger than CHUNK_SIZE */
void dppopNarray(char *x, unsigned int n) {
  unsigned int tailSize = n%CHUNK_SIZE ;
  char *xin = x+n-tailSize ;
  if (tailSize>0) dppopN(xin, tailSize) ;
  while(xin>x) {
    xin -= CHUNK_SIZE ;
    dppopN(xin, CHUNK_SIZE) ;
  }
}
/* LOOKs a large number "n" of consecutive chars,
 * to a location starting at address "x".
 * This "n"-sized array is cut into pieces no larger than CHUNK_SIZE */
void dplookNarray(char *x, unsigned int n) {
  unsigned int tailSize = n%CHUNK_SIZE ;
  char *xin = x+n-tailSize ;
  if (tailSize>0) dplookN(xin, tailSize) ;
  while(xin>x) {
    xin -= CHUNK_SIZE ;
    dplookN(xin, CHUNK_SIZE) ;
  }
}

/********** Exported PUSH/POP/LOOK functions : ************/

void dppushcharacter_(char *x) {
  dppushN(x,1) ;
}
void dppopcharacter_(char *x) {
  dppopN(x,1) ;
}
void dplookcharacter_(char *x) {
  dplookN(x,1) ;
}

void dppushboolean_(char *x) {
  dppushN(x,4) ;
}
void dppopboolean_(char *x) {
  dppopN(x,4) ;
}
void dplookboolean_(char *x) {
  dplookN(x,4) ;
}

void dppushinteger4_(char *x) {
  dppushN(x,4) ;
}
void dppopinteger4_(char *x) {
  dppopN(x,4) ;
}
void dplookinteger4_(char *x) {
  dplookN(x,4) ;
}

void dppushinteger8_(char *x) {
  dppushN(x,8) ;
}
void dppopinteger8_(char *x) {
  dppopN(x,8) ;
}
void dplookinteger8_(char *x) {
  dplookN(x,8) ;
}

void dppushinteger16_(char *x) {
  dppushN(x,16) ;
}
void dppopinteger16_(char *x) {
  dppopN(x,16) ;
}
void dplookinteger16_(char *x) {
  dplookN(x,16) ;
}

void dppushreal4_(char *x) {
  dppushN(x,4) ;
}
void dppopreal4_(char *x) {
  dppopN(x,4) ;
}
void dplookreal4_(char *x) {
  dplookN(x,4) ;
}

void dppushreal8_(char *x) {
  dppushN(x,8) ;
}
void dppopreal8_(char *x) {
  dppopN(x,8) ;
}
void dplookreal8_(char *x) {
  dplookN(x,8) ;
}

void dppushreal16_(char *x) {
  dppushN(x,16) ;
}
void dppopreal16_(char *x) {
  dppopN(x,16) ;
}
void dplookreal16_(char *x) {
  dplookN(x,16) ;
}

void dppushreal32_(char *x) {
  dppushN(x,32) ;
}
void dppopreal32_(char *x) {
  dppopN(x,32) ;
}
void dplookreal32_(char *x) {
  dplookN(x,32) ;
}

void dppushcomplex4_(char *x) {
  dppushN(x,4) ;
}
void dppopcomplex4_(char *x) {
  dppopN(x,4) ;
}
void dplookcomplex4_(char *x) {
  dplookN(x,4) ;
}

void dppushcomplex8_(char *x) {
  dppushN(x,8) ;
}
void dppopcomplex8_(char *x) {
  dppopN(x,8) ;
}
void dplookcomplex8_(char *x) {
  dplookN(x,8) ;
}

void dppushcomplex16_(char *x) {
  dppushN(x,16) ;
}
void dppopcomplex16_(char *x) {
  dppopN(x,16) ;
}
void dplookcomplex16_(char *x) {
  dplookN(x,16) ;
}

void dppushcomplex32_(char *x) {
  dppushN(x,32) ;
}
void dppopcomplex32_(char *x) {
  dppopN(x,32) ;
}
void dplookcomplex32_(char *x) {
  dplookN(x,32) ;
}

/******************** The same for arrays: ****************/

void dppushcharacterarray_(char *x, unsigned int *n) {
  dppushNarray(x,*n) ;
}
void dppopcharacterarray_(char *x, unsigned int *n) {
  dppopNarray(x,*n) ;
}
void dplookcharacterarray_(char *x, unsigned int *n) {
  dplookNarray(x,*n) ;
}

void dppushbooleanarray_(char *x, unsigned int *n) {
  dppushNarray(x,(*n*4)) ;
}
void dppopbooleanarray_(char *x, unsigned int *n) {
  dppopNarray(x,(*n*4)) ;
}
void dplookbooleanarray_(char *x, unsigned int *n) {
  dplookNarray(x,(*n*4)) ;
}

void dppushinteger4array_(char *x, unsigned int *n) {
  dppushNarray(x,(*n*4)) ;
}
void dppopinteger4array_(char *x, unsigned int *n) {
  dppopNarray(x,(*n*4)) ;
}
void dplookinteger4array_(char *x, unsigned int *n) {
  dplookNarray(x,(*n*4)) ;
}

void dppushinteger8array_(char *x, unsigned int *n) {
  dppushNarray(x,(*n*8)) ;
}
void dppopinteger8array_(char *x, unsigned int *n) {
  dppopNarray(x,(*n*8)) ;
}
void dplookinteger8array_(char *x, unsigned int *n) {
  dplookNarray(x,(*n*8)) ;
}

void dppushinteger16array_(char *x, unsigned int *n) {
  dppushNarray(x,(*n*16)) ;
}
void dppopinteger16array_(char *x, unsigned int *n) {
  dppopNarray(x,(*n*16)) ;
}
void dplookinteger16array_(char *x, unsigned int *n) {
  dplookNarray(x,(*n*16)) ;
}

void dppushreal4array_(char *x, unsigned int *n) {
  dppushNarray(x,(*n*4)) ;
}
void dppopreal4array_(char *x, unsigned int *n) {
  dppopNarray(x,(*n*4)) ;
}
void dplookreal4array_(char *x, unsigned int *n) {
  dplookNarray(x,(*n*4)) ;
}

void dppushreal8array_(char *x, unsigned int *n) {
  dppushNarray(x,(*n*8)) ;
}
void dppopreal8array_(char *x, unsigned int *n) {
  dppopNarray(x,(*n*8)) ;
}
void dplookreal8array_(char *x, unsigned int *n) {
  dplookNarray(x,(*n*8)) ;
}

void dppushreal16array_(char *x, unsigned int *n) {
  dppushNarray(x,(*n*16)) ;
}
void dppopreal16array_(char *x, unsigned int *n) {
  dppopNarray(x,(*n*16)) ;
}
void dplookreal16array_(char *x, unsigned int *n) {
  dplookNarray(x,(*n*16)) ;
}

void dppushreal32array_(char *x, unsigned int *n) {
  dppushNarray(x,(*n*32)) ;
}
void dppopreal32array_(char *x, unsigned int *n) {
  dppopNarray(x,(*n*32)) ;
}
void dplookreal32array_(char *x, unsigned int *n) {
  dplookNarray(x,(*n*32)) ;
}

void dppushcomplex4array_(char *x, unsigned int *n) {
  dppushNarray(x,(*n*4)) ;
}
void dppopcomplex4array_(char *x, unsigned int *n) {
  dppopNarray(x,(*n*4)) ;
}
void dplookcomplex4array_(char *x, unsigned int *n) {
  dplookNarray(x,(*n*4)) ;
}

void dppushcomplex8array_(char *x, unsigned int *n) {
  dppushNarray(x,(*n*8)) ;
}
void dppopcomplex8array_(char *x, unsigned int *n) {
  dppopNarray(x,(*n*8)) ;
}
void dplookcomplex8array_(char *x, unsigned int *n) {
  dplookNarray(x,(*n*8)) ;
}

void dppushcomplex16array_(char *x, unsigned int *n) {
  dppushNarray(x,(*n*16)) ;
}
void dppopcomplex16array_(char *x, unsigned int *n) {
  dppopNarray(x,(*n*16)) ;
}
void dplookcomplex16array_(char *x, unsigned int *n) {
  dplookNarray(x,(*n*16)) ;
}

void dppushcomplex32array_(char *x, unsigned int *n) {
  dppushNarray(x,(*n*32)) ;
}
void dppopcomplex32array_(char *x, unsigned int *n) {
  dppopNarray(x,(*n*32)) ;
}
void dplookcomplex32array_(char *x, unsigned int *n) {
  dplookNarray(x,(*n*32)) ;
}

/************* Debug displays of the state of the stack: ***********/

void dpprinttopplace_() {
    DoubleChainedBlock *stack = dpcurStack ;
    int nbBlocks = 0 ;
    while(stack) {
	stack = stack->prev ;
	nbBlocks++ ;
    }
    printf("DP Stack  top: %i+%li\n",nbBlocks,dpcurStackTop - dpcurStackBottom) ;
}

void dpprintlookingplace_() {
    if (dplookStack == NULL)
	dpprinttopplace_() ;
    else {
	DoubleChainedBlock *stack = dplookStack ;
	unsigned int nbBlocks = 0 ;
	while(stack) {
	    stack = stack->prev ;
	    nbBlocks++ ;
	}
	printf("DP Stack look: %i+%li\n",nbBlocks,dplookStackTop - dplookStackBottom) ;
    }
}

/**** Experiment for a numeric sum that does the sum from smallest to largest ******/

static double sumlist[3000] ;
static int sumindex = 0 ;

void dpsumr8reset_() {
  int i ;
  sumindex = 0 ;
  for (i=0 ; i<3000 ; i++)
    sumlist[i] = 0.0 ;
}

void dpsumr8_(double *x) {
  double absnumber ;
  int decali ;
  int inserti = 0 ;
  double number = *x ;
  if (number>0.0 || number<0.0) {
    if (number>0.0)
      absnumber = number ;
    else
      absnumber = -number ;
    while (inserti<sumindex &&
           absnumber > ((sumlist[inserti]>0)?sumlist[inserti]:-sumlist[inserti]))
      inserti++ ;
    if (sumindex>=3000) {
      printf("Help! sum list full!\n") ;
      exit(0) ;
    } else {
      for (decali=sumindex;decali>inserti;decali--)
        sumlist[decali] = sumlist[decali-1] ;
      sumlist[inserti]=number ;
      sumindex++ ;
    }
  }
}

void dpsumr8compute_() {
  int i ;
  double sum = 0.0 ;
  for (i=0 ; i<sumindex ; i++) {
    printf(" %25.20e",sumlist[i]) ;
    sum += sumlist[i] ;
  }
  printf("\n\n Sum value: %25.20e\n", sum) ;
}
