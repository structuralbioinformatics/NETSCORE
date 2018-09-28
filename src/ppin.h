#include        <stdio.h>
#include        <math.h>
#include        <ctype.h>
#include        <string.h>
#include        <sys/types.h>
#include        <sys/wait.h>
#include        <malloc.h>
#include        <stdlib.h>
#include        <unistd.h>
#include        <errno.h>
#include        <time.h>

#define         LOCAL   70
#define         MAXS    1024
#define         MAXP    65000
#define         MAXI    8000
#define         MAXID   5000
#define         MAXIDE  50
#define         MAXD    100
#define         MINSTR  1.0
#define         LOW     1.0e-100
#define         UPPER   1.0e+100
#define         MAXSCR  1.0
#define         MAXSIGMA 10
#define         MINSCR  0.0
#define         MINLINK 3  
#define         MIN     1.0e-5
#define         MAXWN   14500
#define         MAXWE   25000
#define         MAXPHE  3 
#define         MAXSAMPLE 45
#define         MAXXG   6000
#define         MAXXSG  8000
#define         MINXLINK 0.5
#define         MAXXDEV  0.5
#define         MINXCOMM 0.5
#define         MINXDEV  0.5
#define         MAXPARSCR 10
#define         MAXPARLNK  5 
#define         MAXPARRND  10
#define         MAXPARCLUS 20
#define         MAXNAME   25 
#define         MAXNAMELARGE   80 
#define         MAXERROR    5
#define         MAXMETHOD  25 
#define         MINXSIM     0.25
#define         MAXSTPG     20
#define         MAXSTPS     10
#define         MAXEST      350
#define         MINSHORT    0.01
#define         MAXGUILD    4
#define         NEXP        50
#define         NGUMB       50
#define         NGAUSS      50
#define         MAXHPDF     200
#define         MINHPDF     5
#define         MAXXPDF     50

typedef struct abundance
{
 float   local[LOCAL];
 float   rms;
 float   cell;
 float   score;
} abundance;

typedef struct node
{
 char        name1[MAXNAME];
 char        name2[MAXNAME];
 abundance   copy;
 int         interact[MAXI];
 int         degree;
} node;

typedef struct edge
{
 int     i;
 int     j;
 node    a;
 node    b;
 float   association;
} edge;

typedef struct gauss
{
 float  average;
 float  rmsd;
} gauss;

typedef struct gradient
{
 int    j;
 float  size;
} gradient;

typedef struct array
{
 float  sample[MAXSAMPLE];
 int    defined[MAXSAMPLE];
} array;

typedef struct probe
{
 array  phenotype[MAXPHE];
 char   name[MAXNAME];
} probe;

typedef struct express
{
 probe  fragment[MAXEST];
 int    split;
} express;

