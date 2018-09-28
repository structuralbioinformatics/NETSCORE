#include "ppin.h"

void PrintShortList()
{
    printf("\nPARAMETERS\n");
    printf("\
    \n\t -i \t File with interactions and affinities \
    \n\t -c \t File with cell abundances and seed scores\
    \n\t -l \t Local Distribution of proteins \
    \n\t -x \t Flag to exclude network analysis and do only phenotype scores without network\
    \n\t -o \t Output File \
    \n\t -h \t Print help \
    \n\t -hf \t Print full list of options \
    \n\n> Parameters for Analysis of  the network: \
    \n  (NetScores are  scores of Nodes and Edges modified by the network)\
    \n\t -mth\t Method of network-base prioritization:\
    \n\t     \t   0 NetShort (default) \
    \n\t     \t   1 NetScore \
    \n\t     \t   2 NetZscore with Static Random Networks \
    \n\t     \t   3 NetZscore with Dynamic Random Networks \
    \n\t     \t   4 NetCombo Static \
    \n\t     \t   5 NetCombo Dynamic \
    \n\t -n \t Number of iterations (if n=0 the flow of information is null)\
    \n\t -ns\t Number of iterations only for NetScore when using NetCombo (default is n)\
    \n\t -nz\t Number of iterations only for NetZscore (static or dynamic) when using NetCombo (default is n)\
    \n\t -nh\t Number of iterations only for NetShort when using NetCombo (default is n)\
    \n\t -dxi\t Use indirect links that can substitute direct neighbours\
    \n\t     \t   0 adds all scores of indirect links (default)\
    \n\t     \t  >0 add highest scores of indirect links but not all \
    \n\t     \t    (use it for large # of seeds or when using PhenoScores)\
    \n\n> Parameters for analysis of expression data: \
    \n  (PhenoScores are scores of nodes according to phenotype relevance)  \
    \n\t -xphe\t phenotype file with the name of the files of samples and expression data \
    \n\t -xnot\t file with array annotations for protein-nodes \
    \n\t -xphs\t phenotype filename used as reference to calculate the scores \
    \n\t -xiter\t number of iterations to modify phenoscores with netscores (default=0)\
    \n\t -xcsg\t flag to skip singleton genes (unlinked on the clustering) to calculate phenoscore\
    \n\n");

    exit(1);
    //obsolete - \n\t -xana\t  flag to performe comparative analyses (default off) - \

}
