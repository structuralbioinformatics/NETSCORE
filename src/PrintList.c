#include "ppin.h"

void PrintList()
{

    printf("\nBASIC PARAMETERS\n");
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
    \n\n");
    printf("============================\n");
    printf("\nADDITIONAL PARAMETERS\n");
    printf("\
    \n\t -stop\t Flag to stop and check gene-gene correlation values before clustering (it requires flag -xprn)\
    \n\n> Parameters for Analysis of  the network: \
    \n  (NetScores are  scores of Nodes and Edges modified by the network)\
    \n\t -t \t Threshold of minimum input interaction-score  (default 0) \
    \n\t -z \t Threshold of minimum output interaction-score (default 0)  \
    \n\t -nr\t Number of random Networks to calculate Z-score (default 100) \
    \n\t -r \t Random type \
    \n\t     \t   0 Full random \
    \n\t     \t   1 Conserving network topology (default)   \
    \n\t -zp\t Population to calculate Zscore \
    \n\t     \t   0 Score vs NR random networks for all nodes (default)\
    \n\t     \t   1 Score vs population of NR random networks specific for each node \
    \n\t -ncz\t Netcombo Z-score type \
    \n\t     \t   0 Use scores of all nodes (default)\
    \n\t     \t   1 Use only the scores of nodes affected or modified \
    \n\t -mpv\t Calculation of Z-score p-values (default is 2)\
    \n\t     \t   0 Neglect calculation \
    \n\t     \t   1 Use #NGUMB Gumbell functions to model the Cumulative Distribution Function (CDF) of random scores \
    \n\t     \t   2 use #NGUMB Gumbell functions to model the Density of Probability Function (PDF) of random scores\
    \n\t     \t   3 Use #NGAUSS Gaussian functions to model the Cumulative Distribution Function (CDF) of random scores \
    \n\t     \t   4 use #NGAUSS Gaussian functions to model the Density of Probability Function (PDF) of random scores\
    \n\t -npv\t number of repetitions to model random scores (default is 10 using theoretical rms 0.01) \
    \n\t -nfpv1\t starting number of functions (exponential or gaussian) to fit (default 1)\
    \n\t -nfpv2\t maximum number of functions (exponential or gaussian) to fit (check also NGUMB and NGAUSS)(default 3)\
    \n\t -nd\t number of iterations to recalculate dual nodes/edges \
    \n\t -mx\t maximum  increase of node score (default MAXSCR=1)\
    \n\t -mxe\t maximum  increase of edge score (default 0)\
    \n\t -ms\t maximum  number of standard deviaton for increasing score (default MAXSIGMA=10)\
    \n\t -mn\t minimum  to increase node score (default MINSCR=0)\
    \n\t -mne\t minimum  to increase edge score (default MINSCR=0)\
    \n\t -mnd\t minimum  to transfer the increase of node-score to dual edge-score (default MINSCR=0)\
    \n\t -mnde\t minimum  to transfer the increase of edge-score to node-score (default MINSCR=0)\
    \n\t -mnst\t minimum number of nodes grouped under the same degree (default MINSTR=0)\
    \n\t -mnste\t minimum number of edges grouped under the same degree (default MINSTR=0)\
    \n\t -dn\t node-degree to weight nodes  (default MINLINK=3)\
    \n\t -de\t edge-degree to weight edges (default MINLINK=3) \
    \n\t -dxn\t degree to weight node score with dual transfer(default 0, means no increase) \
    \n\t -dxe\t degree to weight edge score with dual transfer(default 0, means no increase) \
    \n\t -dxi\t Use indirect links that can substitute direct neighbours\
    \n\t     \t   0 adds all scores of indirect links (default)\
    \n\t     \t  >0 add highest scores of indirect links but not all (use it for large # of seeds or when using PhenoScores)\
    \n\t -dnn\t weight increase of scores in nodes and edges by its degree \
    \n\t     \t   0 no normalization (default)\
    \n\t     \t   1 normalizes increase by k (for netscore and netzcore)\
    \n\t -e \t tolerance to stop the iteration in netscores (default 1.0e-10)\
    \n\n> parameters for analysis of expression data: \
    \n  (phenoscores are scores of nodes according to phenotype relevance)  \
    \n\t -xcp [0..5]\t clustering protocol \
    \n\t          0: density search clustering\
    \n\t          1: markov cluster (mcl default)\
    \n\t          2: hierarchical clustering\
    \n\t          3: k-means clustering\
    \n\t          4: som clustering\
    \n\t          5: read clusters from external file\
    \n\t\t* density search clustering:\
    \n\t\t -xcl\t minimum score to link a putative new member to a cluster \
    \n\t\t -xca\t maximum deviation of the mean produced by the inclusion of a new member \
    \n\t\t -xcr\t minimum ratio of common elements with respect to the smallest group for merging two clusters\
    \n\t\t -xcm\t maximum ratio of deviation of the means when merging two clusters \
    \n\t\t -xcma\t minimum mean cut-off to stop cluster growth \
    \n\t\t* mcl clustering: \
    \n\t\t -xci\t inflation \
    \n\t\t -xcw\t threshold of distance similarity to generate graph (default 0.2 or 5.0 if xng flag is ON)\
    \n\t\t* hierarchical clustering: \
    \n\t\t -xhm [0..3] specifies which hierarchical clustering method to use \
    \n\t\t        0: pairwise complete-linkage (m)\
    \n\t\t        1: pairwise single-linkage (s)\
    \n\t\t        2: pairwise centroid-linkage (c)\
    \n\t\t        3: pairwise average-linkage (a)\
    \n\t\t        (default: 0) \
    \n\t\t -xcn\t number of clusters to use \
    \n\t\t* k-means clustering: \
    \n\t\t -xcn\t number of clusters to use \
    \n\t\t -xkr\t number of times the k-means clustering algorithm is run \
    \n\t\t* som clustering: \
    \n\t\t -xxn\t specifies the horizontal dimension of the som grid (default: 2) \
    \n\t\t -xyn\t specifies the vertical dimension of the som grid (default: 1) \
    \n\t\t* read external file with clusters \
    \n\t\t -xcf\t file with clusters \
    \n\t -xst [0..10] statistical test to assign genes to a phenotype\
    \n\t     \t   0 kolmogorov-smirnov (default)\
    \n\t     \t   1 student t-test of populations with unequal variances\
    \n\t     \t   2 student t-test of populations with the same variances\
    \n\t     \t   3 ratio: |median(phenotype)-median(no-phenotype)|/|max-min|    \
    \n\t     \t   4 chi square test (recommended for large number of samples or large clusters) \
    \n\t     \t   5 cramer's v (correlation of expression vs phenotype)  \
    \n\t     \t   6 contingency coefficient c (correlation of expression vs phenotype) \
    \n\t     \t   7 symmetric mutual information (correlation of expression vs phenotype) (p-value by chi-square)\
    \n\t     \t   8 1 - cramer's v (expression phenotype_1 vs expression phenotype_2, correspondence of samples is required)  \
    \n\t     \t   9 1 - contingency coefficient c (expression phenotype_1 vs expression phenotype_2, correspondence of samples is required) \
    \n\t     \t  10 1 - symmetric mutual information (expression phenotype_1 vs expression phenotype_2, correspondence of samples is required)(p-value by chi-square) \
    \n\t -xsts [0..3] phenoscore strength (default is 2)\
    \n\t     \t   0 strength= 1 -  p-value \
    \n\t     \t   1 strength= statistic hypothesis (ks, t, cramer v, ccc, mi) \
    \n\t     \t   2 strength= (1 -  p-value) x (statistic hypothesis)\
    \n\t     \t   3 strength= delta-dirac(p-value < 0.05) x (statistic hypothesis) \
    \n\t     \t   4 strength= delta-dirac(p-value < 0.05) \
    \n\t -xsto [0..2] use of seed-scores (default is 2)\
    \n\t     \t   0 score= score_seed * strength \
    \n\t     \t   1 score= (score_seed + strength) \
    \n\t     \t   2 score= score_seed * (1 + strength) \
    \n\t -xstg [0..3] phenoscores calculation\
    \n\t     \t   0 compare one to all phenotypes and use individual genes: each gene has a score (default)\
    \n\t     \t   1 compare one to all phenotypes and use cluster genes: all genes of a cluster have the same score\
    \n\t     \t   2 use maximum p-value of one to one comparisons by individual genes\
    \n\t     \t   3 use maximum p-value of one to one comparisons by cluster genes\
    \n\t     \t   4 RankProd: compare one to all phenotypes using individual genes\
    \n\t -xrkpp\t P-value threshold to select seeds (default 0.05)\
    \n\t -xrkpf\t Fold threshold to select seeds with overexpression using RankProd (symmetric threshold is chosen for underexpression, default 1.5)\
    \n\t -xrkpt\t Top ranking genes (up/down regulated) according to the lowest RankProd values (default 100, if 0 then use all)\
    \n\t -xsign [0..3] use of positive/negative distance measure for gene clustering \
    \n\t     \t   0 only positive distance is used to compare/cluster genes (default)\
    \n\t     \t   1 clustering uses absolute values of distance if only if netscore has positive increase\
    \n\t     \t     (do not use with option -xstg=1)\
    \n\t     \t   2 activity (expression zscore) is positive in the largest set of the reference-phenotype \
    \n\t     \t     (it forces merging over & under expression when clustering after first iteration)\
    \n\t     \t   3 always use absolute values of distance \
    \n\t -xedge [0..6] edge scoring based on expression \
    \n\t     \t   0 no modification (default)\
    \n\t     \t   1 similar expression: weighted by phenoscore (kolmogorov-smirnov) of similar expression distribution\
    \n\t     \t   2 similar expression: weighted by phenoscore (student t-test unequal variances) of similar expression distribution\
    \n\t     \t   3 similar expression: weighted by phenoscore (student t-test same variances) of similar expression distribution\
    \n\t     \t   4 correlated expression: weighted by absolute value of distance measure of expression (pce) \
    \n\t     \t   5 positive correlation: weighted by (1+pce)/2 \
    \n\t     \t   6 anti-correlation: weighted by (1-pce)/2 \
    \n\t -xsim [0..2]  matrix of similarities between protein nodes calculated with their est fragments (probes)\
    \n\t     \t   0 average of probe-probe similarities (default)\
    \n\t     \t   1 worst of probe-probe similarities (close to null)\
    \n\t     \t   2 best of probe-probe similarities (close to 1 or -1)\
    \n\t -xmx\t top percentage of scores from  netscore/netshort to modify distance similarities and recluster (sefault 50) \
    \n\t -xtol\t tolerance to stop the iteration in phenoscores (default 1.0e-8)\
    \n\t -xcorr [1..11]     specifies the distance measure (correlation) for gene clustering\
    \n\t          1: uncentered correlation\
    \n\t          2: pearson correlation\
    \n\t          3: uncentered correlation, absolute value\
    \n\t          4: pearson correlation, absolute value\
    \n\t          5: spearman's rank correlation\
    \n\t          6: kendall's tau\
    \n\t          7: euclidean distance\
    \n\t          8: city-block distance\
    \n\t          9: symmetric mutual information \
    \n\t         10: cramer's v         \
    \n\t         11: contingency coefficient c      \
    \n\t          (default: 2)\
    \n\t -xpath [0..2] weighting the distance measure between genes with the shortest-path distance\
    \n\t          0: no weighting (default)\
    \n\t          1: weight is: power(wp,(1 - shortest_path)), default is wp=1\
    \n\t          2: weight is: 0 (if shortest_path > wp) and 1 otherwise, default is wp=100\
    \n\t -xwp\t parameter wp if xpath is not null \
    \n\t -xcg [0..2] specifies whether to center each row (gene) in the data \
    \n\t          0: no centering (default) \
    \n\t          1: subtract the mean of each row (a) \
    \n\t          2: subtract the median of each row (m) \
    \n\t -xncr [0..3]  method to remove the backround noise of the correlation between genes  \
    \n\t          0: no modification (default) \
    \n\t          1: method clr by sqrt(z1*z1 + z2*z2) \
    \n\t          2: method apc \
    \n\t          3: method asc \
    \n\t -xlog\t flag to specify to log-transform the data before clustering \
    \n\t -xng\t flag to specify to normalize each row (gene) in the data \
    \n\t -xcsg\t flag to skip singleton genes (unlinked on the clustering) to calculate phenoscore \
    \n\t -xprn\t flag to verbose output with information on clustering\
    \n\t -xana\t flag to analyze the clusters on the fly along iterations\
    \n\n");

    exit(1);
    //obsolete - \n\t -xana\t  flag to performe comparative analyses (default off) - \

}
