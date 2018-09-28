# The program requires stacksize unlimited
# bash: ulimit -s unlimited
# csh:  limit  stacksize unlimited
C    = gcc
OPT=	-g -O3  -fexpensive-optimizations 
OBJS = netscore.o ReadInput.o ReadInteractions.o ReadNode.o allocate.o AssignNode.o \
ModifyDual.o ModifyNodeScore.o  ModifyNodeZScoreD.o ModifyNodeZScoreS.o \
ModifyScore.o ModifyZScoreD.o ModifyZScoreS.o ShortModify.o ShortModifyNode.o \
DualEdge.o DualNode.o Zscore.o InitZscore.o ZscorePopulation.o InitZscorePopulation.o \
InitScore.o InitScorePopulation.o InitZscoreS.o InitZscoreSPopulation.o \
ScorePopulation.o NetZScoreMotorS.o NetZScoreMotorD.o  ZscoreSPopulation.o\
ran0.o ran3.o CumulativeMatrix.o RandomSwapNodes.o RandomScores.o ScoreDistribution.o EvalueNode.o \
mrqmin.o mrqcof.o covsrt.o gaussj.o fgauss.o fsexp.o fgumbel.o ffgauss.o ffgumbel.o \
fdgumbel.o ffdgumbel.o qromo.o qtrap.o trapzd.o polint.o midexp.o midpnt.o midinf.o InitParameters.o\
NetScoreMotor.o NetShortMotor.o  NetPrioritization.o NetComboMotorD.o NetComboMotorS.o ZscoreS.o Score.o \
RandomCumulativeMatrix.o Normal.o InitMatrix.o Transform.o EdgeDistribution.o CmpMatrix.o XPhenoScore.o\
mod.o DualNodeSCR.o gdistribute.o LocalProduct.o ReadPhenotypeExpress.o  XAssignPhenoScore.o XAssignModifyPhenoScore.o \
PearsonScore.o XFilterMatrix.o  XGeneZscore.o XSampleZscore.o XGenePhenotype.o XGenePhenoScore.o XRankProd.o \
kstwo.o tutest.o ttest.o Dmedian.o sort.o avevar.o betai.o probks.o gammln.o betacf.o PrintList.o PrintShortList.o PrintNetScore.o \
XModifyPhenoScore.o Dijkstra.o Graph.o XAnalyzePhenoNet.o PrintGaussian.o PearsonScoreMatrix.o FloydWarshall.o \
PearsonScoreMatrix2i.o XAnalyzeCluster.o chstwo.o gammq.o gcf.o gser.o BinChi.o BinCramer.o BinMI.o Bin2Cramer.o Bin2MI.o\
WGraph.o WDijkstra.o DummyInteractions.o PrintCorrelation.o WFloydWarshall.o \
Xcluster.o NormalizeCLR.o NormalizeAPC.o NormalizeASC.o SelectNbins.o SelectNbins2.o cntab1.o cntab2.o Clustering.o ReadClusterFile.o


LIBS =  -lm -fstack-protector-all 
  
.SUFFIXES: .c .o
.c.o:
	$(C) -c -O3  -fexpensive-optimizations $*.c 

all: $(OBJS) 
	$(C) -o netscore  $(OBJS)  $(LIBS) $(OPT)
clean:
	\rm *.o 

