#include "ppin.h"

void   RandomScores(edges,size2,interaction,protein,nodescr,average_nodescr,rms_nodescr,n_nodescr,edgescr,average_edgescr,rms_edgescr,n_edgescr,maxmin,distribution,ndstr,linker,rnd,tolerance,threshold,thero,iteration,diter)
int       **distribution,*ndstr;
int       *size2;
int       edges;
node      *protein;
edge      *interaction;
int       *n_edgescr,*n_nodescr;
float     tolerance,threshold,thero;
int       *iteration,diter;
float     maxmin[];
int       rnd[],linker[];
float    *nodescr,*edgescr,*average_nodescr,*average_edgescr,*rms_nodescr,*rms_edgescr;
{
 node      *rndprotein, *avector();
 int        i,j,n,size,intact,iteration_global, iteration_netscore, iteration_netzscore, iteration_netshort,nonswap; 
 int       **swap, **imatrix();     
 float     *error;
 int       *rnd2,*ivector();
 float     *NetScoreMotor(),*NetShortMotor(),*NetZScoreMotorS(),*NetZScoreMotorD(),*NetComboMotorD(),*NetComboMotorS();
 float     *fvector();
 void       free_fvector(),free_ivector(), free_avector(),free_imatrix(),free_matrix();

   size=size2[0];
   intact=size2[1];

   error=fvector(0,MAXERROR);
   swap=imatrix(0,size,0,1);
   rndprotein=avector(0,size);
   rnd2=ivector(0,MAXPARRND);
   for (i=0;i<MAXPARRND;i++)rnd2[i]=rnd[i];
   rnd2[0]=1;
   printf("Random Swap Nodes for E-value calculation\n");
   RandomSwapNodes(size2,interaction,protein,distribution,ndstr,MAXID,rnd2,swap);
   nonswap=0;
   for (i=0;i<size;i++){ 
     rndprotein[i]=protein[i];
     rndprotein[i].copy.score=protein[swap[i][0]].copy.score;
     if (protein[i].copy.score>0.5 && abs(rndprotein[i].copy.score-protein[i].copy.score) < 0.1) nonswap++;
     }
   if (nonswap>0) printf("Total of unchanged scores > 0.5 : %d \n",nonswap);

   switch(rnd[4]){
    case 0:
     iteration_netshort=iteration[0];
     error=NetShortMotor(edges,size2,interaction,rndprotein,nodescr,average_nodescr,rms_nodescr,n_nodescr,edgescr,average_edgescr,rms_edgescr,n_edgescr,maxmin,distribution,ndstr,linker,rnd,tolerance,threshold,thero,iteration_netshort,diter);
     break;
    case 1:
     iteration_netscore=iteration[0];
     error=NetScoreMotor(edges,size2,interaction,rndprotein,nodescr,average_nodescr,rms_nodescr,n_nodescr,edgescr,average_edgescr,rms_edgescr,n_edgescr,maxmin,distribution,ndstr,linker,rnd,tolerance,threshold,thero,iteration_netscore,diter);
     break;
    case 2:
     iteration_netzscore=iteration[0];
     error=NetZScoreMotorS(edges,size2,interaction,rndprotein,nodescr,average_nodescr,rms_nodescr,n_nodescr,edgescr,average_edgescr,rms_edgescr,n_edgescr,maxmin,distribution,ndstr,linker,rnd,tolerance,threshold,thero,iteration_netzscore,diter);
     break;
    case 3:
     iteration_netzscore=iteration[0];
     error=NetZScoreMotorD(edges,size2,interaction,rndprotein,nodescr,average_nodescr,rms_nodescr,n_nodescr,edgescr,average_edgescr,rms_edgescr,n_edgescr,maxmin,distribution,ndstr,linker,rnd,tolerance,threshold,thero,iteration_netzscore,diter);
     break;
    case 4:
     error=NetComboMotorS(edges,size2,interaction,rndprotein,nodescr,average_nodescr,rms_nodescr,n_nodescr,edgescr,average_edgescr,rms_edgescr,n_edgescr,maxmin,distribution,ndstr,linker,rnd,tolerance,threshold,thero,iteration,diter);
     break;
    case 5:
     error=NetComboMotorD(edges,size2,interaction,rndprotein,nodescr,average_nodescr,rms_nodescr,n_nodescr,edgescr,average_edgescr,rms_edgescr,n_edgescr,maxmin,distribution,ndstr,linker,rnd,tolerance,threshold,thero,iteration,diter);
     break;
    default:
     iteration_netshort=iteration[0];
     error=NetShortMotor(edges,size2,interaction,rndprotein,nodescr,average_nodescr,rms_nodescr,n_nodescr,edgescr,average_edgescr,rms_edgescr,n_edgescr,maxmin,distribution,ndstr,linker,rnd,tolerance,threshold,thero,iteration_netshort,diter);
   }
   free_fvector(error,0,MAXERROR);
   free_avector(rndprotein,0,size);
   free_ivector(rnd2,0,MAXPARRND);
   free_imatrix(swap,0,size,0,1);
   
}
