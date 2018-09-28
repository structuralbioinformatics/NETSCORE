#include "ppin.h"

void PrintNetScore(OUT,size,edges,protein,interaction,nodepval,nodescr,average_nodescr,rms_nodescr,edgescr,average_edgescr,rms_edgescr,maxmin)
FILE  *OUT;
int   size,edges;
node      *protein;
edge      *interaction;
float *nodepval,*nodescr,*edgescr,*average_nodescr,*average_edgescr,*rms_nodescr,*rms_edgescr,*maxmin;
{

 int k,iout;

  fprintf(OUT,"\nFINAL SCORES \n");
  fprintf(OUT,"NODES\n");
  for(k=0;k<size;k++){
     fprintf(OUT,"%10d\t%10s [%10.5f] => [%10.5f]\tE-value %8.3e\n",k,protein[k].name1,protein[k].copy.score,(nodescr[k]+protein[k].copy.score),nodepval[k]);
   }
  if (maxmin[3]>0 ){
  fprintf(OUT,"EDGES\n");
  for(k=0;k<edges;k++){
     fprintf(OUT,"%10d\t%10s interact %10s [%10.5f] => [%10.5f]\n",k,interaction[k].a.name1,interaction[k].b.name1,interaction[k].association,(edgescr[k]+interaction[k].association));
   }}


  fprintf(OUT,"\nAVERAGE OF SCORES AND FLUCTUATION\n");
  fprintf(OUT,"NODES\n");
  for(k=0;k<size;k++){
    fprintf(OUT,"%10d\t%10s [%10.5f] => [%10.5f] +/- [%10.5f]\n",k,protein[k].name1,protein[k].copy.score,(average_nodescr[k]+protein[k].copy.score),rms_nodescr[k]);
  }
  if (maxmin[3]>0 ){
  fprintf(OUT,"EDGES\n");
  for(k=0;k<edges;k++){
     fprintf(OUT,"%10d\t%10s interact %10s [%10.5f] => [%10.5f] +/- [%10.5f]\n",k,interaction[k].a.name1,interaction[k].b.name1,interaction[k].association,(average_edgescr[k]+interaction[k].association),rms_edgescr[k]);
  }}
  iout=fflush(OUT);

}
