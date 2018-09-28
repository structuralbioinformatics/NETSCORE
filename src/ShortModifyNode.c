#include "ppin.h"
#define eps 1.0e-6

float  ShortModifyNode(iter,dnode,dedge,size,edges,linker,scr,escr,maxmin)
int       iter,size,edges,linker;
float    *scr,*escr,*maxmin;
node     *dnode;
edge     *dedge;
{
  int   i,j,k,kk;
  float error,dscr,scri,qu,qu2,squ,max,min,sigma,weight;
  float **weight_adjacent,**dist;
  float *fvector(),**matrix();
  void  free_ivector(),free_fvector(),free_matrix();
  void  WGraph(), WDijkstra(),WFloydWarshall();
  clock_t   launch, done;
  float  diff;

  max=maxmin[0];
  min=maxmin[1];
  sigma=maxmin[2];

//  printf("Subroutine ModifyNode\n");
  weight_adjacent=matrix(0,size,0,size);
  dist=matrix(0,size,0,size);
  error=0.0;
  if (size==0){printf("Error in ModifyNode\n");return error;}
  if (linker==0){linker=MINLINK;} 
  printf("Do Graph\n");
  launch = clock();
  WGraph(weight_adjacent,dnode,dedge,scr,escr,edges,size);
  done = clock();
  diff = 1000*(done - launch) / CLOCKS_PER_SEC;
  printf("Recalculating the Graph with scored distance. Time %e ms\n",diff);
 // for (i=0;i<size;i++){WDijkstra(weight_adjacent,size,dist[i],i); }
  for (i=0;i<size;i++) memcpy(dist[i],weight_adjacent[i],size*sizeof(float));
  printf("Do shortest path\n");
  launch = clock();
  WFloydWarshall(dist,size);
  done = clock();
  diff = (done - launch) / CLOCKS_PER_SEC;
  printf("Running FloydWarshall with weighted edges. Time %f sec\n",diff);
  dscr=0.0;
  qu=qu2=0.0;
  for (i=0;i<size;i++){
       scri=0;
       for (j=0;j<size;j++){
       if (i!=j){
         if (dist[i][j]>eps)scri+= 1./dist[i][j];
       }}
       dscr=scri;
       qu+=dscr;
       qu2+=dscr*dscr;
     }
  printf("\n");
  qu=qu/size;
  qu2=qu2/size;
  squ=qu2-qu*qu;
  if (squ>eps)squ=sqrt(squ);else squ=1.0;

  dscr=0.0;
  weight=1./(squ*exp(iter*log((float)linker)));
  //weight=1./(squ*pow(linker,(iter)));
  for (i=0;i<size;i++){
       scri=0;
       for (j=0;j<size;j++){
       if (i!=j){
         if (dist[i][j]>eps)scri+= 1.0/dist[i][j];
       }}
       scri=scri*weight;
       if (scri>sigma*squ && sigma>eps ) scri=sigma*squ;
       if (scri>max && max>0 ) scri=max;
       if (scri<min) scri=0.0;
       //if (scri>10) printf("ITER %d LINKER %d POW %f SQU %f SCRI %f\n",iter, linker,pow(linker,(iter)),squ,scri);
       dscr= fabs(scr[i]-scri);
       error+=dscr;
       scr[i]=scr[i]+scri;
     }


  free_matrix(weight_adjacent,0,size,0,size);
  free_matrix(dist,0,size,0,size);


  return error;
}

#undef eps
