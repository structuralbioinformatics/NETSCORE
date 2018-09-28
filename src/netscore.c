#include "ppin.h"

main (int argc, char *argv[])
{

 int       i,j,ii,jj,k,n,m,iout,edges,size,*size2,intact,*iteration,diter,flag_dist,help;
 int       *n_edgescr,*n_nodescr,**shortest_path;
 int       number_phe,number_samples[MAXPHE];
 char      fint[MAXS],fcopy[MAXS],flocal[MAXS],fout[MAXS],fref[MAXS],fphe[MAXPHE][MAXS],input_U[MAXNAMELARGE],fnot[MAXS],fclus[MAXS];
 FILE      *OUT,*INP;
 float     tolerance,threshold,thero,xtolerance;
 float     rho,mind,xerror;
 float    *error,*nodepval,*nodescr,*edgescr,*average_nodescr,*average_edgescr,*rms_nodescr,*rms_edgescr;
 float    **xsimilar;
 int      **xfilter,**xgroup,*xdgroup,xngroup,**xgroup2,*xdgroup2,xngroup2,**adjacent;
 node      *protein;
 edge      *interaction;
 int       linker[MAXPARLNK];
 int       rnd[MAXPARRND];
 int      **distribution,*ndstr;
 node     *dnode;
 edge     *dedge;
 express  *xprobe;
 express  *xZscore;
 float      maxmin[MAXPARSCR];
 float      clusparam[MAXPARCLUS]; 
 float     **xPvalue,**xES,**xPscore;
 int       *xmethod,xPhenoRef,xstep,xiteration,xana;
 char       dist;
 void      Dijkstra();
 void      FloydWarshall();
 clock_t   launch, done, start_launch, end_done;
 double     diff,total_time;

 float     *NetPrioritization(),**XSimilarityMatrix(),XModifyPhenoScore();
 int       *ReadInput(),**XFilterMatrix(),ReadInteractions();
 float    **matrix(),*fvector();
 int       *ivector(),**imatrix(),**Graph();
 express   *xvector();
 edge      *evector();
 node      *avector();
 gradient  *gvector(),**gmatrix();
 void      free_fvector(),free_gvector(),free_ivector(),free_evector(),free_avector(),free_xvector(),free_gmatrix(),free_matrix(),free_imatrix();
 void      ReadPhenotypeExpress(),XPhenoScore(),PrintList(),PrintShortList(),PrintNetScore(),PrintCorrelation();

 help=0;

  for (i=0;i<argc;i++){
     if (strcmp(argv[i], "-h") == 0 )  {help=1;}          
     if (strcmp(argv[i], "-hf") == 0 ) {help=2;}          
  }
  if (argc<=1){help=1;}          
  if (help==2){PrintList();}
  if (help==1){PrintShortList();}

  size2=ivector(0,2);
  error=fvector(0,MAXERROR);
  xmethod=ivector(0,MAXMETHOD);
  iteration=ivector(0,MAXGUILD);

  launch = clock();
  start_launch = clock();
  size2=ReadInput(argc,argv,fint,fcopy,flocal,fout,fphe,&xPhenoRef,&tolerance,&threshold,iteration,&diter,&number_phe,number_samples,&thero,linker,rnd,maxmin,clusparam,xmethod,&xtolerance,&xiteration,&xana,fnot,fclus);
  done = clock();
  diff = 1000*(done - launch) / CLOCKS_PER_SEC;
  printf("Reading Input. Time %e ms\n",diff);
  switch(xmethod[5])
  { case 7: dist='e';break; 
    case 8: dist='b';break;
    case 2: dist='c';break;
    case 4: dist='a';break;
    case 1: dist='u';break;
    case 3: dist='x';break;
    case 5: dist='s';break;
    case 6: dist='k';break;
    case 9: dist='m';break;
    case 10: dist='v';break;
    case 11: dist='g';break;
    default: dist='c';
  }

  size=size2[0];
  intact=size2[1];
  printf("\nNumber of Nodes: %5d Number of Interactions: %10d\n",size,intact);
  if (size>=MAXWN)printf("Warning: MAXWN < Number of Nodes (It can produce future errors in Subroutine CumulativeMatrix)\n");

  protein=avector(0,size);
  interaction=evector(0,intact);
  nodepval=fvector(0,size);
  nodescr=fvector(0,size);
  edgescr=fvector(0,intact);
  average_nodescr=fvector(0,size);
  average_edgescr=fvector(0,intact);
  rms_nodescr=fvector(0,size);
  rms_edgescr=fvector(0,intact);
  n_nodescr=ivector(0,size);
  n_edgescr=ivector(0,intact);
  distribution=imatrix(0,MAXD,0,MAXID);
  ndstr=ivector(0,MAXD);
  xprobe=xvector(0,size);
  xZscore=xvector(0,size);
  xsimilar=matrix(0,size,0,size);
  xfilter=imatrix(0,size,0,size);
  xgroup=imatrix(0,MAXXG,0,MAXXSG);
  xdgroup=ivector(0,MAXXG);
  xgroup2=imatrix(0,MAXXG,0,MAXXSG);
  xdgroup2=ivector(0,MAXXG);
  xPvalue=matrix(0,size,0,MAXPHE);
  xES=matrix(0,size,0,MAXPHE);
  xPscore=matrix(0,size,0,MAXPHE);
  adjacent=imatrix(0,size,0,size);
  shortest_path=imatrix(0,size,0,size);

  for (i=0;i<MAXD;i++){for (j=0;j<MAXID;j++){ distribution[i][j]=0;}}
  for (i=0;i<MAXXG;i++){xdgroup[i]=0;xdgroup2[i]=0;}
  for (i=0;i<MAXXG;i++){for (j=0;j<MAXXSG;j++){xgroup[i][j]=0;xgroup2[i][j]=0;}}
  for (i=0;i<size;i++){for (j=0;j<size;j++){xsimilar[i][j]=0.0;xfilter[i][j]=0;}}
  launch = clock();
  if (xmethod[14]==0){
      edges=ReadInteractions(fint,fcopy,flocal,threshold,interaction,protein,size,distribution,ndstr,maxmin);
  }else{
      edges=DummyInteractions(fcopy,flocal,protein,size);
      xmethod[3]=0;
  }
  done = clock();
  diff = 1000*(done - launch) / CLOCKS_PER_SEC;
  printf("Reading Interactions. Time %e ms\n",diff);

  if (edges > intact) {printf("ERROR: Intact (%d) < Edges (%d) \n"); exit(1);}

  OUT=fopen(fout,"w");

  if (number_phe>0){
     if (xmethod[14]==0){
      launch = clock();
      adjacent=Graph(interaction,edges,size);
      done = clock();
      diff = 1000*(done - launch) / CLOCKS_PER_SEC;
      printf("Generate Graph Time %e ms\n",diff);
      //for (i=0;i<size;i++){ Dijkstra(adjacent,size,shortest_path[i],i); }
      launch = clock();
      for (i=0;i<size;i++){memcpy(shortest_path[i],adjacent[i],size*sizeof(int));}
      FloydWarshall(shortest_path,size);
      done = clock();
      diff = (done - launch) / CLOCKS_PER_SEC;
      printf("FloydWarshall execution Time %f sec\n",diff);
     }else{
      for (i=0;i<size;i++){for (j=0;j<size;j++){adjacent[i][j]=1;shortest_path[i][j]=1;}}
      for (i=0;i<size;i++){adjacent[i][i]=0;shortest_path[i][i]=0;}
     }
     for (i=0;i<size;i++){
       for (ii=0;ii<MAXEST;ii++){ 
        for (j=0;j<number_phe;j++){ 
           for (k=0;k<number_samples[j];k++){xprobe[i].fragment[ii].phenotype[j].sample[k]=0.0;xprobe[i].fragment[ii].phenotype[j].defined[k]=0;}}}}
     launch = clock();
     ReadPhenotypeExpress(fphe,fnot,size,number_phe,number_samples,protein,xprobe,xmethod);
     if (xana==1){XAnalyzePhenoNet(OUT,protein,interaction,size,edges,threshold,xprobe,number_phe,number_samples,xPhenoRef,shortest_path,clusparam,xmethod,fclus);}
     done = clock();
     diff = 1000*(done - launch) / CLOCKS_PER_SEC;
     printf("Done with ReadPhenotypeExpress. Time %e ms\n",diff);
     launch = clock();
     xsimilar=XSimilarityMatrix(protein,xprobe,size,number_phe,number_samples,xmethod,clusparam,shortest_path,dist);
     printf("Done with SimilarityMatrix. \n");
     xfilter=XFilterMatrix(protein,interaction,size,threshold);
     done = clock();
     diff = 1000*(done - launch) / CLOCKS_PER_SEC;
     printf("Done with SimilarityMatrix & FilterMatrix. Time %e ms\n",diff);
     if (xmethod[15]==1) PrintCorrelation(xsimilar,size); 
     launch = clock();
     XPhenoScore(protein,interaction,xPhenoRef,fphe,xsimilar,xfilter,xES,xPvalue,xZscore,xPscore,xprobe,xgroup,&xngroup,xdgroup,xmethod,number_phe,number_samples,size,edges,clusparam,fclus);
     done = clock();
     diff = 1000*(done - launch) / CLOCKS_PER_SEC;
     printf("Done with XPhenoScore. Time %e ms\n",diff);
     }
  if (xmethod[14]==0){
     launch = clock();
     error=NetPrioritization(edges,size2,interaction,protein,nodepval,nodescr,average_nodescr,rms_nodescr,n_nodescr,edgescr,average_edgescr,rms_edgescr,n_edgescr,maxmin,distribution,ndstr,linker,rnd,tolerance,threshold,thero,iteration,diter);
     done = clock();
     diff = (done - launch) / CLOCKS_PER_SEC;
     printf("Done with NetPrioritization. Time %f sec\n",diff);
    } 
  PrintNetScore(OUT,size,edges,protein,interaction,nodepval,nodescr,average_nodescr,rms_nodescr,edgescr,average_edgescr,rms_edgescr,maxmin);

  if (number_phe>0 && xiteration>0 && xmethod[14]==0){
     xerror=xtolerance + 1.0;
     xstep=1;
     while ( xerror >= xtolerance && xstep < xiteration){
         printf("\nMODIFICATION OF THE PHENOSCORE\n");
         launch = clock();
         xfilter=XFilterMatrix(protein,interaction,size,threshold);
         xerror=XModifyPhenoScore(xstep,protein,interaction,xPhenoRef,fphe,fclus,xsimilar,xfilter,xES,xPvalue,xZscore,xPscore,xprobe,xgroup,&xngroup,xdgroup,xgroup2,&xngroup2,xdgroup2,xmethod,number_phe,number_samples,size,edges,clusparam,average_nodescr,rms_nodescr,average_edgescr,rms_edgescr);
         done = clock();
         diff = 1000*(done - launch) / CLOCKS_PER_SEC;
         printf("Done with FilterMatrix & ModifyPhenoScore. Time %e ms\n",diff);
         launch = clock();
         error=NetPrioritization(edges,size2,interaction,protein,nodepval,nodescr,average_nodescr,rms_nodescr,n_nodescr,edgescr,average_edgescr,rms_edgescr,n_edgescr,maxmin,distribution,ndstr,linker,rnd,tolerance,threshold,thero,iteration,diter);
         done = clock();
         diff = (done - launch) / CLOCKS_PER_SEC;
         printf("Done with NetPrioritization. Time %f sec\n",diff);
         PrintNetScore(OUT,size,edges,protein,interaction,nodepval,nodescr,average_nodescr,rms_nodescr,edgescr,average_edgescr,rms_edgescr,maxmin);
         xstep++;
         printf("\nModify PhenoScore: Step %d up to %d (Error %e vs limit %e)\n",xstep,xiteration,xerror,xtolerance);
         }
     }
  
  end_done = clock();
  total_time = (end_done - start_launch) / CLOCKS_PER_SEC;
  printf("Total Time of execution %f sec\n",total_time);
   
  fclose(OUT);

  

  free_avector(protein,0,size);
  free_evector(interaction,0,intact);
  free_fvector(error,0,MAXERROR);
  free_fvector(nodepval,0,size);
  free_fvector(nodescr,0,size);
  free_fvector(edgescr,0,intact);
  free_fvector(average_nodescr,0,size);
  free_fvector(average_edgescr,0,intact);
  free_fvector(rms_nodescr,0,size);
  free_fvector(rms_edgescr,0,intact);
  free_ivector(n_nodescr,0,size);
  free_ivector(n_edgescr,0,intact);
  free_ivector(size2,0,2);
  free_ivector(ndstr,0,MAXD);
  free_imatrix(distribution,0,MAXD,0,MAXID);
  free_imatrix(xgroup,0,MAXXG,0,MAXXSG);
  free_ivector(xdgroup,0,MAXXG);
  free_imatrix(xgroup2,0,MAXXG,0,MAXXSG);
  free_ivector(xdgroup2,0,MAXXG);
  free_imatrix(xfilter,0,size,0,size);
  free_imatrix(adjacent,0,size,0,size);
  free_matrix(xsimilar,0,size,0,size);
  free_xvector(xprobe,0,size);
  free_xvector(xZscore,0,size);
  free_matrix(xPvalue,0,size,0,MAXPHE);
  free_matrix(xES,0,size,0,MAXPHE);
  free_matrix(xPscore,0,size,0,MAXPHE);
  free_ivector(xmethod,0,MAXMETHOD);
  free_ivector(iteration,0,MAXGUILD);
  free_imatrix(shortest_path,0,size,0,size);

  exit(1);

}
