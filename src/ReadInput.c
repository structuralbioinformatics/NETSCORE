#include "ppin.h"
int *ReadInput(argc,argv,fint,fcopy,flocal,fout,fphe,nref,tolerance,threshold,iteration,diter,number_phe,number_samples,thero,linker,rnd,maxmin,clusparam,xmethod,xtolerance,xiteration,xana,fnot,fclus)
int     argc;
char   *argv[];
char    fint[MAXS],fcopy[MAXS],flocal[MAXS],fout[MAXS],fphe[MAXPHE][MAXS],fnot[MAXS],fclus[MAXS];
float   *threshold,*tolerance,*xtolerance,*thero,maxmin[],clusparam[];
int    *iteration,*xiteration,*xana,*diter,*xmethod,*nref,*number_phe,number_samples[],rnd[],linker[];
{

  FILE     *INP,*EXP;
  char      a[MAXNAME],b[MAXNAME],protein[MAXP][MAXNAME],buffer[MAXS];
  char      array[MAXS],fref[MAXS];
  float     affinity,max,min,maxe,mine,mind,minde,minst,minste,sigma,xcl,xca,xcr,xcm,xcma,xci,xcw,xmx,wpath,c,r,p,min_tpvalue,method_foldrkp, method_pmaxrkp;
  int       iteration_global, iteration_netscore, iteration_netzscore, iteration_netshort;
  int       nph,size,i,j,k,skipa,skipb,*n,linker_edge,linker_node,linker_edge_ext,linker_node_ext,linker_direct,linker_degree,nrand,nrande, method_toprkp, method_pvalue, method_fpvalue1, method_fpvalue2, method_rpvalue,method_netcombo, method_guild, rtype,nzscr,method_statistic,method_strength,method_seed,method_genecluster, method_sign,method_edge,method_sim,method_norm,method_corr,method_path,method_cluster,method_hcluster,method_kcluster,method_ncluster,method_sxcluster,method_sycluster,method_center;
  int      *ivector();
  express  *probe;
  express  *xvector();
  void      nrerror();
  void      free_xvector();

/* Reading input line */

  memset(fint,'\0',MAXS); 
  memset(fout,'\0',MAXS); 
  memset(fcopy,'\0',MAXS); 
  memset(flocal,'\0',MAXS); 
  memset(fnot,'\0',MAXS); 
  memset(fclus,'\0',MAXS); 
  memset(fref,'\0',MAXS); 
  memset(buffer,'\0',MAXS); 

  if (MAXGUILD < 4) nrerror("Check MAXGUILD < 4 to define ITERATIONS on methods");
  iteration[0]=iteration_global=1;
  iteration[1]=iteration_netscore=0;
  iteration[2]=iteration_netzscore=0;
  iteration[3]=iteration_netshort=0;
  
  *xiteration=0;
  *xana=0;
  *diter=0;
  *thero=0.1;
  *tolerance=1.0e-10;
  *xtolerance=1.0e-8;
  *threshold=0.0;
  *nref=0;
   nph=0;
   linker[0]=MINLINK;
   linker[1]=MINLINK;
   linker[2]=0;
   linker[3]=0;
   linker[4]=0;   /* indirect links*/
   linker[5]=0;
   rnd[0]=100;
   rnd[1]=1;
   rnd[2]=0;
   rnd[3]=10;
   rnd[4]=0;
   rnd[5]=0;
   rnd[6]=2;
   rnd[7]=1;
   rnd[8]=3;
   maxmin[0]=MAXSCR;
   maxmin[1]=MINSCR;
   maxmin[2]=MAXSIGMA;
   maxmin[3]=0;
   maxmin[4]=MINSCR;
   maxmin[5]=MINSCR;
   maxmin[6]=MINSCR;
   maxmin[7]=MINSTR;
   maxmin[8]=MINSTR;
   maxmin[9]=MINSCR;
   clusparam[0]=MINXLINK;
   clusparam[1]=MAXXDEV;
   clusparam[2]=MINXCOMM;
   clusparam[3]=MINXDEV; 
   clusparam[4]=0;
   clusparam[5]=1;
   clusparam[6]=2.0;
   clusparam[7]=0.0;
   clusparam[8]=10.0;
   clusparam[9]=1.0;
   clusparam[10]=2.0;
   clusparam[11]=1.0;
   clusparam[12]=0.0; 
   clusparam[13]=50.0; 
  

   for(i=0;i<MAXMETHOD;i++){xmethod[i]=0;}
   xmethod[11]=1;
   xmethod[13]=2;
   xmethod[16]=1;
   xmethod[17]=2;
   xmethod[19]=1500;
   xmethod[20]=100;
   xmethod[21]=50000;

  strcpy(fout,"NetScore.out");
  printf("Executing NETSCORE EXPRESS: \n");
  for (i=0;i<argc;i++){
      if (strcmp(argv[i], "-i") == 0) {	/*File with Interactions  */
         strcpy(fint, argv[i+1]);
         printf("\tINPUT -i: %s\n",fint);
      }
      if (strcmp(argv[i], "-o") == 0) {	/*File with output  */
         strcpy(fout, argv[i+1]);
         printf("\tINPUT -o: %s\n",fout);
      }
      if (strcmp(argv[i], "-c") == 0) {	/*File with Copy number  */
         strcpy(fcopy, argv[i+1]);
         printf("\tINPUT -c: %s\n",fcopy);
      }
      if (strcmp(argv[i], "-l") == 0) {	/*File with Localization  */
         strcpy(flocal, argv[i+1]);
         printf("\tINPUT -l: %s\n",flocal);
      }
      if (strcmp(argv[i], "-dn") == 0) {	/*Linker degree to predict nodes  */
         sscanf(argv[i+1], "%d", &linker_node);
         linker[0]=linker_node;
         printf("\tINPUT -dn: %d\n",linker_node);
      }
      if (strcmp(argv[i], "-de") == 0) {	/*Linker degree to predict edges  */
         sscanf(argv[i+1], "%d", &linker_edge);
         linker[1]=linker_edge;
         printf("\tINPUT -de: %d\n",linker_edge);
      }
      if (strcmp(argv[i], "-dxn") == 0) {	/*Linker degree to increase nodes score  */
         sscanf(argv[i+1], "%d", &linker_node_ext);
         linker[2]=linker_node_ext;
         printf("\tINPUT -dxn: %d\n",linker_node_ext);
      }
      if (strcmp(argv[i], "-dxe") == 0) {	/*Linker degree to increase edges score  */
         sscanf(argv[i+1], "%d", &linker_edge_ext);
         linker[3]=linker_edge_ext;
         printf("\tINPUT -dxe: %d\n",linker_edge_ext);
      }
      if (strcmp(argv[i], "-dxi") == 0) {	/*Indirect links can substitute direct neighbors  */
         sscanf(argv[i+1], "%d", &linker_direct);
         linker[4]=linker_direct;
         printf("\tINPUT -dxi: %d\n",linker_direct);
      }
      if (strcmp(argv[i], "-dnn") == 0) {	/*Weight score by degree */
         sscanf(argv[i+1], "%d", &linker_degree);
         linker[5]=linker_degree;
         printf("\tINPUT -dnn: %d\n",linker_degree);
      }
      if (strcmp(argv[i], "-nr") == 0) {	/*Number of random iterations for Z-score */
         sscanf(argv[i+1], "%d", &nrand);
         rnd[0]=nrand;
         printf("\tINPUT -nr: %d\n",nrand);
      }
      if (strcmp(argv[i], "-r") == 0) {	/*Number of random Type  for Z-score*/
         sscanf(argv[i+1], "%d", &rtype);
         if (rtype>1||rtype<0){rtype=1;}
         rnd[1]=rtype;
         printf("\tINPUT -r: %d\n",rtype);
      }
      if (strcmp(argv[i], "-zp") == 0) {	/*Population to calculate Z-score */
         sscanf(argv[i+1], "%d", &nzscr);
         rnd[2]=nzscr;
         printf("\tINPUT -zp: %d\n",nzscr);
      }

      if (strcmp(argv[i], "-npv") == 0) {	/*  Number of functions to model random Scores */
         sscanf(argv[i+1], "%d", &method_rpvalue);
         rnd[3]=method_rpvalue;
         printf("\tINPUT -npv: %d\n",method_rpvalue);
      }


      if (strcmp(argv[i], "-mth") == 0) {	/* Selection of method of Prioritization*/
         sscanf(argv[i+1], "%d", &method_guild);
         rnd[4]=method_guild;
         printf("\tINPUT -mth: %d\n",method_guild);
      }
      if (strcmp(argv[i], "-ncz") == 0) {	/* Selection of method of NetCombo Z-scoring of NetScore NetZscore and NetShort */
         sscanf(argv[i+1], "%d", &method_netcombo);
         rnd[5]=method_netcombo;
         printf("\tINPUT -ncz: %d\n",method_netcombo);
      }
      if (strcmp(argv[i], "-mpv") == 0) {	/*  Type of functions to model random Scores */
         sscanf(argv[i+1], "%d", &method_pvalue);
         if (method_pvalue<5) rnd[6]=method_pvalue;
         printf("\tINPUT -mpv: %d\n",rnd[6]);
      }

      if (strcmp(argv[i], "-nfpv1") == 0) {	/*  Number of functions to model random Scores */
         sscanf(argv[i+1], "%d", &method_fpvalue1);
         rnd[7]=method_fpvalue1;
         printf("\tINPUT -nfpv1: %d\n",method_fpvalue1);
      }
      if (strcmp(argv[i], "-nfpv2") == 0) {	/*  Number of functions to model random Scores */
         sscanf(argv[i+1], "%d", &method_fpvalue2);
         rnd[8]=method_fpvalue2;
         printf("\tINPUT -nfpv2: %d\n",method_fpvalue2);
      }
     // if (strcmp(argv[i], "-nts") == 0) {	/*  Minimum score to consider Seeds */
     //    sscanf(argv[i+1], "%f", &min_tpvalue);
     //    maxmin[9]=min_tpvalue;
     //    printf("\tINPUT -npv: %d\n",min_tpvalue);
     // }
      if (strcmp(argv[i], "-n") == 0) {	/*Limit number of iterations */
         sscanf(argv[i+1], "%d", &iteration_global);
         iteration[0]=iteration_global;
         printf("\tINPUT -n: %d\n",iteration_global);
      }
      if (strcmp(argv[i], "-ns") == 0) {	/*Limit number of iterations */
         sscanf(argv[i+1], "%d", &iteration_netscore);
         iteration[1]=iteration_netscore;
         printf("\tINPUT -ns: %d\n",iteration_netscore);
      }
      if (strcmp(argv[i], "-nz") == 0) {	/*Limit number of iterations */
         sscanf(argv[i+1], "%d", &iteration_netzscore);
         iteration[2]=iteration_netzscore;
         printf("\tINPUT -nz: %d\n",iteration_netscore);
      }
      if (strcmp(argv[i], "-nh") == 0) {	/*Limit number of iterations */
         sscanf(argv[i+1], "%d", &iteration_netshort);
         iteration[3]=iteration_netshort;
         printf("\tINPUT -nh: %d\n",iteration_netshort);
      }
      if (strcmp(argv[i], "-nd") == 0) {	/*Updating Dual Espace on number of iterations */
         sscanf(argv[i+1], "%d", diter);
         printf("\tINPUT -nd: %d\n",*diter);
      }
      if (strcmp(argv[i], "-mx") == 0) {	/*MAX Score transfer */
         sscanf(argv[i+1], "%f", &max);
         maxmin[0]=max;
         printf("\tINPUT -mx: %e\n",max);
      }
      if (strcmp(argv[i], "-mxe") == 0) {	/*MAX Edge Score transfer */
         sscanf(argv[i+1], "%f", &maxe);
         maxmin[3]=maxe;
         printf("\tINPUT -mxe: %e\n",maxe);
      }
      if (strcmp(argv[i], "-mn") == 0) {	/*MIN Score transfer */
         sscanf(argv[i+1], "%f", &min);
         maxmin[1]=min;
         printf("\tINPUT -mn: %e\n",min);
      }
      if (strcmp(argv[i], "-mne") == 0) {	/*MIN Edge Score transfer */
         sscanf(argv[i+1], "%f", &mine);
         maxmin[4]=mine;
      }
      if (strcmp(argv[i], "-mnd") == 0) {	/*MIN Edge Score transfer */
         sscanf(argv[i+1], "%f", &mind);
         maxmin[5]=mind;
         printf("\tINPUT -mnd: %e\n",mind);
      }
      if (strcmp(argv[i], "-mnde") == 0) {	/*MIN Edge Score transfer */
         sscanf(argv[i+1], "%f", &minde);
         maxmin[6]=minde;
         printf("\tINPUT -mnde: %e\n",minde);
      }
      if (strcmp(argv[i], "-mnst") == 0) {	/*MIN Edge Score transfer */
         sscanf(argv[i+1], "%f", &minst);
         maxmin[7]=minst;
         printf("\tINPUT -mnst: %e\n",minst);
      }
      if (strcmp(argv[i], "-mnste") == 0) {	/*MIN Edge Score transfer */
         sscanf(argv[i+1], "%f", &minste);
         maxmin[8]=minste;
         printf("\tINPUT -mnste: %e\n",minste);
      }
      if (strcmp(argv[i], "-ms") == 0) {	/*Max #times of sigma in Score transfer */
         sscanf(argv[i+1], "%f", &sigma);
         maxmin[2]=sigma;
         printf("\tINPUT -ms: %e\n",sigma);
      }
      if (strcmp(argv[i], "-e") == 0) {	/*Tolerance  to stop iterations */
         sscanf(argv[i+1], "%f", tolerance);
         printf("\tINPUT -e: %e\n",tolerance);
      }
      if (strcmp(argv[i], "-z") == 0) {	/*Score  threshold */
         sscanf(argv[i+1], "%f", thero);
         printf("\tINPUT -z: %e\n",*thero);
      }
      if (strcmp(argv[i], "-t") == 0) {	/*Score  threshold*/
         sscanf(argv[i+1], "%f", threshold);
         printf("\tINPUT -t: %e\n",*threshold);
      }
      if (strcmp(argv[i], "-xphe") == 0) {	/*File with samples of one phenotype  */
         strcpy(fphe[nph], argv[i+1]);
         nph++;
         if (nph>MAXPHE)nrerror("Too many phenotypes, increase MAXPH");
      }
      if (strcmp(argv[i], "-xphs") == 0) {	/*File with samples of one phenotype  */
         strcpy(fref, argv[i+1]);
      }
      if (strcmp(argv[i], "-xst") == 0) {	/* Statistical test to assign genes to a phenotype */
         sscanf(argv[i+1], "%d", &method_statistic);
         xmethod[1]=method_statistic;
         printf("\tINPUT -xst: %d\n",method_statistic);
      }
      if (strcmp(argv[i], "-xsign") == 0) {	/* Use of negative correlation */
         sscanf(argv[i+1], "%d", &method_sign);
         xmethod[2]=method_sign;
         printf("\tINPUT -xsign: %d\n",method_sign);
      }
      if (strcmp(argv[i], "-xedge") == 0) {	/* Edge scoring based on expression */
         sscanf(argv[i+1], "%d", &method_edge);
         xmethod[3]=method_edge;
         printf("\tINPUT -xedge: %d\n",method_edge);
      }
      if (strcmp(argv[i], "-xsim") == 0) {	/* Edge scoring based on expression */
         sscanf(argv[i+1], "%d", &method_sim);
         xmethod[4]=method_sim;
         printf("\tINPUT -xsim: %d\n",method_sim);
      }
      if (strcmp(argv[i], "-xiter") == 0) {	/*Limit number of iterations in PhenoScore */
         sscanf(argv[i+1], "%d", xiteration);
         printf("\tINPUT -xiter: %d\n",*xiteration);
      }      
      if (strcmp(argv[i], "-xana") == 0) {	/* Flag to performe comparative analyses */
         *xana=1;
         xmethod[18]=1;
         printf("\tINPUT -xana active\n");
      }
      if (strcmp(argv[i], "-xtol") == 0) {	/*Tolerance  to stop iterations in PhenoScore*/
         sscanf(argv[i+1], "%f", xtolerance);
         printf("\tINPUT -xtol: %e\n",*xtolerance);
      }
      if (strcmp(argv[i], "-xstg") == 0) {	/* Statistical scores assigned to cluster-genes or individual genes*/
         sscanf(argv[i+1], "%d", &method_genecluster);
         xmethod[0]=method_genecluster;
         printf("\tINPUT -xstg: %d\n",method_genecluster);
      }
      if (strcmp(argv[i], "-xrkpf") == 0) {	/*Fold threshold to select seeds with RankProd */
         sscanf(argv[i+1], "%f", &method_foldrkp);
         xmethod[19]=(int)(method_foldrkp*1000);
         printf("\tINPUT -xrkpf: %f\n",method_foldrkp);
      }
      if (strcmp(argv[i], "-xrkpt") == 0) {	/*Fold threshold to select seeds with RankProd */
         sscanf(argv[i+1], "%d", &method_toprkp);
         xmethod[20]=method_toprkp;
         printf("\tINPUT -xrkpt: %d\n",method_toprkp);
      }
      if (strcmp(argv[i], "-xrkpp") == 0) {	/*P-value threshold to select seeds with RankProd */
         sscanf(argv[i+1], "%f", &method_pmaxrkp);
         xmethod[21]=(int)(method_pmaxrkp*1000000);
         printf("\tINPUT -xrkpp: %f\n",method_pmaxrkp);
      }
      if (strcmp(argv[i], "-xnot") == 0) {	/* File with Array annotations for protein-nodes  */
         strcpy(fnot, argv[i+1]);
         printf("\tINPUT -xnot: %s\n",fnot);
      }
      if (strcmp(argv[i], "-xcorr") == 0) {	/* Specifies the distance measure  */
         sscanf(argv[i+1], "%d", &method_corr);
         xmethod[5]=method_corr;
         printf("\tINPUT -xcorr: %d\n",method_corr);
      }
      if (strcmp(argv[i], "-xpath") == 0) {	/* Weighting the distance measure between genes with the shortest-path distance  */
         sscanf(argv[i+1], "%d", &method_path);
         xmethod[6]=method_path;
         printf("\tINPUT -xpath: %d\n",method_path);
      }
      if (strcmp(argv[i], "-xprn") == 0) {	/* Flag to print matrix of distances  and  inputs */
         xmethod[7]=1;
         printf("\tINPUT -xprn active\n");
      }
      if (strcmp(argv[i], "-xlog") == 0) {	/* Flag to transform data to log */
         xmethod[8]=1;
         printf("\tINPUT -xlog active\n");
      }
      if (strcmp(argv[i], "-xng") == 0) {	/* Flag normalize  genes */
         xmethod[9]=1;
         printf("\tINPUT -xng active\n");
      }
      if (strcmp(argv[i], "-xcg") == 0) {	/* Center genes */
         sscanf(argv[i+1], "%d", &method_center);
         xmethod[10]=method_center;
         printf("\tINPUT -xcg: %d\n",method_center);
      }
      if (strcmp(argv[i], "-xcp") == 0) {	/* Method of clustering */
         sscanf(argv[i+1], "%d", &method_cluster);
         xmethod[11]=method_cluster;
         printf("\tINPUT -xcp: %d\n",method_cluster);
      }
      if (strcmp(argv[i], "-xncr") == 0) {	/* Flag normalize  genes */
         sscanf(argv[i+1], "%d", &method_norm);
         xmethod[12]=method_norm;
         printf("\tINPUT -xncr: %d\n",method_norm);
      }
      if (strcmp(argv[i], "-xsts") == 0) {	/* Statistical test to assign genes to a phenotype */
         sscanf(argv[i+1], "%d", &method_strength);
         xmethod[13]=method_strength;
         printf("\tINPUT -xsts: %d\n",method_strength);
      }
      if (strcmp(argv[i], "-x") == 0) {	/* Flag to neglect interactions */
         xmethod[14]=1;
         printf("\tINPUT -x active (interactions are not read)\n");
      }
      if (strcmp(argv[i], "-stop") == 0) {	/* Flag to neglect interactions */
         xmethod[15]=1;
         printf("\tINPUT -stop active \n");
      }
      if (strcmp(argv[i], "-xcsg") == 0) {	/* Flag to neglect interactions */
         xmethod[16]=0;
         printf("\tINPUT -xcsg active \n");
      }
      if (strcmp(argv[i], "-xsto") == 0) {	/* Statistical test to assign genes to a phenotype */
         sscanf(argv[i+1], "%d", &method_seed);
         xmethod[17]=method_seed;
         printf("\tINPUT -xsto: %d\n",method_seed);
      }
      if (strcmp(argv[i], "-xwp") == 0) {	/* Parameter to weight correlartion by shortest path*/
         sscanf(argv[i+1], "%f", &wpath);
         clusparam[5]=wpath;
         printf("\tINPUT -xwp: %e\n",wpath);
      }
      if (strcmp(argv[i], "-xcl") == 0) {	/* Density Search Clustering: minimum score to link a putative new member to a group */
         sscanf(argv[i+1], "%f", &xcl);
         clusparam[0]=xcl;
         printf("\tINPUT -xcl: %e\n",xcl);
         }
      if (strcmp(argv[i], "-xca") == 0) {	/* Density Search Clustering: Maximum deviation of the mean produced by the inclusion of a new member */
         sscanf(argv[i+1], "%f", &xca);
         clusparam[1]=xca;
         printf("\tINPUT -xca: %e\n",xca);
         }
      if (strcmp(argv[i], "-xcr") == 0) {	/* Density Search Clustering: Minimum ratio of common elements for merging clusters(0 for null intersection)  */
         sscanf(argv[i+1], "%f", &xcr);
         clusparam[2]=xcr;
         printf("\tINPUT -xcr: %e\n",xcr);
         }
      if (strcmp(argv[i], "-xcm") == 0) {	/* Density Search Clustering: Minimum mean deviation when merging clusters */
         sscanf(argv[i+1], "%f", &xcm);
         clusparam[3]=xcm;
         printf("\tINPUT -xcm: %e\n",xcm);
         }
      if (strcmp(argv[i], "-xcma") == 0) {	/* Density Search Clustering: Minimum mean cut-off expansion */
         sscanf(argv[i+1], "%f", &xcma);
         clusparam[4]=xcma;
         printf("\tINPUT -xcma: %e\n",xcma);
         }
      if (strcmp(argv[i], "-xci") == 0) {	/* MCL  Clustering Inflation*/
         sscanf(argv[i+1], "%f", &xci);
         clusparam[6]=xci;
         printf("\tINPUT -xci: %e\n",xci);
         }
      if (strcmp(argv[i], "-xcw") == 0) {	/* MCL  Clustering Threshold of distance similarity to generate graph*/
         sscanf(argv[i+1], "%f", &xcw);
         clusparam[12]=xcw;
         printf("\tINPUT -xcw: %e\n",xcw);
         }
      if (strcmp(argv[i], "-xhm") == 0) {	/* Method of Hierarchical clustering */
         sscanf(argv[i+1], "%d", &method_hcluster);
         clusparam[7]=(float) method_hcluster;
         printf("\tINPUT -xhm: %d\n",method_hcluster);
         }
      if (strcmp(argv[i], "-xcn") == 0) {	/* Number of expected clusters (hierarchical|k-means) */
         sscanf(argv[i+1], "%d", &method_ncluster);
         clusparam[8]=(float) method_ncluster;
         printf("\tINPUT -xcn: %d\n",method_ncluster);
         }
      if (strcmp(argv[i], "-xkr") == 0) {	/* Number of times the k-means clustering algorithm is run  */
         sscanf(argv[i+1], "%d", &method_kcluster);
         clusparam[9]=(float) method_kcluster;
         printf("\tINPUT -xkr: %d\n",method_kcluster);
         }
      if (strcmp(argv[i], "-xxn") == 0) {	/* Specifies the horizontal dimension of the SOM grid (default: 2) */
         sscanf(argv[i+1], "%d", &method_sxcluster);
         clusparam[10]=(float) method_sxcluster;
         printf("\tINPUT -xxn: %d\n",method_sxcluster);
         }
      if (strcmp(argv[i], "-xyn") == 0) {	/* Specifies the horizontal dimension of the SOM grid (default: 2) */
         sscanf(argv[i+1], "%d", &method_sycluster);
         clusparam[11]=(float) method_sycluster;
         printf("\tINPUT -xyn: %d\n",method_sycluster);
         }
      if (strcmp(argv[i], "-xmx") == 0) {	/* Specifies the horizontal dimension of the SOM grid (default: 2) */
         sscanf(argv[i+1], "%f", &xmx);
         clusparam[13]=xmx;
         printf("\tINPUT -xmx: %f\n",xmx);
         }
      if (strcmp(argv[i], "-xcf") == 0) {	/* File with clusters  */
         strcpy(fclus, argv[i+1]);
         printf("\tINPUT -xcf: %s\n",fclus);
      }
         



  } 

  

  *number_phe=nph;

   if (clusparam[12]<=0){ if (xmethod[12]==1) {clusparam[12]=5.0;}else{clusparam[12]=0.2;}}
   if (xmethod[6]==2 && clusparam[5]==1){clusparam[5]=100;}

  printf("Options METHOD:\t\t\t");
  for (k=0;k<MAXMETHOD;k++){printf("[%d] %d ;",k,xmethod[k]);}
   printf("\n");
  printf("Options CLUSTER:\t\t");
  for (k=0;k<MAXPARCLUS;k++){printf("[%d] %f ;",k,clusparam[k]);}
   printf("\n");

  printf("Check Files\n");
  printf("File with interactions:    \t%s\n",fint); 
  printf("File with protein score:   \t%s\n",fcopy); 
  printf("File with protein location:\t%s\n",flocal);
  printf("File with annotation:      \t%s\n",fnot);
  printf("Output file:               \t%s\n",fout); 

  for (k=0;k<nph;k++){
   INP=fopen(fphe[k],"r");
   if (!INP){nrerror("No phenotype data");}else{printf("File with Phenotype:       \t%s\n",fphe[k]);}
   if(strcmp(fref,fphe[k])==0){*nref=k;printf("Reference Phenotype:       \t%s\n",fref);}
   number_samples[k]=0;
   while(!feof(INP)) {
    memset(buffer,'\0',MAXS);
    fgets(buffer,MAXS,INP);
    sscanf(buffer,"%s",array);
    EXP=fopen(array,"r");
    if (!EXP){nrerror("No sample data");}else{number_samples[k]++;}
    fclose(EXP);
    if (number_samples[k]>MAXSAMPLE)nrerror("Too many samples per phenotype, increase MAXSAMPLE");
   }
   number_samples[k]=number_samples[k]-1;
   fclose(INP);
  }


  if (xmethod[14]==0){
  INP=fopen(fint,"r");
  if (!INP){nrerror("No Interaction data");}
  size=0;
  i=0;
  while(!feof(INP)) {
    memset(buffer,'\0',MAXS);
    fgets(buffer,MAXS,INP);
    sscanf(buffer,"%s%s%f",a,b,&affinity);
    skipa=0;
    skipb=0;
    for (j=0;j<size;j++){
      if (!strcmp(protein[j],a)){skipa++;}
      if (!strcmp(protein[j],b)){skipb++;}
      if (skipa>0 && skipb>0){break;}
    }
    if (size>MAXP){nrerror("Error: increase MAXP\n");}
    if (skipa == 0 && affinity > *threshold ){strcpy(protein[size],a);size++;}
    if (skipb == 0 && affinity > *threshold ){strcpy(protein[size],b);size++;}
    if (affinity > *threshold ) i++;
  }
  fclose(INP);
  }else{
  INP=fopen(fcopy,"r");
  if (!INP){nrerror("No Protein data");}
  size=0;
  i=0;
  while(!feof(INP)) {
    memset(buffer,'\0',MAXS);
    fgets(buffer,MAXS,INP);
    sscanf(buffer,"%s%f%f%f",a,&c,&r,&p);
    strcpy(protein[size],a);
    size++;
    if (size>MAXP){nrerror("Error: increase MAXP\n");}
  }
  fclose(INP);
  }
  n=ivector(0,2);
  n[0]=size;
  n[1]=i;
  return n;
}

