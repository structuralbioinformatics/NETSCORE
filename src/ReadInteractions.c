#include "ppin.h"
int     ReadInteractions(fint,fcopy,flocal,threshold,interaction,protein,size,distribution,ndstr,maxmin)
char    fint[MAXS],fcopy[MAXS],flocal[MAXS];
float   threshold;
edge   *interaction;
node   *protein;
int     size,*ndstr,**distribution;
float   maxmin[];
{

  FILE       *INP;
  char        a[MAXNAME],b[MAXNAME],buffer[MAXS],proteome[MAXP][MAXNAME];
  int         i,ii,j,jj,kk,k,proteosize,n,nn,m,skipa,skipb,skipc, dim, max_degree,min_dstr,top10;
  float       affinity,max;
  int         *ini;
  abundance   copy[MAXP];
  int         ReadNode();
  abundance   AssignNode();
  void        nrerror();
  void      free_ivector();
  int       *ivector();
  void       gdistribute();

  printf("\nReading interactions\n");

  proteosize=ReadNode(fcopy,flocal,copy,proteome);


  INP=fopen(fint,"r");
  ii=0;
  n=0;
  if (!INP){nrerror("INTERACTIONS file not found\n");}
  while(!feof(INP)) {
   memset(buffer,'\0',MAXS);
   memset(a,'\0',MAXNAME);
   memset(b,'\0',MAXNAME);
   affinity=0;
   fgets(buffer,MAXS,INP);
   sscanf(buffer,"%s%s%f",a,b,&affinity);
   if (affinity > threshold){
    skipa=0;
    skipb=0;
    skipc=0;
    for (j=0;j<n;j++){
      if (!strcmp(protein[j].name1,a) || !strcmp(protein[j].name2,a)){skipa++;interaction[ii].i=j;}
      if (!strcmp(protein[j].name1,b) || !strcmp(protein[j].name2,b)){skipb++;interaction[ii].j=j;}
      if (skipa>0 && skipb>0){break;}
    }
    if (skipa>0 && skipb>0){
      for (k=0;k<ii;k++){
        if (  (interaction[k].i==interaction[ii].i && interaction[k].j==interaction[ii].j)
            ||(interaction[k].i==interaction[ii].j && interaction[k].j==interaction[ii].i)){
          skipc++;break;
        }
      }
    }
    if (skipc == 0){
     if (skipa == 0){
       protein[n].copy=AssignNode(copy,proteome,proteosize,a);
       strcpy(protein[n].name1,a);
       strcpy(protein[n].name2,a);
       interaction[ii].i=n;
       interaction[ii].a=protein[n];
       n++;
     }
     if (skipb == 0 && strcmp(b,a)){
       protein[n].copy=AssignNode(copy,proteome,proteosize,b);
       strcpy(protein[n].name1,b);
       strcpy(protein[n].name2,b);
       interaction[ii].j=n;
       interaction[ii].b=protein[n];
       n++;
     }
     if (n>size*size){nrerror("Error: number of nodes larger than expected");}
     interaction[ii].association=affinity;
     ii++;
    }
   }
  }
  fclose(INP);
  if (MAXI<MAXD){nrerror("Maximum Degree (MAXD) larger than MAXI");}
  for (k=0;k<MAXD;k++){ndstr[k]=0;}
  max_degree=0;
  printf("\nDegree Distribution\n");
  for (k=0;k<size;k++){
    protein[k].degree=0;
    for (j=0;j<ii;j++){
        if (interaction[j].i == k || interaction[j].j == k){
           protein[k].interact[protein[k].degree]=j;
           protein[k].degree++;
           if (protein[k].degree>MAXI){nrerror("Degree larger than MAXI");}
        }
    }
    if (protein[k].degree<MAXD-1){
       if ((protein[k].degree+1) > max_degree){max_degree=protein[k].degree+1;}
       if (ndstr[protein[k].degree]<MAXID){dim=ndstr[protein[k].degree];}else{dim=MAXID-1;printf("NDSTR[%d]=%d > MAXID\n",protein[k].degree,ndstr[protein[k].degree]);}
       distribution[protein[k].degree][dim]=k;
       //printf("Primary Name %s Secondary %s Distribution[%d][%d] = %d\n",protein[k].name1,protein[k].name2,protein[k].degree,ndstr[protein[k].degree],k);
       ndstr[protein[k].degree]+=1;
    }else{
       max_degree=MAXD;
       if (ndstr[MAXD-1]<MAXID){dim=ndstr[MAXD-1];}else{dim=MAXID-1;printf("NDSTR[%d]=%d > MAXID\n",MAXD-1,ndstr[MAXD-1]);}
       distribution[MAXD-1][dim]=k;
       //printf("Primary Name %s Secondary %s Distribution[%d][%d] = %d\n",protein[k].name1,protein[k].name2,MAXD-1,ndstr[MAXD-1],k);
       ndstr[MAXD-1]+=1;
    }
  }

  
  for (k=0;k<size;k++){
   if (protein[k].degree<MAXD-1){
     if (ndstr[protein[k].degree]>MAXID){
     char* errorf=malloc(MAXS*sizeof(char));
     sprintf(errorf,"Error distribution: Limiting  MAXID was surpassed in degree %d (NDSTR[%d]=%d)\n",k,protein[k].degree,ndstr[protein[k].degree]);
     nrerror(errorf);
     }
   }else{
     if (ndstr[MAXD-1]>MAXID){
     char* errorf=malloc(MAXS*sizeof(char));
     sprintf(errorf,"Error distribution: Limiting  MAXID was surpassed in degree %d (NDSTR[%d]=%d)\n",k,MAXD-1,ndstr[MAXD-1]);
     nrerror(errorf);
     }
   }
  }

  printf("Distribution of Nodes by Degrees (printing first 10 nodes)\n");
  for   (k=0;k<max_degree;k++) {
  //for   (j=0;j<ndstr[k];j++){
  top10=ndstr[k];
  if (ndstr[k]>10) top10=10;
  for   (j=0;j<top10;j++){
    printf("Distribution[%d][%d] = %d (Primary Name %s Secondary %s) \n",k,j,distribution[k][j],protein[distribution[k][j]].name1,protein[distribution[k][j]].name2);
  }}
  

  min_dstr=(int)maxmin[7];
  gdistribute(distribution,ndstr,max_degree,MAXID,min_dstr);
  fflush(stdout);
  printf("Re-Distribution of Nodes by Degrees (printing first 10 nodes)\n");
  for   (k=0;k<max_degree;k++) {
  //for   (j=0;j<ndstr[k];j++){
  top10=ndstr[k];
  if (ndstr[k]>10) top10=10;
  for   (j=0;j<top10;j++){
    printf("Re-Distribution[%d][%d] = %d (Primary Name %s Secondary %s) \n",k,j,distribution[k][j],protein[distribution[k][j]].name1,protein[distribution[k][j]].name2);
  }}
  
  fflush(stdout);


  max=0.0;
  for (j=0;j<ii;j++){
    interaction[j].a=protein[interaction[j].i];
    interaction[j].b=protein[interaction[j].j];
    if (max<=interaction[j].association){max=interaction[j].association;}
  }

  printf("MAXIMUM affinity %e\n",max);
  for (j=0;j<ii;j++){ interaction[j].association=interaction[j].association/max;}
  

  return ii;
 
}
  

