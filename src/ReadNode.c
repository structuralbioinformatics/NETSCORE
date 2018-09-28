#include "ppin.h"
int ReadNode(fcopy,flocal,copy,proteome)
char    fcopy[MAXS],flocal[MAXS],proteome[MAXP][MAXNAME];
abundance copy[MAXP];
{

  FILE     *INP;
  char      a[MAXNAME],loc[MAXNAME],buffer[MAXS],plocal[MAXP][MAXNAME];
  int       i,ii,j,k,n,size;
  float     c,r,p;
  int       local[LOCAL];   
  abundance clocal[MAXP];
  void      nrerror();
  
  size=0;
  INP=fopen(fcopy,"r");
  if (!INP) {nrerror("COPY-NUMBER file not found\n");}
  while(!feof(INP)) {
    memset(buffer,'\0',MAXS);
    fgets(buffer,MAXS,INP);
    sscanf(buffer,"%s%f%f%f",a,&c,&r,&p);
    strcpy(proteome[size],a);
    copy[size].cell=c;
    copy[size].rms=r;
    copy[size].score=p;
    size++;
    if (size>MAXP){nrerror("Too many proteins in COPY-NUMBER: Increase MAXP\n");}
  }
  fclose(INP);

  k=0;
  INP=fopen(flocal,"r");
  if (INP) {
   while(!feof(INP)) {
    memset(buffer,'\0',MAXS);
    fgets(buffer,MAXS,INP);
    sscanf(buffer,"%s",a);
    strcpy(plocal[k],a);
    for (j=0;j<LOCAL;j++){clocal[k].local[j]=0.0;}
    n=0;
    for (j=0;j<LOCAL;j++){ 
        memset(loc,'\0',MAXNAME);
        sscanf(buffer+MAXNAME+j*10,"%s",loc);
        if (strlen(loc)>0){
         sscanf(loc,"%d",&local[n]);
         n++;
        } 
    }
    for (j=0;j<n;j++){clocal[k].local[local[j]]=1.0/(float)n;}
    k++;
    if (k>MAXP){nrerror("Too many proteins in LOCALIZATION: Increase MAXP\n");}
   }
   fclose(INP);
  }else{
   for (i=0;i<size;i++){
     strcpy(plocal[i],proteome[i]);
     for (j=0;j<LOCAL;j++){clocal[i].local[j]=0.0;}
     local[0]=1;
     for (j=0;j<1;j++){clocal[i].local[local[j]]=1.0;}
   }
   k=size;
  }

  for (i=0;i<size;i++){
    for (ii=0;ii<k;ii++){
        if (!strcmp(proteome[i],plocal[ii])){
           for (j=0;j<LOCAL;j++){copy[i].local[j]=clocal[ii].local[j]*copy[i].cell;}
        }
    }
  }

  return size;
}
