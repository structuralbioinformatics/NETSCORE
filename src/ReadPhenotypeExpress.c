#include "ppin.h"

void ReadPhenotypeExpress(fphe,fnot,size,number_phe,number_samples,protein,xprobe,xmethod)
char      fphe[MAXPHE][MAXS],fnot[MAXS];
int       size,number_phe,*number_samples,*xmethod;
express   *xprobe;
node      *protein;
{

 char      array[MAXS],buffer[MAXS],a[MAXNAME],b[MAXNAME];
 int       i,j,jj,n,s,count;
 float     exp;
 FILE      *INP,*EXP;
 void      nrerror();
 

  printf("Read Expression files\n"); 
  count=0;
  for (j=0;j<size;j++){xprobe[j].split=0;}
  INP=fopen(fnot,"r");
  if (!INP){nrerror("Annotation File Not Found");}else{printf("\nReading Annotation File: %s\n",fnot);}
  while(!feof(INP)) {
    memset(buffer,'\0',MAXS);
    memset(a,'\0',MAXNAME);
    memset(b,'\0',MAXNAME);
    fgets(buffer,MAXS,INP);
    sscanf(buffer,"%s%s",a,b);
    for (j=0;j<size;j++){
      if (!strcmp(protein[j].name1,a) || !strcmp(protein[j].name2,a)){
       count++;
       strncpy(xprobe[j].fragment[xprobe[j].split].name,b,MAXNAME);
       if (xmethod[7]==1)printf("Node[%d]\t%s\tProbe[%d]\t%s\n",j,protein[j].name1,xprobe[j].split,xprobe[j].fragment[xprobe[j].split].name);
       xprobe[j].split += 1;
       if (xprobe[j].split >= MAXEST) {
	       printf("Node[%d]\t%s\tProbe[%d]\t%s\n",j,protein[j].name1,xprobe[j].split,xprobe[j].fragment[xprobe[j].split].name);
	       printf("Too many fragments for protein %d \n",xprobe[j].split);

	       nrerror("Error: increase MAXEST");
       }
       break;
       }
      }
    }
  fclose(INP);


  if (count==0){nrerror("Annotation file not able to find protein codes\n");}

  for (n=0;n<number_phe;n++){
   count=0;
   INP=fopen(fphe[n],"r");
   if (!INP){nrerror("No phenotype data");}else{printf("\nReading Phenotype[%d]: %s\n",n,fphe[n]);}
   s=0;
   while(!feof(INP)) {
     if (s==number_samples[n]){break;}
     memset(buffer,'\0',MAXS);
     fgets(buffer,MAXS,INP);
     sscanf(buffer,"%s",array);
     printf(" => Reading Sample[%d]: %s\n",s+1,array);
     EXP=fopen(array,"r");
     while(!feof(EXP)){
      memset(buffer,'\0',MAXS);
      memset(a,'\0',MAXNAME);
      fgets(buffer,MAXS,EXP);
      sscanf(buffer,"%s%f",a,&exp);
      for (j=0;j<size;j++){
      for (jj=0;jj<xprobe[j].split;jj++){
       if (!strcmp(xprobe[j].fragment[jj].name,a) ){
	count++;
        xprobe[j].fragment[jj].phenotype[n].sample[s]=exp;
        xprobe[j].fragment[jj].phenotype[n].defined[s]=1;
        if (xmethod[7]==1)printf("\tNode[%d]Probe[%d][%d][%d]: \t%s  \tExpression: %10.3e\tMask: %5d\n",j,jj,n,s,xprobe[j].fragment[jj].name,xprobe[j].fragment[jj].phenotype[n].sample[s],xprobe[j].fragment[jj].phenotype[n].defined[s]);
        break;
       }
      }}
     }
     fclose(EXP);
     s++;
   }
   if (count==0){nrerror("All files of phenotype were enable to find array codes\n");}
   fclose(INP);
   if (s!=number_samples[n])nrerror("Inconsistent number of samples of phenoptype");
 }

}


