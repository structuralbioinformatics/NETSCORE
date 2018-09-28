#include "ppin.h"
abundance AssignNode(copy,proteome,size,name)
abundance copy[MAXP];
char      proteome[MAXP][MAXNAME];
char      name[MAXNAME];
int       size;
{

  int       i,j,skp;
  abundance cp,Cnull;
 
  Cnull.cell=0.0;
  Cnull.rms=0.0;
  Cnull.score=0.0;
  for (j=0;j<LOCAL;j++){ Cnull.local[j]=0.0;}

  skp=0;
  for (i=0;i<size;i++){
    if (!strcmp(proteome[i],name)){
       cp=copy[i];
       skp=1;
       break;
    }
  }
  if (skp==0){cp=Cnull;}

  return   cp;

}
