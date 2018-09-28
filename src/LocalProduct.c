#include "ppin.h"
float LocalProduct(i,j,dnode)
int   i,j;
node *dnode;
{
 float x;
 int   loc;

  x=0.0;
  
  for (loc=0;loc<LOCAL;loc++){
     //x+= (dnode[i].copy.local[loc]*dnode[i].copy.cell) * (dnode[j].copy.local[loc]*dnode[j].copy.cell);
     //printf("Local Nodes I %d and J %d test LOC %d X= %f = %f x %f \n",i,j,loc,x,dnode[i].copy.local[loc],dnode[j].copy.local[loc]);
     x+= (dnode[i].copy.local[loc]) * (dnode[j].copy.local[loc]);
  }

  if (x>0)x=1.0;
  
 return x;
}
