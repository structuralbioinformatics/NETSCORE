#include "ppin.h"
float  Score(d_scr,n_scr,inode,maxmin)
int      *n_scr;
gradient **d_scr;
float     maxmin[];
int      inode;
{
 float x,y,z;
 float Transform();

 y=z=x=0.0; 

 y=Transform(inode,maxmin,n_scr,d_scr);
 if (y<=0){return x;}else{z=y;}
 
 return z;

}


