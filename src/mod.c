int mod(i,j)
int i,j;
{
 int   k,m;
 float r;
 r=(float)(i)/((float)(j));
 m=(int)(r);
 if (j*m == i ){k=1;}else{k=0;}

 return k;
}
