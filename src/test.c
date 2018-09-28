#include        <stdio.h>
#include        <math.h>
#include        <ctype.h>
#include        <string.h>
#include        <sys/types.h>
#include        <sys/wait.h>
#include        <malloc.h>
#include        <openssl/des.h>
#include        <stdlib.h>
#include        <unistd.h>
#include        <errno.h>

main (int argc, char *argv[])
{
 void    InitParameters();
 int     i,k,rnd[10],*lista,ndata,ma,mfit;
 float   *x,*y,*param;
 float     *fvector();
 int     *ivector();
 void    free_fvector(),free_ivector();


    ndata=24;
    ma=50;
    mfit=9;
    rnd[6]=1;
    x=fvector(0,ndata+1);
    y=fvector(0,ndata+1);
    lista=ivector(0,ma+1);
    param=fvector(0,ma+1);

    x[1]=-1.459139; y[1]=7.157560e-03;
    x[2]=-1.202878; y[2]=4.497664e-02;	
    x[3]=-0.946618; y[3]=4.508246e-01; 
    x[4]=-0.690357; y[4]=2.064308e-01 ;
    x[5]=-0.434096; y[5]=2.652290e-02 ;
    x[6]=-0.177835; y[6]=3.080316e-02 ;
    x[7]=0.078426;  y[7]=2.736799e+00;
    x[8]=0.334687;  y[8]=1.689626e-01;
    x[9]=0.590948;  y[9]=2.329539e-02;
    x[10]=0.847208; y[10]=2.350564e-02;
    x[11]=1.103468; y[11]=4.090723e-02;
    x[12]=1.359729; y[12]=5.936255e-02;
    x[13]=1.615990; y[13]=6.855285e-02;
    x[14]=1.872251; y[14]=5.051991e-03;
    x[15]=5.716163; y[15]=7.016586e-05;
    x[16]=5.972425; y[16]=7.016586e-05;
    x[17]=6.228685; y[17]=0.000000e+00;
    x[18]=6.484947; y[18]=0.000000e+00;
    x[19]=6.741207; y[19]=0.000000e+00;
    x[20]=6.997468; y[20]=7.017217e-05;
    x[21]=7.253728; y[21]=7.017217e-05;
    x[22]=7.509990; y[22]=0.000000;
    x[23]=7.766250; y[23]=0.000000e+00;

   InitParameters(rnd,x,y,ndata,param,ma,lista,mfit);

   for (k=0;k<=mfit;k++){printf("Parameter[%5d] =\t%e\n",k,param[k]);}
   for (k=0;k<=mfit;k++){printf("Lista[%5d] =\t%d\n",k,lista[k]);}

   free_fvector(x,0,ndata+1);
   free_fvector(y,0,ndata+1);
   free_ivector(lista,0,ma+1);
   free_fvector(param,0,ma+1);

} 


/*
main (int argc, char *argv[])
{
 void shuffle();
 int  i,j,list[15],list2[15],init,size,idum;
 float ran3(),test;
 
 for (i=0;i<15;i++){list2[i]=list[i]=i+1;}
 size=15;
 for (j=0;j<2;j++){
 idum=size*(j+1);
 printf("Before Init=%d Size= %d Idum= %d\n",init,size,idum);
 test=ran3(&idum)*1000*size;
 init= (int)test * size;
 printf("Before2 Init=%d Size= %d Test=%f \n",init,size, test);
 shuffle(list2,size,init);
 for (i=0;i<15;i++){printf("List[%d]= %d  List2[%d]=%d\n",i,list[i],i,list2[i]);}
 }

}
*/

/*
#define FUNC(x) ((*funk)(x,a,na))

main (int argc, char *argv[])
{
 void (*funcs)(), fgauss();
 float(*funk)(), ffgauss();
 float x,v[5],a[4],ymod,y,dyda[4];
 int   na,i;
  na=3;
  funcs=&fgauss;
  funk=&ffgauss;
  v[0]=v[5]=0;
  v[1]=v[2]=0.5;
  v[3]=1;
  a[0]=0;
  a[1]=a[2]=a[3]=1.0;
  ymod=0.0;
  for (i=0;i<6;i++){
   x=v[i];
   (*funcs)(x,a,&ymod,dyda,na);
    printf("y=%f \n",ymod);
    ymod=FUNC(x);
    printf("y2=%f \n",ymod);
   
  }
  printf("a[1]=%f a[2]=%f a[3]=%f \n",a[1],a[2],a[3]);

}

#undef FUNC
*/

void fgauss(x,a,y,dyda,na)
float x,a[],*y,dyda[];
int na;
{
	int i;
	float fac,ex,arg;

        printf("entre gauss\n");
	*y=0.0;
	for (i=1;i<=na-1;i+=3)
	{
		arg=(x-a[i+1])/a[i+2];
		ex=exp(-arg*arg);
		fac=a[i]*ex*2.0*arg;
		*y += a[i]*ex;
		dyda[i]=ex;
		dyda[i+1]=fac/a[i+2];
		dyda[i+2]=fac*arg/a[i+2];
	}
}
/*
void shuffle(theArr,size,init)
int *theArr;
int size;
int init;
{
	   int temporary, randomNum, last, kdum,idum;
           float ran3(),test;
           printf("Init %d Size %d\n",init,size);
	   for (last = size; last > 1; last--)
	   {
              idum   = last*init;
              test   = ran3(&idum);
              kdum   = (int) (test  * last);
              printf("Idum %d last %d Kdum %d test %f \n",idum,last,kdum,test);
              test   = ran3(&kdum);
	      randomNum =  (int) (test * last) % last;
              printf("Idum %d last %d Kdum %d randomNum %d test %f\n",idum,last,kdum, randomNum,test );
	      temporary = theArr[randomNum];
	      theArr[randomNum] = theArr[last - 1];
	      theArr[last - 1] = temporary;
	   }
}// end shuffleElements( )
*/

#include <math.h>

float ffgauss(x,a,na)
float x,a[];
int na;
{
	int i;
	float y,ex,arg;

	y=0.0;
	for (i=1;i<=na-1;i+=3)
	{
		arg=(x-a[i+1])/a[i+2];
		ex=exp(-arg*arg);
		y += a[i]*ex;
	}

        return y;

}
#include "ppin.h"

void    InitParameters(rnd,x,y,ndata,a,ma,lista,mfit)
float *x,*y,*a;
int   *rnd,ndata,ma,*lista,mfit;
{

 float *xx,*yy,*wk1,*wk2,*fvector();
 int   nn,j,i,k,jmax,*max,*ivector();
 void  derivative(),free_fvector(),free_ivector(),sort2();
 
 nn =ndata*2;
 xx=fvector(0,nn+1);
 yy=fvector(0,nn+1);
 wk1=fvector(0,nn+1);
 wk2=fvector(0,nn+1);
 max=ivector(0,MAXXPDF);

 derivative(x,y,xx,yy,ndata,nn);

 j=1;
 k=1;
 for (i=2;i<=nn;i+=2){
   if (yy[i]>0 && yy[i+2]<0){
     jmax=j;
     max[j]=k+1;
     j++;
     if (j>=MAXXPDF) {printf("Too many maxima on function, increase MAXXPDF\n");break; }
   }
   printf("%d] Y[%f]=%f dY[%d]=%f  dY[%d]=%f \n",k,x[k],y[k],i-1,yy[i-1],i,yy[i]);
   k++;
 }
 for (i=1;i<=jmax;i++){
  wk1[i]=y[max[i]];
  wk2[i]=(float)max[i];
  printf("W1[%d]=%f W2[%d]=%f X=%f \n",i,y[max[i]],i,wk2[i],x[max[i]]);
 }

 if (jmax>1) sort2(jmax,wk1,wk2);
 k=1;
 for (i=jmax;i>=1;i--) {max[k]=(int)wk2[i];k++;}
 if (jmax>mfit/3) jmax=mfit/3;
 
 for (k=0;k<=ma;k++)lista[k]=k;
 a[0]=0.0;
 for (k=1;k<=ma;k+=3){ a[k]=0.0; a[k+1]=0.0; a[k+2]=1.0;}

 if (rnd[6]==1 || rnd[6]==3){
  for (k=1;k<=mfit;k+=3){ a[k]=1.0;}
 }

 k=1;

 if (rnd[6]==2){
   for (i=1;i<=jmax;i++){
     a[k]=1.0;
     a[k+1]=x[max[i]];
     a[k+2]= exp(-exp(0))/y[max[i]];
     k+=3;
   }
 } 

 if (rnd[6]==4){
   for (i=1;i<=jmax;i++){
     a[k]=1.0;
     a[k+1]=x[max[i]];
     k+=3;
   }
 }
 
 

 free_fvector(xx,0,nn+1);
 free_fvector(yy,0,nn+1);
 free_fvector(wk1,0,nn+1);
 free_fvector(wk2,0,nn+1);
 free_ivector(max,0,MAXXPDF);

}

void  derivative(x,y,xx,yy,n,m)
float   *x,*y,*xx,*yy;
int      n,m;
{

 int   i,j;
 float dd,ddn,dy,dx,dyn,dxn;

 j=1;
 yy[0]=0;
 xx[0]=0;
 for (i=1;i<n;i++){
  dd=ddn=0;
  dy=y[i+1]-y[i];
  dx=x[i+1]-x[i];
  if (i<n-1){ 
   dyn=y[i+2]-y[i+1];
   dxn=x[i+2]-x[i+1];
   if (dxn>0) ddn=dyn/dxn;
  }
  if (dx>0) dd=dy/dx;
  xx[j]=x[i];
  if (i==1) yy[j]=dd;
  j++;
  xx[j]=(x[i+1]+x[i])/2;
  yy[j]=dd;
  j++;
  xx[j]=x[i+1];
  if (i==(n-1)){ 
    yy[j]=dd;
  }else{
    yy[j]= (dd+ddn)/2;
  }
  if (j>m)break;
 } 


}

void sort2(n,ra,rb)
int n;
float *ra,*rb;
{
	int l,j,ir,i;
	float rrb,rra;

	l=(n >> 1)+1;
	ir=n;
	for (;;) {
		if (l > 1) {
			rra=ra[--l];
			rrb=rb[l];
		} else {
			rra=ra[ir];
			rrb=rb[ir];
			ra[ir]=ra[1];
			rb[ir]=rb[1];
			if (--ir == 1) {
				ra[1]=rra;
				rb[1]=rrb;
				return;
			}
		}
		i=l;
		j=l << 1;
		while (j <= ir)	{
			if (j < ir && ra[j] < ra[j+1]) ++j;
			if (rra < ra[j]) {
				ra[i]=ra[j];
				rb[i]=rb[j];
				j += (i=j);
			}
			else j=ir+1;
		}
		ra[i]=rra;
		rb[i]=rrb;
	}
}


