#include "ppin.h"
#define EPS 1.0e-12

int    InitParameters(rnd,x,y,ndata,a,ma,lista,mfit)
float *x,*y,*a;
int   *rnd,ndata,ma,*lista,mfit;
{

 float aa,dd,arg,g,ex,*xx,*yy,*wk1,*wk2,*alpha,*beta,*beta2,*zeta,*fvector();
 int   nn,j,i,i0,i1,j1,k,jmax,jj,m,*max,*ivector();
 void  derivative(),free_fvector(),free_ivector(),sort2();


// Initialize Standard parameters
 for (k=0;k<=ma;k++)lista[k]=k;
 a[0]=0.0;
 for (k=1;k<=ma;k+=3){ a[k]=0; a[k+1]=0.0; a[k+2]=1.0;}

 if (rnd[6]==1 || rnd[6]==3){
  for (k=1;k<=mfit;k+=3){ a[k]=1.0;}
  return;
 }

// when using PDF instead of CDF function
 
 nn =ndata*2;
 xx=fvector(0,nn+1);
 yy=fvector(0,nn+1);
 wk1=fvector(0,nn+1);
 wk2=fvector(0,nn+1);
 max=ivector(0,MAXXPDF);
 alpha=fvector(0,mfit);
 beta=fvector(0,mfit);
 beta2=fvector(0,mfit);
 zeta=fvector(0,mfit);

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
   //printf("%d] Y[%f]=%f dY[%d]=%f  dY[%d]=%f \n",k,x[k],y[k],i-1,yy[i-1],i,yy[i]);
   k++;
 }
 for (i=1;i<=jmax;i++){
  wk1[i]=y[max[i]];
  wk2[i]=(float)max[i];
  //printf("W1[%d]=%f W2[%d]=%f X=%f \n",i,y[max[i]],i,wk2[i],x[max[i]]);
 }

 if (jmax>1) sort2(jmax,wk1,wk2);
 k=1;
 for (i=jmax;i>=1;i--) {max[k]=(int)wk2[i];k++;}

 jj=0; 
 printf("\nMaxima found on numerical function of PDF:\n");
 for (i=1;i<=jmax;i++){
  printf("Maximum Bin[%5d]\tX=%f\t PDF= %e\n",max[i],x[max[i]],y[max[i]]);
  if (y[max[i]]<(y[max[1]]/1000.0) && jj==0) {jj=i-1;printf("Selected Maximum Number of functions: %d\n",jj);}
  //else{printf("Y[%d]=%f Y[%d]/1000=%f/1000 JJ=%d \n",max[i],y[max[i]],max[1],y[max[1]],jj);}
 }
 if (jj==0)jj=jmax;
 if (jj>mfit/3) jj=mfit/3;
 printf("Limitted number of non-null starting functions: %d\n\n",jj);

 k=1;

 aa=0.0;
 for (i=1;i<=jj;i++){aa+=fabs(y[max[i]]);}
 if (aa==0)aa=1.0;

 if (rnd[6]==2){
   for (i=1;i<=jj;i++){
     alpha[i]= x[max[i]];
     beta[i]= exp(-exp(0))/(jj*y[max[i]]);
   }
   for (i=1;i<=jj;i++){
     zeta[i]=0.0; 
     for (j=1;j<=jj;j++){
        if (j!=i){
            arg      =  (alpha[i]-alpha[j])/beta[j];
     	    ex       =  exp(-arg);
            g        =  exp(-arg-ex)/beta[j];
            zeta[i] +=   g/jj;
            }
        }
   }
   dd=1.0;
   m=0;
   printf("Iteration to calculate initial values of beta in Gumbel\n");
   while (dd>EPS && m<10){
     for (i=1;i<=jj;i++){
       if (y[max[i]]>zeta[i]) {beta2[i]= exp(-exp(0))/(jj*(y[max[i]]-zeta[i]));}
       else {beta2[i]=beta[i];}
       }
     dd=0;
     for (i=1;i<=jj;i++)dd+=(beta[i]-beta2[i])*(beta[i]-beta2[i]);
     dd=sqrt(dd);
     for (i=1;i<=jj;i++){printf("Beta[%d]=%f Beta2[%d]=%f Zeta[%d]=%f X[%d]=%f Y[%d]=%f\n",i,beta[i],i,beta2[i],i,zeta[i],max[i],x[max[i]],max[i],y[max[i]]);}
     printf("RMSD (beta) = %e\n\n",dd);
     for (i=1;i<=jj;i++)beta[i]=beta2[i];
     for (i=1;i<=jj;i++){
       zeta[i]=0.0; 
       for (j=1;j<=jj;j++){
          if (j!=i){
            arg      =  (alpha[i]-alpha[j])/beta[j];
     	    ex       =  exp(-arg);
            g        =  exp(-arg-ex)/beta[j];
            zeta[i] +=   g/jj;
            }
          }
     }
     m++;
   }
   for (i=1;i<=jj;i++){
     //a[k]=y[max[i]]/exp(-exp(0));
     a[k]=1.0/jj;
     //a[k]=fabs(y[max[i]])/aa;
     a[k+1]=alpha[i];
     a[k+2]=beta[i];
     k+=3;
   }
   for (i=jj+1;i<=jmax;i++){
     a[k]=EPS/jmax;
     a[k+1]= x[max[i]];
     a[k+2]= exp(-exp(0))/y[max[i]];
     k+=3;
   }
 } 

 if (rnd[6]==4){
   for (i=1;i<=jj;i++){
     alpha[i]= y[max[i]];
     i0=max[i];
     i1=max[i]+1;
     j1=max[i]-1;
     beta[i] = 0.5*(1.0/sqrt(log(y[i0]/(EPS+y[i1])))+1.0/sqrt(log(y[i0]/(EPS+y[j1]))));
   }
   for (i=1;i<=jj;i++){
     //a[k]=fabs(y[max[i]])/aa;
     a[k]=alpha[i];
     a[k+1]=x[max[i]];
     a[k+2]=beta[i];
     printf("Starting parameters A[%d]=%f Xo[%d]=%f Sigma[%d]=%f\n",i,alpha[i],i,x[max[i]],i,beta[i]);
     k+=3;
   }
 }
 
 

 free_fvector(xx,0,nn+1);
 free_fvector(yy,0,nn+1);
 free_fvector(wk1,0,nn+1);
 free_fvector(wk2,0,nn+1);
 free_fvector(alpha,0,mfit);
 free_fvector(beta,0,mfit);
 free_fvector(beta2,0,mfit);
 free_fvector(zeta,0,mfit);
 free_ivector(max,0,MAXXPDF);

 return jmax;

}

#undef EPS

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


