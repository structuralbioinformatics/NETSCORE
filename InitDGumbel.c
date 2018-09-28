#include "ppin.h"
#define EPS 1.0e-12

float InitDGumbel(weight,alpha,beta,jj,mfit,x,y,max)
float *weight,*alpha,*beta,*x,*y;
int    jj,mfit,*max;
{
 float  *beta2,*xhi,*xhi2,*zeta,*gamma,*delta,*fvector();
 float   dd,dd1,ddd,ddg,dg,dddg,dg1,e,ex,arg,g,zz;
 int     i,j,k,m,mm,mmm,mmmm,i0,i1;
 void    free_fvector();

  beta2=fvector(0,mfit);
  zeta=fvector(0,mfit);
  gamma=fvector(0,mfit);
  delta=fvector(0,mfit);
  xhi=fvector(0,mfit);
  xhi2=fvector(0,mfit);

  e=exp(exp(0));
// Initialize variables beta xhi alpha
  for (i=1;i<=jj;i++){
     i0=max[i];
     i1=i0+1;
     alpha[i]= x[i0];
     xhi[i]= e*y[i0];
     if (y[i1]<EPS){
       x[i1]=(x[i0]+x[i1])/2;
       y[i1]=y[i0]/2;
       }
     beta[i]=(x[i1]-x[i0])/(2+log(y[i0]/y[i1]));
     if (beta[i]<EPS)beta[i]=EPS;
   }
//Initialize zeta
   for (i=1;i<=jj;i++){
     zeta[i]=0.0;
     for (j=1;j<=jj;j++){
        if (j!=i && beta[j]>EPS){
            arg      =  (alpha[i]-alpha[j])/beta[j];
            ex       =  exp(-arg);
            g        =  exp(-arg-ex);
            zeta[i] +=  xhi[j] * g;
        //printf ("Test Zeta Xhi[%d]=%e Beta[%d]=%e Zeta[%d]=%e Ymax=%e ARG=%e EX=%e G=%e\n",i,xhi[i],i,beta[i],i,zeta[i],y[max[i]],arg,ex,g);
            }
        }
   }

//for (i=1;i<=jj;i++){printf ("Xhi[%d]=%e Beta[%d]=%e Zeta[%d]=%e Ymax=%e \n",i,xhi[i],i,beta[i],i,zeta[i],y[max[i]]);}
// Iteration to calculate beta, xhi, and zeta
   dd1=ddd=dd=1.0;
   mmmm=m=0;
   printf("Iteration to calculate initial values of beta in Gumbel\n");
   while (dd>EPS && m<100){
     //Calculate gamma
     for (i=1;i<=jj;i++){
       gamma[i]=0;
       for (j=1;j<=jj;j++){
       if (j!=i && beta[j]>EPS){
          gamma[i]+=beta[j]*xhi[j];
          }}
       printf("Test on gamma[%d]=%e + %e\n",i,gamma[i],beta[i]*e*y[max[i]]);
       }
     //Initialize xhi2 
     for (i=1;i<=jj;i++){xhi2[i]=xhi[i];}
     // Initialize beta2
     for (i=1;i<=jj;i++){beta2[i]=beta[i];}
     //ReCalculate gamma with beta2 and xhi2
     mmm=mm=0;
     dddg=ddg=1;
     dg1=dg=0;
     while (mm<50 && ddg>EPS){
      // recalculate zeta
      for (i=1;i<=jj;i++){
       zeta[i]=0.0;
       for (j=1;j<=jj;j++){
          if (j!=i&& beta2[j]>EPS){
            arg      =  (alpha[i]-alpha[j])/beta2[j];
            ex       =  exp(-arg);
            g        =  exp(-arg-ex);
            zeta[i] +=  xhi2[j] * g;
            }
          }
      }
      //recalculate xhi as xhi2 
      for (i=1;i<=jj;i++){
       if(y[max[i]]>zeta[i]) {xhi2[i]=e*(y[max[i]]-zeta[i]);}
       else                  {xhi2[i]=xhi[i];}
       //printf ("Xhi2[%d]=%e Beta2[%d]=%e Zeta[%d]=%e Ymax=%e \n",i,xhi2[i],i,beta2[i],i,zeta[i],y[max[i]]);
       }
      for (i=1;i<=jj;i++){
       gamma[i]=0;
       for (j=1;j<=jj;j++){
       if (j!=i && beta2[j]>EPS){
          //printf("Adding on gamma[%d]=%e + %e\n",i,gamma[i],beta2[j]*xhi2[j]);  
          gamma[i]+=beta2[j]*xhi2[j];
          }}
       if (gamma[i]>1) {
           printf("Second(iterative) Test on gamma[%d]=%e + %e ",i,gamma[i],beta2[i]*e*y[max[i]]);  
           gamma[i]=1;
           printf("Reset to gamma[%d]=%e + %e\n",i,gamma[i],beta2[i]*e*y[max[i]]) ;
          }else{printf("Second(iterative) Test on gamma[%d]=%e + %e\n",i,gamma[i],beta2[i]*e*y[max[i]]);}
       }
     // Recalculate beta2
      for (i=1;i<=jj;i++){
       i0=max[i];
       if (gamma[i]>1 && zeta[i]>y[i0]){ beta2[i]= (1-gamma[i])/(e*(y[i0]-zeta[i]));}else{beta2[i]=beta[i];}
       //printf("New Beta2[%d]=%e Zeta[%d]=%e Ymax=%e Gamma[%d]=%e\n",i,beta2[i],i,zeta[i],y[max[i]],i,gamma[i]);
       if (beta2[i]<EPS) beta2[i]=beta[i];
       }
     // Error criteria
      dg1=dg;
      dg=0;
      for (i=1;i<=jj;i++){
        dg+=beta2[i]*xhi2[i];
       }
      dddg=fabs(dg1-dg);
      ddg=sqrt((1-dg)*(1-dg));
      if (dddg<EPS) {mmm++;}else{mmm=0;}
      printf("Iteration Gamma Test[%d]= %e Error= %e Difference=%e Repeat=%d\n",mm,dg,ddg,dddg,mmm);
      if (mmm>3)mm=50;
      mm++;
     }
     // Calculate delta
     for (i=1;i<=jj;i++){
       i0=max[i];
       i1=i0+1;
       arg=(x[i1]-alpha[i])/beta2[i];
       delta[i]=exp(-arg);
       } 
     // Recalculate beta2
     for (i=1;i<=jj;i++){
       i0=max[i];
       i1=i0+1;
       beta2[i]=(x[i1]-x[i0])/(1+log(y[i0]/y[i1])+delta[i]);
       if (beta2[i]<EPS)beta2[i]=EPS;
       }    
     dd1=dd;
     dd=0;
     for (i=1;i<=jj;i++)dd+=(beta[i]-beta2[i])*(beta[i]-beta2[i]);
     for (i=1;i<=jj;i++)dd+=(xhi[i]-xhi2[i])*(xhi[i]-xhi2[i]);
     dd=sqrt(dd);
     ddd=fabs(dd1-dd);
     if (ddd<EPS){mmmm++;}else{mmmm=0;}
     for (i=1;i<=jj;i++){ printf("Beta[%d]=%e Beta2[%d]=%e Xhi[%d]=%e Xhi2[%d]=%e Zeta[%d]=%e Delta[%d]=%e Gamma[%d]=%e X[%d]=%f Y[%d]=%f\n", i,beta[i],i,beta2[i],i,xhi[i],i,xhi2[i],i,zeta[i],i,delta[i],i,gamma[i],max[i],x[max[i]],max[i],y[max[i]]);} 
     printf("RMSD (beta & xhi) = %e  Difference=%e Repeat=%d\n\n",dd,ddd,mmmm);
     for (i=1;i<=jj;i++)beta[i]=beta2[i];
     for (i=1;i<=jj;i++)xhi[i]=xhi2[i];
     // recalculate zeta
     for (i=1;i<=jj;i++){
       zeta[i]=0.0;
       for (j=1;j<=jj;j++){
          if (j!=i&& beta[j]>EPS){
            arg      =  (alpha[i]-alpha[j])/beta[j];
            ex       =  exp(-arg);
            g        =  exp(-arg-ex);
            zeta[i] +=  xhi[j] * g;
            }
          }
     }
     if (mmmm>3)m=100;
     m++;
   }
   
 for (i=1;i<=jj;i++) weight[i]=xhi[i]*beta[i];

 free_fvector(delta,0,mfit);
 free_fvector(zeta,0,mfit);
 free_fvector(beta2,0,mfit);
 free_fvector(xhi,0,mfit);
 free_fvector(xhi2,0,mfit);
 free_fvector(gamma,0,mfit);

 return dd;

}
#undef EPS
