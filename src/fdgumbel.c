#include <math.h>

#define EPS 1.0e-6

void fdgumbel(x,a,y,dyda,na)
float x,a[],*y,dyda[];
int na;
{
	int i;
	float ex,arg,aa,alpha,beta,fac,g,ff,ai,w;
 	*y=0.0;
        aa=0;
        for (i=1;i<=na-1;i+=3){aa+=sqrt(a[i]*a[i]);}
        if (aa==0)aa=1;
        if (na<=3){aa=1;a[1]=1;}
        for (i=1;i<=na-1;i+=3)
	{
                if (a[i+2]<EPS)a[i+2]=EPS;
                alpha  =  a[i+1];
                beta   =  a[i+2];
		arg    =  (x-alpha)/beta;
		ex     =  exp(-arg);
                g      =  exp(-arg-ex)/beta;
                ai     =  sqrt(a[i]*a[i]); 
                w      =  ai/aa;
                if (ai==0)ai=1;
                ff     =  (a[i]/ai)*(1-w);
                fac    =  g*(1-ex);
		*y    +=  w*g;
		dyda[i]   = g*ff/aa;
		dyda[i+1] = w*fac/beta;
                dyda[i+2] = -w*(g - fac*arg)/beta;
	}

        //if (na<=3)dyda[1]=0;

}
#undef EPS
