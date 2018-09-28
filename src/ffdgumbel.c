#include <math.h>

float ffdgumbel(x,a,na)
float x,a[];
int na;
{
	int i;
	float y,ex,arg,aa,alpha,beta,fac,g,w,ai;
 	y=0.0;
        aa=0;
        for (i=1;i<=na-1;i+=3){aa+=sqrt(a[i]*a[i]);}
        if (na<=3){aa=1;a[1]=1;}       
        for (i=1;i<=na-1;i+=3)
	{
                alpha =  a[i+1];
                beta  =  a[i+2];
		arg   =  (x-alpha)/beta;
		ex    =  exp(-arg);
                g     =  exp(-arg-ex)/beta;
                ai    =  sqrt(a[i]*a[i]);
                w     =  ai/aa;
                fac   =  g*(1-ex);
		y    +=  w*g;
	}

        return y;

}
