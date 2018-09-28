#include <math.h>

void fgumbel(x,a,y,dyda,na)
float x,a[],*y,dyda[];
int na;
{
	int i;
	float ex,arg,aa,g,w,ai,ff;
 	*y=0.0;
        aa=0;
        for (i=1;i<=na-1;i+=3){aa+=sqrt(a[i]*a[i]);}
        if (na<=3) {a[1]=1.0;aa=1;}
        for (i=1;i<=na-1;i+=3)
	{
		arg =  (x-a[i+1])/a[i+2];
		ex  =  exp(-arg);
                g   =  exp(-ex);
                ai  =  sqrt(a[i]*a[i]);
                w   =  ai/aa;
                ff  =  (a[i]/ai)*(1-w);
                *y +=  w*g;
                dyda[i]   = g*ff/aa;
		dyda[i+1] = -w*g*ex/a[i+2];
                dyda[i+2] = w*g*ex/(a[i+2]*a[i+2]);
	}
        //if (na<=3) {dyda[1]=1.0e-12;}

}
