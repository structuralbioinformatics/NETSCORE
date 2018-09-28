#include <math.h>

float ffgumbel(x,a,na)
float x,a[];
int na;
{
	int i;
	float y,ex,arg,aa,g,w,ai;
 	y=0.0;
        aa=0;
        for (i=1;i<=na-1;i+=3){aa+=sqrt(a[i]*a[i]);}
        if (na<=3){aa=1;a[1]=1;}
        for (i=1;i<=na-1;i+=3)
	{
		arg =  (x-a[i+1])/a[i+2];
		ex  =  exp(-arg);
                g   =  exp(-ex);
                ai  =  sqrt(a[i]*a[i]);
                w   =  ai/aa;
		y  +=  w*g;
	}

        return  y;

}
