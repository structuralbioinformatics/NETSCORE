#include <math.h>
#include <stdio.h>

#define FUNC(x) ((*func)(x,p,na))

float trapzd(func,a,b,n,p,na)
float a,b,*p;
float (*func)();	/* ANSI: float (*func)(float); */
int n,na;
{
	float x,tnm,sum,del;
	static float s;
	static int it;
	int j;

        if (a==b) return 0.0;
	if (n == 1) {
		it=1;
		return (s=0.5*(b-a)*(FUNC(a)+FUNC(b)));
	} else {
		tnm=it;
		del=(b-a)/tnm;
		x=a+0.5*del;
//		for (sum=0.0,j=1;j<=it;j++,x+=del) {sum += FUNC(x);printf("TRAPZD[%d] till %d B=%f A=%f Sum= %f Func(%e)=%e DEL=%e\n",j,it,b,a,sum,x,FUNC(x),del);}
		for (sum=0.0,j=1;j<=it;j++,x+=del) sum += FUNC(x);
		it *= 2;
		s=0.5*(s+(b-a)*sum/tnm);
		return s;
	}
}

#undef FUNC

