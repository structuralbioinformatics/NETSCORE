#include <math.h>

void polint(xa,ya,n,x,y,dy)
float xa[],ya[],x,*y,*dy;
int n;
{
	int i,m,ns=1;
	float den,dif,dift,ho,hp,w;
	float *c,*d,*fvector();
	void nrerror(),free_fvector();

	dif=fabs(x-xa[1]);
	c=fvector(1,n);
	d=fvector(1,n);
	for (i=1;i<=n;i++) {
		if ( (dift=fabs(x-xa[i])) < dif) {
			ns=i;
			dif=dift;
		}
		c[i]=ya[i];
		d[i]=ya[i];
	}
	*y=ya[ns--];
	for (m=1;m<n;m++) {
		for (i=1;i<=n-m;i++) {
			ho=xa[i]-x;
			hp=xa[i+m]-x;
			w=c[i+1]-d[i];
			if ( (den=ho-hp) == 0.0) nrerror("Error in routine POLINT");
			den=w/den;
			d[i]=hp*den;
			c[i]=ho*den;
		}
		*y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
	}
	free_fvector(d,1,n);
	free_fvector(c,1,n);
}
