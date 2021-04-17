#include<math.h>
#include<stdio.h>
#define FOR(i) for(int i=0;i<n;i++)

void rkstep23(
	int n,
	void f(int n,double t,double y[],double dydt[]),
	double t, double y[], double h,
	double yh[], double dy[]
){
	double k0[n]; f(n,t,y,k0);
	double y1[n]; FOR(i) y1[i]=y[i]+(0.5 *h)*k0[i];
	double k1[n]; f(n,t+0.5*h,y1,k1);
	double y2[n]; FOR(i) y2[i]=y[i]+(0.75*h)*k1[i];
	double k2[n]; f(n,t+0.75,y2,k2);
	double ka[n]; FOR(i) ka[i]=(2*k0[i]+3*k1[i]+4*k2[i])/9;
	FOR(i) yh[i]=y[i]+h*ka[i];
	FOR(i) dy[i]=(ka[i]-k1[i])*h;
}

#define PRINT(x) fprintf(stderr,"%9.3g ",x)
#define TRACE(t,y) PRINT(t);FOR(i)PRINT(y[i]);fprintf(stderr,"\n")

void driver(
	int  n, /* y[n] */
	void f(int n,double t,double*y,double*dydt), /* dy/dt=f(t,y) */
	double a,              /* the start-point a */
	double*y,                    /* y(a) -> y(b) */
	double b,              /* the end-point of the integration */
	double h,                    /* initial step-size */
	double acc,            /* absolute accuracy goal */
	double eps             /* relative accuracy goal */
){
	double t=a;
	TRACE(t,y);
	while(t<b){
		if(t+h>b)h=b-t;
		double yh[n],dy[n];
		rkstep23(n,f,t,y,h,yh,dy);
		double sum=0; FOR(i)sum+=y[i]*y[i];
		double norm_y=sqrt(sum);
		sum=0; FOR(i)sum+=dy[i]*dy[i];
		double err=sqrt(sum);
		double tol=(acc+eps*norm_y)*sqrt(h/(b-a));
		if(err<tol){
			t=t+h;
			FOR(i) y[i]=yh[i];
			TRACE(t,y);
		}
		if(err>0) h*=0.95*pow(tol/err,0.25);
		else h*=2;
	}//while
}
