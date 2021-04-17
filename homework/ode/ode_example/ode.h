void driver(
	int  n, /* y[n] */
	void f(int n,double t,double*y,double*dydt), /* dy/dt=f(t,y) */
	double a,              /* the start-point a */
	double*y,                    /* y(a) -> y(b) */
	double b,              /* the end-point of the integration */
	double h,                    /* initial step-size */
	double acc,            /* absolute accuracy goal */
	double eps             /* relative accuracy goal */
);
