#ifndef _KOMPLEX_H_
#define _KOMPLEX_H_

struct komplex {double re; double im;};
typedef struct komplex komplex;

void 	komplex_print	(const char *s, komplex z);
void 	komplex_set	    (komplex *z, double x, double y);
komplex komplex_new 	(double x, double y);
komplex komplex_add	    (komplex a, komplex b);
komplex komplex_sub	    (komplex a, komplex b);

int     komplex_equal   (komplex a, komplex b, double acc, double eps);

komplex komplex_mul     (komplex a, komplex b);
komplex komplex_div     (komplex a, komplex b);

komplex komplex_conj    (komplex z);
double komplex_abs      (komplex z);
komplex komplex_exp     (komplex z);
komplex komplex_sin     (komplex z);
komplex komplex_cos     (komplex z);
komplex komplex_sqrt    (komplex z);

#endif
