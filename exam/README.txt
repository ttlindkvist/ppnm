Examination project 13 - Two sided Jacobi algorithm for SVD

Author: Thomas Toft Lindkvist
AUID 643642 - Student number 201905635
35 mod 22 = 13

------------ DESCRIPTION ------------
As explored in this class, the elementary Jacobi eigenvalue algorithm only works on symmetric matrices.
But suppose we wanted to find the eigenvalues of a non-symmetric matrix, or the singular values of a 
non-square matrix, could we generalize the known algorithm? It turns out the answer is yes! 
The two sided Jacobi algorithm is such a generalization. Any tall matrix A can, via the implemented algorithm, 
be decomposed into A=U*D*VT, where U and V are orthogonal matrices and D is a diagonal matrix containing the singular values of A.

A short report with the mathematics involved is found in report.pdf

Convergence: The criterion for convergence is when two successive sweeps produce the same diagonalized matrix (within machine precision).

Computational cost: the algorithm cost is O(n^3) for a n*n matrix - loop over upper triangle is O(n^2) and the update is O(n). 
This is depicted in the timing.png plot, where the implemented algorithm is timed across a range of square matrices 
- and as seen the fit shows a n^x tendency, with x ≈ 3.

Versus GSL: As seen from the timing plot - this implementation is MUCH slower than the one-sided algorithm implemented in GSL.
Around 8-11 times slower in the usecases depicted here.


--------- OUTPUT FROM PROGRAM (see also file output.out)  -----------

EXAMINATION PROJECT 13
AUTHOR: Thomas Toft Lindkvist
AUID643642 - Student number: 201905635
35 mod 22 = 13
Two-sided Jacobi alg. for SVD

Starting with 7x5 matrix: 
 0.0334699  0.3299642  0.6906357  0.4224867  0.2062651 
 0.2501285  0.6365585  0.8636224  0.3016556  0.0249239 
 0.3649927  0.7653810  0.3178603  0.1357284  0.1067453 
 0.7606716  0.0816723  0.5509562  0.5646519  0.7432705 
 0.9817601  0.2188498  0.4543959  0.5186540  0.5738417 
 0.5552310  0.7285463  0.6288739  0.7961281  0.5837158 
 0.2747664  0.8295981  0.9136800  0.9654021  0.2520848 
SVD done with 6 sweeps
Check UDV^T=A: 
 0.0334699  0.3299642  0.6906357  0.4224867  0.2062651 
 0.2501285  0.6365585  0.8636224  0.3016556  0.0249239 
 0.3649927  0.7653810  0.3178603  0.1357284  0.1067453 
 0.7606716  0.0816723  0.5509562  0.5646519  0.7432705 
 0.9817601  0.2188498  0.4543959  0.5186540  0.5738417 
 0.5552310  0.7285463  0.6288739  0.7961281  0.5837158 
 0.2747664  0.8295981  0.9136800  0.9654021  0.2520848 
Matrix D: 
 0.2047131 -0.0000000 -0.0000000 -0.0000000  0.0000000 
-0.0000000  0.4226660 -0.0000000  0.0000000  0.0000000 
 0.0000000 -0.0000000  0.5919544  0.0000000  0.0000000 
-0.0000000  0.0000000 -0.0000000  1.1340104  0.0000000 
-0.0000000 -0.0000000  0.0000000  0.0000000  3.0799361 
Matrix V: 
-0.4054365 -0.2399272  0.5001257  0.6038095  0.4041591 
 0.2041337  0.3063950  0.6122151 -0.5433797  0.4408885 
 0.0806004 -0.6978795 -0.3410564 -0.3054001  0.5448661 
-0.4429836  0.5693436 -0.4878287  0.0038263  0.4915506 
 0.7689082  0.1933110 -0.1441204  0.4968589  0.3221361 
Matrix U: 
 0.3951709 -0.2566993 -0.4267653 -0.2344823  0.2628066 
-0.0797416 -0.6887586  0.1174316 -0.3924792  0.3274774 
 0.2727241  0.0544643  0.7789711 -0.2107779  0.2465174 
 0.3617286 -0.1817502 -0.2365868  0.5450749  0.3768349 
-0.5142120 -0.1878265  0.2268678  0.5486792  0.3833392 
 0.3441380  0.5139727  0.0620510  0.0356160  0.4765146 
-0.4994100  0.3525206 -0.2932449 -0.3835712  0.4968910 
Checking orthonormality
Matrix V^T*V: 
 1.0000000  0.0000000 -0.0000000 -0.0000000 -0.0000000 
 0.0000000  1.0000000  0.0000000 -0.0000000 -0.0000000 
-0.0000000  0.0000000  1.0000000  0.0000000 -0.0000000 
-0.0000000 -0.0000000  0.0000000  1.0000000  0.0000000 
-0.0000000 -0.0000000 -0.0000000  0.0000000  1.0000000 
Matrix U^T*V: 
 1.0000000  0.0000000 -0.0000000 -0.0000000 -0.0000000 
 0.0000000  1.0000000 -0.0000000  0.0000000  0.0000000 
-0.0000000 -0.0000000  1.0000000 -0.0000000 -0.0000000 
-0.0000000  0.0000000 -0.0000000  1.0000000 -0.0000000 
-0.0000000  0.0000000 -0.0000000 -0.0000000  1.0000000 


----- TESTING -----
Testing SVD algorithm on multiple random tall matrices, A, size NxM (N>M).
M in [50, 250] and N = M + rand([50, 250]) + 1
Testing is done by checking orthonormality of U and V, that D is diagonal, and U*D*V^T = A
Test 1 (size of A is 155 x 130) with 9 sweeps:  	Success (code 0)
Test 2 (size of A is 211 x 150) with 9 sweeps:  	Success (code 0)
Test 3 (size of A is 246 x 156) with 9 sweeps:  	Success (code 0)
Test 4 (size of A is 286 x 123) with 9 sweeps:  	Success (code 0)
Test 5 (size of A is 226 x 171) with 9 sweeps:  	Success (code 0)
Test 6 (size of A is 324 x 135) with 9 sweeps:  	Success (code 0)
Test 7 (size of A is 387 x 222) with 9 sweeps:  	Success (code 0)
Test 8 (size of A is 182 x 92) with 9 sweeps:  	    Success (code 0)
Test 9 (size of A is 88 x 66) with 8 sweeps:  	    Success (code 0)
Test 10 (size of A is 268 x 83) with 9 sweeps:  	Success (code 0)