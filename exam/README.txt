Examination project 13

Author: Thomas Toft Lindkvist
AUID 643642 - Student number 201905635
35 mod 22 = 13

Two sided Jacobi algorithm for SVD

TODO:
- Testing on large matrices
- Short report in LaTeX
- Timing and graphs
- Comparison with GSL
- Perhaps some use-case?

------------ DESCRIPTION ------------
As explored in this class, the elementary Jacobi eigenvalue algorithm only works on symmetric matrices.
But suppose we wanted to find the eigenvalues of a non-symmetric matrix, or the singular values of a non-square matrix, could we generalize the known algorithm?
It turns out the answer is yes! The two sided Jacobi algorithm is such a generalization.
Any tall matrix A can, via the implemented algorithm, be decomposed into A=U*D*VT, where U and V are orthogonal matrices and D is a diagonal matrix containing the singular values of A.

A short report with the mathematics involved is found here (see report.pdf (TODO))

Convergence: The criterion for convergence is when two successive sweeps produce the same diagonalized matrix (within machine precision). 
