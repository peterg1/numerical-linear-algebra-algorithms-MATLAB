function [approximation,residual] = GMRES(A,b,n)
%{
GMRES gives an approximate solution to Ax = b
 
A is a square matrix

b is the b in Ax = b

n denotes to what iteration we carry GMRES out to. n should be less than or
equal to the size of A.

approximation is an approximation to x in AX = b.

A more sophisticated version of this algprithm is avialable in this respository
under GMRES_Givens.m. It uses Givens rotations to give an apporixmation and residual
at each iteration with little overhead (maybe infact fastor depending on whether or not MATLAB's built in QR factorization 
makes use of the HEssenberg structure of H) and allows for the user to stop the algorithm when the residual
is under an acceptable threshold as apposed to running for n steps. The algorithm presented
here though requires much less code.

This algorithm uses a function called Arnoldi, which is avaialbe in this repository, as a subroutine.

As for most iterative methods, the computation time is domiated by some
matrix vector product. As iterative algorithms are usually applied to large
sparse matrices, the comptuation time of the matrix vector product depends on
how the matrix structure is utilized in computation.

The convergence of GMRES depends heavily on specrum of A and the algorithm
will preform poorly for an arbitrary matrix A. The algorithm will converge
extremely fast for a matrix like A = 2*eye(200) + .5*randn(200)/sqrt(200)
(this matrix will have its eigenvalues evenely distributed around a disk in
 the complex plane of radius 1/2 centers at 2.
radius 1/2 

It can be shown that for the GMRES algorithm, the residue at step n is given by Pn(A)b, where
Pn is the unique polynomial of degree n with 0 as its constant coefficent that minimizes ||Pn(A)b|| in the 2-norm.
%}
[Q,H] = Arnoldi(A,b,n); 
Q = Q(:,1:n);
e = eye(n+1,1);
if n == length(A)
    e = eye(n,1);
end

c = norm(b)*e;
[Y,T] = qr(H,0);
X = linsolve(T,Y'*c);   
residual = norm(H*X-c);
X = Q*X;
approximation = X;
