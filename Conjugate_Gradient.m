function [approximation, residual] = CG(A,b,n) 
%{ 
Conjugate Gradient Iteration -

approximates Ax = b for n iterations, A is a symmetric positive definite matrix
alternitively, n may be replaces with some convergence condition

ith column of approximation is the approxiamtion of x at ith step

ith entry of residual is the residual at ith step.

As for most iterative methods, the computation time is domiated by some
matrix vector product. As iterative algorithms are usually applied to large
sparse matrices, the comptuation time of the matrix vector product depends on
how the matrix structure is utilized in computation. Here the matrix-vector product is A*p in the line r = r - a*A*p;
A is an m by m matric and p is a legnth m vector, so A*p is 2m^2 FLOPS if matrix stucture is not utilized.

convergence - discussion of convergence is similar to that given for GMRES.
also, see Trefethen and Bau section 38.
%}
x = 0;
r = b;
p = b;
[~,q] = size(A);
approximation = zeros(q,n);
residual = zeros(n,1);
for i = 1:n;
    rr = r; %need to store r_(i-1) for later use
    a = (r'*r)/(p'*A*p);
    x = x + a*p;
    approximation(1:q,i) = x;
    r = r - a*A*p;
    residual(i) = norm(r);
    B = (r'*r)/(rr'*rr);
    p = r + B*p;
end
