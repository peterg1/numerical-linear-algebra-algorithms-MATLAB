function [approximation, residue] = BCG(A,b,n) 
%{
Biconjugate gradient iteration for approximating x in Ax = b for n steps of
the iteration. Alternativley, n may be replaced with some convergence
condition.

A is a squre matrix. The algorithm will converge poorly for an arbitrary
square matrix though, a good precondiyioner must be used.
Disscussion of convergence is similar to that of the GMRES 
and congugate gradient iterations. A good preconditioner must be used.

the ith column of approxiamtion stores the approxiamtion of x at the ith step 

the ith entry of residue stores the residue of the ith approximation.

As for most iterative methods, the computation time is domiated by some
matrix vector product. Here it is computing r = r - a*A*p; and s = s - a'*A*q;

The algorithm will preform poorly on a random matrix without
preconditioners applied to it. See trefethen and Bau section 39 for
convergence.
%}
[e,~] = size(A);
q = rand(e,1); 
x = zeros(e,1);
p = b;
r = b;
s = q;
approximation = zeros(e,n);
residue = zeros(n,1);
for i = 1:n
    a = (s'*r)/(q'*A*p);
    x = x + a*p;
    approximation(1:e,i) = x;
    rr = r; % r_(i-1) will be needed later
    r = r - a*A*p;
    residue(i) = norm(r);  
    ss = s; % s_(i-1) will be needed later
    s = s - a'*A*q;
    B = (s'*r)/(ss'*rr);
    p = r + B*p;
    q = s + conj(B)*q;
end

