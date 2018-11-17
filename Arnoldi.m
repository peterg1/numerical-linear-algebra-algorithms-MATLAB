function [Q,H] = Arnoldi(A,b,n)
%{
A is a square matrix

b is the starting vector (for eigenvalue problems b should be a random
vector, for solving systems of equations it will be specially chosen). 

n is the number of iterations to be preformed. Note that n must be less
than or equal to the number of columns of A.

H is the upper n+1 by n upper-left section of a Hessenberg form of A. The characteristic
polynomial, p, of H will be the unique monic polynomail such that
||p(A)b|| is minimized, where ||.|| is the 2-norm. The eigenvlaues of H
will approximate the larger eigenvlaues of A. The larger the n the better
the approximation and the more eigenvalues it will accuratly apporximate.

The columns of Q is an orthonormal basis for the Krylov subspace A, Ab, (A^2)b,...
A^(n-1)b.

As for most iterative methods, the computation time is domiated by some
matrix vector product. As iterative algorithms are usually applied to large
sparse matrices, the comptuation time of the matrix vector product depends on
how the matrix structure is utilized in computation.

Convergence: 
The convergence rate of the eigenvlaues of H produced by the Arnoldi iteration is not fully
understood and is apart of current research in numerical linear algebra (see
Treffethen and Bau chapter 34). 
%}
[~,m] = size(A);
Q = zeros(m); 
H = zeros(n); 
Q(1:m,1) = b/norm(b);
for i = 1:n
    v = A*Q(1:m,i);
    for j = 1:i
        H(j,i) = Q(1:m,j)'*v;
        v = v - H(j,i)*Q(1:m,j);
    end
    if norm(v) == 0 %if v is 0 then the Krylov subspace has now become an invarient
                     %subspace under A at this step. In this case each
                     %eigenvalue of H will be an eigenvalue of A and the
                     %solution x to Ax = b will lie in the generated Krylov
                     %subspace.
        return;
    end
    H(i+1,i) = norm(v);
    Q(1:m,i+1) = v/H(i+1,i);
end
H = H(1:n,1:n);
