function [W, R] = householder_factorization(A)
%{
QR factorization using householder reflectors (backwards stable unlike classic gram schmidt qr factorization)
Input is a matrix A, output is W which stores the vectors used for the householder reflectors used
(1st column = first vector, 2nd column = 2nd vector...) and R is upper triangular.
If v_1 = W(:1), v_2 = W(:2),... and Q_1 is the first reflector, i.e. Q_1 = I - 2(v_1v_1*)/(v_1*v_1), 
Q_2 is the second reflector, i.e, I - 2(v_2v_2*)/(v_2*v_2), ..., (the v_i's are normalized
so the term v_i*v_i is actually not needed), then Q* = Q_nQ_(n-1)...Q_2Q_1, 
Q = Q_1Q_2...Q_n (note that reflectors are hermitian and orthogonal), and A = QR. Note the formulas
Q* = Q_nQ_(n-1)...Q_2Q_1 and  Q = Q_1Q_2...Q_n can often be used in place of Q and Q* and so Q and
Q* should not be formed directly. 

FLOP count - the FLOP count is dominated by applying the householder reflector inside the for loop:
A(k:m,k:n) = A(k:m,k:n) - 2*W(k:m, k)*W(k:m, k)'*A(k:m,k:n), which for a householder vector of length g = m-k+1 does, 
per column, 2g-1 flops for the dot product as well as g scalar multiplications and g subtractions for approximatly
4g FLOPS, so the FLOP count is approximatly 4 flops per entery. This gives us a total of approximatly: If m > n,
(sum from k = 1 to n)(4*(m-k)(n-k)) ~ 4(mn^2 - (mn^2+mn)/2 - (n^3 + n^2)/2 + (2n^3 + 3n^2 + n)/6 which after 
dropping the lower order terms gives us 2mn^2 - 2/3n^3 FLOPS. If n >= m then this count would become 
4/3m^3.
%}

[m,n] = size(A); 
if n >= m
    stop = m-1;
else
    stop = n;
end
W = zeros(m,stop);
for k = 1:stop
    e = eye(m-k+1,1); 
    x = A(k:m,k); 
    W(k:m, k) = sign(x(1))*norm(x)*e + x; % stores the vector used to make the householder reflector in W
    W(k:m, k) = W(k:m, k)/norm(W(k:m, k));
    A(k:m,k:n) = A(k:m,k:n) - 2*W(k:m, k)*W(k:m, k)'*A(k:m,k:n); % applying the reflector to induce 0's below the diagonal in column k
end
R = A(1:m,1:n);
end


