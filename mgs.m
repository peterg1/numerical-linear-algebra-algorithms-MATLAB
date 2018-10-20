function [Q,R] = mgs(X)
%{ 
QR factorization via modified gram schmidt(backwards stable, unlike classic gram schmidt). Input X
is a matrix, output Q is a matrix with orthonormal columns, R is upper triangular, X = QR.

FLOP count - the number of FLOPS is dominated by the inner for loop. For an m by n matrix the
dot product dot(Q(:,i),Q(:,j)) contains m multiplications and m-1 additions, and Q(:,j) - R(i,j)*Q(:,i) 
contains m multiplications and m subtractions for approximatly 4m FLOPS. The FLOP count is the approximatly 
(sum from i = 1 to n)[(sum j = i+1 to n) 4m] ~ (sum from i = 1 to n) i*4m ~ 2mn^2 FLOPS.
%}

[m,n] = size(X); 
Q = X; 
R = zeros (n,n);
for i = 1:n
    R(i,i) = norm(Q(:,i)); 
    Q(:,i) = Q(:,i)/R(i,i); 
    for j = i+1:n
        R(i,j) = dot(Q(:,i),Q(:,j)); 
        Q(:,j) = Q(:,j) - R(i,j)*Q(:,i);
    end
end
 
end
