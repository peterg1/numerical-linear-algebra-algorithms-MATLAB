function [L,U,P] = Gauss(A) 
%{
Gaussian elimination with partial pivoting. The input A is an m by n matrix. If k = min(m,n) then the output is
is an m by m permutation matrix P,
an m by c lower triangular matrix L, and a c by n upper triangular matrix U such that PA = LU

FLOP count - To the leading order, the flop count is the same as gaussian elimination without pivoting (finding
all of the pivots takes ~ 1/2(m^2) comparisons). The FLOP count is dominated by preforming the elimination step 
after pivoting:

for k = 1: z
...
    for j = k+1:m  
        L(j,k) = U(j,k)/U(k,k);
        U(j,k:m) = U(j,k:m) - L(j,k)*U(k,k:m);

the line U(j,k:m) = U(j,k:m) - L(j,k)*U(k,k:m) involves m-k+1
multiplications and m-k+1 subtractions for a total of 2 flops per entry, so the flop count from here is: 
sum from i = 1 to z (sum from j = i+1 to m)(2*(n-i)) = sum from i = 1 to z (2*(m-i)(n-i)) = 2*[zmn - m(z^2-z)/2 -
n(z^2-z)/2 + (2z^3+3z^2+z)/6 ~ 2*[zmn - (1/2)mz^2 - (1/2)nz^2 + (1/3)z^3] where z ~ min(m,n). 
For the case where m = n this reduces to (2/3)m^3 FLOPS.
%}
[m,n] = size(A);
L = eye(m);
P = eye(m);
U = A;
if m <= n
    z = m-1;
else
    z = n;
end

for k = 1:z  
    x = U(k,k);
    pivot_row = k;
    for i = k+1:m % find the largest element in column k with row # > k
        if abs(x) < abs(U(i,k))
            x = U(i,k);
            pivot_row = i;
        end
    end     
    
    %swaps pivot_row with row k in U
    u = U(pivot_row,k:n);             
    U(pivot_row,k:n) = U(k,k:n);       
    U(k,k:n) = u;     
    
    %swaps pivot_row with row k in P
    q = P(pivot_row,:);             
    P(pivot_row,:) = P(k,:);        
    P(k,:) = q;
    
    %swaps pivot_row with k in L
    l = L(pivot_row,1:k-1);         
    L(pivot_row,1:k-1) = L(k,1:k-1);
    L(k,1:k-1) = l;           
    
    %preforms Gaussian elimination after pivoting
    for j = k+1:m  
        L(j,k) = U(j,k)/U(k,k);
        U(j,k:n) = U(j,k:n) - L(j,k)*U(k,k:n);
    end
end
U = U(1:min(m,n),:);
L = L(:,1:min(m,n));
end