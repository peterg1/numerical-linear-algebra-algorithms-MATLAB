function [H,R] = HR_to_Hessenberg(A) 
%{
Input is a square matrix A, and the output is a Hessenberg matrix H, tridiagonal if A is hermitian, and a matrix R used to
store the householder reflectors (1st column holds the 1st reflector, 2nd column hold holds the 2nd reflector...).
This algorithm applies householder reflectors to reduce the matrix A to
Hessenberg form via a unitary similarity transformation, i.e. A=QHQ*, where
if v_i is the ith householde reflector, then Q_i = I -
2(v_iv_i*)/(v_i*v_i), then Q = Q_1Q_2...Q_(m-2).

FLOP count - The flop count is dominated by applying the reflector to the
left and right of H:

for k = 1:m-2;
...
    H(k+1:m,k:m) = H(k+1:m,k:m) - 2*R(k+1:m,k)*(R(k+1:m,k)'*H(k+1:m,k:m)); 
    H(1:m,k+1:m) = H(1:m,k+1:m) - 2*(H(1:m,k+1:m)*R(k+1:m,k))*R(k+1:m,k)';

Both lines preform the same number of FLOPS so we will just count the FLOPS
of one line and double. R(k+1:m,k)'*H(k+1:m,k:m) takes ~ 2(m-k)^2 FLOPS,
then the outter prouct with this resulting row vector and 2*R(k+1:m,k)
takes ~ (m-k)^2 multiplications, then finally subtracting this result from H(k+1:m,k:m) takes ~ (m-k)^2 subtractions.
this gives us a total of ~ 5(m-k)^2 FLOPS. So the FLOP count is 2*(sum from
k = 1 to m-2 of 5(m-k)^2 ~ 10m^3.
%}
H = A; 
[m,n] = size(A);
R = zeros(m); 
for k = 1:m-2;
    x = H(k+1:m,k);
    e = eye(m-k,1);
    R(k+1:m,k) = sign(x(1))*norm(x)*e + x; % stores kth reflector
    R(k+1:m,k) = R(k+1:m,k)/norm(R(k+1:m,k)); % normalizes kth reflector
    H(k+1:m,k:m) = H(k+1:m,k:m) - (2*R(k+1:m,k))*(R(k+1:m,k)'*H(k+1:m,k:m)); %uses reflector to kill elements
    H(1:m,k+1:m) = H(1:m,k+1:m) - (2*(H(1:m,k+1:m)*R(k+1:m,k)))*R(k+1:m,k)'; %multiplies by conjugate transpose(it is hermetian)
end
end