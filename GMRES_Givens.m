function [approximation,residual] = GMRES(A,b,n)
%{

A is a square matrix

b is the b in Ax = b

n denotes to what iteration we carry GMRES out to. n should be less than or
equal to the size of A.

approximation is an approximation to x in AX = b.

GMRES gives approximate solution to Ax = b up to n steps. Givens rotations
are used to approxiamte the answer at each iteration. Since Givens rotations are used to update the QR
factorization at each step given the QR factorization at the previous step, each iteration takes O(n) time. 
Alterantively (how it should be used),
the algorithm can be modified to run until the residual is below a given
threshold (for this case a good precondition must be chosen for the matrix or else the algorithm will
most likely run until n = length(A)).

The convergence of GMRES depends heavily on specrum of A and the algorithm
will preform poorly for an arbitrary matrix A, and a good preconditioner must be used 
to get a good convergence rate. The algorithm will converge
extremely fast for a matrix like A = 2*eye(200) + .5*randn(200)/sqrt(200)
- this matrix will have its eigenvalues evenly distributed around a disk in
 the complex plane of radius 1/2 centered at 2. 

%}
[~,m] = size(A);
Q = zeros(m,n); %will store the orthogonal bases for the Krylov subspace
H = zeros(m,n); % will store the Hessenberg matrix generatex by the arnoldi itteration
residual = zeros(n,1); 
approximation = zeros(m,n);
Q(1:m,1) = b/norm(b);
givens = zeros(2*n+2,2); %stores the 2 by 2 givens rotations, 2i-1:2i,1:2 contains ith givens rotation.
for i = 1:2;  %Due Arnoldi iteration for first 2 steps
    v = A*Q(1:m,i);
    for j = 1:i;
        H(j,i) = Q(1:m,j)'*v;
        v = v - H(j,i)*Q(1:m,j);
    end
    H(i+1,i) = norm(v);
    Q(1:m,i+1) = v/H(i+1,i);
end
Qn = zeros(n+1,n); %This is the Q in the QR factorization.
R = zeros(n); %R in QR. QR factors will be updated each time we go through the Arnoldi Iteration


[G,~] = planerot(H(1:2,1)); %MATLAB's built in Givens/plane rotation
givens(1:2,:) = G;
Qn(1:2,1:2) = G';

R(1:2,1:2) = G*H(1:2,1:2);
[G,y] = planerot([R(2,2),H(3,2)]');
givens(3:4,1:2) = G;
R(2:3,2) = y;
Qn(1:3,1:3) = [[Qn(1,1),Qn(2,1),0]',[Qn(1,2),Qn(2,2),0]',[0,0,1]']*[[1,0,0]',[0,G(1,1:2)]',[0,G(2,1:2)]'];
P = eye(n+1,n);
P(1:3,1:3) = Qn(1:3,1:3);
for i = 3:n;
    v = A*Q(1:m,i); 
    for j = 1:i            %do next step of Arnoldi Iteration
        H(j,i) = Q(1:m,j)'*v;
        v = v - H(j,i)*Q(1:m,j);
    end
    if norm(v) == 0 %if the norm of v is 0 then the Krylov subspace generated by b has now become an invarient
                     %subspace under A at this step. This implies x can be
                     %solved for explicitly at this point and so our last
                     %estimate for x would be our best.
        return;
    end
    H(i+1,i) = norm(v);
    R(1:i+1,i) = H(1:i+1,i);
    for j = 1:i-1 %update R with givens roation. Takes O(i)/O(m) time (linear)
        R(j:j+1,i) = givens((2*j)-1:2*j,1:2)*R(j:j+1,i);
    end
    [G,y] = planerot(R(i:i+1,i));
    givens(2*i-1:2*i,1:2) = G;
    R(i:i+1,i) = y;
    for j = 1:i  %update Qn with givens roation. Takes O(i)/O(m) time (linear)
        Qn(j,i:i+1) = [Qn(j,i)*G(1,1),Qn(j,i)*G(2,1)]';
    end
    Qn(i+1,i:i+1) = [G(1,2),G(2,2)];
    e = eye(i+1,1);
    c = norm(b)*e;
    X = linsolve(R(1:i,1:i),Qn(1:i+1,1:i)'*c);  %gets the appoximation
    X = Q(1:m,1:i)*X;
    approximation(:,i) = X;
    residual(i) = norm(A*X - b); 
    Q(1:m,i+1) = v/H(i+1,i);
end
end

