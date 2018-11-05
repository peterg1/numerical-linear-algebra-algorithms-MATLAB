function [D] = QR_shift(A,tolerance)
%{
This algorithm utilizes two functions called find_submatrix and
HR_to_Hessenberg, which can be found in the same repository as this
algorithm.

A is a symmetric matrix

tolerance is the threshold that subdiagonal elements need to be under before they are set to 0. 

D is backwards stable approximation to the Shur form (diagonalization of A since A is symmetric).

Convergence rate depends on the relative spacing of the eigenvalues, which is why an eigenvalue approxiamtion shift is 
used (to bring one of the eigenvalues close to 0). Full discussion of the convergence rate is too complex to provide here. 
See Trefethen and Bau section 28/29.

Apprixmating the real Shur form for the non-symemtric case while using shifts and minimizing computation becomes
more nuanced, see QR algorithm with implicit double shift.

%}



%This algoithm can be simplified if we assume deflation will only occur at the m,m-1 and m-1, m spot, as will typically occur 
%when the Rayleigh quotient or wilkinson shift is chosen

[D,~] = HR_to_Hessenberg(A); %A is first reduced to Hessenberg form to reduce computation.
a = 1;
c = 1;
while a ~= 0; % a == 0 iff A is diagonalized so we keep doing QR algorith until this occurs
    if c == 1; %this section should only occur if there were off diagonal elements replaced by 0
        [a,b] = find_submatrix(D); %This function finds the first block of a hermetian tridiagonal matrix A of maximal 
                                   %size that has non-zero entries on its subdiagonal. This step is know as deflation. When 
                                   %the Rayleigh quotient shift is used it will usually 
                                   %only result in the removal of the last row and column. This isn't assumed here but if it 
                                   %is a shorter algorithm can be given.
        if a == 0;
            return
        end
        B = D(a:b,a:b); 
         
    end
    [m,n] = size(B);
    u = B(m,m); %This is the prescribed Rayleigh quotient shift. A more sophisticated shift may be calculated here instead.
    [Q,R] = qr(B-u*eye(m)); %gets QR factorization of B-u*I with built in qr algorithm (strategies for qr factorization
                            %of symmetric tridiagonal matrices can be used to greatly speed up the decomposition
                            %which qr algorithms in this repository do not use, hopefully MATLAB's built in qr algorithm does)
    B = R*Q + u*eye(m); % sets B to the next matrix to be used in this algorithm
    c = 0;
    for i = 1:n-1; %checks if any of the off diagonal elements of B are under the tolerance level
        if abs(B(i,i+1)) < tolerance;
            B(i,i+1) = 0;
            B(i+1,i) = 0;
            c = 1; %lets us know if any off diagonal element was set to 0 so we can deflate
        end
    end
    D(a:b,a:b) = B; 
end
