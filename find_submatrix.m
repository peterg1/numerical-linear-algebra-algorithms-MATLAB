function [a,b] = find_submatrix(A) 
%{
This function finds the first block of a hermetian tridiagonal matrix A of maximal size that has non-zero
entries on its subdiagonal. This is used as a subroutine for the QR
algorithms with shifts for hermitain matrices. It returns a,b, a>b which are the rows where the block 
starts and ends.
%}
[m,n] = size(A); 
a = 0;
b = m;
for k = 1:m-1;   
    if A(k,k+1) ~= 0; 
        a = k;
        break
    end
end
if a ==0;
    return;
end
if a == m-1;
    return;
else
    for j = a+1:m-1;
        if A(j,j+1) == 0;
            b = j;
            return;
        end
    end
end
end

