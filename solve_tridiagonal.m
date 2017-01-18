function [ m ] = solve_tridiagonal(n, A, r)
%Crout Factorization for Tridiagonal Linear Systems 
%Solve Am=r, with dim(A)=n
%   A has to be tridiagonal

%Horizontically concatenate A and r
A=horzcat(A,r);

%Declare Helpers
L=[];
U=[];
Z=[];

%Step 1
L(1,1)=A(1,1);
U(1,2)=A(1,2)/L(1,1);
Z(1)=A(1,n+1)/L(1,1);

%Step 2
for i=[2:1:n-1]
L(i,i-1)=A(i,i-1);
L(i,i)=A(i,i)-L(i,i-1)*U(i-1,i);
U(i,i+1)=A(i,i+1)/L(i,i);
Z(i)=( A(i,n+1)-L(i,i-1)*Z(i-1))/L(i,i);

end 

%Step 3
L(n,n-1)=A(n,n-1);
L(n,n)=A(n,n)-L(n,n-1)*U(n-1,n);
Z(n)= ( A(n,n+1)- L(n,n-1)*Z(n-1))/L(n,n);

%display(L);

%display(U);

%display(Z);

%Step 4
m(n)=Z(n);

%Step 5
for i=[(n-1):-1:1]
    m(i)=Z(i)-U(i,i+1)*m(i+1);
end



end

