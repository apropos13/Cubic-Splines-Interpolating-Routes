%{
The script after loading the data.mat file it then prompts the user to 
enter whether he wants to use Natural Boundary Condition (by typing "N") 
or the Clamped Boundary Condition by typing "C").

It then sets up the matrix A, solves the system  to find the value of the 
vector m and then calculates the coefficients. The splines are then 
represented as symbolic expressions and are stored in the arrays s_x and s_y. 

%}

prompt = 'Type N for "Natural Boundary" or C for "Clamped Bundary" ';
response = input(prompt,'s'); %Store selected method by the user 

N=0;
if response=='N' 
    N=1;
    display('Natural Boundary Selected');
    
elseif response=='C'
    C=1;
    display('Clamped Boundary Selected')
    
else %if user types someting else it throws exception 
    error('Error. You must type either N or C (case sensitive)')
end


%import the data 
filename='data.mat'; 
data=importdata(filename);

time = data(:,1);
x= data (:,2);
y=data (:,3);
h=[];
r_x=[];
r_y=[];
n=size(x,1); %find size of first dimension (i.e. 32x1-->32)

%calculate h(i)
for j=[1:1:n-1] %indexing in matlab starts at 1
    
    h(j)=time(j+1)-time(j);

end

%FINDS r in Am=r
r_x=find_r(h, n, x); 
r_y=find_r(h, n, y);

%setup A
A=zeros(n);
for i=[2:1:n-1]
    
    A(i,i+1)=h(i)/ (h(i-1)+h(i)); %miu
    A(i-1,i-1)=2; %diagonal
    A(i,i-1)= 1- A(i,i+1); %lamda
    
    
end

%Fill in the last and second to last value of diagonal
A(n,n)=2;
A(n-1,n-1)=2;

%Natural Boundary 
if N==1
    A(1,2)=0; %miu_0
    A(n,n-1)=0; %lamda_n-1
    
    r_x(1)=0;
    r_x(n)=0;
    
    r_y(1)=0;
    r_y(n)=0;
else %Clamped Boundary
    A(1,2)=1; %miu_0
    A(n,n-1)=1; %lamda_0
    
    
    %use 3point endpoint formula for r_x(1) with h>0
   
    fx0_prime=(1/(2*h(1)))*(-3*x(1)+4*x(2)-x(3));
    r_x(1)=6/(h(1)*h(1)) * (x(2)-x(1)) - 6/h(1) * fx0_prime;
    %use 3point endpoint formula for r_x(n) with h<0
    fxn_prime=(1/(-2*h(1)))*(-3*x(n)+4*x(n-1)-x(n-2));
    r_x(n)=( (6*fxn_prime)/h(n-1)) - (6*(x(n)-x(n-1))/(h(n-1)*h(n-1))); 
    
  
    %use 3point endpoint formula for r_y(1) with h>0
    fy0_prime=(1/(2*h(1)))*(-3*y(1)+4*y(2)-y(3));
    r_y(1)=6/(h(1)*h(1)) * (y(2)-y(1)) - 6/h(1) * fy0_prime;
    %use 3point endpoint formula for r_y(n) with h<0
    fyn_prime=(1/(-2*h(1)))*(-3*y(n)+4*y(n-1)-y(n-2));
    r_y(n)=( (6*fyn_prime)/h(n-1)) - (6*(y(n)-y(n-1))/(h(n-1)*h(n-1))); 
    
end

%Solves Am=r for  x-data and then y-data
m_x=solve_tridiagonal(n, A, r_x);

m_y=solve_tridiagonal(n, A, r_y);


%Up to here we have m0....mn
%Use this to find bi,ci,di.

b_x=[];
c_x=[];
d_x=[];

b_y=[];
c_y=[];
d_y=[];
   
%Coefficients for x-splines
for i=[1:1:n-1]
    b_x(i)=  ( ( x(i+1)- x(i) ) / h(i) ) - ( (h(i)*(m_x(i+1)+2*m_x(i)))/6);
    c_x(i)=m_x(i)/2;
    d_x(i)=( m_x(i+1)-m_x(i) ) / (6*h(i)) ;
    
    
end

%Coefficient for y-splines
for i=[1:1:n-1]
    b_y(i)=( ( y(i+1)- y(i) ) / h(i)) - (h(i)*(m_y(i+1)+2*m_y(i)))/6;
    c_y(i)=m_y(i)/2;
    d_y(i)=( m_y(i+1)-m_y(i) ) / (6*h(i));
end


%define symbolic variable t
syms t

%now define s_x,s_y to store the splines
s_x=sym([]);
s_y=sym([]);

for i=[1:1:n-1]
    s=x(i)+b_x(i)*(t-time(i))+c_x(i)*(t-time(i)).^2+ d_x(i)*(t-time(i)).^3;
    s_x(i)=s;
    
    p=y(i)+b_y(i)*(t-time(i))+c_y(i)*(t-time(i)).^2+ d_y(i)*(t-time(i)).^3;
    s_y(i)=p;
    
end

%Just a sanity check 
if size(s_x,2)~=n-1
    error('System error')
end

if size(s_y,2)~=n-1
    error('System error')
end


%Plot s_x(i) vs x
for i=[1:1:n-1] 
    ezplot(s_x(i), [time(i),time(i+1)]);
    xlim([0 7])
    ylim([-0.8,0.8])
    hold on
end
title('S^x vs t')
xlabel('x(t)')
ylabel('S^x (t)')
hold off

figure
%Plot s_y(i) vs y
for i=[1:1:n-1]
    ezplot(s_y(i), [time(i),time(i+1)]);
    xlim([0 7])
    ylim([-1,1.5])
    hold on
end
title('S^y vs t')
xlabel('y (t)')
ylabel('S^y (t)')
hold off

%Finds parametric equations, produces plots for tracking problem
solve_tracking(n,s_x,s_y,time); 

%Finds velocity, produces quiver plots
solve_velocity(n,h,s_x,s_y,x,y,time);

