%{
n=number of measurements
h=array of time differences 
s_x=array of splines in x
s_y=array of splines in y
time=array of time data points 

Finds velocities by using It calculates the velocity using two methods: 
one by differentiating the x-splines and y-splines 
and second by using the numerical way of three point  midpoint method. 

%}
function [  ] = solve_velocity( n,h,s_x, s_y,x,y, time)

%First CALCULATE VELOCITY USING SPLINES
syms t
u_x=sym([]);
u_y=sym([]);

u_xNum=[];
u_yNum=[];
s_xNum=[];
s_yNum=[];

for i=[1:1:n-1]
    
    %FIND EXPRESSION FOR DERIVATIVE FOR EACH SPLINE
    u1=diff(s_x(i),t);
    u_x(i)=u1;
    
    u2=diff(s_y(i),t);
    u_y(i)=u2;
    
    
end


%CALCULATE FOR EACH TIME t the x,y coord and the corresponding velocity
%In this loop we dont add the nth point. Use s(n) to calculate it
for i=[1:1:n-1]
    
    
    t=time(i);
    
    s_xNum(i)=subs(s_x(i));
    s_yNum(i)=subs(s_y(i));
    
    u_xNum(i)=subs(u_x(i));
    u_yNum(i)=subs(u_y(i));
    
    
end
%FIX THE nth ENTRY
t=time(n);
s_xNum(n)=subs(s_x(n-1));
s_yNum(n)=subs(s_y(n-1));

u_xNum(n)=subs(u_x(n-1));
u_yNum(n)=subs(u_y(n-1));


figure 
%2D PLOT to plot the dashed line
for i=[1:1:n-1 ]
    s=ezplot(s_x(i),s_y(i), [time(i),time(i+1)]);
    set(s,'linestyle','--');
    hold on
end
title({'Velocity Field Using Cubic Splines';'Dotted appears the actual curve'})
xlabel('x(t)')
ylabel('y(t)')

%plot vecotr field 
quiver(s_xNum,s_yNum,u_xNum,u_yNum)

%-----------APPROXIMATE USING THREE POINT MID POINT METHOD-----------
v_x=[];
v_y=[];



%Treat endpoints first
%Use three point endpoint formula for x0,y0 and xn,yn
c=1/(2*h(1));
v_x(1)=c*(-3*x(1)+4*x(2)-x(3));
v_y(1)=c*(-3*y(1)+4*y(2)-y(3));

c=-1/(2*h(1)); %Now we are using -h for xn,yn
v_x(n)=c*(-3*x(n)+4*x(n-1)-x(n-2));
v_y(n)=c*(-3*y(n)+4*y(n-1)-y(n-2));

%Fill in intermidiate values v1....vn-1
for i=[2:1:n-1]
    v_x(i)=(1/(2*h(i)))* ( x(i+1)-x(i-1));
    v_y(i)=(1/(2*h(i)))* ( y(i+1)-y(i-1));
end



figure 
%2D PLOT
for i=[1:1:n-1 ]
    l=ezplot(s_x(i),s_y(i), [time(i),time(i+1)]);
    set(l,'linestyle','--');
    hold on
end
title({'Velocity Field Using Three Point Mid Point Method';'Dotted appears the actual curve'})
xlabel('x(t)')
ylabel('y(t)')
%produce velocity field
quiver(s_xNum,s_yNum,v_x,v_y);


% Last plot: Superimpose the two fields to compare them. The original curve
% is not included
figure 
quiver(s_xNum,s_yNum,u_xNum,u_yNum);
hold on
quiver(s_xNum,s_yNum,v_x,v_y);
title('The two Vector Fields superimposed');
legend('Spline Differentiation','Three point Midpoint Method')

%{
%UNCOMMENT TO FIND AVERAGE DIFFERENCE OF THE TWO VELOCITIES
u_xT=0;
v_xT=0;
for i=[1:1:n ]
    u_xT=u_xT+u_xNum(i);
    v_xT=v_xT+v_x(i);
end
 
Average_ux=u_xT/(n);
Average_vx=v_xT/(n);
difference=abs(Average_ux-Average_vx)
%}
end 