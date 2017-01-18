%TRACKING PROBLEM 
%{
This function does not return anything but plots the solution to the 
tracking problem as specified in Section 2.3 of the project. 
This routine produces three plots s_x vs t , s_y  vs t and s_y  vs s_x .

%}
function [  ] = solve_tracking( n,s_x, s_y, time)

figure 

%2D PLOT
for i=[1:1:n-1 ]
    s=ezplot(s_x(i),s_y(i), [time(i),time(i+1)]);
    hold on
end
title('S^y (t) vs S^x (t)')
xlabel('S^x (t)')
ylabel('S^y (t)')
hold off


end