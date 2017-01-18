%caclulate r_x and r_y in Am=r_x and A_m=r_y
%If the coordinate is x then it returns r_x 
%while if coordinate is y it returns r_y.

function [ r ] = find_r( h, n, coord)

for i=[2:1:n-1] 
    
    k1= 6/(h(i)+h(i-1));
    k2=(coord(i+1)-coord(i))/ h(i);
    k3=(coord(i) - coord(i-1)) / h(i-1);
    
    r(i,1)= k1 * (k2-k3);

end

end
