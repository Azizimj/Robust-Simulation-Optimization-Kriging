function [f] = simulation(x,e)

x_low=[-5 -5];% x lower bound of decision variable size=1*d
x_high=[5 5]; % x upper bound of decision variable size=1*d
e_high=[5 5]; e_low=[-5 -5];
rang_x=x_high-x_low; mid_x=(x_high+x_low)/2;
rang_e=e_high-e_low; mid_e= (e_high+e_low)/2;

%converting of experiment points
x=x.*(rang_x/2)+mid_x;
e=e.*(rang_e/2)+mid_e;

f= 5*(x(1)^2+x(2)^2)^2 + (e(1)^2+e(2)^2)^2 + x(1)*(-e(1)+e(2)+5) + x(2)*(e(1)-e(2)+3);


end