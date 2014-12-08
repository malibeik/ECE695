function x = root( x0,max_iter,epsilon,delta )
%ROOT Summary of this function goes here
%   Detailed explanation goes here

x = x0;
fx = sqrt(3)*lamqDP*cos(x)+sqrt(3)*lamdDP*sin(x)+2*idc*(LqDP-LdDP)*sin(2*x);

for n=1:max_iter
    
    fp = -sqrt(3)*lamqDP*sin(x)+sqrt(3)*lamdDP*cos(x)+4*idc*(LqDP-LdDP)*cos(2*x);
    if fp<delta
        disp('derivative too small (Newton Method)')
        return;
    end
    
    d = fx/fp;
    x = x - d;
    
    fx = sqrt(3)*lamqDP*cos(x)+sqrt(3)*lamdDP*sin(x)+2*idc*(LqDP-LdDP)*sin(2*x);
    if abs(d)<epsilon
        disp('convergence');
        return;
    end
    
end

