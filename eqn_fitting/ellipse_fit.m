function y=ellipse_fit(x,xdata)
% the function ellipse contains the general equation of an ellipse, as
% y=F(x), in cartesian coordinates.  this function is designed to be used
% to fit the boundaries of 2-phase regions (specifically L0+Ld) that appear
% elliptical.

% x=[a b c1 c2 theta]
%rot_matrix=[cos(x(5)) sin(x(5));-sin(x(5)) cos(x(5))];
part1=((x(2).^2)./(x(1).^2));
part2=((xdata-x(3)).^2);
y=x(4)+sqrt((x(2).^2)-(part1.*part2));
%y=nonrot_y*rot_matrix;