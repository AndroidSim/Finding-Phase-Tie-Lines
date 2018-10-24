function [yesorno]=isintersect(line_pts,boundary_pts)
% the function isintersect determines whether the line line_pts intersects the boundary boundary_pts
[nbpts,bxy]=size(boundary_pts);
if (bxy ~= 2)
	error('num cols of boundary points not equal to 2 (x and y)')
end
[nlpts,lxy]=size(boundary_pts);
if (lxy ~= 2)
	error('num cols of line points not equal to 2 (x and y)')
end

[xint,yint]=polyxpoly(line_pts(:,1),line_pts(:,2),boundary_pts(:,1),boundary_pts(:,2),'unique');
if (isempty(xint) | isempty(yint))
    yesorno=0;
elseif (size(xint,1) == 1 & size(yint,1) == 1)
    yesorno=0;
else
    yesorno=1;
end