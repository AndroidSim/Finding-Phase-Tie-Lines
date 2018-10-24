function [nquad]=whichquad(vector)
% the function whichquad determines which quadrant a vector lies in, returning 1 for
% quadrant 1, 2 for quadrant 2, 3 for quadrant 3, and 4 for quadrant 4

if (vector(1) > 0 & vector(2) > 0)
    nquad=1;
elseif (vector(1) < 0 & vector(2) > 0)
    nquad=2;
elseif (vector(1) < 0 & vector(2) < 0)
    nquad=3;
else % vector(1) > 0 & vector(2) < 0
    nquad=4;
end