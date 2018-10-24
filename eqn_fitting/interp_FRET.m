function FRET=interp_FRET(a,alpha,data)
% the function interp_FRET interpolates the FRET value at a from the FRET data in data at alpha
FRET=interp1(alpha,data,a);