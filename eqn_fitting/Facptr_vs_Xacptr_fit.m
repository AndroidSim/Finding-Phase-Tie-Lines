function Facptr=Facptr_vs_Xacptr_fit(x,mf)
% the exponential factor in Jeff Buboltz's FRET equation desribes how the fluorescense of the acceptor
% and thus fret varies with mole fraction of the acceptor.  This is the equation of that factor.

Facptr=x(1)*(1-exp(-x(2)*mf));