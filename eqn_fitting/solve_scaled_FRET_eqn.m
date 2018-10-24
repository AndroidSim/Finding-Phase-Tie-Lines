function F=solve_scaled_FRET_eqn(x,aint,amax,fmax,Fap1,Fap2,xdonr,xaceptr,expconsta,expconstb)

% a=alpha;
ai=aint;
am=amax;
Fa=Fap1;
Fb=Fap2;
C1=expconsta;
C2=expconstb;
Xo=xdonr;
Xi=xaceptr;
% x(1) = Kpo = kpdonr;
% x(2) = Kpi = kpaceptr;

F=[ai.*Fa.*(1-exp(-C1.*Xi./(ai+(1-ai)./x(2))))./(1-exp(-C1.*Xi))./(ai+(1-ai)./x(1))+(1-ai).*Fb.*(1-exp(-C2.*Xi./(ai.*x(2)+1-ai)))./(1-exp(-C2.*Xi))./(ai.*x(1)+1-ai)-ai.*Fa-(1-ai).*Fb;
    (am.*Fa.*(1-exp(-C1.*Xi./(am+(1-am)./x(2))))./(1-exp(-C1.*Xi))./(am+(1-am)./x(1))+(1-am).*Fb.*(1-exp(-C2.*Xi./(am.*x(2)+1-am)))./(1-exp(-C2.*Xi))./(am.*x(1)+1-am)-am.*Fa-(1-am).*Fb)-fmax];