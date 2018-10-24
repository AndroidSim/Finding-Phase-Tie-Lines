function fit=lsqfit_scaled_FRET_eqn(x,data,alpha,Fap1,Fap2,xdonr,xaceptr,expconsta,expconstb)

a=alpha;
Fa=Fap1;
Fb=Fap2;
C1=expconsta;
C2=expconstb;
Xo=xdonr;
Xi=xaceptr;
% x(1) = Kpo = kpdonr;
% x(2) = Kpi = kpaceptr;

F=a.*Fa.*(1-exp(-C1.*Xi./(a+(1-a)./x(2))))./(1-exp(-C1.*Xi))./(a+(1-a)./x(1))+(1-a).*Fb.*(1-exp(-C2.*Xi./(a.*x(2)+1-a)))./(1-exp(-C2.*Xi))./(a.*x(1)+1-a)-a.*Fa-(1-a).*Fb;
                    