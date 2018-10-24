function [dF]=derv_FRET_eqn(alpha,Fap1,Fap2,xdonr,xaceptr,kpdonr,kpaceptr,expconsta,expconstb,order)

a=alpha;
Fa=Fap1;
Fb=Fap2;
C1=expconsta;
C2=expconstb;
Xo=xdonr;
Xi=xaceptr;
Kpo=kpdonr;
Kpi=kpaceptr;

if (isequal(order,1))
    dF=Fa.*(1-exp(-C1.*Xi./(a+(1-a)./Kpi)))./(1-exp(-C1.*Xi))./(a+(1-a)./Kpo)-a.*Fa.*C1.*Xi./(a+(1-a)./Kpi).^2.*(1-1./Kpi).*exp(-C1.*Xi./(a+(1-a)./Kpi))./(1-exp(-C1.*Xi))./(a+(1-a)./Kpo)-a.*Fa.*(1-exp(-C1.*Xi./(a+(1-a)./Kpi)))./(1-exp(-C1.*Xi))./(a+(1-a)./Kpo).^2.*(1-1./Kpo)-Fb.*(1-exp(-C2.*Xi./(a.*Kpi+1-a)))./(1-exp(-C2.*Xi))./(a.*Kpo+1-a)-(1-a).*Fb.*C2.*Xi./(a.*Kpi+1-a).^2.*(Kpi-1).*exp(-C2.*Xi./(a.*Kpi+1-a))./(1-exp(-C2.*Xi))./(a.*Kpo+1-a)-(1-a).*Fb.*(1-exp(-C2.*Xi./(a.*Kpi+1-a)))./(1-exp(-C2.*Xi))./(a.*Kpo+1-a).^2.*(Kpo-1);
end
if (isequal(order,2))
    dF=-2.*Fa.*C1.*Xi./(a+(1-a)./Kpi).^2.*(1-1./Kpi).*exp(-C1.*Xi./(a+(1-a)./Kpi))./(1-exp(-C1.*Xi))./(a+(1-a)./Kpo)-2.*Fa.*(1-exp(-C1.*Xi./(a+(1-a)./Kpi)))./(1-exp(-C1.*Xi))./(a+(1-a)./Kpo).^2.*(1-1./Kpo)+2.*a.*Fa.*C1.*Xi./(a+(1-a)./Kpi).^3.*(1-1./Kpi).^2.*exp(-C1.*Xi./(a+(1-a)./Kpi))./(1-exp(-C1.*Xi))./(a+(1-a)./Kpo)-a.*Fa.*C1.^2.*Xi.^2./(a+(1-a)./Kpi).^4.*(1-1./Kpi).^2.*exp(-C1.*Xi./(a+(1-a)./Kpi))./(1-exp(-C1.*Xi))./(a+(1-a)./Kpo)+2.*a.*Fa.*C1.*Xi./(a+(1-a)./Kpi).^2.*(1-1./Kpi).*exp(-C1.*Xi./(a+(1-a)./Kpi))./(1-exp(-C1.*Xi))./(a+(1-a)./Kpo).^2.*(1-1./Kpo)+2.*a.*Fa.*(1-exp(-C1.*Xi./(a+(1-a)./Kpi)))./(1-exp(-C1.*Xi))./(a+(1-a)./Kpo).^3.*(1-1./Kpo).^2+2.*Fb.*C2.*Xi./(a.*Kpi+1-a).^2.*(Kpi-1).*exp(-C2.*Xi./(a.*Kpi+1-a))./(1-exp(-C2.*Xi))./(a.*Kpo+1-a)+2.*Fb.*(1-exp(-C2.*Xi./(a.*Kpi+1-a)))./(1-exp(-C2.*Xi))./(a.*Kpo+1-a).^2.*(Kpo-1)+2.*(1-a).*Fb.*C2.*Xi./(a.*Kpi+1-a).^3.*(Kpi-1).^2.*exp(-C2.*Xi./(a.*Kpi+1-a))./(1-exp(-C2.*Xi))./(a.*Kpo+1-a)-(1-a).*Fb.*C2.^2.*Xi.^2./(a.*Kpi+1-a).^4.*(Kpi-1).^2.*exp(-C2.*Xi./(a.*Kpi+1-a))./(1-exp(-C2.*Xi))./(a.*Kpo+1-a)+2.*(1-a).*Fb.*C2.*Xi./(a.*Kpi+1-a).^2.*(Kpi-1).*exp(-C2.*Xi./(a.*Kpi+1-a))./(1-exp(-C2.*Xi))./(a.*Kpo+1-a).^2.*(Kpo-1)+2.*(1-a).*Fb.*(1-exp(-C2.*Xi./(a.*Kpi+1-a)))./(1-exp(-C2.*Xi))./(a.*Kpo+1-a).^3.*(Kpo-1).^2;
end