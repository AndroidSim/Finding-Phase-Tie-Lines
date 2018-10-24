function FRET=FRET_eqn_kp_fit(x,alpha,xdonr,xaceptr,expconsta,expconstb,Fap1,Fap2)
% the equation developed by Jeff Buboltz that describes the shape of
% a FRET curve along a thermodynamic tie line.
% x=[kpdonr kpaceptr];

beta=1-alpha;
%xdonr=0.001;
%xaceptr=0.00025;
%expconsta=238;
%expconstb=508;

part1a=alpha;
part1b=Fap1;
part1c=((1-exp(-(expconsta).*xaceptr./(alpha+(beta./x(2)))))./(1-exp(-(expconsta).*xaceptr)));
part1d=(xdonr./(alpha+(beta./x(1)))./xdonr);
part2a=beta;
part2b=Fap2;
part2c=((1-exp(-(expconstb).*xaceptr./((alpha.*x(2))+beta)))./(1-exp(-(expconstb).*xaceptr)));
part2d=(xdonr./((alpha.*x(1))+beta)./xdonr);

FRET=(part1a.*part1b.*part1c.*part1d)+(part2a.*part2b.*part2c.*part2d);