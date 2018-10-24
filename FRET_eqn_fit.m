function FRET=FRET_eqn_fit(x,alpha)
% the equation developed by Jeff Buboltz that describes the shape of
% a FRET curve along a thermodynamic tie line.
% x=[kpdonr kpaceptr Fap1 Fap2];

beta=1-alpha;
xdonr=0.001;
xaceptr=0.00025;

part1a=alpha;
part1b=x(3);
part1c=((1-exp(-xaceptr./(alpha+(beta./x(2)))))./(1-exp(-xaceptr)));
part1d=(xdonr./(alpha+(beta./x(1)))./xdonr);
part2a=(-beta);
part2b=x(4);
part2c=((1-exp(-xaceptr./((alpha.*x(2))+beta)))./(1-exp(-xaceptr)));
part2d=(xdonr./((alpha.*x(1))+beta)./xdonr);

FRET=(part1a.*part1b.*part1c.*part1d)+(part2a.*part2b.*part2c.*part2d);
