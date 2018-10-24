function FRET=FRET_eqn(alpha,Fap1,Fap2,xdonr,xaceptr,kpdonr,kpaceptr)
% the equation developed by Jeff Buboltz that describes the shape of
% a FRET curve along a thermodynamic tie line.

beta=1-alpha;
part1a=alpha;
part1b=Fap1;
part1c=((1-exp(-xaceptr./(alpha+(beta./kpaceptr))))./(1-exp(-xaceptr)));
part1d=((xdonr./(alpha+(beta./kpdonr)))./xdonr);
part2a=(alpha-1);
part2b=Fap2;
part2c=((1-exp(-xaceptr./((alpha.*kpaceptr)+beta)))./(1-exp(-xaceptr)));
part2d=((xdonr./((alpha.*kpdonr)+beta))./xdonr);

FRET=(part1a.*part1b.*part1c.*part1d)+(part2a.*part2b.*part2c.*part2d);
