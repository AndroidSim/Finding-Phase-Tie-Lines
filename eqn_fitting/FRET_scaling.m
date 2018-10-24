function [FRET_scaled]=FRET_scaling(FRET_data,alpha,Faa,Fab,xdonr,xaceptr,kpdonr,kpaceptr,expconsta,expconstb)

beta=1-alpha;

if (Faa > Fab)
    part1a=alpha;
    part1b=Faa-Fab;
    part1c=((1-exp(-(expconsta).*xaceptr./(alpha+(beta./kpaceptr))))./(1-exp(-(expconsta).*xaceptr)));
    part1d=((xdonr./(alpha+(beta./kpdonr)))./xdonr);
    FRET_scaled=FRET_data-(part1a.*part1b.*part1c.*part1d)
else
    part2a=beta;
    part2b=Fab-Faa;
    part2c=((1-exp(-(expconstb).*xaceptr./((alpha.*kpaceptr)+beta)))./(1-exp(-(expconstb).*xaceptr)));
    part2d=((xdonr./((alpha.*kpdonr)+beta))./xdonr);
    FRET_scaled=FRET_data-(part2a.*part2b.*part2c.*part2d);
end