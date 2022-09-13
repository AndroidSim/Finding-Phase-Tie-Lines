function angles = slope_profile(dxO,dxA,a)

iu = [0:0.01:1]';
[Olength,u] = chord_length(dxO);
[Alength,v] = chord_length(dxA);
bC = interp1(u,dxO,iu);
aC = interp1(v,dxA,tieline_config_fxn(a,iu));
angles = point2angle(bC,aC).*(180/pi);

return