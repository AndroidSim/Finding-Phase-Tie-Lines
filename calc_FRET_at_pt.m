function FRET = calc_FRET_at_pt(exp_surface,alpha,xt,Kp_fxn_type,gamma3,gamma4,alpha_pt,beta_pt,Xd,Xa,Ca,Cb)
% calc_FRET_at_pt calculates the model FRET at a point (data point) given at variables alpha and xt

% get Kps
[Kpd,Kpa]=Kp_fxns(xt,Kp_fxn_type,gamma3,gamma4);
% get fluorescence values at alpha phase and beta phase boundary points
Fax=[alpha_pt(1); beta_pt(1)];
Fay=[alpha_pt(2); beta_pt(2)];
Fa=interp2(exp_surface',Fax,Fay);
% evaluate model FRET equation  
FRET=feval(@FRET_eqn,alpha,Fa(1),Fa(2),Xd,Xa,Kpd,Kpa,Ca,Cb);