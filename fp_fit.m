function [kp,chisq,fpc] = fp_fit(x,y,aorb)
% x = fraction of phase alpha or fraction of phase beta
% y = fraction of probe in alpha or fraction of probe in beta
% aorb = 'a' for alpha or 'b' for beta
%
% searching for kp to best fit fraction of probe data
% for grid search: kp = [0.1:0.02:0.9 1:0.2:10]';
% fb = linspace(0,1);
% [fpa,fpb] = fraction_probe(kp,fb);
% t = linspace(0,1,length(xdata))';
% kp0 = 1;
% kp1 = 1;
% options = optimset('display','iter');
% [x,resnorm] = lsqcurvefit(@fp_fit_fxn,[kp0 kp1]',x,y,[kp(1) kp(1)]',[kp(end) kp(end)]',options);

options = optimset('display','iter');
[kp,chisq] = lsqcurvefit(@fp_fit_fxn,1,x,y,0.1,10,options,aorb);
fpc = feval(@fp_fit_fxn,kp,x,aorb);