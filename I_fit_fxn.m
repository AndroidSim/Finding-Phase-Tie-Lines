function v = fp_fit_fxn(x,xdata,kp)
% xdata = [1 fa@2_phase 0]
% fpalpha = (kp.*fa)./(fb+kp.*fa);
% fpbeta = fb./(fb+kp.*fa);

v = x(1)*((kp.*xdata)./((1-xdata)+kp.*xdata)) + x(2)*((1-xdata)./((1-xdata)+kp.*xdata));

return