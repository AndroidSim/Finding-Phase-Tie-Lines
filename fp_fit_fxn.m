function v = fp_fit_fxn(x,xdata,aorb)
% x = kp;
% xdata = fa
% fpa = (fa./(fa+(kp*fb)));
% fpb = ((kp*fb)./(fa+(kp*fb)));

if aorb == 'a'
    v = (xdata./(xdata+(x*(1-xdata))));
end

if aorb == 'b'
    v = ((x*xdata)./((1-xdata)+(x*xdata)));
end

return