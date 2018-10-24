function F=exp_fit(x,xdata)

F=x(1).*exp(-x(2).*xdata+x(3));