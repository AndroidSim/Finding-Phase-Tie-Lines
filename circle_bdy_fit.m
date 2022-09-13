function [best_p,best_chisq] = circle_bdy_fit(xdata,ydata)

%lb = [0;0;0;0;.1];
%ub = [2*pi;2*pi;2*pi;2*pi;10];
%lb = [pi;0;pi/2;3*pi/2;.1];
%ub = [3*pi/2;pi/2;pi;2*pi;10];
%options = optimset('display','iter','maxfunevals',2000);
%[p,chisq] = lsqcurvefit(@circle_bdy_fit_fxn,x0,xdata,ydata,lb,ub,options);
best_p = zeros(5,5);
best_chisq = zeros(5,1);
nfilled = 1;
traja = [0:pi/6:2*pi];
trajb = [0:pi/6:2*pi];
cptc = [0:pi/6:2*pi];
cptd = [0:pi/6:2*pi];
kp = [0.4:0.1:0.9 1:1:4];
for a = 1:length(traja)
    for b = 1:length(trajb)
        for c = 1:length(cptc)
            for d = 1:length(cptd)
                for p = 1:length(kp)
                    if traja(a) == trajb(b) | (any(traja(a) == [0 2*pi]) & any(trajb(b) == [0 2*pi]))
                        continue;
                    end
                    
                    if cptc(c) == cptd(d) | (any(cptc(c) == [0 2*pi]) & any(cptd(d) == [0 2*pi]))
                        continue;
                    end
                    
                    ytheor = feval(@circle_bdy_fit_fxn,[traja(a);trajb(b);cptc(c);cptd(d);kp(p)],xdata);
                    chisq = norm(ytheor-ydata);
                    
                    % update stack of best fit parameters
                    if (isempty(find(best_chisq == 0)))
                        [best_chisq,index] = sort(best_chisq);
                        best_p = best_p(index,:);
                        i = find(chisq <= best_chisq);
                        if (~isempty(i))
                            i = max(i);
                            best_chisq(i) = chisq;
                            best_p(i,:) = [traja(a) trajb(b) cptc(c) cptd(d) kp(p)];
                        end
                    else
                        best_chisq(nfilled) = chisq;
                        best_p(nfilled,:) = [traja(a) trajb(b) cptc(c) cptd(d) kp(p)];
                        nfilled = nfilled+1;
                    end
                end
            end
        end
    end
end

return