function [best_fit,fretcalcalfa,fretcalcbeta,fretdata]=find_tielines(surface,bounds,ntpts,Fap1,Fap2,xdonr,xaceptr,kpdonr,kpaceptr)
% the function find_tielines is an iterative method to find the thermodynamic
% tie lines of a 2-phase coexistence region given the data points of the trajectory, 
% and the boundaries of the coexistence region

% find sizes of passed arrays
[nx,ny]=size(surface);
[nbpts,bxy]=size(bounds);

for i=1:nbpts
% find test trajectories from a starting boundary point
    [trajs,alpha,beta]=test_trajs(bounds,bounds(i,:),ntpts);

% find min least squares error of fret data along each trajectory and fret calculated
% with fret equation
    [ntrajs,ntpts,txy]=size(trajs);
    for j=1:ntrajs
        fretdata=interp2(surface',trajs(j,:,1),trajs(j,:,2));
        fretcalcalfa=feval(@FRET_eqn,alpha(j,:),Fap1,Fap2,xdonr,xaceptr,kpdonr,kpaceptr);
%       fretcalcalfa=feval(@FRET_eqn,alpha(j,:),beta(j,:),Fap1,Fap2,xdonr,xaceptr,kpdonr,kpaceptr);
        fretcalcbeta=feval(@FRET_eqn,beta(j,:),Fap1,Fap2,xdonr,xaceptr,kpdonr,kpaceptr);
%       fretcalcalfa=feval(@FRET_eqn,beta(j,:),alpha(j,:),Fap1,Fap2,xdonr,xaceptr,kpdonr,kpaceptr);
        Error_Vector_alfa=fretcalcalfa-fretdata;
        Error_Vector_beta=fretcalcbeta-fretdata;
%       sse_alfa(j)=sum(Error_Vector_alfa.^2);
%       sse_beta(j)=sum(Error_Vector_beta.^2);
        sse_alfa=sum(Error_Vector_alfa.^2);
        sse_beta=sum(Error_Vector_beta.^2);
	if (sse_alfa <= sse_beta)
		least(j)=sse_alfa;
    	else
    		least(j)=sse_beta;
    	end
    end
    lowest=find(min(least)==least);
    best_fit(i,:,1)=trajs(lowest,:,1);
    best_fit(i,:,2)=trajs(lowest,:,2);
%   if (min(sse_alfa) <= min(sse_beta))
%	least=find(min(sse_alfa)==sse_alfa);
%   else
%  	least=find(min(sse_beta)==sse_beta);
%   end
%   best_fit(i,:,1)=trajs(least,:,1);
%   best_fit(i,:,2)=trajs(least,:,2);
end