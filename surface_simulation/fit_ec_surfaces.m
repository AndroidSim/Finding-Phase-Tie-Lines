function [metric]=fit_ec_surfaces(expr_surface,calc_surface)

% fit_ec_surfaces fits the calculated surface calc_surface to the experimental surface exp_surface
% exp_surface = [x_coords y_coords fret_value], calc_surface = [x_coords y_coords fret_value]

% data_surface=FdiI_surface;
% data_surface(isnan(p2_best_fit_surface))=sqrt(NaN);
% [ei,ej]=find(~isnan(data_surface));
% [ci,cj]=find(~isnan(p2_best_fit_surface));
% expr_surface=[ei ej data_surface(~isnan(data_surface))];
% calc_surface=[ci cj p2_best_fit_surface(~isnan(p2_best_fit_surface))];

%Error_Vector=calc_surface-exp_surface;
%sse(i,j)=sum((Error_Vector.^2)./exp_surface); % length(Error_Vector)

% the following code implements a shape fitting algorithm given in gunsel+tekalp

% create distance matrix D and scaled distance matrix H for both experimental and calculated surface  
nespts=length(expr_surface);
ncspts=length(calc_surface);
if (~isequal(nespts,ncspts))
    error('num experimental surface pts not equal to num calculated surface pts');
else
    nspts=nespts;
end
eD=zeros(nspts,nspts);
cD=zeros(nspts,nspts);
eH=zeros(nspts,nspts);
cH=zeros(nspts,nspts);

for i=1:nspts
    for j=1:nspts
        eD(i,j)=norm(expr_surface(i,:)-expr_surface(j,:));
        cD(i,j)=norm(calc_surface(i,:)-calc_surface(j,:));
    end
end
e_sigma=sum(sum(eD))/nspts; % both expr and calc sigmas = mean distance between surface pts
c_sigma=sum(sum(cD))/nspts;

eH=exp(1/2*(eD/e_sigma).^2);
cH=exp(1/2*(cD/c_sigma).^2);

% solve eigenvalue and eigenvalue equation of matrix H for both experimental and calculated surfaces
[eV,eE]=eig(eH);
[cV,cE]=eig(cH);

% compute mismatch matrix Z 
for i=1:nspts
    for j=1:nspts
        Z(i,j)=norm(eV(i,:)-cV(j,:));
    end
end

% find matching feature points (num matching pts = K)
Kmax=nspts;
K=0;
n=1;
for i=1:nspts
    [minvi,minj]=min(Z(i,:));
    [minvj,mini]=min(Z(:,minj));
    if (isequal(i,mini) & isequal(minvi,minvj))
        K=K+1;
        simv(n)=minvi; % or minvj
        n=n+1;
    end
end

% compute similarity metric
if (isequal(K,Kmax))
    metric=(norm(simv).^2)/K; % or just norm(simv) ?
elseif (K < Kmax)
    metric=((norm(simv).^2)+2*(Kmax-K))/Kmax; % or just norm(simv)+2(Kmax-K) ?
end