function [best_gammas,data_surface] = param_grid_search(exp_surface,bd_pts,bd_rep,Kp_fxn_type,Xd,Xa,Ca,Cb,nsave)

% param_grid_search performs a search on a grid covering parameter space.  the parameters belong to the
% FRET model used to fit an experimental FRET surface

% get x,y coords of data points within phase boundary
minx=min(bd_pts(:,1));
maxx=max(bd_pts(:,1));
miny=min(bd_pts(:,2));
maxy=max(bd_pts(:,2));
x=[round(minx):round(maxx)];
y=[round(miny):round(maxy)];
y_within=cell(length(x),1);
coords=cell(length(x),1);
b_pts=bd_pts;
b_pts(end,:)=[];
xdata=zeros(1,2);
for i=1:length(x)
    xx=linspace(x(i),x(i),length(y));
    [xi,yi]=polyxpoly(xx,y,b_pts(:,1),b_pts(:,2),'unique');
    if (length(yi) == 2)
        y_within{i}=find(yi(1) < [1:size(exp_surface,2)] & [1:size(exp_surface,2)] < yi(2));
    elseif (length(yi) == 1)
        continue
        %y_within{i}=yi;
    else
        continue
        %error('either 0 or > 2 intersection points')
    end
    matrix=zeros(length(y_within{i}),2);
    xv=linspace(x(i),x(i),length(y_within{i}));
    matrix(:,1)=xv';
    matrix(:,2)=y_within{i}';
    xdata=cat(1,xdata,matrix);
    clear matrix
end
xdata(1,:)=[];
xin=xdata;
%xdata(1:2:end,:)=[];
%xdata(1:2:end,:)=[];
%xdata(1:2:end,:)=[];
%xdata(1:2:end,:)=[];
%keyboard
%xdata=xdata(round(end/2),:);
%ydata=exp_surface(xdata(1),xdata(2));
%keyboard
ydata=zeros(length(xdata),1);
for i=1:length(xdata)
    ydata(i)=exp_surface(xdata(i,1),xdata(i,2));
end

% perform grid search
% gamma1 0->1, gamma2 0->1, gamma3 0.1->10, gamma4 0.1->10
gammas=zeros(1,4);
best_gammas=zeros(nsave,4);
best_fit_values=zeros(nsave,1);
nfilled=1;
for g1=0:0.1:1
    for g2=g1:0.1:1
        for g3=0.1:0.1:10
            for g4=0.1:0.1:10
                % calculated FRET at each xdata according to model
                %gammas=[g1 g2 g3 g4];
                gammas=[0.2 0.8 1 1];
                calc_FRET=FRET_theory_calc(gammas,xdata,exp_surface,bd_pts,bd_rep,Kp_fxn_type,Xd,Xa,Ca,Cb);
                % least squares fit to experimental FRET
                fit=sum((calc_FRET-ydata).^2);
                
                % update stack of best fit parameters
                if (isempty(find(0 == best_fit_values)))
                    [best_fit_values,index]=sort(best_fit_values);
                    best_gammas=best_gammas(index,:);
                    i=find(fit <= best_fit_values);
                    if (~isempty(i))
                        i=max(i);
                        best_fit_values(i)=fit;
                        best_gammas(i,:)=gammas;
                    end
                else
                    best_fit_values(nfilled)=fit;
                    best_gammas(nfilled,:)=gammas;
                    nfilled=nfilled+1;
                end
            end
        end
    end
end

[best_fit_values,index]=sort(best_fit_values);
best_gammas=best_gammas(index,:);

nt=30;
draw_tie_lines(bd_pts,bd_rep,nt,gammas(1),gammas(2))
calc_surface=draw_surface(exp_surface,xdata,calc_FRET);

data_surface=exp_surface;
[nx,ny]=size(data_surface);
index=zeros([nx ny]);
for i=1:size(xin,1)
    index(xin(i,1),xin(i,2))=1;
end
index=logical(index);
data_surface(~index)=sqrt(NaN); % isnan(calc_surface)