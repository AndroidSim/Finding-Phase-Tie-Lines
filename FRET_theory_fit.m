function [gammas,lsq,varargout] = FRET_theory_fit(exp_surface,bd_pts,bd_rep,Kp_fxn_type,Xd,Xa,Ca,Cb,init_pt,varargin)

% FRET_theory_fit is the fxn that fits the experimental FRET surface data (exp_surface) within the phase boundary
% given in bd_pts with the theoretical FRET model developed by Jeff Buboltz
% varargout = data_surface,calc_surface

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
ydata=zeros(length(xdata),1);
for i=1:length(xdata)
    ydata(i)=exp_surface(xdata(i,1),xdata(i,2));
end

% evoke minimization routine
options=optimset('MaxFunEvals',1000,'TolFun',1e-10,'TolX',1e-10,'MaxIter',1000,'LargeScale','on','LevenbergMarquardt','on','Display','iter','Jacobian','off');
if (length(init_pt) == 2 & length(varargin) == 3)
    if (varargin{1} == 1)
        x0=init_pt;
        lb=[0.1 0.1]';
        ub=[10 10]';
        gammas = lsqcurvefit(@FRET_theory_calc,x0,xdata,ydata,lb,ub,options,exp_surface,bd_pts,bd_rep,Kp_fxn_type,Xd,Xa,Ca,Cb,varargin{1},varargin{2},varargin{3});
        g_params=[varargin{2};varargin{3};gammas];
        ycalc=feval(@FRET_theory_calc,g_params,xdata,exp_surface,bd_pts,bd_rep,Kp_fxn_type,Xd,Xa,Ca,Cb);
        % display output
        nt=30;
        draw_tie_lines(bd_pts,bd_rep,nt,g_params(1),g_params(2));
        calc_surface=draw_surface(exp_surface,bd_pts,g_params,bd_rep,Kp_fxn_type,Xd,Xa,Ca,Cb);
    elseif (varargin{1} == 2)
        x0=init_pt;
        lb=[0 0]';
        ub=[1 1]';
        gammas = lsqcurvefit(@FRET_theory_calc,x0,xdata,ydata,lb,ub,options,exp_surface,bd_pts,bd_rep,Kp_fxn_type,Xd,Xa,Ca,Cb,varargin{1},varargin{2},varargin{3});
        g_params=[gammas;varargin{2};varargin{3}];
        ycalc=feval(@FRET_theory_calc,g_params,xdata,exp_surface,bd_pts,bd_rep,Kp_fxn_type,Xd,Xa,Ca,Cb);
        % display output
        nt=30;
        draw_tie_lines(bd_pts,bd_rep,nt,g_params(1),g_params(2));
        calc_surface=draw_surface(exp_surface,bd_pts,g_params,bd_rep,Kp_fxn_type,Xd,Xa,Ca,Cb);
    end
else
    x0=init_pt';
    gammas=x0;
    lb=[0 0 0.1 0.1]';
    ub=[1 1 10 10]';
    gammas = lsqcurvefit(@new_FRET_theory_calc,x0,xdata,ydata,lb,ub,options,exp_surface,bd_pts,bd_rep,Kp_fxn_type,Xd,Xa,Ca,Cb); 
    ycalc=feval(@new_FRET_theory_calc,gammas,xdata,exp_surface,bd_pts,bd_rep,Kp_fxn_type,Xd,Xa,Ca,Cb);
    % display output
    %nt=30;
    draw_tie_lines(bd_pts,bd_rep,nt,gammas(1),gammas(2));
    calc_surface=draw_surface(exp_surface,bd_pts,gammas,bd_rep,Kp_fxn_type,Xd,Xa,Ca,Cb);
end 
lsq=sum((ycalc-ydata).^2);

if (nargout == 3 | nargout == 4)
    varargout{1}=calc_surface;
end
if (nargout == 4)
    data_surface=exp_surface;
    [nx,ny]=size(data_surface);
    index=zeros([nx ny]);
    for i=1:size(xin,1)
        index(xin(i,1),xin(i,2))=1;
    end
    index=logical(index);
    data_surface(~index)=sqrt(NaN); % isnan(calc_surface)
    varargout{2}=data_surface;
end