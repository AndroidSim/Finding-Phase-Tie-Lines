function calc_surface = draw_surface(exp_surface,bd_pts,gammas,bd_rep,Kp_fxn_type,Xd,Xa,Ca,Cb,varargin)
% draw_surface draws the best fit surface from tie line determinization routine

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
    else
        continue
    end
    matrix=zeros(length(y_within{i}),2);
    xv=linspace(x(i),x(i),length(y_within{i}));
    matrix(:,1)=xv';
    matrix(:,2)=y_within{i}';
    xdata=cat(1,xdata,matrix);
    clear matrix
end
xdata(1,:)=[];

if (~isempty(varargin))
    inc=varargin{1};
    for i=1:inc
        xdata(1:2:end,:)=[];
    end
    calc_FRET=feval(@FRET_theory_calc,gammas,xdata,exp_surface,bd_pts,bd_rep,Kp_fxn_type,Xd,Xa,Ca,Cb);
    calc_surface=calc_FRET;
    return
end
calc_FRET=feval(@FRET_theory_calc,gammas,xdata,exp_surface,bd_pts,bd_rep,Kp_fxn_type,Xd,Xa,Ca,Cb);
calc_surface=zeros(size(exp_surface));
for i=1:size(xdata,1) % length(xdata)
    calc_surface(xdata(i,1),xdata(i,2))=calc_FRET(i);
end
calc_surface(calc_surface == 0)=sqrt(NaN);
calc_surface(isnan(calc_surface))=exp_surface(isnan(calc_surface));
figure
plot_surface(calc_surface);