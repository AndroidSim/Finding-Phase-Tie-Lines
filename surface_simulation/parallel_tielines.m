function [tielines,alpha_beta_pts]=parallel_tielines(boundary_pts,npts,varargin) %varagrin = degree_inc or 'slope',slope
% the function parallel_tielines generates many sets of parallel tie lines with
% different slopes (given by degree_inc) that traverse the 2-phase region with 
% the boundary given by the points in boundary_pts.

[nbpts,bxy]=size(boundary_pts);
if (bxy ~= 2) 
	error('num cols of boundary points not equal to 2 (x and y)')
end    
xmin=min(boundary_pts(:,1))-10;
xmax=max(boundary_pts(:,1))+10;
ymin=min(boundary_pts(:,2))-10;
ymax=max(boundary_pts(:,2))+10;
div=(ymax-ymin)/(xmax-xmin);

if (nargin == 4)
    if (~strcmp(varargin{1},'slope')) %varargin{1} ~= 'slope'
        error('slope passed but not recognized')
    else
        slope=varargin{2};
        if (~isnumeric(slope))
            error('slope passed is not a number')
        end
    end
    % generate set of parallel lines
    parlines=gen_par_lines(xmin,xmax,ymin,ymax,slope,div,npts);
    [npls,nppts,pxy]=size(parlines);
    % with 'unique' only 2 intersection points per parallel line
%    [x_intersect,y_intersect]=polyxpoly(parlines(:,:,1)',parlines(:,:,2)',boundary_pts(:,1),boundary_pts(:,2),'unique');
    iint=1;
    for i=1:npls
        [xint,yint]=polyxpoly(parlines(i,:,1),parlines(i,:,2),boundary_pts(:,1),boundary_pts(:,2),'unique');
        if (~isempty(xint) & ~isempty(yint))
            if (size(xint,1) == 2 & size(yint,1) == 2)
                x_intersect(:,iint)=xint;
                y_intersect(:,iint)=yint;
                alpha_pts(iint,1)=xint(1);
                alpha_pts(iint,2)=yint(1);
                beta_pts(iint,1)=xint(2);
                beta_pts(iint,2)=yint(2);
                iint=iint+1;
            end
        end
    end
    alpha_beta_pts=cat(3,alpha_pts,beta_pts);% in 3rd dim, 1 = alpha, 2 = beta;
    intersect_pts=cat(3,x_intersect,y_intersect);
 %   [pair,ntie,txy]=size(i_intersect_pts);
 %   intersect_pts(1:pair,1:ntie,1:txy,iset)=i_intersect_pts;
 %   iset=iset+1;
    clear x_intersect;
    clear y_intersect;
    
    [pair,ntie,txy]=size(intersect_pts);
    tylines=zeros(ntie,npts,txy);
    for j=1:ntie
        xline=linspace(intersect_pts(1,j,1),intersect_pts(2,j,1),npts);
        yline=linspace(intersect_pts(1,j,2),intersect_pts(2,j,2),npts);
        tylines(j,:,1)=xline;
        tylines(j,:,2)=yline;
    end
    tielines=tylines;
else 
    if (nargin ~= 3)
        error('nargin not equal to 3 if not passing a slope')
    end
    iset=1;
    degree_inc=varargin{1};
    %radian=degree_inc*(pi/180);
    radian_inc=degree_inc*(pi/180);
    %radian=radian_inc;
    for radian=0:radian_inc:pi
        slope=tan(radian);
        if (radian == pi/2)
            slope = Inf;
        end
        % generate set of parallel lines
        parlines=gen_par_lines(xmin,xmax,ymin,ymax,slope,div,npts);
        [npls,nppts,pxy]=size(parlines);
        % with 'unique' only 2 intersection points per parallel line
%       [x_intersect,y_intersect]=polyxpoly(parlines(:,:,1)',parlines(:,:,2)',boundary_pts(:,1),boundary_pts(:,2),'unique');
        iint=1;
        for i=1:npls
            [xint,yint]=polyxpoly(parlines(i,:,1),parlines(i,:,2),boundary_pts(:,1),boundary_pts(:,2),'unique');
            if (~isempty(xint) & ~isempty(yint))
                if (size(xint,1) == 2 & size(yint,1) == 2)
                    x_intersect(:,iint)=xint;
                    y_intersect(:,iint)=yint;
                    alpha_pts(iint,1)=xint(1);
                    alpha_pts(iint,2)=yint(1);
                    beta_pts(iint,1)=xint(2);
                    beta_pts(iint,2)=yint(2);
                    iint=iint+1;
                end
            end
        end
        alpha_beta_pts{iset}=cat(3,alpha_pts,beta_pts);% in 3rd dim, 1 = alpha, 2 = beta;
        i_intersect_pts=cat(3,x_intersect,y_intersect);
 %      [pair,ntie,txy]=size(i_intersect_pts);
        intersect_pts{iset}=i_intersect_pts;
 %      intersect_pts(1:pair,1:ntie,1:txy,iset)=i_intersect_pts;
        iset=iset+1;
        clear x_intersect;
        clear y_intersect;
    end
    
    nsets=size(intersect_pts,2);
    for i=1:nsets
        int_pts=intersect_pts{i};
        [pair,ntie,txy]=size(int_pts);
        tylines=zeros(ntie,npts,txy);
        for j=1:ntie
            xline=linspace(int_pts(1,j,1),int_pts(2,j,1),npts);
            yline=linspace(int_pts(1,j,2),int_pts(2,j,2),npts);
            tylines(j,:,1)=xline;
            tylines(j,:,2)=yline;
        end
        tielines{i}=tylines;
    end
end

