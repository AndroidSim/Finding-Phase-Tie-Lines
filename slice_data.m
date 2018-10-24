function [linedata,xline,yline]=slice_data(surface,p1,p2,nps)
%slice_data is a function to pick out a linear trajectory of points from 
%a surface of FRET data

%surface is a matrix with row indices as x conc (Xplsat) and 
%column indices as y conc (Xchol)
[nx,ny]=size(surface);
%change coord system to nx-1 and ny-1
nx=nx-1;
ny=ny-1;
p1=p1-1;
p2=p2-1;
%find the eqn of the line going through points p1 and p2
%line of data will contain nps points
m=(p2(2)-p1(2))/(p2(1)-p1(1));
m
%point-slope form for eqn of a line: y=m*(x-p1(1))+p1(2);
%find x intercept
xint=(0+(m*p1(1))-p1(2))/m;
%find y intercept
yint=m*(0-p1(1))+p1(2);
xint
yint
%initialize enter and exit points
enter=0;
exit=0;
%find range x and y coords. of points on the line that enter and 
%exit surface
if (isnan(m) | isinf(m))
	enter(1)=p1(1);
	enter(2)=0;
	exit(1)=p1(1);
	exit(2)=ny;
	enter
	exit
elseif (m==0)
	enter(1)=0;
	enter(2)=yint;
	exit(1)=nx;
	exit(2)=yint;
	enter
	exit
else
	if (xint <= nx & xint > 0 & yint <= ny & yint > 0)
    		enter(1)=0;
    		enter(2)=yint;
		exit(1)=xint;
		exit(2)=0;
		enter
		exit
	else
    		if (xint >= nx & yint >= ny)
        		exit(1)=nx;
        		exit(2)=m*(exit(1)-p1(1))+p1(2);
			enter(2)=ny;
			enter(1)=(enter(1)+m*p1(1)-p1(2))/m;
			enter
			exit
    		end
		if (xint < 0 & yint < ny)
			enter(1)=0;
			enter(2)=yint;
			point1(1)=nx;
			point1(2)=m*(point1(1)-p1(1))+p1(2);
			point2(2)=ny;
			point2(1)=(point2(2)+m*p1(1)-p1(2))/m;
			point1
			point2
			if (point1(2) < ny & point2(1) > nx)
				exit(1)=point1(1);
				exit(2)=point1(2);
			end
			if (point1(2) > ny & point2(1) < nx)
				exit(1)=point2(1);
				exit(2)=point2(2);
			end
			enter
			exit
		end
    		if (yint >= ny & xint <= nx)
        		enter(2)=ny;
        		enter(1)=(enter(2)+m*p1(1)-p1(2))/m;
			exit(1)=xint;
			exit(2)=0;
			enter
			exit
		end
		if (xint >= nx & yint <= ny)
			enter(2)=yint;
			enter(1)=0;
			exit(1)=nx;
			exit(2)=m*(exit(1)-p1(1))+p1(2);
			enter
			exit
		end
		if (yint < 0 & xint < nx)
			enter(1)=xint;
			enter(2)=0;
			point1(1)=nx;
			point1(2)=m*(point1(1)-p1(1))+p1(2);
			point2(2)=ny;
			point2(1)=(point2(2)+m*p1(1)-p1(2))/m;
			if (point1(2) < ny & point2(1) > nx)
				exit(1)=point1(1);
				exit(2)=point1(2);
			end
			if (point1(2) > ny & point2(1) < nx)
				exit(1)=point2(1);
				exit(2)=point2(2);
			end
			enter
			exit
		end
		if (xint == 0 & yint == 0)
			enter(1)=0;
			enter(2)=0;
			point1(1)=nx;
			point1(2)=m*(point1(1)-p1(1))+p1(2);
			point2(2)=ny;
			point2(1)=(point2(2)+m*p1(1)-p1(2))/m;
			if (point1(2) <= ny & point2(1) >= nx)
				exit(1)=point1(1);
				exit(2)=point1(2);
			end
			if (point1(2) >= ny & point2(1) <= nx)
				exit(1)=point2(1);
				exit(2)=point2(2);
			end
			enter
			exit
		end
	end
end
%enter
%exit
%find the coordinates of the points on surface surrounding each point
%on the line given by the enter and exit points
%change back to original coordinates
nx=nx+1;
ny=ny+1;
enter=enter+1;
exit=exit+1;
x=[1:nx];
y=[1:ny];
m=(exit(2)-enter(2))/(exit(1)-enter(1));
if (isnan(m) | isinf(m))
	xline=enter(1);
	yline=linspace(enter(2),exit(2),nps);	
elseif (m==0)
	xline=linspace(enter(1),exit(1),nps);
	yline=enter(2);
else
	xline=linspace(enter(1),exit(1),nps);
%	yline=m*(xline-enter(1))+enter(2);
	yline=linspace(enter(2),exit(2),nps);
end
%ndata=nx*ny;
%resurface=zeros(ndata,3);
%resurface(:,1)=repmat(x',(ndata/nx),1);
%matrix=repmat(y,nx,1);
%matrixtwo=reshape(matrix,ndata,1);
%resurface(:,2)=matrixtwo(:);
%%for i=1:2:x
%%	matrix2=[matrix(:,i);matrix(:,i+1)];
%%	matrix2=cat(1,matrix(:,i),matrix(:,i+1));
%%end
%resurface(:,3)=reshape(surface,ndata,1);
%matrix=reshape(surface,ndata,1);
%resurface(:,3)=matrix(:);
%linedata=griddata(resurface(:,1),resurface(:,2),resurface(:,3),xline,yline);
linedata=interp2(surface',xline,yline);