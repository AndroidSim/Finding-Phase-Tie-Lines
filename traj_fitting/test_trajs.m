function [trajs,alpha,beta]=test_trajs(bounds,startpt,ntpts)
% the function test_trajs finds the trajectories to test the fit of the
% FRET equation to determine tie lines.  The trajectories are the lines
% from a starting point on the boundary of the 2-phase region to every
% other point on the boundary

% find sizes of passed arrays
[npts,xy]=size(startpt);
[nbpts,bxy]=size(bounds);

% find the fraction of each phase along trajectory with every boundary point of 
% 2-phase region 

% check to make sure dims of arrays are right
if (xy ~= 2 | npts ~= 1)
	error('> 1 starting point or num cols not equal to 2(x and y)');
	break;
end
if (bxy ~= 2)
	error('num cols of boundary points not equal to 2 (x and y)');
	break;
end
if (~all(ismember(startpt,bounds,'rows')))
	error('start point does not lie on the boundary');
	break;
end

% find distances between every boundary point (symmetric matrix)
dbounds=zeros(nbpts,nbpts);
for i=1:nbpts
	for j=1:nbpts
		distbounds(i,j)=sqrt(((bounds(i,1)-bounds(j,1)).^2)+((bounds(i,2)-bounds(j,2)).^2));
	end
end

% find the linear trajectories that go from the startpt on the boundary to every other
% point on the boundary
[repsp,isp,ib]=intersect(startpt,bounds,'rows');
if (repsp ~= startpt)
	error('startpt not a pt on the boundary');
	break;
end
trajs=zeros(nbpts-1,ntpts,2);
ntrajs=nbpts-1; 
repbs=bounds;
repbs(ib,:)=[];
for i=1:ntrajs
	trajs(i,:,1)=linspace(bounds(ib,1),repbs(i,1),ntpts);
    trajs(i,:,2)=linspace(bounds(ib,2),repbs(i,2),ntpts);
end

% find the two distances between each trajectory point and the end points of the 
% trajectory.  the start point is the common end point of all trajectories
disttpts=zeros(ntrajs,ntpts);
for i=1:ntrajs
    for j=1:ntpts
	    disttpts(i,j)=sqrt(((trajs(i,j,1)-bounds(ib,1)).^2)+((trajs(i,j,2)-bounds(ib,2)).^2));
    end
end
	
% find fraction of phase: alpha and beta. amount of alpha defined as fraction from startpt
alpha=zeros(ntrajs,ntpts);
beta=zeros(ntrajs,ntpts);
repdb=distbounds;
repdb(ib,:)=[];
for i=1:ntrajs
    beta(i,1:ntpts)=disttpts(i,1:ntpts)./repdb(i,ib);
end
alpha=1-beta;

% find distances between each trajectory point and every boundary point
%dalpha=zeros(ntpts,nbpts);
%for i=1:ntpts
%	for j=1:nbpts
%		disttraj(i,j)=sqrt(((traj(i,1)-bounds(j,1))^2)+((traj(i,2)-bounds(j,2))^2));
%	end
%end

% i expect only 2 endpoints (i.e. places where the trajectory intersects phase boundaries)
%[endpts,it,ib]=intersect(traj,bounds,'rows');
%endpts=zeros(2);
%n=1;
%for i=1:ntpts
%	for j=nbpts
%		if (all(traj(i,:)==bounds(j,:)))
%			endpts(n,:)=bounds(traj(i,:)==bounds(j,:));
%			n=n+1;
%		end
%	end
%end
%if (~any(ismember(traj,bounds,'rows')))
%	endpts=interp1(x,Y,xi)
%end	
	