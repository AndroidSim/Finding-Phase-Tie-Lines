function [bd_length,bd_intervals] = boundary(bd_pts,bd_rep)
% the boundary function (right now) calculates the perimeter of the phase boundary given 
% by the boundary points in bd_pts.  bd_rep is how the boundary is represented, linear or cubic spline

[nbpts,bxy]=size(bd_pts);
if (~isequal(bxy,2))
	error('number of columns of boundary point matrix is not equal to 2 (Xs and Xc coords)')
end

switch bd_rep
    case 'linear'
        bd_length=sum(sqrt(diff(bd_pts(:,1)).^2+diff(bd_pts(:,2)).^2));
        bd_intervals(1,:)=[0 0];
        %sum_dist(1)=0;
        %sum_n_dist(1)=0;
        for i=1:nbpts-1
            pt2pt_dist=norm(bd_pts(i+1,:)-bd_pts(i,:));
            n_p2p_dist=pt2pt_dist./bd_length;
            bd_intervals(i+1,1)=bd_intervals(i,1)+n_p2p_dist;
            bd_intervals(i+1,2)=bd_intervals(i,2)+pt2pt_dist;
            %sum_n_dist(i+1)=sum_n_dist(i)+n_p2p_dist;
            %sum_dist(i+1)=sum_dist(i)+pt2pt_dist;
        end   
    case 'cubic'
%           code
    otherwise
        % error('no boundary representation selected')
        % default is linear
        [bd_length,bd_intervals]=boundary(bd_pts,'linear')
end
return