function [bdy_length,bdy_intervals] = bdy_fxn(bdy_pts,bdy_rep)
% [bdy_length,bdy_intervals] = bdy_fxn(bdy_pts,bdy_rep)
% the boundary function (right now) calculates the perimeter of the phase boundary given 
% by the boundary points in bd_pts.  bd_rep is how the boundary is represented, linear or cubic spline

[nbpts,bxy] = size(bdy_pts);

if bxy > 3 | bxy < 2
	error('the coordinate dimension of the boundary pts must be 3 (ternary) or 2 (cartesian)');
else
    if bxy == 3
        bdy_pts = tern2cart(bdy_pts,1);
        tern = true;
        cart = false;
    end
    
    if bxy == 2
        tern = false;
        cart = true;
    end
end

switch bdy_rep
    case 'linear'
        % bdy_length = sum(sqrt(diff(bdy_pts(:,1)).^2+diff(bdy_pts(:,2)).^2));
        bdy_length = sum(sqrt(sum(diff(bdy_pts).^2,2)));
        bdy_intervals(1,:) = [0 0];
        
        for i=1:nbpts-1
            pt2pt_dist = norm(bdy_pts(i+1,:)-bdy_pts(i,:));
            n_p2p_dist = pt2pt_dist./bdy_length;
            bdy_intervals(i+1,1) = bdy_intervals(i,1)+n_p2p_dist;
            bdy_intervals(i+1,2) = bdy_intervals(i,2)+pt2pt_dist;
        end
        
    case 'spline'
%           code
    otherwise
        % default is linear
        [bdy_length,bdy_intervals] = bdy_fxn(bdy_pts,'linear')
end

return