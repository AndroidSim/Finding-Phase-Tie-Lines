function bpt = get_bdy_pt(bdy_pts,bdy_rep,s,total_dist,bdy_intervals)

% s == the fraction of total arc length
% total_dist == boundary perimeter
% output arg == [Xs Xc], the point giving the distance == frac*total arc length

[nbpts,bxy] = size(bdy_pts);

if ~isequal(bxy,2)
	error('number of columns of boundary point matrix is not equal to 2 [X Y]')
end

if s < 0 | s > 2
    error('s must be between 0 and 2')
end

if s > 1
    s = s-1;
end

if isequal(s,0) 
    bpt = bdy_pts(1,:);
    return
end

if isequal(s,1) 
    bpt = bdy_pts(end,:);
    return
end

switch bdy_rep
    case 'linear'
        jl = 0;
        ju = nbpts+1;
        ascnd = (bdy_intervals(end,1) > bdy_intervals(1,1));
        
        while (ju-jl > 1)
            jm = round((ju+jl)./2);
            if (s > bdy_intervals(jm,1) == ascnd) % s > sum_n_dist(jm)  ==  ascnd
                jl = jm;
            else
                ju = jm;
            end
            j = jl;
        end
        
        % m  =  slope of line to next boundary pt
        deltax = bdy_pts(j+1,1)-bdy_pts(j,1);
        deltay = bdy_pts(j+1,2)-bdy_pts(j,2);
        warning off
        m = deltay./deltax;
        warning on
        % d  =  distance from i boundary pt to the point of interest
        d = s*total_dist-bdy_intervals(j,2);
        x_inc = (d./(sqrt(1+(m.^2))));
        y_inc = ((m.*d)./(sqrt(1+(m.^2))));
        % determine direction of changing x and y
        if (deltax < 0 & deltay < 0)
            x_inc = -x_inc;
            y_inc = -y_inc;
        elseif (deltax < 0 & deltay == 0)
            x_inc = -x_inc;
            y_inc = 0;
        elseif (deltax < 0 & deltay > 0)
            x_inc = -x_inc;
            y_inc = -y_inc;
        elseif (deltax == 0 & deltay < 0)
            x_inc = 0;
            y_inc = -d;
        elseif (deltax == 0 & deltay > 0)
            x_inc = 0;
            y_inc = d;
        elseif (deltax == 0 & deltay == 0)%this should never happen
            x_inc = 0;
            y_inc = 0;
        end
        
        if (j == 1)
            Xs = bdy_pts(1,1)+x_inc;
            Xc = bdy_pts(1,2)+y_inc;
        else
            % x coord of point of interest
            Xs = bdy_pts(j,1)+x_inc;
            % y coord of point of interest
            Xc = bdy_pts(j,2)+y_inc;
        end
        
        bpt = [Xs Xc]; 
        
        if (isequal(s,bdy_intervals(j,1))) % sum_n_dist(j)
            bpt = bdy_pts(j,:);
            return
        end
    case 'spline'
%           code
    otherwise
        % default is linear
        bpt = get_bdy_pt(bdy_pts,'linear',s,total_dist)
end
return