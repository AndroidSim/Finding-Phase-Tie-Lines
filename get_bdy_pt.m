function bpt = get_bdy_pt(bd_pts,bd_rep,s,total_dist,bd_intervals)

% s == the fraction of total arc length
% total_dist == boundary perimeter
% output arg == [Xs Xc], the point giving the distance == frac*total arc length

[nbpts,bxy]=size(bd_pts);
if (~isequal(bxy,2))
	error('number of columns of boundary point matrix is not equal to 2 (Xs and Xc coords)')
end
if (s < 0 | s > 2)
    error('s must be between 0 and 2')
end
if (s > 1)
    s=s-1;
end
if (iseqwithn(s,0,10^-4)) % s == 0
    bpt=bd_pts(1,:);
    return
end
if (iseqwithn(s,1,10^-4)) % s == 1
    bpt=bd_pts(end,:);
    return
end
switch bd_rep
    case 'linear'
        jl=0;
        ju=nbpts+1;
        %ascnd=(sum_n_dist(end) > sum_n_dist(1));
        ascnd=(bd_intervals(end,1) > bd_intervals(1,1));
        while (ju-jl > 1)
            jm=round((ju+jl)./2);
            if (s > bd_intervals(jm,1) == ascnd) % s > sum_n_dist(jm) == ascnd
                jl=jm;
            else
                ju=jm;
            end
            j=jl;
        end
        % m = slope of line to next boundary pt
        deltax=bd_pts(j+1,1)-bd_pts(j,1);
        deltay=bd_pts(j+1,2)-bd_pts(j,2);
        %warning off
        m=deltay./deltax;
        %warning on
        % d = distance from i boundary pt to the point of interest
        %d=s*total_dist-sum_dist(j);
        d=s*total_dist-bd_intervals(j,2);
        x_inc=(d./(sqrt(1+(m.^2))));
        y_inc=((m.*d)./(sqrt(1+(m.^2))));
        % determine direction of changing x and y
        %keyboard
        if (deltax < 0 & deltay < 0)
            x_inc=-x_inc;
            y_inc=-y_inc;
        elseif (deltax < 0 & deltay == 0)
            x_inc=-x_inc;
            y_inc=0;
        elseif (deltax < 0 & deltay > 0)
            x_inc=-x_inc;
            y_inc=-y_inc;
        elseif (deltax == 0 & deltay < 0)
            x_inc=0;
            y_inc=-d;
        elseif (deltax == 0 & deltay > 0)
            x_inc=0;
            y_inc=d;
        elseif (deltax == 0 & deltay == 0)%this should never happen
            x_inc=0;
            y_inc=0;
        end
        if (j == 1)
            Xs=bd_pts(1,1)+x_inc;
            Xc=bd_pts(1,2)+y_inc;
        else
            % x coord of point of interest
            Xs=bd_pts(j,1)+x_inc;
            % y coord of point of interest
            Xc=bd_pts(j,2)+y_inc;
        end
        bpt=[Xs Xc];  
        if (iseqwithn(s,bd_intervals(j,1),10^-4)) % sum_n_dist(j)
            bpt=bd_pts(j,:);
            return
        end
    case 'cubic'
%           code
    otherwise
        % error('no boundary representation selected')
        % default is linear
        bpt=get_bd_pt(bd_pts,'linear',s,total_dist)
end
return