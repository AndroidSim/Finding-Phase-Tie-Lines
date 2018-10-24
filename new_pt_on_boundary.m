function bpt = new_pt_on_boundary(bd_pts,bd_rep,s,origin,end_pt)

% s == the fraction of total arc length
% origin == the origin on boundary from which the fraction refers == [Xso Xco]
% end_pt == pt where the distance from orgin to end_pt = total arc length

% output arg == [Xs Xc], the point giving the distance == frac*total arc length

[nbpts,bxy]=size(bd_pts);
if (~isequal(bxy,2))
	error('number of columns of boundary point matrix is not equal to 2 (Xs and Xc coords)')
end
if (~isequal(size(origin,1),1) & ~isequal(size(origin,2),2)) % or ~isequal(length(origin),prod(origin))
    error('origin must be a 1x2 vector')
end
if (s < 0 | s > 2)
    error('s must be between 0 and 2')
end
if (s > 1)
    s=s-1;
end
if (s == 0)
    bpt=origin;
    return
end
if (s == 1)
    bpt=end_pt;
    return
end
switch bd_rep
    case 'linear'
        [total_dist,dist_intervals]=bd_arc_length(bd_pts,bd_rep,bd_pts(1,:),bd_pts(end,:));
        sum_dist(1)=0;
        sum_n_dist(1)=0;
        for i=1:nbpts-1
            pt2pt_dist=norm(bd_pts(i+1,:)-bd_pts(i,:));
            n_p2p_dist=pt2pt_dist./total_dist;
            sum_n_dist(i+1)=sum_n_dist(i)+n_p2p_dist;
            sum_dist(i+1)=sum_dist(i)+pt2pt_dist;
        end
        
        jl=0;
        ju=nbpts+1;
        ascnd=(sum_n_dist(end) > sum_n_dist(1));
        while (ju-jl > 1)
            jm=round((ju+jl)./2);
            if (s > sum_n_dist(jm) == ascnd)
                jl=jm;
            else
                ju=jm;
            end
            j=jl;
        end
        % m = slope of line to next boundary pt
        deltax=bd_pts(j+1,1)-bd_pts(j,1);
        deltay=bd_pts(j+1,2)-bd_pts(j,2);
        warning off
        m=deltay./deltax;
        warning on
        % d = distance from i boundary pt to the point of interest
        d=s*total_dist-sum_dist(j);
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
            Xs=origin(1)+x_inc;
            Xc=origin(2)+y_inc;
        else
            % x coord of point of interest
            Xs=bd_pts(j,1)+x_inc;
            % y coord of point of interest
            Xc=bd_pts(j,2)+y_inc;
        end
        bpt=[Xs Xc];  
        if (iseqwithn(s,sum_n_dist(j)))
            bpt=bd_pts(j,:);
            return
        end
    case 'cubic'
%           code
    otherwise
        % error('no boundary representation selected')
        % default is linear
        bpt=new_pt_on_boundary(bd_pts,'linear',s,origin,end_pt)
end
return