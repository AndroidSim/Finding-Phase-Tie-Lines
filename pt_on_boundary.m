function bpt = pt_on_boundary(bd_pts,bd_rep,s,origin,end_pt)

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
        %keyboard
        [total_dist,dist_intervals]=bd_arc_length(bd_pts,bd_rep,origin,end_pt);
        dist_2_travel=s*total_dist;
        oonbd=0;
        eonbd=0;
        if (all(origin == end_pt))
            if (all(origin == bd_pts(1,:)))
                io=1;
                ie=length(bd_pts);
            else
                bpt=origin;
                return
            end
        else
            if (any(ismember(unique(bd_pts,'rows'),origin,'rows')))
                io=find(ismember(bd_pts,origin,'rows'));
                if (length(io) == 2)
                    io=io(1);
                end
                if (any(ismember(unique(bd_pts,'rows'),end_pt,'rows')))
                    ie=find(ismember(bd_pts,end_pt,'rows'));
                    if (length(ie) == 2)
                        ie=ie(2);
                    end
                    if (ie < io)
                        error('index of end point must be > index of origin')
                    end
                else
                    for i=1:nbpts-1
                        warning off
                        m=(bd_pts(i+1,2)-bd_pts(i,2))/(bd_pts(i+1,1)-bd_pts(i,1));
                        warning on
                        if ((iseqwithn(end_pt(2),m*(end_pt(1)-bd_pts(i,1))+bd_pts(i,2)) | (isinf(m) & bd_pts(i,1) == end_pt(1))) & isinbox(bd_pts(i,:),bd_pts(i+1,:),end_pt)) % end_pt(2) == m*(end_pt(1)-bd_pts(i,1))+bd_pts(i,2)
                            eonbd=1;
                            ie=i;
                            if (io > i)
                                error('index of end point must be > index of origin')
                            else
                                break
                            end      
                        end
                    end
                    if (eonbd ~= 1)
                        error('end point does not lie on the boundary')
                    end
                end
            else
                if (any(ismember(unique(bd_pts,'rows'),end_pt,'rows')))
                    ie=find(ismember(bd_pts,end_pt,'rows'));
                    if (length(ie) == 2)
                        ie=ie(2);
                    end
                    for i=1:nbpts-1
                        warning off
                        m=(bd_pts(i+1,2)-bd_pts(i,2))/(bd_pts(i+1,1)-bd_pts(i,1));
                        warning on
                        if ((iseqwithn(origin(2),m*(origin(1)-bd_pts(i,1))+bd_pts(i,2)) | (isinf(m) & bd_pts(i,1) == origin(1))) & isinbox(bd_pts(i,:),bd_pts(i+1,:),origin)) % origin(2) == m*(origin(1)-bd_pts(i,1))+bd_pts(i,2)
                            oonbd=1;
                            io=i;
                            if (ie < i)
                                error('index of end point must be > index of origin')
                            else
                                break
                            end
                        end
                    end
                    if (oonbd ~= 1)
                        error('origin does not lie on the boundary')
                    end
                else
                    for i=1:nbpts-1
                        warning off
                        m=(bd_pts(i+1,2)-bd_pts(i,2))/(bd_pts(i+1,1)-bd_pts(i,1));
                        warning on
                        if ((iseqwithn(origin(2),m*(origin(1)-bd_pts(i,1))+bd_pts(i,2)) | (isinf(m) & bd_pts(i,1) == origin(1))) & isinbox(bd_pts(i,:),bd_pts(i+1,:),origin)) % origin(2) == m*(origin(1)-bd_pts(i,1))+bd_pts(i,2)
                            oonbd=1;
                            io=i;
                            for j=1:nbpts-1
                                warning off
                                m=(bd_pts(j+1,2)-bd_pts(j,2))/(bd_pts(j+1,1)-bd_pts(j,1));
                                warning on
                                if ((iseqwithn(end_pt(2),m*(end_pt(1)-bd_pts(j,1))+bd_pts(j,2)) | (isinf(m) & bd_pts(i,1) == end_pt(1))) & isinbox(bd_pts(j,:),bd_pts(j+1,:),end_pt)) % end_pt(2) == m*(end_pt(1)-bd_pts(j,1))+bd_pts(j,2)
                                    eonbd=1;
                                    ie=j;
                                    if (i > j+1)
                                        error('index of end point must be > index of origin')
                                    else
                                        break
                                    end
                                end
                            end
                            if (eonbd ~= 1)
                                error('end point does not lie on the boundary')
                            end
                        end
                    end
                    if (oonbd ~= 1)
                        error('origin does not lie on the boundary')
                    end
                end
            end    
        end
       
        avail_dist=dist_intervals(1);
        dist_2_go=dist_2_travel;
        dist_traveled=0;
        i=io;
        j=i+1;
        for k=1:length(dist_intervals)
            avail_dist=dist_intervals(k);
            if (dist_2_go < avail_dist)
                % m = slope of line to next boundary pt
                deltax=bd_pts(j,1)-bd_pts(i,1);
                deltay=bd_pts(j,2)-bd_pts(i,2);
                warning off
                m=deltay./deltax;
                warning on
                % d = distance from i boundary pt to the point of interest
                d=dist_2_go;
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
                if (i == 1)
                    Xs=origin(1)+x_inc;
                    Xc=origin(2)+y_inc;
                else
                    % x coord of point of interest
                    Xs=bd_pts(i,1)+x_inc;
                    % y coord of point of interest
                    Xc=bd_pts(i,2)+y_inc;
                end
                bpt=[Xs Xc];
                break
            elseif (iseqwithn(dist_2_go,avail_dist))
                bpt=bd_pts(j,:);
                break
            else
                i=i+1;
                j=j+1;
            end
            dist_traveled=dist_traveled+avail_dist;
            dist_2_go=dist_2_travel-dist_traveled;
        end
    case 'cubic'
%           code
    otherwise
        % error('no boundary representation selected')
        % default is linear
        bpt=pt_on_boundary(bd_pts,'linear',s,origin,end_pt)
end
return