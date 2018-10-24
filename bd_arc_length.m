function [arc_length,dist_intervals] = bd_arc_length(bd_pts,bd_rep,a,b)

% bd_arc_length(varargin) calculates the arc length of a section of the phase boundary from a to b.
% a == coords of a [Xs Xc] on boundary == [x y]
% b == coords of b [Xs Xc] on boundary == [x y]                                    

% boundary is approximated as piece-wise linear fxn (bd_rep == 'linear') or cubic spline fxn (bd_rep == 'cubic')

[nbpts,bxy]=size(bd_pts);
if (~isequal(bxy,2))
	error('number of columns of boundary point matrix is not equal to 2 (Xs and Xc coords)')
end
if (~isequal(size(a,1),1) & ~isequal(size(a,2),2))
    error('a must be a 1x2 vector')
end
if (~isequal(size(b,1),1) & ~isequal(size(b,2),2))
    error('b must be a 1x2 vector')
end
bdy_pts=bd_pts;
bdy_pts(end,:)=[];
switch bd_rep
    case 'linear'
        % linear is summing euclidean distances: sqrt(((x(j+1)-x(j)).^2)+((y(j+1)-y(j)).^2))
        % of segments between boundary points
        arc_length=0;
        aonbd=0;
        bonbd=0;
        n=1;
        if (all(a == b))
            if (all(a == bd_pts(1,:)) & all(b == bd_pts(1,:)))
                for i=1:nbpts-1
                    pt2pt_dist=norm(bd_pts(i+1,:)-bd_pts(i,:));
                    arc_length=arc_length+pt2pt_dist;
                    dist_intervals(n)=pt2pt_dist;
                    n=n+1;
                end
            else
                arc_length=0;
                dist_intervals(n)=0;
                return
            end
        else
            if (any(ismember(bdy_pts,a,'rows'))) % unique(bd_pts,'rows')
                ia=find(ismember(bdy_pts,a,'rows'));
                if (any(ismember(bdy_pts,b,'rows')))
                    ib=find(ismember(bdy_pts,b,'rows'));
                    if (ib > ia)
                        for i=ia:ib-1
                            pt2pt_dist=norm(bd_pts(i+1,:)-bd_pts(i,:));             
                            arc_length=arc_length+pt2pt_dist;
                            dist_intervals(n)=pt2pt_dist;
                            n=n+1;
                        end
                    else
                        for i=ib:ia-1
                            pt2pt_dist=norm(bd_pts(i+1,:)-bd_pts(i,:));             
                            arc_length=arc_length+pt2pt_dist;
                            dist_intervals(n)=pt2pt_dist;
                            n=n+1;
                        end
                    end
                else
                    for i=1:nbpts-1
                        warning off
                        m=(bd_pts(i+1,2)-bd_pts(i,2))/(bd_pts(i+1,1)-bd_pts(i,1));
                        warning on
                        if ((iseqwithn(b(2),m*(b(1)-bd_pts(i,1))+bd_pts(i,2)) | (isinf(m) & bd_pts(i,1) == b(1))) & isinbox(bd_pts(i,:),bd_pts(i+1,:),b)) % b(2) == m*(b(1)-bd_pts(i,1))+bd_pts(i,2)
                            bonbd=1;
                            if (ia < i)
                                for j=ia:i-1
                                    pt2pt_dist=norm(bd_pts(j+1,:)-bd_pts(j,:));             
                                    arc_length=arc_length+pt2pt_dist;
                                    dist_intervals(n)=pt2pt_dist;
                                    n=n+1;
                                end
                                arc_length=arc_length+norm(b-bd_pts(i,:));
                                dist_intervals(n)=norm(b-bd_pts(i,:));
                                break
                            elseif (ia > i+1)
                                arc_length=arc_length+norm(bd_pts(i+1,:)-b);
                                dist_intervals(n)=norm(bd_pts(i+1,:)-b);
                                n=n+1;
                                for j=i+1:ia-1
                                    pt2pt_dist=norm(bd_pts(j+1,:)-bd_pts(j,:));             
                                    arc_length=arc_length+pt2pt_dist;
                                    dist_intervals(n)=pt2pt_dist;
                                    n=n+1;
                                end
                                break
                            else
                                arc_length=arc_length+norm(b-bd_pts(ia,:));
                                dist_intervals(n)=norm(b-bd_pts(ia,:));
                                break
                            end      
                        end
                    end
                    if (bonbd ~= 1)
                        error('b does not lie on the boundary')
                    end
                end
            else
                if (any(ismember(bdy_pts,b,'rows')))
                    ib=find(ismember(bdy_pts,b,'rows'));
                    for i=1:nbpts-1
                        warning off
                        m=(bd_pts(i+1,2)-bd_pts(i,2))/(bd_pts(i+1,1)-bd_pts(i,1));
                        warning on
                        if ((iseqwithn(a(2),m*(a(1)-bd_pts(i,1))+bd_pts(i,2)) | (isinf(m) & bd_pts(i,1) == a(1))) & isinbox(bd_pts(i,:),bd_pts(i+1,:),a)) % (a(2) == m*(a(1)-bd_pts(i,1))+bd_pts(i,2)
                            aonbd=1;
                            if (ib < i)
                                for j=ib:i-1
                                    pt2pt_dist=norm(bd_pts(j+1,:)-bd_pts(j,:));             
                                    arc_length=arc_length+pt2pt_dist;
                                    dist_intervals(n)=pt2pt_dist;
                                    n=n+1;
                                end
                                arc_length=arc_length+norm(a-bd_pts(i,:));
                                dist_intervals(n)=norm(a-bd_pts(i,:));
                                break
                            elseif (ib > i+1)
                                arc_length=arc_length+norm(bd_pts(i+1,:)-a);
                                dist_intervals=norm(bd_pts(i+1,:)-a);
                                n=n+1;
                                for j=i+1:ib-1
                                    pt2pt_dist=norm(bd_pts(j+1,:)-bd_pts(j,:));             
                                    arc_length=arc_length+pt2pt_dist;
                                    dist_intervals(n)=pt2pt_dist;
                                    n=n+1;
                                end
                                break
                            else
                                arc_length=arc_length+norm(a-bd_pts(ib,:));
                                dist_intervals(n)=norm(a-bd_pts(ib,:));
                                break
                            end
                        end
                    end
                    if (aonbd ~= 1)
                        error('a does not lie on the boundary')
                    end
                else
                    for i=1:nbpts-1
                        warning off
                        m=(bd_pts(i+1,2)-bd_pts(i,2))/(bd_pts(i+1,1)-bd_pts(i,1));
                        warning on
                        %keyboard
                        if ((iseqwithn(a(2),m*(a(1)-bd_pts(i,1))+bd_pts(i,2)) | (isinf(m) & bd_pts(i,1) == a(1))) & isinbox(bd_pts(i,:),bd_pts(i+1,:),a)) % (a(2) == m*(a(1)-bd_pts(i,1))+bd_pts(i,2)
                            aonbd=1;
                            for j=1:nbpts-1
                                warning off
                                m=(bd_pts(j+1,2)-bd_pts(j,2))/(bd_pts(j+1,1)-bd_pts(j,1));
                                warning on
                                %keyboard
                                if ((iseqwithn(b(2),m*(b(1)-bd_pts(j,1))+bd_pts(j,2)) | (isinf(m) & bd_pts(i,1) == b(1))) & isinbox(bd_pts(j,:),bd_pts(j+1,:),b)) % (b(2) == m*(b(1)-bd_pts(j,1))+bd_pts(j,2)
                                    bonbd=1;
                                    if (i+1 < j)
                                        arc_length=arc_length+norm(a-bd_pts(i+1,:));
                                        dist_intervals(n)=norm(a-bd_pts(i+1,:));
                                        n=n+1;
                                        for k=i+1:j-1
                                            pt2pt_dist=norm(bd_pts(k+1,:)-bd_pts(k,:));             
                                            arc_length=arc_length+pt2pt_dist;
                                            dist_intervals(n)=pt2pt_dist;
                                            n=n+1;
                                        end
                                        arc_length=arc_length+norm(b-bd_pts(j,:));
                                        dist_intervals(n)=norm(b-bd_pts(j,:));
                                        break
                                    elseif (i > j+1)
                                        arc_length=arc_length+norm(bd_pts(j+1,:)-b);
                                        dist_intervals(n)=norm(bd_pts(j+1,:)-b);
                                        n=n+1;
                                        for k=j+1:i
                                            pt2pt_dist=norm(bd_pts(k+1,:)-bd_pts(k,:));             
                                            arc_length=arc_length+pt2pt_dist;
                                            dist_intervals(n)=pt2pt_dist;
                                            n=n+1;
                                        end
                                        arc_length=arc_length+norm(bd_pts(i,:)-a);
                                        dist_intervals(n)=norm(bd_pts(i,:)-a);
                                        break
                                    elseif (i+1 == j)
                                        arc_length=arc_length+norm(a-bd_pts(i+1,:));
                                        dist_intervals(n)=norm(a-bd_pts(i+1,:));
                                        n=n+1;
                                        arc_length=arc_length+norm(b-bd_pts(j,:));
                                        dist_intervals(n)=norm(b-bd_pts(j,:));
                                        break
                                    elseif (i == j+1)
                                        arc_length=arc_length+norm(a-bd_pts(i,:));
                                        dist_intervals(n)=norm(a-bd_pts(i,:));
                                        n=n+1;
                                        arc_length=arc_length+norm(b-bd_pts(j+1,:));
                                        dist_intervals(n)=norm(b-bd_pts(j+1,:));
                                        break
                                    else
                                        arc_length=arc_length+norm(a-b);
                                        dist_intervals(n)=norm(a-b);
                                        break
                                    end
                                end
                            end
                            if (bonbd ~= 1)
                                error('b does not lie on the boundary')
                            end
                        end
                    end
                    if (aonbd ~= 1)
                        error('a does not lie on the boundary')
                    end
                end
            end    
        end     
    case 'cubic'
%           code
    otherwise
        % error('no boundary representation selected')
        % default is linear
        arc_length=bd_arc_length(bd_pts,'linear',a,b);
end
return