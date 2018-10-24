function [tielines,alpha_beta_pts,tangent_pts]=tangent_tielines(boundary_pts,npts,ntielines)
% the function tangent_tielines generates multiple sets of tielines from the tangent lines
% of the choosen critical points
[nbpts,bxy]=size(boundary_pts);
if (bxy ~= 2)
	error('num cols of boundary points not equal to 2 (x and y)')
end
xmin=min(boundary_pts(:,1));
xmax=max(boundary_pts(:,1));
ymin=min(boundary_pts(:,2));
ymax=max(boundary_pts(:,2));

%bs_pts=[10 5;15 5;20 5;25 10;25 15;25 20;20 25;15 25;10 25;5 20;5 15;5 10;10 5]
%plot(points(:,1),points(:,2),'or');
%p = polyfit(points(:,1),points(:,2),2);
%x=0:0.5:4;
%y = polyval(p,x);
%hold on,plot(x,y,'-sg'),hold off
%line2=cat(2,zeros(1,21)',[-10:10]')

% define rotation matrix R for rotating a vector counterclockwise 90 degrees
R=[0 -1;1 0];
% loop over all boundary pts as first critical point
tielines=cell(nbpts,(nbpts-2));
alpha_beta_pts=cell(nbpts,(nbpts-2));
tangent_pts=cell(nbpts,(nbpts-2));
for i=1:nbpts
    % loop over possible critical point pairs
    jset=1;
    for j=1:(nbpts-1)
        if (j == i)
            % critical points not same point
            continue;
        end
        crit_pt1=boundary_pts(i,:);
        crit_pt2=boundary_pts(j,:);
        % picking to points surrounding critical points, assumes first boundary point and last
        % boundary point are the same
        if (i == 1)
            if (j == nbpts)
                b4_crit_pt1=boundary_pts(nbpts-1,:);
                af_crit_pt1=boundary_pts(i+1,:);
                b4_crit_pt2=boundary_pts(j-1,:);
                af_crit_pt2=boundary_pts(2,:);
            else
                b4_crit_pt1=boundary_pts(nbpts-1,:);
                af_crit_pt1=boundary_pts(i+1,:);
                b4_crit_pt2=boundary_pts(j-1,:);
                af_crit_pt2=boundary_pts(j+1,:);
            end
        elseif (i == nbpts)
            if (j == 1)
                b4_crit_pt1=boundary_pts(i-1,:);
                af_crit_pt1=boundary_pts(2,:);
                b4_crit_pt2=boundary_pts(nbpts-1,:);
                af_crit_pt2=boundary_pts(j+1,:);
            else
                b4_crit_pt1=boundary_pts(i-1,:);
                af_crit_pt1=boundary_pts(2,:);
                b4_crit_pt2=boundary_pts(j-1,:);
                af_crit_pt2=boundary_pts(j+1,:);
            end
        else
            if (j == 1)
                b4_crit_pt1=boundary_pts(i-1,:);
                af_crit_pt1=boundary_pts(i+1,:);
                b4_crit_pt2=boundary_pts(nbpts-1,:);
                af_crit_pt2=boundary_pts(j+1,:);
            elseif (j == nbpts)
                b4_crit_pt1=boundary_pts(i-1,:);
                af_crit_pt1=boundary_pts(i+1,:);
                b4_crit_pt2=boundary_pts(j-1,:);
                af_crit_pt2=boundary_pts(2,:);
            else
                b4_crit_pt1=boundary_pts(i-1,:);
                af_crit_pt1=boundary_pts(i+1,:);
                b4_crit_pt2=boundary_pts(j-1,:);
                af_crit_pt2=boundary_pts(j+1,:);
            end
        end
        around_cpt1=cat(1,cat(1,b4_crit_pt1,crit_pt1),af_crit_pt1);
        around_cpt2=cat(1,cat(1,b4_crit_pt2,crit_pt2),af_crit_pt2);
        around_cpts=cat(1,around_cpt1,around_cpt2);
        % find each line that is tangent to a critical point:
        %---orthogonal vector to normal vector, or
        %---fit the critical point and two surrounding points to a quadratic poly
        %       need to rotate point so i can fit with y=p(x)
        icrit=1;
        tan_pts=zeros(npts,2,2);
        for k=2:3:5
            b4_v=around_cpts(k-1,:)-around_cpts(k,:);
            af_v=around_cpts(k+1,:)-around_cpts(k,:);
            %angle_between=acos((dot(b4_v,af_v)./(norm(b4_v).*norm(af_v))));
            %if (~isreal(angle_between))
            %    error('angle between before and after vector imaginary')
            %end
            bisect_v=(norm(b4_v)*af_v+norm(af_v)*b4_v);
            % calculate tangent as vector orthogonal to normal vector, which is the
            % bisecting vector bisect_v
            normal_v=-bisect_v;
            tangent_v=R*normal_v';
            tangent_v=tangent_v';
            % get line that corresponds to vector, line contains critical point
            tan_slope(icrit)=tangent_v(2)./tangent_v(1);
            if (tangent_v == [0 0])
                if (b4_v(1) == 0 & af_v(1) == 0)
                    tan_slope(icrit)=Inf;
                elseif (b4_v(2) == 0 & af_v(2) == 0)
                    tan_slope(icrit)=0;
                end
            end        
            if (tan_slope(icrit) == 0)
                xtan=linspace(xmin,xmax,npts);
                ytan=tan_slope(icrit).*(xtan-around_cpts(k,1))+around_cpts(k,2);
            elseif (isinf(tan_slope(icrit)))
                ytan=linspace(ymin,ymax,npts);
                xtan=ones(1,npts).*(around_cpts(k,1));
            else
                xtan=linspace(xmin,xmax,npts);
                ytan=tan_slope(icrit).*(xtan-around_cpts(k,1))+around_cpts(k,2);
            end
            tan_pts(1:end,1:end,icrit)=cat(2,xtan',ytan');
            icrit=icrit+1;
        end
        tangent_pts{i,jset}=tan_pts;
        %tan_pts=zeros(length(xmin:xmax),2,2);
        %icrit=1;
        %for k=2:3:5
        %    b4_v=around_cpts(k-1,:)-around_cpts(k,:);
        %    af_v=around_cpts(k+1,:)-around_cpts(k,:);
        %    angle_between=acos((dot(b4_v,af_v)./(norm(b4_v).*norm(af_v))));
        %    bisect_v=(norm(b4_v)*af_v+norm(af_v)*b4_v);
            % calculate tangent from the derivative of a quadratic polynomical fit to
            % the critical point and its two surrounding points
        %    if (bisect_v(1) == 0)
                %tangent_v=R*normal_v';
                %tangent_v=tangent_v';
                % get line that corresponds to vector, line contains critical point
                %tan_slope=tangent_v(2)./tangent_v(1);
        %        tan_slope(icrit)=0;
        %        xtan=xmin:xmax;
        %        ytan=tan_slope.*(xtan-around_cpts(k,1))+around_cpts(k,2);
        %        tan_pts(1:end,1:end,icrit)=cat(2,xtan',ytan');
        %    elseif (bisect_v(2) == 0)
                %tangent_v=R*normal_v';
                %tangent_v=tangent_v';
                % get line that corresponds to vector, line contains critical point
                %tan_slope=tangent_v(2)./tangent_v(1);
        %        tan_slope(icrit)=Inf;
        %        xtan=ones(1,length(xmin:xmax)).*around_cpts(k,1);
        %        ytan=linspace(ymin,ymax,length(xtan));
        %        tan_pts(1:end,1:end,icrit)=cat(2,xtan',ytan');
        %    else
        %        if (bisect_v(2) > 0)
                    % rotate so critical point on y-axis (vertical) so to fit with polyfit
        %            angle2subtract=atan(bisect_v(2)/bisect_v(1));
        %            if (angle2subtract > 0)
        %                theta=(pi/2)-angle2subtract;
        %                rotation_matrix=[cos(theta) -sin(theta);sin(theta) cos(theta)];
        %            else
        %                theta=-angle2subtract;
        %                rotation_matrix=[cos(theta) sin(theta);-sin(theta) cos(theta)];
        %            end
        %            r_b4_pt=rotation_matrix*around_cpts(k-1,:)';
        %            r_b4_pt=r_b4_pt';
        %            r_crit_pt=rotation_matrix*around_cpts(k,:)';
        %            r_crit_pt=r_crit_pt';
        %            r_af_pt=rotation_matrix*around_cpts(k+1,:)';
        %            r_af_pt=r_af_pt';
        %            rot_pts=cat(1,cat(1,r_b4_pt,r_crit_pt),r_af_pt);
        %            p=polyfit(rot_pts(:,1),rot_pts(:,2),2);
        %            p_der=polyder(p);
        %            r_tan_slope=polyval(p_der,r_crit_pt(1));
        %            r_xtan=xmin:xmax;
        %            r_ytan=r_tan_slope.*(r_xtan-rot_pts(2,1))+rot_pts(2,2);
        %            r_tan_pts=cat(1,r_xtan,r_ytan);
        %            tan_points=zeros(size(r_tan_pts));
        %            for r=1:length(r_tan_pts)
        %                tan_points(:,r)=rotation_matrix'*r_tan_pts(:,r);
        %            end
        %            tan_pts(1:end,1:end,icrit)=tan_points';
        %            tan_slope(icrit)=(tan_pts(1,2,icrit)-tan_pts(2,2,icrit)./tan_pts(1,1,icrit)-tan_pts(2,1,icrit));
        %        end
        %        if (bisect_v(2) < 0)
                    % rotate so critical point on y-axis (vertical) so to fit with polyfit
        %            angle2subtract=atan(bisect_v(2)/bisect_v(1));
        %            if (angle2subtract > 0)
        %                theta=angle2subtract;
        %                rotation_matrix=[cos(theta) -sin(theta);sin(theta) cos(theta)];
        %            else
        %                theta=(-pi/2)-angle2subtract;
        %                rotation_matrix=[cos(theta) sin(theta);-sin(theta) cos(theta)];
        %            end
        %            r_b4_pt=rotation_matrix*around_cpts(k-1,:)';
        %            r_b4_pt=r_b4_pt';
        %            r_crit_pt=rotation_matrix*around_cpts(k,:)';
        %            r_crit_pt=r_crit_pt';
        %            r_af_pt=rotation_matrix*around_cpts(k+1,:)';
        %            r_af_pt=r_af_pt';
        %            rot_pts=cat(1,cat(1,r_b4_pt,r_crit_pt),r_af_pt);
        %            p=polyfit(rot_pts(:,1),rot_pts(:,2),2);
        %            p_der=polyder(p);
        %            r_tan_slope=polyval(p_der,r_crit_pt(1));
        %            r_xtan=xmin:xmax;
        %            r_ytan=r_tan_slope.*(r_xtan-rot_pts(2,1))+rot_pts(2,2);
        %            r_tan_pts=cat(1,r_xtan,r_ytan);
        %            tan_points=zeros(size(r_tan_pts));
        %            for r=1:length(r_tan_pts)
        %                tan_points(:,r)=rotation_matrix'*r_tan_pts(:,r);
        %            end
        %            tan_pts(1:end,1:end,icrit)=tan_points';
        %            tan_slope(icrit)=(tan_pts(1,2,icrit)-tan_pts(2,2,icrit)./tan_pts(1,1,icrit)-tan_pts(2,1,icrit));
        %        end
        %    end
        %    icrit=icrit+1;
        %end
        %tangent_pts{i,jset}=tan_pts;
        % find intersection point of both tangent lines
        x_int=(tan_slope(2).*crit_pt2(1)-crit_pt2(2)+crit_pt1(2)-tan_slope(1).*crit_pt1(1))./(tan_slope(2)-tan_slope(1));
        y_int=(tan_slope(1).*crit_pt2(2)-tan_slope(2).*crit_pt1(2)+tan_slope(1).*tan_slope(2).*crit_pt1(1)-tan_slope(1).*tan_slope(2).*crit_pt2(1))./(tan_slope(1)-tan_slope(2));
        % there are four possibilities: 1) no intersection points (tangent lines are parallel)
        %                                   -so generate tielines by "parallel_tielines" function
        %                               2) infinite number intersection points (tangent lines are the same)
        %                               3) intersection point = one of the critical points
        %                               4) one unique intersection ~= to either critical point
        % find vectors starting at intersection point and ending at critical points and angle between them
        if (isnan(x_int) | isnan(y_int))
            if (tan_slope(1) == 0 & tan_slope(2) == 0)
                if (crit_pt2(2) == tan_slope(1).*(crit_pt2(1)-crit_pt1(1))+crit_pt1(2))
                    % tangent line of crit_pt1 == tangent line of crit_pt2
                    % pick test vector to rotate to find where boundary is
                    if (crit_pt1(1) > crit_pt2(1))
                        x1=crit_pt1(1)+1;
                        y1=crit_pt1(2);
                        x2=crit_pt2(1)-1;
                        y2=crit_pt2(2);
                        v1=[x1 y1]-[crit_pt1(1) crit_pt1(2)];
                        v2=[x2 y2]-[crit_pt2(1) crit_pt2(2)];
                    else
                        x1=crit_pt1(1)-1;
                        y1=crit_pt1(2);
                        x2=crit_pt2(1)+1;
                        y2=crit_pt2(2);
                        v1=[x1 y1]-[crit_pt1(1) crit_pt1(2)];
                        v2=[x2 y2]-[crit_pt2(1) crit_pt2(2)];
                    end
                    x_int=(crit_pt1(1)+crit_pt2(1))./2;
                    y_int=(crit_pt1(2)+crit_pt2(2))./2;
                    test_v=v1;
                    angle_between=acos((dot(v1,v2)./(norm(v1).*norm(v2)))); %should be pi
                    % divide angle into ntielines number of smaller angles
                    angle_inc=angle_between./(ntielines+1);
                    % rotate test vector 90 degrees counterclockwise and test for intersection with boundary
                    rot_v=R*test_v';
                    test_line=vec2line(rot_v,crit_pt1,xmin,xmax,ymin,ymax,npts);
                    if (isintersect(test_line,boundary_pts))
                        rot_mat=[cos(angle_inc) -sin(angle_inc);sin(angle_inc) cos(angle_inc)];
                        rot_vec=v1';
                    else
                        % rotate test vector 90 degrees clockwise and test for intersection with boundary
                        rot_v=R'*test_v';
                        test_line=vec2line(rot_v,crit_pt1,xmin,xmax,ymin,ymax,npts);
                        if (isintersect(test_line,boundary_pts))
                            rot_mat=[cos(angle_inc) sin(angle_inc);-sin(angle_inc) cos(angle_inc)];
                            rot_vec=v1';
                        else
                            error('no boundary intersection for tan_slopes = 0 and same line')
                        end
                    end
                else
                    % tangent line of crit_pt1 is parallel to tangent line of crit_pt2
                    [tielines{i,jset},alpha_beta_pts{i,jset}]=parallel_tielines(boundary_pts,npts,'slope',0);
                    continue;
                end
            elseif (isinf(tan_slope(1)) & tan_slope(2) == 0)
                % pick test vector to rotate to find where boundary is
                % intersection point exists but not determined above, so find intersection of lines 
                % traversing from xmin-10 -> xmax+10 or ymin-10 -> ymax+10
                x1ine1=ones(1,npts).*crit_pt1(1);
                yline1=linspace((ymin-10),(ymax+10),npts);
                xline2=linspace((xmin-10),(xmax+10),npts);
                yline2=ones(1,npts).*crit_pt2(2);
                [xint,yint]=polyxpoly(xline1',yline1',xline2',yline2','unique')
                if (isempty(xint) | isempty(yint) | size(xint,1) > 1 | size(yint,1) > 1)
                    error('no intersection, tan_slope 1=Inf and tan_slope 2=0')
                end
                x_int=xint;
                y_int=yint;
                v1=[crit_pt1(1) crit_pt1(2)]-[xint yint];
                v2=[crit_pt2(1) crit_pt2(2)]-[xint yint];
                test_v=v1;
                angle_between=acos((dot(v1,v2)./(norm(v1).*norm(v2)))); %should be pi/2
                % divide angle into ntielines number of smaller angles
                angle_inc=angle_between./(ntielines+1);
                % rotate test vector 45 degrees counterclockwise and test for intersection with boundary
                test_rot=[cos(pi/4) -sin(pi/4);sin(pi/4) cos(pi/4)];
                rot_v=test_rot*test_v';
                test_line=vec2line(rot_v,crit_pt1,xmin,xmax,ymin,ymax,npts);
                if (isintersect(test_line,boundary_pts))
                    rot_mat=[cos(angle_inc) -sin(angle_inc);sin(angle_inc) cos(angle_inc)];
                    rot_vec=v1';
                else
                    % rotate test vector 45 degrees clockwise and test for intersection with boundary
                    rot_v=test_rot'*test_v';
                    test_line=vec2line(rot_v,crit_pt1,xmin,xmax,ymin,ymax,npts);
                    if (isintersect(test_line,boundary_pts))
                        rot_mat=[cos(angle_inc) sin(angle_inc);-sin(angle_inc) cos(angle_inc)];
                        rot_vec=v1';
                    end
                end
                % rotate test vector 45 degrees counterclockwise and test for intersection with boundary
                test_v=-v1;
                rot_v=test_rot*test_v';
                test_line=vec2line(rot_v,crit_pt1,xmin,xmax,ymin,ymax,npts);
                if (isintersect(test_line,boundary_pts))
                    rot_mat=[cos(angle_inc) -sin(angle_inc);sin(angle_inc) cos(angle_inc)];
                    rot_vec=-v1';
                else
                    % rotate test vector 45 degrees clockwise and test for intersection with boundary
                    rot_v=test_rot'*test_v';
                    test_line=vec2line(rot_v,crit_pt1,xmin,xmax,ymin,ymax,npts);
                    if (isintersect(test_line,boundary_pts))
                        rot_mat=[cos(angle_inc) sin(angle_inc);-sin(angle_inc) cos(angle_inc)];
                        rot_vec=-v1';
                    else
                        error('no boundary intersection for tan_slope1 = Inf and tan_slope2 = 0')
                    end
                end
            elseif (tan_slope(1) == 0 & isinf(tan_slope(2)))
                % pick test vector to rotate to find where boundary is
                % intersection point exists but not determined above, so find intersection of lines 
                % traversing from xmin-10 -> xmax+10 or ymin-10 -> ymax+10
                x1ine1=linspace((xmin-10),(xmax+10),npts);
                yline1=ones(1,npts).*crit_pt1(2);
                xline2=ones(1,npts).*crit_pt2(1);
                yline2=linspace((ymin-10),(ymax+10),npts);
                [xint,yint]=polyxpoly(xline1',yline1',xline2',yline2','unique')
                if (isempty(xint) | isempty(yint) | size(xint,1) > 1 | size(yint,1) > 1)
                    error('no intersection, tan_slope 1=Inf and tan_slope 2=0')
                end
                x_int=xint;
                y_int=yint;
                v1=[crit_pt1(1) crit_pt1(2)]-[xint yint];
                v2=[crit_pt2(1) crit_pt2(2)]-[xint yint];
                test_v=v2;
                angle_between=acos((dot(v1,v2)./(norm(v1).*norm(v2)))); %should be pi/2
                % divide angle into ntielines number of smaller angles
                angle_inc=angle_between./(ntielines+1);
                % rotate test vector 45 degrees counterclockwise and test for intersection with boundary
                test_rot=[cos(pi/4) -sin(pi/4);sin(pi/4) cos(pi/4)];
                rot_v=test_rot*test_v';
                test_line=vec2line(rot_v,crit_pt2,xmin,xmax,ymin,ymax,npts);
                if (isintersect(test_line,boundary_pts))
                    rot_mat=[cos(angle_inc) -sin(angle_inc);sin(angle_inc) cos(angle_inc)];
                    rot_vec=v2';
                else
                    % rotate test vector 45 degrees clockwise and test for intersection with boundary
                    rot_v=test_rot'*test_v';
                    test_line=vec2line(rot_v,crit_pt2,xmin,xmax,ymin,ymax,npts);
                    if (isintersect(test_line,boundary_pts))
                        rot_mat=[cos(angle_inc) sin(angle_inc);-sin(angle_inc) cos(angle_inc)];
                        rot_vec=v2';
                    end
                end
                % rotate test vector 45 degrees counterclockwise and test for intersection with boundary
                test_v=-v2;
                rot_v=test_rot*test_v';
                test_line=vec2line(rot_v,crit_pt2,xmin,xmax,ymin,ymax,npts);
                if (isintersect(test_line,boundary_pts))
                    rot_mat=[cos(angle_inc) -sin(angle_inc);sin(angle_inc) cos(angle_inc)];
                    rot_vec=-v2';
                else
                    % rotate test vector 45 degrees clockwise and test for intersection with boundary
                    rot_v=test_rot'*test_v';
                    test_line=vec2line(rot_v,crit_pt2,xmin,xmax,ymin,ymax,npts);
                    if (isintersect(test_line,boundary_pts))
                        rot_mat=[cos(angle_inc) sin(angle_inc);-sin(angle_inc) cos(angle_inc)];
                        rot_vec=-v2';
                    else
                        error('no boundary intersection for tan_slope1 = 0 and tan_slope2 = Inf')
                    end
                end 
            elseif (isinf(tan_slope(1)) & isinf(tan_slope(2)))
                if (crit_pt2(2) == tan_slope(1).*(crit_pt2(1)-crit_pt1(1))+crit_pt1(2))
                    % tangent line of crit_pt1 == tangent line of crit_pt2
                    % pick test vector to rotate to find where boundary is
                    if (crit_pt1(2) > crit_pt2(2))
                        x1=crit_pt1(1);
                        y1=crit_pt1(2)+1;
                        x2=crit_pt2(1);
                        y2=crit_pt2(2)-1;
                        v1=[x1 y1]-[crit_pt1(1) crit_pt1(2)];
                        v2=[x2 y2]-[crit_pt2(1) crit_pt2(2)];
                    else
                        x1=crit_pt1(1);
                        y1=crit_pt1(2)-1;
                        x2=crit_pt2(1);
                        y2=crit_pt2(2)+1;
                        v1=[x1 y1]-[crit_pt1(1) crit_pt1(2)];
                        v2=[x2 y2]-[crit_pt2(1) crit_pt2(2)];
                    end
                    x_int=(crit_pt1(1)+crit_pt2(1))./2;
                    y_int=(crit_pt1(2)+crit_pt2(2))./2;
                    test_v=v1;
                    angle_between=acos((dot(v1,v2)./(norm(v1).*norm(v2)))); %should be pi
                    % divide angle into ntielines number of smaller angles
                    angle_inc=angle_between./(ntielines+1);
                    % rotate test vector 90 degrees counterclockwise and test for intersection with boundary
                    rot_v=R*test_v';
                    test_line=vec2line(rot_v,crit_pt1,xmin,xmax,ymin,ymax,npts);
                    if (isintersect(test_line,boundary_pts))
                        rot_mat=[cos(angle_inc) -sin(angle_inc);sin(angle_inc) cos(angle_inc)];
                        rot_vec=v1';
                    else
                        % rotate test vector 90 degrees clockwise and test for intersection with boundary
                        rot_v=R'*test_v';
                        test_line=vec2line(rot_v,crit_pt1,xmin,xmax,ymin,ymax,npts);
                        if (isintersect(test_line,boundary_pts))
                            rot_mat=[cos(angle_inc) sin(angle_inc);-sin(angle_inc) cos(angle_inc)];
                            rot_vec=v1';
                        else
                            error('no boundary intersection for tan_slopes = Inf and same line')
                        end
                    end
                else
                    % tangent line of crit_pt1 is parallel to tangent line of crit_pt2
                    [tielines{i,jset},alpha_beta_pts{i,jset}]=parallel_tielines(boundary_pts,npts,'slope',Inf);
                    continue;
                end
            end
        elseif (isinf(x_int) & isinf(y_int))
            % happens when tan_slope(1) == tan_slope(2)
            % check
            if (tan_slope(1) ~= tan_slope(2))
                error('tangent slope 1 ~= tangent slope 2, eventhough x_int = y_int = Inf')
            end
            if (crit_pt2(2) == tan_slope(1).*(crit_pt2(1)-crit_pt1(1))+crit_pt1(2))
                % tangent line of crit_pt1 == tangent line of crit_pt2
                % pick test vector to rotate to find where boundary is
                if (crit_pt1(1) > crit_pt2(1))
                    x1=crit_pt1(1)+1;
                    y1=tan_slope(1).*(x1-crit_pt1(1))+crit_pt1(2);
                    x2=crit_pt2(1)-1;
                    y2=tan_slope(2).*(x2-crit_pt2(1))+crit_pt2(2);
                    v1=[x1 y1]-[crit_pt1(1) crit_pt1(2)];
                    v2=[x2 y2]-[crit_pt2(1) crit_pt2(2)];
                else
                    x1=crit_pt1(1)-1;
                    y1=tan_slope(1).*(x1-crit_pt1(1))+crit_pt1(2);
                    x2=crit_pt2(1)+1;
                    y2=tan_slope(2).*(x2-crit_pt2(1))+crit_pt2(2);
                    v1=[x1 y1]-[crit_pt1(1) crit_pt1(2)];
                    v2=[x2 y2]-[crit_pt2(1) crit_pt2(2)];
                end
                x_int=(crit_pt1(1)+crit_pt2(1))./2;
                y_int=(crit_pt1(2)+crit_pt2(2))./2;
                test_v=v1;
                angle_between=acos((dot(v1,v2)./(norm(v1).*norm(v2)))); %should be pi
                % divide angle into ntielines number of smaller angles
                angle_inc=angle_between./(ntielines+1);
                % rotate test vector 90 degrees counterclockwise and test for intersection with boundary
                rot_v=R*test_v';
                test_line=vec2line(rot_v,crit_pt1,xmin,xmax,ymin,ymax,npts);
                if (isintersect(test_line,boundary_pts))
                    rot_mat=[cos(angle_inc) -sin(angle_inc);sin(angle_inc) cos(angle_inc)];
                    rot_vec=v1';
                else
                    % rotate test vector 90 degrees clockwise and test for intersection with boundary
                    rot_v=R'*test_v';
                    test_line=vec2line(rot_v,crit_pt1,xmin,xmax,ymin,ymax,npts);
                    if (isintersect(test_line,boundary_pts))
                        rot_mat=[cos(angle_inc) sin(angle_inc);-sin(angle_inc) cos(angle_inc)];
                        rot_vec=v1';
                    else
                        error('no boundary intersection for tan_slopes = m and same line')
                    end
                end
            else
                % tangent line of crit_pt1 is parallel to tangent line of crit_pt2
                [tielines{i,jset},alpha_beta_pts{i,jset}]=parallel_tielines(boundary_pts,npts,'slope',tan_slope(1));
                continue;
            end
        elseif (crit_pt1 == [x_int y_int])
            if (tan_slope(1) == 0)
                x=crit_pt1(1)+1;
                y=tan_slope(1).*(x-crit_pt1(1))+crit_pt1(2);
            elseif (tan_slope(1) == Inf)
                y=crit_pt1(2)+1;
                x=crit_pt1(1);
            else
                x=crit_pt1(1)+1;
                y=tan_slope(1).*(x-crit_pt1(1))+crit_pt1(2);
            end
            v1=[x y]-[x_int y_int];
            to_cpt2=crit_pt2-[x_int y_int];
            test_v=v1;
            bisect_v=(norm(test_v)*to_cpt2+norm(to_cpt2)*test_v);
            %test_rot=[cos(angle_between/2) -sin(angle_between/2);sin(angle_between/2) cos(angle_between/2)];
            %rot_v=test_rot*test_v';
            test_line=vec2line(bisect_v,crit_pt1,xmin,xmax,ymin,ymax,npts);
            if (isintersect(test_line,boundary_pts))
                angle_between=acos((dot(test_v,to_cpt2)./(norm(test_v).*norm(to_cpt2))));
                angle_inc=angle_between./(ntielines+1);
                theta=rot_angle((test_v./norm(test_v)),(bisect_v./norm(bisect_v)),'asin');
                if (theta > 0)
                    rot_mat=[cos(angle_inc) -sin(angle_inc);sin(angle_inc) cos(angle_inc)];
                else
                    rot_mat=[cos(angle_inc) sin(angle_inc);-sin(angle_inc) cos(angle_inc)]';
                end
                rot_vec=v1';
            else
                test_v=-v1;
                bisect_v=(norm(test_v)*to_cpt2+norm(to_cpt2)*test_v);
                test_line=vec2line(bisect_v,crit_pt1,xmin,xmax,ymin,ymax,npts);
                if (isintersect(test_line,boundary_pts))
                    angle_between=acos((dot(test_v,to_cpt2)./(norm(test_v).*norm(to_cpt2))));
                    angle_inc=angle_between./(ntielines+1);
                    %keyboard
                    theta=rot_angle((test_v./norm(test_v)),(bisect_v./norm(bisect_v)),'asin');
                    if (theta > 0)
                        rot_mat=[cos(angle_inc) -sin(angle_inc);sin(angle_inc) cos(angle_inc)];
                    else
                        rot_mat=[cos(angle_inc) sin(angle_inc);-sin(angle_inc) cos(angle_inc)]';
                    end
                    rot_vec=-v1';
                else
                    error('no boundary intersection for crit_pt1 = [x_int y_int]')
                end
            end
        elseif (crit_pt2 == [x_int y_int])
            if (tan_slope(2) == 0)
                x=crit_pt2(1)+1;
                y=tan_slope(2).*(x-crit_pt2(1))+crit_pt2(2);
            elseif (tan_slope(2) == Inf)
                y=crit_pt2(2)+1;
                x=crit_pt2(1);
            else
                x=crit_pt2(1)+1;
                y=tan_slope(2).*(x-crit_pt2(1))+crit_pt2(2);
            end
            v2=[x y]-[x_int y_int];
            to_cpt1=crit_pt1-[x_int y_int];
            test_v=v2;
            bisect_v=(norm(test_v)*to_cpt2+norm(to_cpt2)*test_v);
            %test_rot=[cos(angle_between/2) -sin(angle_between/2);sin(angle_between/2) cos(angle_between/2)];
            %rot_v=test_rot*test_v';
            test_line=vec2line(bisect_v,crit_pt2,xmin,xmax,ymin,ymax,npts);
            if (isintersect(test_line,boundary_pts))
                angle_between=acos((dot(test_v,to_cpt1)./(norm(test_v).*norm(to_cpt1))));
                angle_inc=angle_between./(ntielines+1);
                theta=rot_angle((test_v./norm(test_v)),(bisect_v./norm(bisect_v)),'asin');
                if (theta > 0)
                    rot_mat=[cos(angle_inc) -sin(angle_inc);sin(angle_inc) cos(angle_inc)];
                else
                    rot_mat=[cos(angle_inc) sin(angle_inc);-sin(angle_inc) cos(angle_inc)]';
                end
                rot_vec=v2';
            else
                test_v=-v2;
                bisect_v=(norm(test_v)*to_cpt1+norm(to_cpt1)*test_v);
                test_line=vec2line(bisect_v,crit_pt2,xmin,xmax,ymin,ymax,npts);
                if (isintersect(test_line,boundary_pts))
                    angle_between=acos((dot(test_v,to_cpt1)./(norm(test_v).*norm(to_cpt1))));
                    angle_inc=angle_between./(ntielines+1);
                    theta=rot_angle((test_v./norm(test_v)),(bisect_v./norm(bisect_v)),'asin');
                    if (theta > 0)
                        rot_mat=[cos(angle_inc) -sin(angle_inc);sin(angle_inc) cos(angle_inc)];
                    else
                        rot_mat=[cos(angle_inc) sin(angle_inc);-sin(angle_inc) cos(angle_inc)]';
                    end
                    rot_vec=-v2';
                else
                    error('no boundary intersection for crit_pt2 = [x_int y_int]')
                end
            end
        else
            to_cpt1=crit_pt1-[x_int y_int];
            to_cpt2=crit_pt2-[x_int y_int];
            angle_between=acos((dot(to_cpt1,to_cpt2)./(norm(to_cpt1).*norm(to_cpt2))));
            % divide angle into ntielines number of smaller angles
            angle_inc=angle_between./(ntielines+1);
            % find vector to rotate to other vector so that the rotated vectors span 2-phase region
            tan_angles(1)=atan(to_cpt1(2)./to_cpt1(1));
            tan_angles(2)=atan(to_cpt2(2)./to_cpt2(1));
            % circ_angles are between 0 and 2*pi
            for k=1:2
                if (tan_angles(k) > 0)
                    if (to_cpt_vs(k,2) > 0)
                        circ_angle(k)=tan_angles(k);
                    else
                        circ_angle(k)=pi+tan_angles(k);
                        %circ_angle(k)=-(pi/2)-tan_angles(k);
                    end
                elseif (tan_angles(k) < 0)
                    if (to_cpt_vs(k,2) > 0)
                        circ_angle(k)=pi+tan_angles(k);
                        %circ_angle(k)=(pi/2)-tan_angles(k);
                    else
                        circ_angle(k)=2*pi+tan_angles(k);
                        %circ_angle(k)=tan_angles(k);
                    end
                else
                    circ_angle(k)=tan_angles(k);
                end
            end
            nquad_1=whichquad(to_cpt1);
            nquad_2=whichquad(to_cpt2);
            if (nquad_2 == nquad_1)
                theta=rot_angle((to_cpt1./norm(to_cpt1)),(to_cpt2./norm(to_cpt2)),'asin');
                if (theta > 0)
                    rot_mat=[cos(angle_inc) -sin(angle_inc);sin(angle_inc) cos(angle_inc)];
                    rot_vec=to_cpt1';
                else
                    rot_mat=[cos(angle_inc) -sin(angle_inc);sin(angle_inc) cos(angle_inc)]';
                    rot_vec=to_cpt2';
                end
            end
          
                
                
            if (circ_angle(1) < pi/2 & circ_angle(2) > (3/2)*pi)
                if (2*pi-circ_angle(2)+circ_angle(1) ~= angle_between)
                    error('to_cpt1 in 1st quad,to_cpt2 in 4th quad, but angle between them not angle_between')
                end
                rot_vec=to_cpt_vs(2,:)';
                rot_mat=[cos(angle_inc) -sin(angle_inc);sin(angle_inc) cos(angle_inc)];
            elseif (circ_angle(2) < pi/2 & circ_angle(1) > (3/2)*pi)
                if (2*pi-circ_angle(1)+circ_angle(2) ~= angle_between)
                    error('to_cpt2 in 1st quad,to_cpt1 in 4th quad, but angle between them not angle_between')
                end
                rot_vec=to_cpt_vs(1,:)';
                rot_mat=[cos(angle_inc) -sin(angle_inc);sin(angle_inc) cos(angle_inc)];
            else
                if ((circ_angle(1)-circ_angle(2)) == (angle_between-4.4409e-016)) %-4.4409e-016
                    rot_vec=to_cpt_vs(2,:)';
                    rot_mat=[cos(angle_inc) -sin(angle_inc);sin(angle_inc) cos(angle_inc)]; 
                elseif ((circ_angle(2)-circ_angle(1)) == (angle_between-4.4409e-016)) %-4.4409e-016
                    rot_vec=to_cpt_vs(1,:)';
                    rot_mat=[cos(angle_inc) -sin(angle_inc);sin(angle_inc) cos(angle_inc)];
                else
                    keyboard
                    error('angle between not equal to difference in angles of to_cpt_vs')
                end
            end
        end
        
        % rotate vector either counterclockwise or clockwise so that the line of the rotated vector intersects
        % the boundary (should only intersect twice):  the first intersection point is on the alpha phase
        % boundary and the second intersection point is on the beta phase boundary
        iint=1;
        while (iint <= ntielines) %for k=1:ntielines
            %keyboard
            rot_vec=rot_mat*rot_vec;
            line_pts=vec2line(rot_vec,[x_int y_int],(xmin-10),(xmax+10),(ymin-10),(ymax+10),npts);
            %if (any(isnan(line_pts)))
                %keyboard
                %end
            [xint,yint]=polyxpoly(line_pts(:,1)',line_pts(:,2)',boundary_pts(:,1),boundary_pts(:,2),'unique');
            if (isempty(xint) & isempty(yint))
                error('rotated vector does not intersect boundary')
            end
            if (size(xint,1) == 2 & size(yint,1) == 2)
                dist_1=sqrt(((x_int-xint(1)).^2)+((y_int-yint(1)).^2));
                dist_2=sqrt(((x_int-xint(2)).^2)+((y_int-yint(2)).^2));
                if (dist_1 < dist_2)
                    alpha_pts(iint,1)=xint(1);
                    alpha_pts(iint,2)=yint(1);
                    beta_pts(iint,1)=xint(2);
                    beta_pts(iint,2)=yint(2);
                    iint=iint+1;
                else
                    alpha_pts(iint,1)=xint(2);
                    alpha_pts(iint,2)=yint(2);
                    beta_pts(iint,1)=xint(1);
                    beta_pts(iint,2)=yint(1);
                    iint=iint+1;
                end
            end  
        end
        %keyboard
        na_points=size(alpha_pts,1);
        nb_points=size(beta_pts,1);
        if (na_points ~= ntielines)
            error('number alpha points not equal to number tie lines')
        end
        if (nb_points ~= ntielines)
            error('number beta points not equal to number tie lines')
        end
        % connect each alpha pt with appropriate beta point with a line
        % this line is a tie line, and the collection of these lines is the tie line set
        % for these two critical points
        ab_pts=zeros(ntielines,bxy,bxy);
        ab_pts=cat(3,alpha_pts,beta_pts);% in 3rd dim, 1 = alpha, 2 = beta
        alpha_beta_pts{i,jset}=ab_pts; 
        clear alpha_pts;
        clear beta_pts;
        ntie=size(ab_pts,1);
        tylines=zeros(ntie,npts,2);
        for k=1:ntie
            xline=linspace(ab_pts(k,1,1),ab_pts(k,1,2),npts);
            yline=linspace(ab_pts(k,2,1),ab_pts(k,2,2),npts);
            tylines(k,:,1)=xline;
            tylines(k,:,2)=yline;
        end
        tielines{i,jset}=tylines;
        jset=jset+1;
    end % for j=1:nbpts
end