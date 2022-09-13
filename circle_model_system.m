function varargout = circle_model_system(boundary,nnodes,varargin)
% circle_model_system is a function to analyze thermodynamic tie line
% structure of a phase boundary modeled as the unit circle.

% for a 2-phase system, there are only three possibilities:
% 1)  2 critical points:  2-phase region surrounded by 1-phase region.
% 2)  1 critical point/1 ending tie line:  1-phase region near critical point
% and a 3-phase invariant region near ending tie line.
% 3)  2 ending tie lines:  3-phase regions near both ending tie lines.

% argument boundary is a 1x2 input structure containing 1 of the above 3 possibilites
% for the phase boundary with fields:
% type = either 'pt' or 'line'.
% coords = if type == 'pt', then coords = scalar angle.
%          if type == 'line', then coords = vector of angles.
%
% (angle between -pi:pi specifies a unique point on the circle).

% checking input arguments
if isstruct(boundary)
    if size(boundary) ~= 2
        error('boundary input argument must be a 1x2 structure array');
    else
        if strcmp(boundary(1).type,'pt') & strcmp(boundary(2).type,'pt')
            if isequal(boundary(1).coords,boundary(2).coords)
                error('critical points cannot be the same point');
            end
        elseif strcmp(boundary(1).type,'pt') & strcmp(boundary(2).type,'line')
            if any(boundary(1).coords == boundary(2).coords)
                error('a critical point cannot lie on an ending tie line');
            end
        elseif strcmp(boundary(1).type,'line') & strcmp(boundary(2).type,'pt')
            if any(boundary(1).coords == boundary(2).coords)
                error('a critical point cannot lie on an ending tie line');
            end
        elseif strcmp(boundary(2).type,'line') & strcmp(boundary(1).type,'line')
            if any(boundary(1).coords(1) == boundary(2).coords) | any(boundary(1).coords(2) == boundary(2).coords)
                error('ending tie lines cannot share any points');
            end
        end
    end
else
    error('boundary input argument must be a 1x2 structure array');
end

% create circle phase boundary
circle_center = cart2tcart([0.5 0.3]);
circle_radius = 0.2;
circle_angles = linspace(-pi,pi)'; % default = 100 angles
circle_pts = angle2pt(circle_angles,circle_radius,circle_center);
figure;
hold on;
ternary_plot(cart2tern(circle_pts,1),'-k','linewidth',3);
% plot(circle_pts(:,1),circle_pts(:,2),'-k','linewidth',4);
if strcmp(boundary(1).type,'pt')
    cpt1 = angle2pt(boundary(1).coords,circle_radius,circle_center);
%     ternary_plot(cart2tern(cpt1,1),'dm','markersize',10,'markerfacecolor','m');
%     plot(cpt1(1),cpt1(2),'dm','markersize',10,'markerfacecolor','m');
elseif strcmp(boundary(1).type,'line')
    etl1a = angle2pt(boundary(1).coords(1),circle_radius,circle_center);
    etl1b = angle2pt(boundary(1).coords(2),circle_radius,circle_center);
%     ternary_plot(cart2tern([etl1a;etl1b],1),'-md','markersize',5,'markerfacecolor','m','linewidth',2);
%     plot([etl1a(1) etl1b(1)],[etl1a(2) etl1b(2)],'-md','markersize',5,'markerfacecolor','m','linewidth',2);
end
if strcmp(boundary(2).type,'pt')
    cpt2 = angle2pt(boundary(2).coords,circle_radius,circle_center);
%     ternary_plot(cart2tern(cpt2,1),'dm','markersize',10,'markerfacecolor','m');
%     plot(cpt2(1),cpt2(2),'dm','markersize',10,'markerfacecolor','m');
elseif strcmp(boundary(2).type,'line')
    etl2a = angle2pt(boundary(2).coords(1),circle_radius,circle_center);
    etl2b = angle2pt(boundary(2).coords(2),circle_radius,circle_center);
%     ternary_plot(cart2tern([etl2a;etl2b]),'-md','markersize',5,'markerfacecolor','m','linewidth',2);
%     plot([etl2a(1) etl2b(1)],[etl2a(2) etl2b(2)],'-md','markersize',5,'markerfacecolor','m','linewidth',2);
end

% specify data points 
x = repmat([0:0.05:1]',length([0:0.05:1]),1);
y = [];
for i = [0:0.05:1]
    y = cat(1,y,repmat(i,length([0:0.05:1]),1));
end
data_pts = cart2tcart([x y]);
i = find((data_pts(:,1)-circle_center(1)).^2 + (data_pts(:,2)-circle_center(2)).^2 < circle_radius^2);
data_pts = cart2tcart([x(i) y(i)]);
clear x y;
% x = data_pts(:,1);
% y = data_pts(:,2);
% data_pts = data_pts(1:2:end,:);
% data_pts = data_pts(1:2:end,:);
ternary_plot(cart2tern(data_pts,1),'.c','markersize',15);

% tie lines are specified using a curve from one point on the boundary or
% ending tie line to another point on the boundary or ending tie line.
% this curve is divided into a certain (random) number of nodal points
% where at each of these nodal points a slope is randomly chosen but with
% no intersection constraints.  this nodal point and slope define a line and
% the nodal point is chosen as the midpoint of the line segment where the
% end points of the line segment is where the line intersects the
% boundary.

if strcmp(boundary(1).type,'pt')
    if strcmp(boundary(2).type,'pt')
        mc1 = angle2tan_slope(boundary(1).coords);
        mc2 = angle2tan_slope(boundary(2).coords);
        if mc1 == Inf
            lina = [cpt1;cpt1(1) 1];
		elseif mc1 > 0
            lina = [cpt1;1 mc1*(1-cpt1(1))+cpt1(2)];
		elseif mc1 == 0
            lina = [cpt1;1 cpt1(2)];
		elseif mc1 < 0 
            lina = [cpt1;-1 mc1*(-1-cpt1(1))+cpt1(2)];
		end       
		if mc2 == Inf
            linb = [cpt2;cpt2(1) 1];
		elseif mc2 > 0
            linb = [cpt2;1 mc2*(1-cpt2(1))+cpt2(2)];
		elseif mc2 == 0
            linb = [cpt2;1 cpt2(2)];
		elseif mc2 < 0 
            linb = [cpt2;-1 mc2*(-1-cpt2(1))+cpt2(2)];
		end
    elseif strcmp(boundary(2).type,'line')
        mc1 = angle2tan_slope(boundary(1).coords);
        if mc1 == Inf
            lina = [cpt1;cpt1(1) 1];
		elseif mc1 > 0
            lina = [cpt1;1 mc1*(1-cpt1(1))+cpt1(2)];
		elseif mc1 == 0
            lina = [cpt1;1 cpt1(2)];
		elseif mc1 < 0 
            lina = [cpt1;-1 mc1*(-1-cpt1(1))+cpt1(2)];
		end
        linb = [etl2a;etl2b];
    end   
elseif strcmp(boundary(1).type,'line')
    if strcmp(boundary(2).type,'pt')
        lina = [etl1a;etl1b];
        mc2 = angle2tan_slope(boundary(2).coords);
        if mc2 == Inf
            linb = [cpt2;cpt2(1) 1];
		elseif mc2 > 0
            linb = [cpt2;1 mc2*(1-cpt2(1))+cpt2(2)];
		elseif mc2 == 0
            linb = [cpt2;1 cpt2(2)];
		elseif mc2 < 0 
            linb = [cpt2;-1 mc2*(-1-cpt2(1))+cpt2(2)];
		end
    elseif strcmp(boundary(2).type,'line')
        lina = [etl1a;etl1b];
        linb = [etl2a;etl2b];
    end   
end
ma = angle2slope(pt2angle(diff(lina)));
mb = angle2slope(pt2angle(diff(linb)));  

if isequal(nnodes,0)
    % uniform tie lines, slopes of tie lines changes monotonically
    if round(ma*1000)/1000 == round(mb*1000)/1000 % no intersection, lines parallel
        linform = 'pt-m';
        m = repmat(ma,ndp,1);
    else
        pt = line_line_int(lina,linb);
        linform = '2pts';  
	end
	% loop over data points and calculate properties
	for d = 1:ndp
        % [p1,p2] = circle_line_int([pt;x y]); % p1 = alpha point and p2 = beta point
        switch linform
            case 'pt-m'
                int_pts = circle_line_int(circle_center,circle_radius,linform,m(d),data_pts(d,:));
            case '2pts'
                int_pts = circle_line_int(circle_center,circle_radius,linform,[pt;data_pts(d,:)]);
            otherwise
                error('invalid line form');
        end
        
        if isempty(int_pts) % no intersection
            error('tie line does not intersect boundary');
        elseif size(int_pts,1) == 1 % tangent intersection
            p1 = int_pts;
            p2 = pt;
        elseif size(int_pts,1) == 2 % 2 intersections
            apt(d,:) = int_pts(1,:); % alpha pt
            bpt(d,:) = int_pts(2,:); % beta pt
        end
	end
elseif nnodes > 0
    % choose random point on boundary
    ra = -pi + (2*pi)*rand(1);
    rpt = angle2pt(ra,circle_radius,circle_center);
    rm = angle2tan_slope(ra);
    nodept(1,:) = rpt;
    napt(1,:) = rpt;
    nbpt(1,:) = [rpt(1)+0.01 rm*(0.01)+rpt(2)];
    nodem(1) = rm;
    
    % travel a random distance along vector normal to the tangent line at
    % the random point.  the normal vector should point inward.
    % this new point is first node.
    int_pts = circle_line_int(circle_center,circle_radius,'pt-m',-1/rm,rpt);
    if size(int_pts,1) == 2 % 2 intersections
        if round(int_pts(1,:).*1000)./1000 == round(rpt.*1000)./1000 % int_pts(1,:) == rpt 
            pt = int_pts(2,:);    
        else
            pt = int_pts(1,:);
        end
        nv(1,:) = pt-rpt;
        d2pt = norm(nv(1,:));
        rd = rand(1)*d2pt;
        nodept(2,:) = nodept(1,:) + rd*(nv(1,:)/d2pt);
    else
        error('no intersection for first node')
    end 
    
    % choose a random slope so that the tie line does not coincide with
    % previous normal vector. 
    anv = pt2angle(nv(1,:));
    rnum = rand(1);
    while rnum == 0 | rnum == 1
        rnum = rand(1);
    end
    nodem(2) = angle2slope(anv + rnum*(pi-anv));
    
    % find alpha and beta pts for the node using the node point and node
    % slope
    int_pts = circle_line_int(circle_center,circle_radius,'pt-m',nodem(2),nodept(2,:));
	if size(int_pts,1) == 2
        if pt2angle(diff(int_pts)) >= 0
            napt(2,:) = int_pts(1,:); % nodal alpha pt
            nbpt(2,:) = int_pts(2,:); % nodal beta pt
        else
            napt(2,:) = int_pts(2,:); % nodal alpha pt
            nbpt(2,:) = int_pts(1,:); % nodal beta pt
        end
    else
        error('no intersection from node')
	end
    
    % debugging
    ternary_plot(cart2tern(nodept,1),'.m','markersize',20);
    ternary_plot(cart2tern(napt(2,:),1),'.r','markersize',20);
    ternary_plot(cart2tern(nbpt(2,:),1),'.b','markersize',20);
    line([napt(2,1);nbpt(2,1)],[napt(2,2);nbpt(2,2)],'color','m');
    
    % loop over data points, find those between node 1 and 2, and
    % determine their alpha and beta points
    ndp = size(data_pts,1);
    done = zeros(ndp,1);
	for d = 1:ndp 
        d2n_1 = pt2lin_dist(nodept(1,:),[nbpt(2,:);napt(2,:)],2);
        d2n = pt2lin_dist(data_pts(d,:),[nbpt(2,:);napt(2,:)],2);
        if isequal(round(d2n*1000)./1000,0) % isequal(d2n,0) % data point lies on nodal line
            apt(d,:) = napt(2,:);
            bpt(d,:) = nbpt(2,:);
        else 
            if sign(d2n) == sign(d2n_1)
                nlin1 = [napt(1,:);nbpt(1,:)];
                nlin2 = [napt(2,:);nbpt(2,:)];
                pt = line_line_int(nlin1,nlin2);
                int_pts = circle_line_int(circle_center,circle_radius,'2pts',[pt;data_pts(d,:)]);
			    if size(int_pts,1) == 2
                    apt(d,:) = int_pts(1,:); % alpha pt
                    bpt(d,:) = int_pts(2,:); % beta pt
                    % debugging
%                     line([apt(d,1);bpt(d,1)],[apt(d,2);bpt(d,2)],'color','m');
                else
                    error('no intersection from node')
			    end
                done(d) = 1;
            end
        end
	end
    % find data points not done
    inotdone = find(~done);
    
    n = 2;
    while n < nnodes+1
		% again travel a random distance along new normal vector of new tie
		% line. this is next node.
		int_pts = circle_line_int(circle_center,circle_radius,'pt-m',-1/nodem(n-1),nodept(n-1,:));
		if size(int_pts,1) == 2 % 2 intersections
            pt1 = int_pts(1,:);
            pt2 = int_pts(2,:);
            nv1 = pt1-nodept(1,:);
            nv2 = pt2-nodept(1,:);
            one = dot(nv(n-1,:),nv1);
            two = dot(nv(n-1,:),nv2);
            if sign(one) > 0 & sign(two) < 0
                nv(n,:) = nv1;
            elseif sign(one) < 0 & sign(two) > 0
                nv(n,:) = nv2;
            else
                error('normal vectors of node point in same direction');
            end
            d2pt = norm(nv(n,:));
            nd = rand(1)*d2pt;
            nodept(n,:) = nodept(n-1,:) + nd*(nv(n,:)/d2pt);
		else
            error('no intersection from normal vector')
		end 
        
		% again choose random slope but this time making sure tie line does
		% not intersect previous tie line.
        % R clockwise rotation = [cos(theta) sin(theta);-sin(theta) cos(theta)]
        theta = pt2angle(nv(n,:));
        av = napt(n-1,:)-nodept(n,:);
        bv = nbpt(n-1,:)-nodept(n,:);
        % rotate frame
        rav = [cos(theta) sin(theta);-sin(theta) cos(theta)]*av';
        rbv = [cos(theta) sin(theta);-sin(theta) cos(theta)]*bv';
		arav = pt2angle(rav);
        arbv = pt2angle(rbv);
        if arav > 0 & arbv < 0
            arbv = arbv+pi;
            if arav > arbv
                ra = arbv + rand(1)*(arav-arbv);
            else
                ra = arav + rand(1)*(arbv-arav);
            end
        elseif arav < 0 & arbv > 0
            arav = arav+pi;
            if arav > arbv
                ra = arbv + rand(1)*(arav-arbv);
            else
                ra = arav + rand(1)*(arbv-arav);
            end
        else
            error('alpha pt and beta pt on same side of horizontal');
        end
        % rotate back to original frame
        nodem(n) = angle2slope(pt2angle([cos(-theta) sin(-theta);-sin(-theta) cos(-theta)]*angle2pt(ra)'));
        
        % find alpha and beta pts for the node using the node point and node
        % slope
        int_pts = circle_line_int(circle_center,circle_radius,'pt-m',nodem(n),nodept(n,:));
		if size(int_pts,1) == 2 
            if pt2angle(diff(int_pts)) >= 0
                napt(2,:) = int_pts(1,:); % nodal alpha pt
                nbpt(2,:) = int_pts(2,:); % nodal beta pt
            else
                napt(2,:) = int_pts(2,:); % nodal alpha pt
                nbpt(2,:) = int_pts(1,:); % nodal beta pt
            end
        else
            error('no intersection from node')
		end
        % loop over data points, find those between node n-1 and n, and
        % determine their alpha and beta points
        ndp = size(inotdone,1);
		for d = 1:ndp
            if pt2angle(nbpt(n,:)-napt(n,:)) >= 0
                d2n_1 = pt2lin_dist(nodept(n-1,:),[nbpt(n,:);napt(n,:)],2);
                d2n = pt2lin_dist(data_pts(inotdone(d),:),[nbpt(n,:);napt(n,:)],2);
            else
                d2n_1 = pt2lin_dist(nodept(n-1,:),[napt(n,:);nbpt(n,:)],2);
                d2n = pt2lin_dist(data_pts(inotdone(d),:),[napt(n,:);nbpt(n,:)],2);
            end
            if isequal(d2n,0) % data point lies on nodal line
                apt(inotdone(d),:) = napt(n,:);
                bpt(inotdone(d),:) = nbpt(n,:);
            else 
                if sign(d2n) == sign(d2n_1)
                    nlin1 = [napt(n-1,:);nbpt(n-1,:)];
                    nlin2 = [napt(n,:);nbpt(n,:)];
                    pt = line_line_int(nlin1,nlin2);
                    int_pts = circle_line_int(circle_center,circle_radius,'2pts',[pt;data_pts(inotdone(d),:)]);
				    if size(int_pts,1) == 2 
                        if pt2angle(diff(int_pts)) >= 0
                            apt(inotdone(d),:) = int_pts(1,:); % alpha pt
                            bpt(inotdone(d),:) = int_pts(2,:); % beta pt
                        else
                            apt(inotdone(d),:) = int_pts(2,:); % alpha pt
                            bpt(inotdone(d),:) = int_pts(1,:); % beta pt  
                        end 
                    else
                        error('no intersection from node')
				    end
                    done(inotdone(d)) = 1;
                end
            end
		end
        % find data points not done
        inotdone = find(~done);
        
		% repeat above procedure until on all nodes are placed
        n = n +1;
    end
    % extend normal vector until intersects boundary.
    % this intersection is ending node.
    int_pts = circle_line_int(circle_center,circle_radius,'pt-m',-1/nodem(end),nodept(end,:));
	if size(int_pts,1) == 2 % 2 intersections
        pt1 = int_pts(1,:);
        pt2 = int_pts(2,:);
        nv1 = pt1-nodept(end,:);
        nv2 = pt2-nodept(end,:);
        one = dot(nv(end,:),nv1);
        two = dot(nv(end,:),nv2);
        if sign(one) > 0 & sign(two) < 0
            nodept(end+1,:) = pt1;   
        elseif sign(one) < 0 & sign(two) > 0
            nodept(end+1,:) = pt2;
        else
            error('normal vectors of end-1 node point in same direction');
        end
        nodem(end+1) = angle2tan_slope(pt2angle(nodept(end,:),circle_center));
        napt(end+1,:) = nodept(end,:);
        nbpt(end+1,:) = [nodept(end,1)+0.01 nodem(end)*(0.01)+nodept(end,2)];
	else
        error('no intersection for end node')
	end
    ternary_plot(cart2tern(nodept,1),'.m','markersize',20);
    
    % loop over data points, find those between node end-1 and end, and
    % determine their alpha and beta points
%     ndp = size(inotdone,1);
	for d = inotdone' 
        d2n_1 = pt2lin_dist(nodept(end,:),[nbpt(end-1,:);napt(end-1,:)],2);
        d2n = pt2lin_dist(data_pts(d,:),[nbpt(end-1,:);napt(end-1,:)],2);
        if isequal(round(d2n*1000)/1000,0) % data point lies on nodal line
            apt(d,:) = napt(end-1,:);
            bpt(d,:) = nbpt(end-1,:);
        else 
            if sign(d2n) == sign(d2n_1)
                nlin1 = [napt(end-1,:);nbpt(end-1,:)];
                nlin2 = [napt(end,:);nbpt(end,:)];
                pt = line_line_int(nlin1,nlin2);
                int_pts = circle_line_int(circle_center,circle_radius,'2pts',[pt;data_pts(d,:)]);
			    if size(int_pts,1) == 2 
                    apt(d,:) = int_pts(1,:); % alpha pt
                    bpt(d,:) = int_pts(2,:); % beta pt
                    % debugging
%                     line([apt(inotdone(d),1);bpt(inotdone(d),1)],[apt(inotdone(d),2);bpt(inotdone(d),2)],'color','m');
                else
                    error('no intersection from node')
			    end
                done(d) = 1;
            end
        end
	end
    inotdone = find(~done);
else
    error('number of tie line nodes must be 0 or positive scalar');
end

% up until now the alpha points and beta points may or may not be
% consistently on one side of the phase boundary or not because the ending
% nodal point was unknown.  the following is to put all the alpha points on
% one side of the phase boundary and the beta points on the other side.
ndp = size(data_pts,1);
beginning = pt2angle(nodept(1,:),circle_center);
ending =pt2angle(nodept(end,:),circle_center);
aapt = pt2angle(apt,circle_center);
abpt = pt2angle(bpt,circle_center);
for d = 1:ndp
    if beginning > ending
        if aapt(d) < beginning & aapt(d) > ending
            temp = bpt(d,:);
            bpt(d,:) = apt(d,:);
            apt(d,:) = temp;
        end
    else
        if aapt(d) < ending & aapt(d) > beginning
            temp = bpt(d,:);
            bpt(d,:) = apt(d,:);
            apt(d,:) = temp;
        end
    end
end
    
% loop over data_pts and calculate the fraction of phase alpha and beta,
% the kp composition function, and the fraction of probe in each phase
for d = 1:ndp
	fb(d,1) = norm(data_pts(d,:)-apt(d,:))/norm(apt(d,:)-bpt(d,:));
	fa(d,1) = 1-fb(d,1);
	xas = apt(d,1);
	xac = apt(d,2);
	xbs = bpt(d,1);
	xbc = bpt(d,2);
	kp(d,1) = (xas + xac)/(xbs + xbc);
	fpb(d,1) = fb(d,1)/(fb(d,1) + kp(d,1)*fa(d,1));
	fpa(d,1) = (kp(d,1)*fa(d,1))/(fb(d,1) + kp(d,1)*fa(d,1));
	% plot tie line
	line([xas;xbs],[xac;xbc],'color','m');
end
% plot alpha and beta points
ternary_plot(cart2tern(apt,1),'.r','markersize',20);
ternary_plot(cart2tern(bpt,1),'.b','markersize',20);

% plot properties (ie, fraction of phase, fraction of probe in
% phases, kp, etc)
x = data_pts(:,1);
y = data_pts(:,2);
i = [0:0.05:1]; 
[xi,yi] = meshgrid(i,i);
zifa = griddata(x,y,fa,xi,yi);
% figure;
% [dx,dy] = gradient(zifa);
% contour(xi,yi,zifa,10);
% hold on
% quiver(xi,yi,dx,dy,2);
% hold off
figure;
ternary_plot(cart2tern(circle_pts,1),'-k','linewidth',3);
hold on
surf(xi,yi,zifa);
hold off
figure;
ternary_plot(cart2tern(circle_pts,1),'-k','linewidth',3);
zikp = griddata(x,y,kp,xi,yi);
hold on
surf(xi,yi,zikp);
hold off
figure;
ternary_plot(cart2tern(circle_pts,1),'-k','linewidth',3);
zifpa = griddata(x,y,fpa,xi,yi);
hold on
surf(xi,yi,zifpa);
hold off
% figure;
% ternary_plot(cart2tern(circle_pts,1),'-k','linewidth',3);
% hold on
% contour(xi,yi,zikp,10);
% hold off

return