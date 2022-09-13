function [fa,fb] = circle_analysis(traj,crit_pts,ntl,kp)
% crit_pts = vector giving angles of critical points
%        or = 2 x 2 matrix specifying critical point locations on circle
%
% traj = 2 x 2 matrix specifying end points of trajectory
%   or = vector giving angles of end points
%
% ntl = number of tie lines

warning off;

if all(size(traj) == 2) & ndims(traj) == 2 % ie if matrix
    traj_apt = traj(1,:);
    ang_traja = pt2angle(traj(1,:));
    traj_bpt = traj(2,:);
    ang_trajb = pt2angle(traj(2,:));
elseif any(size(traj) == 1) & ~all(size(traj) == 1) % ie if vector
    traj_apt = angle2pt(traj(1));
    ang_traja = traj(1);
    traj_bpt = angle2pt(traj(2));
    ang_travb = traj(2);
else
    error('1st argument must be a matrix specifying trajectory end points');
end

if all(size(crit_pts) == 2) & ndims(crit_pts) == 2 % ie if matrix
    cpt1 = crit_pts(1,:);
    ang_cpt1 = pt2angle(crit_pts(1,:));
    cpt2 = crit_pts(2,:);
    ang_cpt2 = pt2angle(crit_pts(2,:));
elseif any(size(crit_pts) == 1) & ~all(size(crit_pts) == 1) % ie if vector
    cpt1 = angle2pt(crit_pts(1));
    ang_cpt1 = crit_pts(1);
    cpt2 = angle2pt(crit_pts(2));
    ang_cpt2 = crit_pts(2);
else
    error('2nd argument must be scalar or matrix specifying tie lines');
end

if cpt1 == cpt2
    error('critical point 1 can not equal critical point 2');
end

ft = linspace(0,1,ntl)'; % parmater ft specifies fraction of trajectory traversed
% parameter x specifies fraction of arc length from critical point 1 
% to critical point 2, which is the alpha phase boundary.
% according to my model p1 = alpha point on tie line, p2 = beta point on
% tie line, p3 = alpha point of the trajectory, and p4 = beta point of the
% trajectory.
% p1 x p2 = sin(theta12), where theta12 = angle between p1 and p2
% also p1 x p2 = (L12)^2 * sqrt((1/L12)^2 - 1/4), where L12 is the length
% of the chord between p1 and p2

% create circle boundary
angles = linspace(0,2*pi)';
circle_pts = angle2pt(angles);
figure;
hold on;
plot(circle_pts(:,1),circle_pts(:,2),'linewidth',5,'color','k');
plot(cpt1(1),cpt1(2),'dm','markersize',10,'markerfacecolor','m')
plot(cpt2(1),cpt2(2),'dm','markersize',10,'markerfacecolor','m')

% find tie lines
% parallel tie lines: find tangent line to critical point, use slope
%   tangent line = [-sin(angle) cos(angle)], slope = cos(angle)/-sin(angle)
%   tangent line in point slope form: y - y(critpt) = slope*(x - x(critpt))

% find fraction of alpha and beta phases
p3 = traj_apt;
p4 = traj_bpt;
p = repmat(p3,length(ft),1)+ft*(p4-p3);
line([p3(1) p4(1)]',[p3(2) p4(2)]','color','g');
mc1 = angle2tan_slope(ang_cpt1);
mc2 = angle2tan_slope(ang_cpt2);

if mc1 == Inf
    lina = [cpt1;cpt1(1) 2];
elseif mc1 > 0
    lina = [cpt1;2 mc1*(2-cpt1(1))+cpt1(2)];
elseif mc1 == 0
    lina = [cpt1;2 cpt1(2)];
elseif mc1 < 0 
    lina = [cpt1;-2 mc1*(-2-cpt1(1))+cpt1(2)];
end
        
if mc2 == Inf
    linb = [cpt2;cpt2(1) 2];
elseif mc2 > 0
    linb = [cpt2;2 mc2*(2-cpt2(1))+cpt2(2)];
elseif mc2 == 0
    linb = [cpt2;2 cpt2(2)];
elseif mc2 < 0 
    linb = [cpt2;-2 mc2*(-2-cpt2(1))+cpt2(2)];
end

% finding fraction of phases
for t = 1:length(ft)
    x = p(t,1);
    y = p(t,2);
    
    % special condition
    if isequal(round([x y].*1000),round(cpt1.*1000)) | isequal(round([x y].*1000),round(cpt2.*1000))
        fa(t,1) = 0.5;
        fb(t,1) = 0.5;
        continue;
    end
    
    % find alpha and beta pts
    if mc1 == mc2
        if mc2 == Inf;
            pt = [x 2];
        elseif mc2 > 0
            pt = [2 mc2*(2-x)+y];
        elseif mc2 == 0
            pt = [2 y];
        elseif mc2 < 0 
            pt = [-2 mc2*(-2-x)+y];
        end
    else
        pt = line_line_int(lina,linb);  
    end
    
    % [p1,p2] = circle_line_int([pt;x y]); % p1 = alpha point and p2 = beta point
    int_pts = circle_line_int([pt;x y]);
    
    if isempty(int_pts) % no intersection
        error('tie line does not intersect boundary');
    elseif size(int_pts,1) == 1 % tangent intersection
        p1 = int_pts;
        p2 = pt;
    elseif size(int_pts,1) == 2 % 2 intersections
        p1 = int_pts(1,:);
        p2 = int_pts(2,:);
    end
    
    % calculate fraction of alpha phase and fraction of beta phase
    fb(t,1) = norm([x y]-p1)/norm(p1-p2);
    fa(t,1) = 1-fb(t,1);
    
    % plot fraction of alpha phase as red line and fraction of beta phase
    % as blue line
    line([x p1(1)]',[y p1(2)]','color','r');
    line([x p2(1)]',[y p2(2)]','color','b');
end

fa = fa./(fa+kp.*fb);
fb = (kp.*fb)./(fa+kp.*fb);

hold off;
warning on;

return