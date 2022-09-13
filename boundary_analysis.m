function [fpa,fpb] = boundary_analysis(bdy,traj,crit,kp,ntl)
% bdy_pts = an nxc matrix specifying the phase boundary
%          
% traj = a 2x3 matrix specifying trajectory end points
%        a 2x1 vector specifying the location of trajectory end points by
%        parameter b (from 0 to 1) = fraction of boundary length
% crit = a px3 matrix specifying critical points
%        a px1 vector specifying the location of critical point by
%        parameter b
% kp = partition coefficient of spin probe
%
% ntl = number of tie lines 

if all(size(bdy) > 1) & ndims(bdy) == 2 % ie if matrix
    [nbp,nbc] = size(bdy); % nbp = number of boundary pts, nc = number of components
    ncv = nbc-1; % ncv = number of compositional variables
    % make boundary convex
    bdy = make_convex_bdy(bdy);
    % declare compositional basis vectors
    % for ternary phase diagram using Xsm and Xchol
    es = [1 0]; % [x y]
    ec = [cos(pi/3) sin(pi/3)]; % [x y]
    bdy = bdy(:,1)*es+bdy(:,3)*ec;
    %b = linspace(0,1,length(bdy))';
    b = linspace(0,1,2*nbp)'; % b = [0:0.02:1]';
    [bdy_length,bdy_intervals] = bdy_fxn(bdy,'linear');
    bdy = interp1(bdy_intervals(:,1),bdy,b,'pchip');
    %b = linspace(0,1,4*nbp)'; % b = [0:0.01:1]';
    %[bdy_length,bdy_intervals] = bdy_fxn(bdy,'linear');
    %bdy = interp1(bdy_intervals(:,1),bdy,b,'spline');
else
    error('1st argument must be an nx3 matrix (Xsm Xdopc Xchol) specifying boundary points');
end

if all(size(traj) > 1) & ndims(traj) == 2 % ie if matrix
    [ntp,ntc] = size(traj); % assuming 2x3 matrix
    % transform to compositional basis
    traj_p = traj(:,1)*es+traj(:,3)*ec;
    bt = bdypt2b(b,bdy,traj_p);
elseif any(size(traj) == 1) & any(size(traj) > 1) % ie if vector
    ntp = length(traj);
    bt = traj;
    traj_p = interp1(b,bdy,bt,'pchip');
else
    error('2nd argument must be a 2x3 matrix or 2x1 vector specifying trajectory end points');
end

for i = 1:ntp-1
    if ismember(traj_p(i,:),traj_p(i+1:end,:),'rows')
        error('trajectory end points can not equal each other');
    end
end

if all(size(crit) > 1) & ndims(crit) == 2 % ie if matrix
    [ncp,ncc] = size(crit); % assuming 2x3 matrix, np = number of phases
    % transform to compositional basis
    crit_p = crit(:,1)*es+crit(:,3)*ec;
    bc = bdypt2b(b,bdy,crit_p);
elseif any(size(crit) == 1) & any(size(crit) > 1) % ie if vector
    ncp = length(crit);
    bc = crit;
    crit_p = interp1(b,bdy,bc,'pchip');
else
    error('3rd argument must be a px3 matrix or px1 vector specifying trajectory end points');
end

for i = 1:ncp-1
    if ismember(crit_p(i,:),crit_p(i+1:end,:),'rows')
        error('critical points can not equal each other');
    end
end

% plot boundary
figure;
ternary_plot(cart2tern(bdy,1),'-k','linewidth',3);
hold on;
ternary_plot(cart2tern(crit_p,1),'dm','markersize',10,'markerfacecolor','m');
ternary_plot(cart2tern(traj_p,1),'-g');

t = linspace(0,1,ntl)';
pt3 = traj_p(1,:);
pt4 = traj_p(2,:);
cpt1 = crit_p(1,:);
cpt2 = crit_p(2,:);
tpt = repmat(pt3,length(t),1)+t*(pt4-pt3);
ptacpt1 = interp1(b,bdy,[bc(1)-0.002:0.001:bc(1)+0.002]','pchip')';
mc1 = ppval(fnder(spline([bc(1)-0.002:0.001:bc(1)+0.002]',ptacpt1)),bc(1));
mc1 = mc1(2)/mc1(1);
ptacpt2 = interp1(b,bdy,[bc(2)-0.002:0.001:bc(2)+0.002]','pchip')';
mc2 = ppval(fnder(spline([bc(2)-0.002:0.001:bc(2)+0.002]',ptacpt2)),bc(2));
mc2 = mc2(2)/mc2(1);

if mc1 == Inf
    lina = [cpt1;cpt1(1) 1];
elseif mc1 > 0
    lina = [cpt1;1 mc1*(1-cpt1(1))+cpt1(2)];
elseif mc1 == 0
    lina = [cpt1;1 cpt1(2)];
elseif mc1 < 0 
    lina = [cpt1;0 mc1*(0-cpt1(1))+cpt1(2)];
end
        
if mc2 == Inf
    linb = [cpt2;cpt2(1) 1];
elseif mc2 > 0
    linb = [cpt2;1 mc2*(1-cpt2(1))+cpt2(2)];
elseif mc2 == 0
    linb = [cpt2;1 cpt2(2)];
elseif mc2 < 0 
    linb = [cpt2;0 mc2*(0-cpt2(1))+cpt2(2)];
end

cpt_m_int = line_line_int(lina,linb);

% find alpha and beta pts
if mc1 == mc2
    if mc2 == Inf;
        pt = [tpt(:,1) ones(length(t),1)];
    elseif mc2 > 0
        pt = [ones(length(t),1) mc2.*(1-tpt(:,1))+tpt(:,2)];
    elseif mc2 == 0
        pt = [ones(length(t),1) tpt(:,2)];
    elseif mc2 < 0 
        pt = [zeros(length(t),1) mc2.*(0-tpt(:,1))+tpt(:,2)];
    end
else
    pt = repmat(cpt_m_int,length(t),1);  
end

% finding fraction of phases
for c = 1:length(t)    
    % special condition
    if isequal(round(tpt(c,:).*1000),round(cpt1.*1000)) | isequal(round(tpt(c,:).*1000),round(cpt2.*1000))
        fa(c,1) = 0.5;
        fb(c,1) = 0.5;
        continue;
    end
    
    int_pts = boundary_line_int(bdy,[pt(c,:);tpt(c,:)]);
    
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
    fb(c,1) = norm(tpt(c,:)-p1)/norm(p1-p2);
    fa(c,1) = 1-fb(c,1);
    
    % plot fraction of alpha phase as red line and fraction of beta phase
    % as blue line
    ternary_plot(cart2tern([tpt(c,:);p1],1),'-b');
    ternary_plot(cart2tern([tpt(c,:);p2],1),'-r');
end

% kp defined as into beta (Lo) phase
fpa = fa./(fa+kp.*fb);
fpb = (kp.*fb)./(fa+kp.*fb);
hold off;

return