function bdy = make_convex_bdy(bdy)
% converts a simple polygon into a convex polygon

if all(size(bdy) > 1) & ndims(bdy) == 2 % ie if matrix
    [nbp,nbc] = size(bdy); % nbp = number of boundary pts, nc = number of components
    % declare compositional basis vectors
    % for ternary phase diagram using Xsm and Xchol
    es = [1 0]; % [x y]
    ec = [cos(pi/3) sin(pi/3)]; % [x y]
    bdy = bdy(:,1)*es+bdy(:,3)*ec;
else
    error('1st argument must be an nx3 matrix (Xsm Xdopc Xchol) specifying boundary points');
end

% rotation matrix R
R = [0 1; -1 0];
notdone = true;

while notdone
    convex = true;
    
    for i = 1:nbp-2
        vi = bdy(i+1,:)-bdy(i,:);
        vi1 = bdy(i+2,:)-bdy(i+1,:);
        x = dot(R*vi',vi1); % perpendicular dot product, same as 2D cross product       
        
        if x >= 0 % ie vi1 points in direction of vi-perp (R*vi')
            % transform bdy(i+1,:)
            % ui = bdy(i,:)-bdy(i+1,:);
            % ui1 = bdy(i+2,:)-bdy(i+1,:);
            % a = 0.75*norm(vi1-vi);
            % bisect_v = norm(ui)*ui1+norm(ui1)*ui;
            % bisect_v = bisect_v./norm(bisect_v);
            % bdy(i+1,:) = bdy(i+1,:)+a*bisect_v;
            bdy(i+1,:) = bdy(i+1,:)+0.75*(vi1-vi);
            convex = false;
        end
    end
    
    if convex
        notdone = false;
    end
end

bdy = cart2tern(bdy,1);

return