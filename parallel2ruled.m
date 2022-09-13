function [u,directrixO,v,directrixA] = parallel2ruled(x,C_bdy,coexpt,boundary)
% [u,directrixO] = parallel2ruled(slope,bdy,coexpt,bdypts,bdyb)
% converts parallel tieline configuration to a ruled tieline configuration
% by finding the two critical points for a specific slope of tielines
% and then determining the two directrices.

switch boundary.config
    case '2critpts'
        angle = pi.*x(1);
        slope = angle2slope(angle);
        b = boundary.parameterization;
        int_pts = boundary_line_int(C_bdy,slope,coexpt);% doesn't matter which coexistence point
        nipt = size(int_pts,1);
        if isempty(int_pts) % no intersection
            error('tie line does not intersect boundary');
        elseif isequal(nipt,1) % tangent intersection
            error('tie line intersects boundary at only 1 point');
        elseif isequal(nipt,2) % 2 intersections
            aC = int_pts(1,:);
            bC = int_pts(2,:);
        elseif nipt > 2 
            error('tie line intersects boundary at > 2 points');
        end
        ab = bdypt2b(aC,b,C_bdy);
        bb = bdypt2b(bC,b,C_bdy);
        
        [s1,section1,s2,section2] = divide_bdy(ab,bb,C_bdy);
        options = optimset('Display','off','NonlEqnAlgorithm','gn');
        if ab < bb
            bi = s2;
            bdyi = section2;
            i2 = fzero(@eqn2zero,[0 1],options);
            cpt2 = b2bdypt(bi,bdyi,i2);
            bi = s1;
            bdyi = section1;
            i1 = fzero(@eqn2zero,[0 1],options);
            cpt1 = b2bdypt(bi,bdyi,i1);
        elseif ab > bb
            bi = s1;
            bdyi = section1;
            i1 = fzero(@eqn2zero,[0 1],options);
            cpt1 = b2bdypt(bi,bdyi,i1);
            bi = s2;
            bdyi = section2;
            i2 = fzero(@eqn2zero,[0 1],options);
            cpt2 = b2bdypt(bi,bdyi,i2);
        else
            error('ab cannot equal bb');
        end

        cpt1b = bdypt2b(cpt1,b,C_bdy);
        cpt2b = bdypt2b(cpt2,b,C_bdy);
        [s1,section1,s2,section2] = divide_bdy(cpt1b,cpt2b,C_bdy);
        if cpt1b <= cpt2b
            u = s1;
            directrixO = section1;
            v = s2;
            directrixA = section2;
        else
            u = s2;
            directrixO = section2;
            v = s1;
            directrixA = section1;
        end
    case '1critpt/1endtl'
        b = boundary.parameterization;
        e1 = x(1);
        e2 = x(2);
        % only one configuration
        ept1 = b2bdypt(b,C_bdy,e1);
        ept2 = b2bdypt(b,C_bdy,e2);
        slope = diff([ept1;ept2]);
        slope = slope(2)/slope(1);
        
        [s1,section1,s2,section2] = divide_bdy(e1,e2,C_bdy);% x(1),x(2) = e1,e2
        [s1length,temp1] = chord_length(section1);
        [s2length,temp2] = chord_length(section2);
        clear temp1 temp2;
        options = optimset('Display','off','NonlEqnAlgorithm','gn');
        if s1length >= s2length
            pt0 = section1(1,:);
            pt1 = section1(end,:);
            int_pts = boundary_line_int(C_bdy,slope,pt0); % or pt1
            nipt = size(int_pts,1);
            if isempty(int_pts) % no intersection
                error('line does not intersect boundary');
            elseif isequal(nipt,1) % tangent intersection
                ptA = int_pts;
                ptB = int_pts;
            elseif isequal(nipt,2) % 2 intersections
                ptA = int_pts(1,:);
                ptB = int_pts(2,:);
            elseif nipt > 2 
                error('line intersects boundary at > 2 points');
            end
            
            if (all(pt0 == ptA) && all(pt1 == ptB)) || (all(pt0 == ptB) && all(pt1 == ptA))
                bi = s1;
                bdyi = section1;
            else
                bdyi = section1(2:end-1,:);
                [bdyilength,bi] = chord_length(bdyi);
                clear bdyilength;
            end
            
            cpbi = fzero(@eqn2zero,[0 1],options);
            [cpti,braket] = b2bdypt(bi,bdyi,cpbi);
            [cps1,braket] = bdypt2b(cpti,s1,section1);
            [cpt,braket] = b2bdypt(s1,section1,cps1);
            bracket = braket{:};
            if isscalar(bracket)
                i = bracket;
                directrixA = flipud(section1(1:i,:));
                directrixO = section1(i:end,:);
            elseif isvector(bracket)
                il = bracket(1);
                iu = bracket(2);
                directrixA = [cpt;flipud(section1(1:il,:))];
                directrixO = [cpt;section1(iu:end,:)];
            else
                error('bracket cannot be a matrix');
            end
        else
            pt0 = section2(1,:);
            pt1 = section2(end,:);
            int_pts = boundary_line_int(C_bdy,slope,pt0); % or pt1
            nipt = size(int_pts,1);
            if isempty(int_pts) % no intersection
                error('line does not intersect boundary');
            elseif isequal(nipt,1) % tangent intersection
                ptA = int_pts;
                ptB = int_pts;
            elseif isequal(nipt,2) % 2 intersections
                ptA = int_pts(1,:);
                ptB = int_pts(2,:);
            elseif nipt > 2 
                error('line intersects boundary at > 2 points');
            end
            
            if (all(pt0 == ptA) && all(pt1 == ptB)) || (all(pt0 == ptB) && all(pt1 == ptA))
                bi = s2;
                bdyi = section2;
            else
                bdyi = section2(2:end-1,:);
                [bdyilength,bi] = chord_length(bdyi);
                clear bdyilength;
            end
            
            cpbi = fzero(@eqn2zero,[0 1],options);
            [cpti,braket] = b2bdypt(bi,bdyi,cpbi);
            [cps2,braket] = bdypt2b(cpti,s2,section2);
            [cpt,braket] = b2bdypt(s2,section2,cps2);
            bracket = braket{:};
            if isscalar(bracket)
                i = bracket;
                directrixA = flipud(section2(1:i,:));
                directrixO = section2(i:end,:);
            elseif isvector(bracket)
                il = bracket(1);
                iu = bracket(2);
                directrixA = [cpt;flipud(section2(1:il,:))];
                directrixO = [cpt;section2(iu:end,:)];
            else
                error('bracket cannot be a matrix');
            end
        end
        [Olength,u] = chord_length(directrixO);
        [Alength,v] = chord_length(directrixA);
    case '2endtls'
    otherwise
        error('invalid boundary configuration');
end

    function y = eqn2zero(x)
        pt = b2bdypt(bi,bdyi,x);
        int_pts = boundary_line_int(C_bdy,slope,pt);
        nipt = size(int_pts,1);
        if isempty(int_pts) % no intersection
            error('line does not intersect boundary');
        elseif isequal(nipt,1) % tangent intersection
            y = 0;
            return;
        elseif isequal(nipt,2) % 2 intersections
            ptA = int_pts(1,:);
            ptB = int_pts(2,:);
        elseif nipt > 2 
            error('line intersects boundary at > 2 points');
        end
        if norm(pt - ptA) < norm(pt - ptB)
            y = norm(ptA - ptB);
        else
            y = -norm(ptA - ptB);
        end
    end  
end