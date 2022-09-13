function [u,directrixO] = tangent2ruled(x,C_bdy,boundary)
% [u,directrixO] = tangent2ruled(slope,bdy,coexpt,bdypts,bdyb)
% converts tangent tieline configuration to a ruled tieline configuration
% by finding the two directrices given the two critical points.

switch boundary.config
    case '2critpts'
        [s1,section1,s2,section2] = divide_bdy(x(1),x(2),C_bdy);
        if x(1) <= x(2)
            u = s1;
            directrixO = section1;
        else
            u = s2;
            directrixO = section2;
        end
    case '1critpt/1endtl'
        b = boundary.parameterization;
        cp = x(1);
        e1 = x(2);
        e2 = x(3);
        if e1 > e2
            % divide boundary into two sections (same procedure
            % when going from tangent to ruled in the 2
            % critical point case above).
            [s1,section1,s2,section2] = divide_bdy(x(2),x(3),C_bdy);% x(2),x(3) = e1,e2
%             [u,directrixO,v,directrixA] = tangent2ruled(e1,e2,C_bdy);
            if cp < e1 & cp > e2 % configuration D in notes
                Cbdy = section1;
                bdy = [Cbdy;Cbdy(1,:)];
                cpt = b2bdypt(b,C_bdy,x(1));
                cps1 = bdypt2b(cpt,s1,section1);
                [cpt,braket] = b2bdypt(s1,section1,cps1);
                bracket = braket{:};
                if isscalar(bracket)
                    i = bracket;
                    directrixO = flipud(section1(1:i,:));
                    directrixA = section1(i:end,:);
                elseif isvector(bracket)
                    il = bracket(1);
                    iu = bracket(2);
                    directrixO = [cpt;flipud(section1(1:il,:))];
                    directrixA = [cpt;section1(iu:end,:)];
                else
                    error('bracket cannot be a matrix');
                end
            elseif cp > e1 | cp < e2 % configuration B in notes
                Cbdy = section2;
                bdy = [Cbdy;Cbdy(1,:)];
                cpt = b2bdypt(b,C_bdy,x(1));
                cps2 = bdypt2b(cpt,s2,section2);
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
            else % cp == e1 | cp == e2
                error('cp cannot equal e1 or e2');
            end
            [Olength,u] = chord_length(directrixO);
            [Alength,v] = chord_length(directrixA);
        elseif e1 < e2
            % divide boundary into two sections (same procedure
            % when going from tangent to ruled in the 2
            % critical point case above).
            [s1,section1,s2,section2] = divide_bdy(x(2),x(3),C_bdy);% x(2),x(3) = e1,e2
%             [u,directrixO,v,directrixA] = tangent2ruled(e1,e2,C_bdy);
            if cp > e1 & cp < e2 % configuration C in notes
                Cbdy = section2;
                bdy = [Cbdy;Cbdy(1,:)];
                cpt = b2bdypt(b,C_bdy,x(1));
                cps2 = bdypt2b(cpt,s2,section2);
                [cpt,braket] = b2bdypt(s2,section2,cps2);
                bracket = braket{:};
                if isscalar(bracket)
                    i = bracket;
                    directrixO = flipud(section2(1:i,:));
                    directrixA = section2(i:end,:);
                elseif isvector(bracket)
                    il = bracket(1);
                    iu = bracket(2);
                    directrixO = [cpt;flipud(section2(1:il,:))];
                    directrixA = [cpt;section2(iu:end,:)];
                else
                    error('bracket cannot be a matrix');
                end
            elseif cp < e1 | cp > e2 % configuration A in notes
                Cbdy = section1;
                bdy = [Cbdy;Cbdy(1,:)];
                cpt = b2bdypt(b,C_bdy,x(1));
                cps1 = bdypt2b(cpt,s1,section1);
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
            else % cp == e1 | cp == e2
                error('cp cannot equal e1 or e2');
            end
            [Olength,u] = chord_length(directrixO);
            [Alength,v] = chord_length(directrixA);
        end
    case '2endtls'
    otherwise
        error('invalid boundary configuration');
end

return