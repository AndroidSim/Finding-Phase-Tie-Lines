function [s1,section1,s2,section2] = divide_bdy(b1,b2,bdy)

% the section 2 is the section of the boundary from b1 to b2
% in increasing order and section 1 is from b1 to b2 in
% decreasing order.
% if b1 <= b2, section1 = directrixO, section2 = directrixA
% if b1 > b2, section1 = directrixA, section2 = directrixO

[blength,b] = chord_length(bdy);
[pt1,braket1] = b2bdypt(b,bdy,b1);
[pt2,braket2] = b2bdypt(b,bdy,b2);
bracket1 = braket1{1}; % braket1{:}
bracket2 = braket2{1}; % braket2{:}
if isscalar(bracket1)
    if isscalar(bracket2)
        ib1 = bracket1;
        ib2 = bracket2;
        if b1 <= b2
            section2 = [bdy(ib1:ib2,:)];
            section1 = [flipud(bdy(1:ib1,:));flipud(bdy(ib2:end-1,:))];
        else
            section2 = [bdy(ib1:end-1,:);bdy(1:ib2,:)];
            section1 = [flipud(bdy(ib2:ib1,:))];
        end  
    elseif isvector(bracket2)
        ib1 = bracket1;
        ib2l = bracket2(1);
        ib2u = bracket2(2);
        if b1 <= b2
            section2 = [bdy(ib1:ib2l,:);pt2];
            section1 = [flipud(bdy(1:ib1,:));flipud(bdy(ib2u:end-1,:));pt2];
        else
            section2 = [bdy(ib1:end-1,:);bdy(1:ib2l,:);pt2];
            section1 = [flipud(bdy(ib2u:ib1,:));pt2];
        end
    else
        error('bracket2 cannot be a matrix');
    end
elseif isvector(bracket1)
    if isscalar(bracket2)
        ib1l = bracket1(1);
        ib1u = bracket1(2);
        ib2 = bracket2;
        if b1 <= b2
            section2 = [pt1;bdy(ib1u:ib2,:)];
            section1 = [pt1;flipud(bdy(1:ib1l,:));flipud(bdy(ib2:end-1,:))];
        else
            section2 = [pt1;bdy(ib1u:end-1,:);bdy(1:ib2,:)];
            section1 = [pt1;flipud(bdy(ib2:ib1l,:))];
        end
    elseif isvector(bracket2)
        ib1l = bracket1(1);
        ib1u = bracket1(2);
        ib2l = bracket2(1);
        ib2u = bracket2(2);
        if b1 <= b2
            section2 = [pt1;bdy(ib1u:ib2l,:);pt2];
            section1 = [pt1;flipud(bdy(1:ib1l,:));flipud(bdy(ib2u:end-1,:));pt2];
        else
            section2 = [pt1;bdy(ib1u:end-1,:);bdy(1:ib2l,:);pt2];
            section1 = [pt1;flipud(bdy(ib2u:ib1l,:));pt2];
        end
    else
        error('bracket2 cannot be a matrix');
    end
else
    error('bracket1 cannot be a matrix');
end
[s1length,s1] = chord_length(section1);
[s2length,s2] = chord_length(section2);
return

% b = boundary.parameterization;
%         bdy = C_bdy;
% %         [u,directrixO,v,directrixA] = divide_bdy(x(1),x(2),C_bdy);
%         % x(1) = the boundary parameter of critical point 1
%         % x(2) = the boundary parameter of critical point 2
%         % the La directrix is the section of the boundary from x(1) to x(2)
%         % in increasing order and the Lo directrix is from x(1) to x(2) in
%         % decreasing order.
%         [cpt1,braket1] = b2bdypt(b,bdy,x(1));
%         [cpt2,braket2] = b2bdypt(b,bdy,x(2));
%         bracket1 = braket1{:};
%         bracket2 = braket2{:};
%         if isscalar(bracket1)
%             if isscalar(bracket2)
%                 ib1 = bracket1;
%                 ib2 = bracket2;
%                 if x(1) <= x(2)
%                     directrixA = [bdy(ib1:ib2,:)];
%                     directrixO = [flipud(bdy(1:ib1,:));flipud(bdy(ib2:end-1,:))];
%                 else
%                     directrixA = [bdy(ib1:end-1,:);bdy(1:ib2,:)];
%                     directrixO = [flipud(bdy(ib2:ib1,:))];
%                 end  
%             elseif isvector(bracket2)
%                 ib1 = bracket1;
%                 ib2l = bracket2(1);
%                 ib2u = bracket2(2);
%                 if x(1) <= x(2)
%                     directrixA = [bdy(ib1:ib2l,:);cpt2];
%                     directrixO = [flipud(bdy(1:ib1,:));flipud(bdy(ib2u:end-1,:));cpt2];
%                 else
%                     directrixA = [bdy(ib1:end-1,:);bdy(1:ib2l,:);cpt2];
%                     directrixO = [flipud(bdy(ib2u:ib1,:));cpt2];
%                 end
%             else
%                 error('bracket2 cannot be a matrix');
%             end
%         elseif isvector(bracket1)
%             if isscalar(bracket2)
%                 ib1l = bracket1(1);
%                 ib1u = bracket1(2);
%                 ib2 = bracket2;
%                 if x(1) <= x(2)
%                     directrixA = [cpt1;bdy(ib1u:ib2,:)];
%                     directrixO = [cpt1;flipud(bdy(1:ib1l,:));flipud(bdy(ib2:end-1,:))];
%                 else
%                     directrixA = [cpt1;bdy(ib1u:end-1,:);bdy(1:ib2,:)];
%                     directrixO = [cpt1;flipud(bdy(ib2:ib1l,:))];
%                 end
%             elseif isvector(bracket2)
%                 ib1l = bracket1(1);
%                 ib1u = bracket1(2);
%                 ib2l = bracket2(1);
%                 ib2u = bracket2(2);
%                 if x(1) <= x(2)
%                     directrixA = [cpt1;bdy(ib1u:ib2l,:);cpt2];
%                     directrixO = [cpt1;flipud(bdy(1:ib1l,:));flipud(bdy(ib2u:end-1,:));cpt2];
%                 else
%                     directrixA = [cpt1;bdy(ib1u:end-1,:);bdy(1:ib2l,:);cpt2];
%                     directrixO = [cpt1;flipud(bdy(ib2u:ib1l,:));cpt2];
%                 end
%             else
%                 error('bracket2 cannot be a matrix');
%             end
%         else
%             error('bracket1 cannot be a matrix');
%         end
% %         sa = gamma1+x*(gamma2-gamma1);
% %         sb = 1+gamma1-x*(1-gamma2+gamma1);
%         [Olength,u] = chord_length(directrixO);
%         [Alength,v] = chord_length(directrixA);
% %                 if any(isnan(directrixO))
% %                     keyboard;
% %                 end
% %                 if any(isnan(directrixA))
% %                     keyboard;
% %                 end