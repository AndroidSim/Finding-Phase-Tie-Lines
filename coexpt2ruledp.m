function [u,v] = coexpt2ruledp(a,coexpt,directrixO,directrixA)
% v = coexpt2ruledp(Ccoex,directrix_O,directrix_A)
% input: coexpt = coexistence point 
%        directrix_O = Lo section of boundary
%        directrix_A = La section of boundary
% output: v = values to zero

[np,nd] = size(coexpt);
[nop,ndo] = size(directrixO);
ndo = ndo-1; % 1st column is u parameters
[nap,nda] = size(directrixA);
nda = nda-1; % 1st column is v(gamma) parameters
if nd ~= ndo | nd ~= nda | ndo ~= nda
    error('dimension of points must agree');
else
%     options = optimset('MaxFunEvals',500,'TolFun',1e-8,'TolX',1e-8,'MaxIter',500,'Display','off','NonlEqnAlgorithm','gn');
    if nd == 2 % cartesian coordinates
        u = 0.5;
        ustart = 0;
        uend = 1;
        oC = interp1(directrixO(:,1),directrixO(:,2:end),u);
        aC = interp1(directrixA(:,1),directrixA(:,2:end),tieline_config_fxn(a,u));
        oCstart = interp1(directrixO(:,1),directrixO(:,2:end),ustart);
        oCend = interp1(directrixO(:,1),directrixO(:,2:end),uend);
        dstart = det([(oC'-aC') (aC'-oCstart')])/norm(oC-aC);
        dend = det([(oC'-aC') (aC'-oCend')])/norm(oC-aC);
        
%         u = fsolve(@eqn2zero2,0.5,options);
%         try
%             u = fzero(@eqn2zero,[0 1],options);
%         catch
%             keyboard;
%         end
        
        done = logical(0);
%         il = 1;
%         iu = length(b);
% start bisection search
        while ~done % iu-il > 1
            oC = interp1(directrixO(:,1),directrixO(:,2:end),u);
            aC = interp1(directrixA(:,1),directrixA(:,2:end),tieline_config_fxn(a,u));
            d = det([(oC'-aC') (aC'-coexpt')])/norm(oC-aC);
            if abs(ustart-uend) < 1e-8 || sign(d) == 0
                done = logical(1);
            elseif sign(d) == sign(dend)
                ustart = u;
            elseif sign(d) == sign(dstart)
                uend = u;
            else
                keyboard;% should not happen!
            end
            u = (ustart+uend)/2;
        end
%             im = round((iu+il)/2);
%             if bi(k) > b(im)
%                 il = im;
%             else
%                 iu = im;
%             end
            % determine new range for binary search (bisection)
    %         data_pt_dist=(((a_pt(2)-b_pt(2)).*data_pt(1))+((b_pt(1)-a_pt(1)).*data_pt(2))+((a_pt(1).*b_pt(2))-(b_pt(1).*a_pt(2))))./(norm(b_pt-a_pt));
    %         if (sign(data_pt_dist) == sign(U)) 
    %             xs=x;
    %         elseif (sign(data_pt_dist) == sign(L))
    %             xe=x;
    %         elseif (iseqwithn(data_pt_dist,0,10^-2)) % sign(data_pt_dist) == 0
    %             done=logical(1);
    %             xt=x;
    %             alpha_pt=a_pt;
    %             beta_pt=b_pt;
    %         else
    %             xs,xe,x
    %             keyboard
    %             error('no convergence');
    %         end
    %         % choose new x based on half-way between new range
    %         x=(xs+xe)./2;
    end
    if nd == 3 % ternary homogeneous coordinates
        u = fsolve(@eqn2zero3,0.5,options); 
%         u = fzero(@eqn2zero,[0 1],options)
    end
end
% if u < 0 
%     u = -u;
% end
% if u > 1
%     u = 1-u;
% end
v = tieline_config_fxn(a,u);
    function y = eqn2zero2(x)
        if x < 0 
            x = -x;
        end
        if x > 1
            x = 1-x;
        end   
        oC = interp1(directrixO(:,1),directrixO(:,2:end),x);
        aC = interp1(directrixA(:,1),directrixA(:,2:end),tieline_config_fxn(a,x));
        if all(aC == oC)
            fo = 1;
            fa = 1;
        else
            fo = norm(coexpt-aC)/norm(aC-oC);
            fa = 1-fo;
        end
        y = coexpt' - [oC.*fo]' - [aC.*fa]';
    end  
    function y = eqn2zero3(x)
        if x < 0 
            x = -x;
        end
        if x > 1
            x = 1-x;
        end 
        oC = interp1(directrixO(:,1),directrixO(:,2:end),x);
        aC = interp1(directrixA(:,1),directrixA(:,2:end),tieline_config_fxn(a,x));
        if all(aC == oC)
            fo = 1;
            fa = 1;
        else
            fo = det([coexpt(3) coexpt(1);aC(3) aC(1)])/det([aC(1) aC(3);oC(1) oC(3)]);
            fa = 1-fo;
        end
        y = coexpt' - [oC.*fo]' - [aC.*fa]';
    end 
    function y = eqn2zero(x)
        oC = interp1(directrixO(:,1),directrixO(:,2:end),x);
        aC = interp1(directrixA(:,1),directrixA(:,2:end),tieline_config_fxn(a,x));
        if x == 0
            y = dstart;
        elseif x == 1
            y = dend;
        else
            y = det([(oC'-aC') (aC'-coexpt')])/norm(oC-aC);
        end
%         if x == 0
%             y = norm(coexpt-aC);
%         elseif x == 1
%             y = -norm(coexpt-aC);
%         else
%             y = det([(oC'-aC') (aC'-coexpt')])/norm(oC-aC);
%         end
    end
end             

% function [alpha,xt,alpha_pt,beta_pt] = get_tieline_coords(bd_pts,bd_rep,bd_length,bd_intervals,gamma1,gamma2,Xsd,Xcd)
% % tieline_pt_coords(args) determines the model variables alpha (distance along tie line) and xt (which tie line) 
% % along tie line between alpha_pt and beta_pt for data point ([Xsd Xcd]) and critical points parameters gamma1 
% % and gamma2 with boundary bd_pts and its representation bd_rep
% data_pt=[Xsd Xcd];
% done=logical(0);
% swapped=logical(0);
% %n=1;
% if (gamma1 > gamma2)
%     swap=gamma1;
%     gamma1=gamma2;
%     gamma2=swap;
%     swapped=logical(1);
% end
% % initialization
% x=0.5;
% xs=0; % critical point 1 given by gamma 1
% xe=1; % critical point 2 given by gamma 2
% if (swapped)
%     sb=gamma1+x*(gamma2-gamma1);
%     sa=1+gamma1-x*(1-gamma2+gamma1);
%     a_pt=get_bd_pt(bd_pts,bd_rep,sa,bd_length,bd_intervals);
%     b_pt=get_bd_pt(bd_pts,bd_rep,sb,bd_length,bd_intervals);
%     L_cp=get_bd_pt(bd_pts,bd_rep,gamma2,bd_length,bd_intervals);
%     U_cp=get_bd_pt(bd_pts,bd_rep,gamma1,bd_length,bd_intervals);
%     U=(((a_pt(2)-b_pt(2)).*U_cp(1))+((b_pt(1)-a_pt(1)).*U_cp(2))+((a_pt(1).*b_pt(2))-(b_pt(1).*a_pt(2))))./(norm(b_pt-a_pt)); 
%     L=(((a_pt(2)-b_pt(2)).*L_cp(1))+((b_pt(1)-a_pt(1)).*L_cp(2))+((a_pt(1).*b_pt(2))-(b_pt(1).*a_pt(2))))./(norm(b_pt-a_pt));
% else % gamma1 <= gamma2
%     sa=gamma1+x*(gamma2-gamma1);
%     sb=1+gamma1-x*(1-gamma2+gamma1);
%     a_pt=get_bd_pt(bd_pts,bd_rep,sa,bd_length,bd_intervals);
%     b_pt=get_bd_pt(bd_pts,bd_rep,sb,bd_length,bd_intervals);
%     L_cp=get_bd_pt(bd_pts,bd_rep,gamma1,bd_length,bd_intervals);
%     U_cp=get_bd_pt(bd_pts,bd_rep,gamma2,bd_length,bd_intervals);
%     U=(((a_pt(2)-b_pt(2)).*U_cp(1))+((b_pt(1)-a_pt(1)).*U_cp(2))+((a_pt(1).*b_pt(2))-(b_pt(1).*a_pt(2))))./(norm(b_pt-a_pt)); 
%     L=(((a_pt(2)-b_pt(2)).*L_cp(1))+((b_pt(1)-a_pt(1)).*L_cp(2))+((a_pt(1).*b_pt(2))-(b_pt(1).*a_pt(2))))./(norm(b_pt-a_pt));
% end
% while (~done) 
%     % find the tie line (xt == value of x from 0 to 1) that passes through [Xsd Xcd]
%     if (swapped) 
%         sb=gamma1+x*(gamma2-gamma1);
%         sa=1+gamma1-x*(1-gamma2+gamma1);
%     else % gamma1 <= gamma2
%         sa=gamma1+x*(gamma2-gamma1);
%         sb=1+gamma1-x*(1-gamma2+gamma1);
%     end
%     a_pt=get_bd_pt(bd_pts,bd_rep,sa,bd_length,bd_intervals);
%     b_pt=get_bd_pt(bd_pts,bd_rep,sb,bd_length,bd_intervals);
%     m=(a_pt(2)-b_pt(2))/(a_pt(1)-b_pt(1));
%     if ((iseqwithn(data_pt(2),m*(data_pt(1)-a_pt(1))+a_pt(2),10^-2)) | (isinf(m) & iseqwithn(a_pt(1),data_pt(1),10^-2))) % a_pts(i,1) == data_pt(1)
%         done=logical(1);
%         xt=x;
%         alpha_pt=a_pt;
%         beta_pt=b_pt;
%         break
%     else
%         % determine new range for binary search (bisection)
%         data_pt_dist=(((a_pt(2)-b_pt(2)).*data_pt(1))+((b_pt(1)-a_pt(1)).*data_pt(2))+((a_pt(1).*b_pt(2))-(b_pt(1).*a_pt(2))))./(norm(b_pt-a_pt));
%         if (sign(data_pt_dist) == sign(U)) 
%             xs=x;
%         elseif (sign(data_pt_dist) == sign(L))
%             xe=x;
%         elseif (iseqwithn(data_pt_dist,0,10^-2)) % sign(data_pt_dist) == 0
%             done=logical(1);
%             xt=x;
%             alpha_pt=a_pt;
%             beta_pt=b_pt;
%         else
%             xs,xe,x
%             keyboard
%             error('no convergence');
%         end
%         %if (sign(data_pt_dist) == sign(L_dist) == sign(U_dist))
%         %    error('data pt has same sign as L and U')
%         %end
%         % choose new x based on half-way between new range
%         x=(xs+xe)./2;
%         %[xs,xe]=tie_bisect(x,ixs,ixe,data_pt,a_pt,b_pt,cp1,cp2);
%     end
%     %n=n+1;
% end
% alpha=norm(data_pt-beta_pt)/norm(alpha_pt-beta_pt);
% %disp(sprintf('number tie line iterations for xt = %f\n',n));
% return