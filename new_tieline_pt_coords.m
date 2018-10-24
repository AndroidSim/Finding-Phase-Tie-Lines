function [alpha,xt,alpha_pt,beta_pt] = new_tieline_pt_coords(bd_pts,bd_rep,gamma1,gamma2,Xsd,Xcd)

% tieline_pt_coords(args) determines the model variables alpha (distance along tie line) and xt (which tie line) 
% along tie line between alpha_pt and beta_pt for data point ([Xsd Xcd]) and critical points parameters gamma1 
% and gamma2 with boundary bd_pts and its representation bd_rep
data_pt=[Xsd Xcd];
done=logical(0);
swapped=logical(0);
%n=1;
if (gamma1 > gamma2)
    swap=gamma1;
    gamma1=gamma2;
    gamma2=swap;
    swapped=logical(1);
end
% initialization
x=0.5;
xs=0; % critical point 1 given by gamma 1
xe=1; % critical point 2 given by gamma 2
if (swapped) 
    L_cp=new_pt_on_boundary(bd_pts,bd_rep,gamma2,bd_pts(1,:),bd_pts(end,:));
    U_cp=new_pt_on_boundary(bd_pts,bd_rep,gamma1,bd_pts(1,:),bd_pts(end,:));
else % gamma1 <= gamma2
    L_cp=new_pt_on_boundary(bd_pts,bd_rep,gamma1,bd_pts(1,:),bd_pts(end,:));
    U_cp=new_pt_on_boundary(bd_pts,bd_rep,gamma2,bd_pts(1,:),bd_pts(end,:));
end
while (~done) 
    % find the tie line (xt == value of x from 0 to 1) that passes through [Xsd Xcd]
    if (swapped) 
        sb=gamma1+x*(gamma2-gamma1);
        sa=1+gamma1-x*(1-gamma2+gamma1);
    else % gamma1 <= gamma2
        sa=gamma1+x*(gamma2-gamma1);
        sb=1+gamma1-x*(1-gamma2+gamma1);
    end
    a_pt=new_pt_on_boundary(bd_pts,bd_rep,sa,bd_pts(1,:),bd_pts(end,:));
    b_pt=new_pt_on_boundary(bd_pts,bd_rep,sb,bd_pts(1,:),bd_pts(end,:));
    warning off
    m=(a_pt(2)-b_pt(2))/(a_pt(1)-b_pt(1));
    warning on
    if ((iseqwithn(data_pt(2),m*(data_pt(1)-a_pt(1))+a_pt(2))) | (isinf(m) & iseqwithn(a_pt(1),data_pt(1)))) % a_pts(i,1) == data_pt(1)
        done=logical(1);
        xt=x;
        alpha_pt=a_pt;
        beta_pt=b_pt;
        break
    end
    if (~done)
        % determine new range for binary search (bisection)
        sxs=gamma1+xs*(gamma2-gamma1);
        sxe=gamma1+xe*(gamma2-gamma1);
        L=new_pt_on_boundary(bd_pts,bd_rep,sxs,bd_pts(1,:),bd_pts(end,:));
        U=new_pt_on_boundary(bd_pts,bd_rep,sxe,bd_pts(1,:),bd_pts(end,:));
        U_dist=(((a_pt(2)-b_pt(2)).*U(1))+((b_pt(1)-a_pt(1)).*U(2))+((a_pt(1).*b_pt(2))-(b_pt(1).*a_pt(2))))./(norm(b_pt-a_pt)); 
        L_dist=(((a_pt(2)-b_pt(2)).*L(1))+((b_pt(1)-a_pt(1)).*L(2))+((a_pt(1).*b_pt(2))-(b_pt(1).*a_pt(2))))./(norm(b_pt-a_pt)); 
        data_pt_dist=(((a_pt(2)-b_pt(2)).*data_pt(1))+((b_pt(1)-a_pt(1)).*data_pt(2))+((a_pt(1).*b_pt(2))-(b_pt(1).*a_pt(2))))./(norm(b_pt-a_pt));
        if ((data_pt_dist < 0 & U_dist < 0) | (data_pt_dist > 0 & U_dist > 0))
            xs=x;
        end
        if ((data_pt_dist < 0 & L_dist < 0) | (data_pt_dist > 0 & L_dist > 0))
            xe=x;
        end
        % choose new x based on half-way between new range
        x=(xs+xe)./2;
        %[xs,xe]=tie_bisect(x,ixs,ixe,data_pt,a_pt,b_pt,cp1,cp2);
    end
    %n=n+1;
end
alpha=norm(data_pt-beta_pt)/norm(alpha_pt-beta_pt);
%disp(sprintf('number tie line iterations for xt = %f\n',n));