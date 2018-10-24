function [alpha,xt,alpha_pt,beta_pt] = tieline_pt_coords(bd_pts,bd_rep,gamma1,gamma2,Xsd,Xcd)

% tieline_pt_coords(args) determines the model variables alpha (distance along tie line) and xt (which tie line) 
% along tie line between alpha_pt and beta_pt for data point ([Xsd Xcd]) and critical points parameters gamma1 
% and gamma2 with boundary bd_pts and its representation bd_rep
data_pt=[Xsd Xcd];
acc=1;
xs=0;
xe=1;
%xparams=[xs xe];
done=logical(0);
swapped=logical(0);
a_pts=zeros(11,2);
b_pts=zeros(11,2);
%n=1;
if (gamma1 > gamma2)
    swap=gamma1;
    gamma1=gamma2;
    gamma2=swap;
    swapped=logical(1);
end
while (~done) 
    % find the two tie lines bracketing the tie line that passes through [Xsd Xcd]
    x=[xs:10^-acc:xe];
    if (swapped) 
        sb=gamma1+x*(gamma2-gamma1);
        sa=1+gamma1-x*(1-gamma2+gamma1);
    else % gamma1 <= gamma2
        sa=gamma1+x*(gamma2-gamma1);
        sb=1+gamma1-x*(1-gamma2+gamma1);
    end
    for i=1:length(x)
        a_pts(i,:)=pt_on_boundary(bd_pts,bd_rep,sa(i),bd_pts(1,:),bd_pts(end,:));
        b_pts(i,:)=pt_on_boundary(bd_pts,bd_rep,sb(i),bd_pts(1,:),bd_pts(end,:));
        warning off
        m=(a_pts(i,2)-b_pts(i,2))/(a_pts(i,1)-b_pts(i,1));
        warning on
        if ((iseqwithn(data_pt(2),m*(data_pt(1)-a_pts(i,1))+a_pts(i,2))) | (isinf(m) & iseqwithn(a_pts(i,1),data_pt(1)))) % a_pts(i,1) == data_pt(1)
            done=logical(1);
            xt=x(i);
            alpha_pt=a_pts(i,:);
            beta_pt=b_pts(i,:);
            break
        end
    end
    if (~done)
        % find tie lines bracketing data point
        [xs,xe]=brackets_test(x,data_pt,a_pts,b_pts);
        acc=acc+1;
    end
    %n=n+1;
end
alpha=norm(data_pt-beta_pt)/norm(alpha_pt-beta_pt);
%disp(sprintf('number tie line iterations for xt = %f\n',n));