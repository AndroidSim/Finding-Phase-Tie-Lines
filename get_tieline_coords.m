function [alpha,xt,alpha_pt,beta_pt] = get_tieline_coords(bd_pts,bd_rep,bd_length,bd_intervals,gamma1,gamma2,Xsd,Xcd)

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
    sb=gamma1+x*(gamma2-gamma1);
    sa=1+gamma1-x*(1-gamma2+gamma1);
    a_pt=get_bd_pt(bd_pts,bd_rep,sa,bd_length,bd_intervals);
    b_pt=get_bd_pt(bd_pts,bd_rep,sb,bd_length,bd_intervals);
    L_cp=get_bd_pt(bd_pts,bd_rep,gamma2,bd_length,bd_intervals);
    U_cp=get_bd_pt(bd_pts,bd_rep,gamma1,bd_length,bd_intervals);
    U=(((a_pt(2)-b_pt(2)).*U_cp(1))+((b_pt(1)-a_pt(1)).*U_cp(2))+((a_pt(1).*b_pt(2))-(b_pt(1).*a_pt(2))))./(norm(b_pt-a_pt)); 
    L=(((a_pt(2)-b_pt(2)).*L_cp(1))+((b_pt(1)-a_pt(1)).*L_cp(2))+((a_pt(1).*b_pt(2))-(b_pt(1).*a_pt(2))))./(norm(b_pt-a_pt));
else % gamma1 <= gamma2
    sa=gamma1+x*(gamma2-gamma1);
    sb=1+gamma1-x*(1-gamma2+gamma1);
    a_pt=get_bd_pt(bd_pts,bd_rep,sa,bd_length,bd_intervals);
    b_pt=get_bd_pt(bd_pts,bd_rep,sb,bd_length,bd_intervals);
    L_cp=get_bd_pt(bd_pts,bd_rep,gamma1,bd_length,bd_intervals);
    U_cp=get_bd_pt(bd_pts,bd_rep,gamma2,bd_length,bd_intervals);
    U=(((a_pt(2)-b_pt(2)).*U_cp(1))+((b_pt(1)-a_pt(1)).*U_cp(2))+((a_pt(1).*b_pt(2))-(b_pt(1).*a_pt(2))))./(norm(b_pt-a_pt)); 
    L=(((a_pt(2)-b_pt(2)).*L_cp(1))+((b_pt(1)-a_pt(1)).*L_cp(2))+((a_pt(1).*b_pt(2))-(b_pt(1).*a_pt(2))))./(norm(b_pt-a_pt));
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
    a_pt=get_bd_pt(bd_pts,bd_rep,sa,bd_length,bd_intervals);
    b_pt=get_bd_pt(bd_pts,bd_rep,sb,bd_length,bd_intervals);
    m=(a_pt(2)-b_pt(2))/(a_pt(1)-b_pt(1));
    if ((iseqwithn(data_pt(2),m*(data_pt(1)-a_pt(1))+a_pt(2),10^-2)) | (isinf(m) & iseqwithn(a_pt(1),data_pt(1),10^-2))) % a_pts(i,1) == data_pt(1)
        done=logical(1);
        xt=x;
        alpha_pt=a_pt;
        beta_pt=b_pt;
        break
    else
        % check to make sure data point is bracketed by upper and lower x boundary
        %if (xs ~= 0 & xe ~= 1)
        %    if (swapped)
        %        sbxs=gamma1+xs*(gamma2-gamma1);
        %        saxs=1+gamma1-xs*(1-gamma2+gamma1);
        %        sbxe=gamma1+xe*(gamma2-gamma1);
        %        saxe=1+gamma1-xe*(1-gamma2+gamma1);
        %        xsa_pt=get_bd_pt(bd_pts,bd_rep,saxs,bd_length,bd_intervals);
        %        xsb_pt=get_bd_pt(bd_pts,bd_rep,sbxs,bd_length,bd_intervals);
        %        xea_pt=get_bd_pt(bd_pts,bd_rep,saxe,bd_length,bd_intervals);
        %        xeb_pt=get_bd_pt(bd_pts,bd_rep,sbxe,bd_length,bd_intervals);
        %        U_dist=(((xea_pt(2)-xeb_pt(2)).*data_pt(1))+((xeb_pt(1)-xea_pt(1)).*data_pt(2))+((xea_pt(1).*xeb_pt(2))-(xeb_pt(1).*xea_pt(2))))./(norm(xeb_pt-xea_pt)); 
        %        L_dist=(((xsa_pt(2)-xsb_pt(2)).*data_pt(1))+((xsb_pt(1)-xsa_pt(1)).*data_pt(2))+((xsa_pt(1).*xsb_pt(2))-(xsb_pt(1).*xsa_pt(2))))./(norm(xsb_pt-xsa_pt));
        %    else % gamma1 <= gamma2
        %        saxs=gamma1+xs*(gamma2-gamma1);
        %        sbxs=1+gamma1-xs*(1-gamma2+gamma1);
        %        saxe=gamma1+xe*(gamma2-gamma1);
        %        sbxe=1+gamma1-xe*(1-gamma2+gamma1);
        %        xsa_pt=get_bd_pt(bd_pts,bd_rep,saxs,bd_length,bd_intervals);
        %        xsb_pt=get_bd_pt(bd_pts,bd_rep,sbxs,bd_length,bd_intervals);
        %        xea_pt=get_bd_pt(bd_pts,bd_rep,saxe,bd_length,bd_intervals);
        %        xeb_pt=get_bd_pt(bd_pts,bd_rep,sbxe,bd_length,bd_intervals);
        %        U_dist=(((xea_pt(2)-xeb_pt(2)).*data_pt(1))+((xeb_pt(1)-xea_pt(1)).*data_pt(2))+((xea_pt(1).*xeb_pt(2))-(xeb_pt(1).*xea_pt(2))))./(norm(xeb_pt-xea_pt)); 
        %        L_dist=(((xsa_pt(2)-xsb_pt(2)).*data_pt(1))+((xsb_pt(1)-xsa_pt(1)).*data_pt(2))+((xsa_pt(1).*xsb_pt(2))-(xsb_pt(1).*xsa_pt(2))))./(norm(xsb_pt-xsa_pt));
        %    end
        %    if (sign(U_dist) == sign(L_dist))
        %        error('data point not bracketed')
        %    end
        %    if (sign(U_dist) == 0 | sign(L_dist) == 0)
        %        error('U and L bounds lie on data point')
        %    end
        %end
        % determine new range for binary search (bisection)
        data_pt_dist=(((a_pt(2)-b_pt(2)).*data_pt(1))+((b_pt(1)-a_pt(1)).*data_pt(2))+((a_pt(1).*b_pt(2))-(b_pt(1).*a_pt(2))))./(norm(b_pt-a_pt));
        if (sign(data_pt_dist) == sign(U)) 
            xs=x;
        elseif (sign(data_pt_dist) == sign(L))
            xe=x;
        elseif (iseqwithn(data_pt_dist,0,10^-2)) % sign(data_pt_dist) == 0
            done=logical(1);
            xt=x;
            alpha_pt=a_pt;
            beta_pt=b_pt;
        else
            xs,xe,x
            keyboard
            error('no convergence');
        end
        %if (sign(data_pt_dist) == sign(L_dist) == sign(U_dist))
        %    error('data pt has same sign as L and U')
        %end
        % choose new x based on half-way between new range
        x=(xs+xe)./2;
        %[xs,xe]=tie_bisect(x,ixs,ixe,data_pt,a_pt,b_pt,cp1,cp2);
    end
    %n=n+1;
end
alpha=norm(data_pt-beta_pt)/norm(alpha_pt-beta_pt);
%disp(sprintf('number tie line iterations for xt = %f\n',n));
return

%[bd_length,bd_intervals]=boundary(bd_pts,'linear');
%for i=1:length(xdata)
%xdata(i,:)
%[alpha,xt,alpha_pt,beta_pt] = get_tieline_coords(bd_pts,'linear',bd_length,bd_intervals,0.3,0.8,xdata(i,1),xdata(i,2));
%xt
%end
%figure
%if (swapped)
%    sbxs=gamma1+xs*(gamma2-gamma1);
%    saxs=1+gamma1-xs*(1-gamma2+gamma1);
%    sbxe=gamma1+xe*(gamma2-gamma1);
%    saxe=1+gamma1-xe*(1-gamma2+gamma1);
%    xsa_pt=get_bd_pt(bd_pts,bd_rep,saxs,bd_length,bd_intervals);
%    xsb_pt=get_bd_pt(bd_pts,bd_rep,sbxs,bd_length,bd_intervals);
%    xea_pt=get_bd_pt(bd_pts,bd_rep,saxe,bd_length,bd_intervals);
%    xeb_pt=get_bd_pt(bd_pts,bd_rep,sbxe,bd_length,bd_intervals);
%else % gamma1 <= gamma2
%    saxs=gamma1+xs*(gamma2-gamma1);
%    sbxs=1+gamma1-xs*(1-gamma2+gamma1);
%    saxe=gamma1+xe*(gamma2-gamma1);
%    sbxe=1+gamma1-xe*(1-gamma2+gamma1);
%    xsa_pt=get_bd_pt(bd_pts,bd_rep,saxs,bd_length,bd_intervals);
%    xsb_pt=get_bd_pt(bd_pts,bd_rep,sbxs,bd_length,bd_intervals);
%    xea_pt=get_bd_pt(bd_pts,bd_rep,saxe,bd_length,bd_intervals);
%    xeb_pt=get_bd_pt(bd_pts,bd_rep,sbxe,bd_length,bd_intervals);
%end
%plot_points(bd_pts)
%xline=linspace(xsa_pt(1),xsb_pt(1),50);
%yline=linspace(xsa_pt(2),xsb_pt(2),50);
%hold on,plot(xline,yline,'-r'),hold off
%xline=linspace(a_pt(1),b_pt(1),50);
%yline=linspace(a_pt(2),b_pt(2),50);
%hold on,plot(xline,yline,'-b'),hold off
%xline=linspace(xea_pt(1),xeb_pt(1),50);
%yline=linspace(xea_pt(2),xeb_pt(2),50);
%hold on,plot(xline,yline,'-g'),hold off
%hold on,plot(data_pt(1),data_pt(2),'dk'),hold off

%g1=[0:0.1:1];
%g2=[0:0.1:1];
%g3=[0.1 0.3 0.5 0.7 0.9 2 4 6 8 10];
%g4=[0.1 0.3 0.5 0.7 0.9 2 4 6 8 10];
%lsq_grid_search=zeros(length(g1),length(g2),length(g3),length(g4));
%for i=1:length(g1)
%for j=1:length(g2)
%for k=1:length(g3)
%for r=1:length(g4)
%if (g1(i) == g2(j) | ((g1(i) == 0 & g2(j) == 1) | (g1(i) == 1 & g2(j) == 0)))
%continue
%end
%[gammas,lsq] = FRET_theory_fit(Pauls_data,bd_pts,'linear','para',0.001,0.001,238,508,[g1(i) g2(j) g3(k) g4(r)]);
%lsq_grid_search(i,j,k,r)=lsq;
%end
%end
%end
%end