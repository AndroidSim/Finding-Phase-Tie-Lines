function tie_lines = draw_tie_lines(bd_pts,bd_rep,nt,gamma1,gamma2)
% draw_tie_lines plots tie lines for relevant model parameters of new tie line determination routine

[bd_length,bd_intervals] = boundary(bd_pts,bd_rep);
a_pts=zeros(nt,2);
b_pts=zeros(nt,2); 
% find the two tie lines bracketing the tie line that passes through [Xsd Xcd]
x=linspace(0,1,nt);
if (gamma1 <= gamma2)
    sa=gamma1+x*(gamma2-gamma1);
    sb=1+gamma1-x*(1-gamma2+gamma1);
    %sa=gamma1+x*(gamma2-gamma1);
    %sb=1+gamma1-x*(1-gamma2+gamma1);
else
    store=gamma1;
    gamma1=gamma2;
    gamma2=store;
    sb=gamma1+x*(gamma2-gamma1);
    sa=1+gamma1-x*(1-gamma2+gamma1);
    %sb=gamma1+x*(gamma2-gamma1);
    %sa=1+gamma1-x*(1-gamma2+gamma1);
end
tie_lines=zeros(nt,50,2);
for i=1:nt
    alpha_pts(i,:)=get_bd_pt(bd_pts,bd_rep,sa(i),bd_length,bd_intervals);
    beta_pts(i,:)=get_bd_pt(bd_pts,bd_rep,sb(i),bd_length,bd_intervals);
    tie_lines(i,:,1)=linspace(alpha_pts(i,1),beta_pts(i,1),50);
    tie_lines(i,:,2)=linspace(alpha_pts(i,2),beta_pts(i,2),50);
end
ab_pts=cat(3,alpha_pts,beta_pts);

plot_points(bd_pts);
hold on;
plot_tielinesandpoints(tie_lines,ab_pts);
hold off;