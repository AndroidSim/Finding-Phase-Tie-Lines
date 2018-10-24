function [alpha,beta]=alphas_betas(tie_line)
% the function alphas_betas calculates the fraction of the disordered phase alpha
% and the fraction of the ordered phase beta for each point along a tie line.  the
% assumption is that the start point of the tie line (first row of vector) begins on 
% the disordered phase boundary the the last point of the tie line (last row of vector)
% ends on the ordered phase boundary
% tie_line=cat(2,tielines(i,:,1)',tielines(i,:,2)');

%[ntpts,xy]=size(tie_line);
ntpts=size(tie_line,1);
%tot_dist=sqrt(((tie_line(1,1)-tie_line(ntpts,1)).^2)+((tie_line(1,2)-tie_line(ntpts,2)).^2));
%tot_dist=sqrt(((tie_line(1,1)-tie_line(end,1)).^2)+((tie_line(1,2)-tie_line(end,2)).^2));
tot_dist=norm(tie_line(1,:)-tie_line(end,:));
for i=1:ntpts
    %dist_pt=sqrt(((tie_line(i,1)-tie_line(1,1)).^2)+((tie_line(i,2)-tie_line(1,2)).^2));
    %beta(i)=dist_pt./tot_dist;
    dist_from_alfa(i)=norm(tie_line(i,:)-tie_line(1,:));
end
beta=dist_from_alfa./tot_dist;
alpha=1-beta;