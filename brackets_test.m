function [xs,xe] = brackets_test(xp,qpt,a_pts,b_pts)
% brackets is a function that determines the two x parameter values of the tie lines that bracket 
% the query point (qpt) 

%for i=1:length(xp)
%    tie_dist(i)=norm(a_pts(i,:)-b_pts(i,:));
%    sum_side(i)=norm(a_pts(i,:)-qpt)+norm(b_pts(i,:)-qpt);
%end
   
%[dist_diff,index]=sort(sum_side-tie_dist);

%i1=index(1);
%i2=index(2);
%x1=xp(i1);
%x2=xp(i2); 
%keyboard
%xs=x1;
%xe=x2;
%if (x1 > x2)
%    xs=x2;
%    xe=x1;
%end
warning off
for i=1:length(xp)-1
    if (xp(1) == 0 & xp(end) == 1) 
        dist_2_cp1=norm(a_pts(1,:)-qpt); % or b_pts(1,:)
        dist_2_cp2=norm(a_pts(end,:)-qpt); % or b_pts(end,:)
        dist_2_qpt=zeros(length(xp)-2,1);
        for j=2:length(xp)-1
            dist_2_qpt(j-1)=(((a_pts(j,2)-b_pts(j,2)).*qpt(1))+((b_pts(j,1)-a_pts(j,1)).*qpt(2))+((a_pts(j,1).*b_pts(j,2))-(b_pts(j,1).*a_pts(j,2))))./(norm(b_pts(j,:)-a_pts(j,:)));
        end
        if (all(dist_2_qpt < 0) | all(dist_2_qpt > 0))
            if (dist_2_cp1 < dist_2_cp2)
                i1=1;
                i2=2;
                break
            else
                i1=length(xp)-1;
                i2=i1+1;
                break
            end
        elseif (i == 1)
            continue
        end
    elseif (xp(1) == 0 & xp(end) ~= 1)
        dist_2_qpt=zeros(length(xp)-1,1);
        for j=2:length(xp)
            dist_2_qpt(j-1)=(((a_pts(j,2)-b_pts(j,2)).*qpt(1))+((b_pts(j,1)-a_pts(j,1)).*qpt(2))+((a_pts(j,1).*b_pts(j,2))-(b_pts(j,1).*a_pts(j,2))))./(norm(b_pts(j,:)-a_pts(j,:)));
        end
        if (all(dist_2_qpt < 0) | all(dist_2_qpt > 0))
            i1=1;
            i2=2;
            break
        elseif (i == 1)
            continue
        end
    elseif (xp(1) ~= 0 & xp(end) == 1)
        dist_2_qpt=zeros(length(xp)-1,1);
        for j=1:length(xp)-1
            dist_2_qpt(j)=(((a_pts(j,2)-b_pts(j,2)).*qpt(1))+((b_pts(j,1)-a_pts(j,1)).*qpt(2))+((a_pts(j,1).*b_pts(j,2))-(b_pts(j,1).*a_pts(j,2))))./(norm(b_pts(j,:)-a_pts(j,:)));
        end
        if (all(dist_2_qpt < 0) | all(dist_2_qpt > 0))
            i1=length(xp)-1;
            i2=i1+1;
            break
        end
    end
    %dist_2_qpt(i)=(((a_pts(i,2)-b_pts(i,2)).*qpt(1))+((b_pts(i,1)-a_pts(i,1)).*qpt(2))+((a_pts(i,1).*b_pts(i,2))-(b_pts(i,1).*a_pts(i,2))))./(norm(b_pts(i,:)-a_pts(i,:)));
    dist_now=(((a_pts(i,2)-b_pts(i,2)).*qpt(1))+((b_pts(i,1)-a_pts(i,1)).*qpt(2))+((a_pts(i,1).*b_pts(i,2))-(b_pts(i,1).*a_pts(i,2))))./(norm(b_pts(i,:)-a_pts(i,:)));
    dist_next=(((a_pts(i+1,2)-b_pts(i+1,2)).*qpt(1))+((b_pts(i+1,1)-a_pts(i+1,1)).*qpt(2))+((a_pts(i+1,1).*b_pts(i+1,2))-(b_pts(i+1,1).*a_pts(i+1,2))))./(norm(b_pts(i+1,:)-a_pts(i+1,:)));
    if ((dist_now > 0 & dist_next < 0) | (dist_now < 0 & dist_next > 0))
        i1=i;
        i2=i+1;
        break
    end
    %if (isnan(dist_2_qpt(i))) % a_pts(i,:) == b_pts(i,:) i.e. critical point
    %    dist_2_qpt(i)=norm(a_pts(i,:)-qpt);
    %end
end
warning on

%[dist_sort,index]=sort(abs(dist_2_qpt));
%i1=index(1);
%i2=index(2);
x1=xp(i1);
x2=xp(i2);
xs=x1;
xe=x2;
if (x1 > x2)
    xs=x2;
    xe=x1;
end

return