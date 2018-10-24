function [p3_tielines,p3_ab_pts]=p3_boundary_tielines(tielines,alpha_beta_pts)
% the function p3_boundary_tielines takes a set of tie lines in tielines and generates a another
% set of tielines were some of the tie lines have been deleting leaving a tie line as a possible
% boundary for an adjacent three-phase region

[ni,nj]=size(tielines);
p3_tielines=cell(ni,nj);
p3_ab_pts=cell(ni,nj);
for i=1:ni
    for j=1:nj
        tylines=tielines{i,j};
        [ntielines,ntpts,txy]=size(tylines);
        if (rem(ntielines,2) ~= 0)
	        n=ntielines-1;
            half_ntie=n/2;
        else
            half_ntie=ntielines/2;
        end
        p3_tielines{i,j}=cell(2,half_ntie);
        p3_ab_pts{i,j}=cell(2,half_ntie);
        for k=1:half_ntie
            tylines=tielines{i,j};
            ab_pts=alpha_beta_pts{i,j};
            tylines(1:k,:,:)=[];
            ab_pts(1:k,:,:)=[];
            p3_tielines{i,j}{1,k}=tylines;
            p3_ab_pts{i,j}{1,k}=ab_pts;
        end
        iset=1;
        for k=1:half_ntie
            tylines=tielines{i,j};
            ab_pts=alpha_beta_pts{i,j};
            tylines(ntielines:-1:ntielines-(k-1),:,:)=[];
            ab_pts(ntielines:-1:ntielines-(k-1),:,:)=[];
            p3_tielines{i,j}{2,iset}=tylines;
            p3_ab_pts{i,j}{2,iset}=ab_pts;
            iset=iset+1;
        end
    end
end    

%figure;
%p3=p3t{2,5};
%plot_tielines(p3);
%ab_pts=p3_ab_pts{2,20}{2,5};
%n=size(ab_pts,1);
%for i=1:n
%hold on,plot(ab_pts(i,1,1),ab_pts(i,2,1),'or'),plot(ab_pts(i,1,2),ab_pts(i,2,2),'sb'),hold off;
%pause;
%end