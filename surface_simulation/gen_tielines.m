function [tielines,alpha_beta_pts]=gen_tielines(boundary_pts,npts,ntielines)
% the function gen_tielines generates non-intersecting tielines traversing a
% 2-phase region defined by the array of points in boundary_pts
[nbpts,bxy]=size(boundary_pts);
if (bxy ~= 2)
	error('num cols of boundary points not equal to 2 (x and y)')
end

% loop over all boundary pts
tielines=cell(nbpts,(nbpts-2));
alpha_beta_pts=cell(nbpts,(nbpts-2));
for i=1:nbpts
    % loop over possible critical point pairs
    for j=2:(nbpts-1)
        crit_pt1=boundary_pts(1,:); 
        crit_pt2=boundary_pts(j,:);
        % calc distance of alpha region (1->j) and beta region (j+1->nbpts-1)
        % calc total distance (approx perimeter as piece-wise linear) and alpha dist
        total_dist=0;
        total_alfa_dist=0;
        total_beta_dist=0;
        ialfa=1;
        ibeta=1;
        for current=1:(nbpts-1)
            next=current+1;
            pt2pt_dist=sqrt(((boundary_pts(next,1)-boundary_pts(current,1)).^2)+((boundary_pts(next,2)-boundary_pts(current,2)).^2));
            if (next <= j)
                alfa_dists(ialfa)=pt2pt_dist;
                total_alfa_dist=total_alfa_dist+pt2pt_dist;
                ialfa=ialfa+1;
            end
            if (next > j & next <= nbpts)
                beta_dists(ibeta)=pt2pt_dist;
                total_beta_dist=total_beta_dist+pt2pt_dist;
                ibeta=ibeta+1;
            end
            total_dist=total_dist+pt2pt_dist;
        end
        % divide the alpha region into number of points = ntielines
        % alpha_pts=zeros(ntielines,2);
        dist_per_apt=total_alfa_dist./(ntielines+1);
        nalfa_pts=0;
        avail_dist=alfa_dists(1);
        dist_2go=dist_per_apt;
        dist_traveled=0;
        ialfa=1;
        jalfa=2;
        while (dist_traveled < (total_alfa_dist-dist_per_apt))  
            if (dist_2go > avail_dist)
                ialfa=ialfa+1;
                jalfa=ialfa+1;
                if (jalfa > j)
                    break;
                end
                dist_2go=dist_2go-avail_dist;
                %dist_traveled=dist_traveled+avail_dist;
                avail_dist=alfa_dists(jalfa-1);      
            elseif (dist_2go < avail_dist)
                %keyboard
                % m = slope of line to next boundary pt
                deltax=boundary_pts(jalfa,1)-boundary_pts(ialfa,1);
                deltay=boundary_pts(jalfa,2)-boundary_pts(ialfa,2);
                warning off
                m=deltay./deltax;
                warning on
                % d = distance from ialfa boundary pt to the alpha point
                d=dist_2go;
                x_inc=(d./(sqrt(1+(m.^2))));
                y_inc=((m.*d)./(sqrt(1+(m.^2))));
                % determine direction of changing x and y
                if (deltax < 0 & deltay < 0)
                    x_inc=-x_inc;
                    y_inc=-y_inc;
                elseif (deltax < 0 & deltay == 0)
                    x_inc=-x_inc;
                    y_inc=0;
                elseif (deltax < 0 & deltay > 0)
                    x_inc=-x_inc;
                    y_inc=-y_inc;
                elseif (deltax == 0 & deltay < 0)
                    x_inc=0;
                    y_inc=-d;
                elseif (deltax == 0 & deltay > 0)
                    x_inc=0;
                    y_inc=d;
                elseif (deltax == 0 & deltay == 0)%this should never happen
                    x_inc=0;
                    y_inc=0;
                end
                % increment number alpha points counter
                nalfa_pts=nalfa_pts+1; 
                if (nalfa_pts == 1 | avail_dist == alfa_dists(jalfa-1))
                    % x coord of alpha point
                    alpha_pts(nalfa_pts,1)=boundary_pts(ialfa,1)+x_inc;
                    % y coord of alpha point
                    alpha_pts(nalfa_pts,2)=boundary_pts(ialfa,2)+y_inc;
                    dist_from_ialfa=sqrt(((alpha_pts(nalfa_pts,1)-boundary_pts(ialfa,1)).^2)+((alpha_pts(nalfa_pts,2)-boundary_pts(ialfa,2)).^2));
                    %if (nalfa_pts == 1)
                    %    dist_from_preva=dist_from_ialfa;
                    %else
                    %    dist_2_ialfa=sqrt(((boundary_pts(ialfa,1)-alpha_pts(nalfa_pts-1,1)).^2)+((boundary_pts(ialfa,2)-alpha_pts(nalfa_pts-1,2)).^2));
                    %    dist_from_preva=dist_2_ialfa+dist_from_ialfa;
                    %end
                    dist_from_preva=dist_from_ialfa;
                else
                    % x coord of alpha point
                    alpha_pts(nalfa_pts,1)=alpha_pts(nalfa_pts-1,1)+x_inc;
                    % y coord of alpha point
                    alpha_pts(nalfa_pts,2)=alpha_pts(nalfa_pts-1,2)+y_inc;
                    dist_from_ialfa=sqrt(((alpha_pts(nalfa_pts,1)-boundary_pts(ialfa,1)).^2)+((alpha_pts(nalfa_pts,2)-boundary_pts(ialfa,2)).^2));
                    dis_from_preva=sqrt(((alpha_pts(nalfa_pts,1)-alpha_pts(nalfa_pts-1,1)).^2)+((alpha_pts(nalfa_pts,2)-alpha_pts(nalfa_pts-1,2)).^2));
                end
                avail_dist=alfa_dists(jalfa-1)-dist_from_ialfa;
                dist_2go=dist_per_apt;
                %dist_traveled=dist_traveled+dist_from_preva;
                dist_traveled=dist_traveled+dist_per_apt;
            else
                %keyboard
                %dist_traveled=dist_traveled+avail_dist;
                dist_traveled=dist_traveled+dist_per_apt;
                nalfa_pts=nalfa_pts+1;
                alpha_pts(nalfa_pts,1)=boundary_pts(jalfa,1);
                alpha_pts(nalfa_pts,2)=boundary_pts(jalfa,2);
                dist_2go=dist_per_apt;
                avail_dist=alfa_dists(jalfa-1);            
                ialfa=ialfa+1;
                jalfa=ialfa+1;
            end
        end
        %keyboard
        % divide the beta region into number of points = ntielines
        dist_per_bpt=total_beta_dist./(ntielines+1);
        nbeta_pts=0;
        % beta_pts=zeros(ntielines,2);
        avail_dist=beta_dists(1);
        dist_2go=dist_per_bpt;
        dist_traveled=0;
        ibeta=j;
        jbeta=j+1;
        while (dist_traveled < (total_beta_dist-dist_per_bpt)) 
            if (dist_2go > avail_dist)
                %keyboard
                ibeta=ibeta+1;
                jbeta=ibeta+1;
                if (jbeta > nbpts)
                    break;
                end
                dist_2go=dist_2go-avail_dist;
                %dist_traveled=dist_traveled+avail_dist;
                avail_dist=beta_dists(jbeta-j);      
            elseif (dist_2go < avail_dist)
                %keyboard
                % m = slope of line to next boundary pt
                deltax=boundary_pts(jbeta,1)-boundary_pts(ibeta,1);
                deltay=boundary_pts(jbeta,2)-boundary_pts(ibeta,2);
                warning off
                m=deltay./deltax;
                warning of
                % d = distance from ibeta boundary pt to the beta point
                d=dist_2go;
                x_inc=(d./(sqrt(1+(m.^2))));
                y_inc=((m.*d)./(sqrt(1+(m.^2))));
                % determine direction of changing x and y
                if (deltax < 0 & deltay < 0)
                    x_inc=-x_inc;
                    y_inc=-y_inc;
                elseif (deltax < 0 & deltay == 0)
                    x_inc=-x_inc;
                    y_inc=0;
                elseif (deltax < 0 & deltay > 0)
                    x_inc=-x_inc;
                    y_inc=-y_inc;
                elseif (deltax == 0 & deltay < 0)
                    x_inc=0;
                    y_inc=-d;
                elseif (deltax == 0 & deltay > 0)
                    x_inc=0;
                    y_inc=d;
                elseif (deltax == 0 & deltay == 0)
                    x_inc=0;
                    y_inc=0;
                end
                % increment number beta points counter
                nbeta_pts=nbeta_pts+1; 
                if (nbeta_pts == 1 | avail_dist == beta_dists(jbeta-j))
                    % x coord of beta point
                    beta_pts(nbeta_pts,1)=boundary_pts(ibeta,1)+x_inc;
                    % y coord of beta point
                    beta_pts(nbeta_pts,2)=boundary_pts(ibeta,2)+y_inc;
                    dist_from_ibeta=sqrt(((beta_pts(nbeta_pts,1)-boundary_pts(ibeta,1)).^2)+((beta_pts(nbeta_pts,2)-boundary_pts(ibeta,2)).^2));
                    %if (nbeta_pts == 1)
                    %    dist_from_prevb=dist_from_ibeta;
                    %else
                    %    dist_2_ibeta=sqrt(((boundary_pts(ibeta,1)-beta_pts(nbeta_pts-1,1)).^2)+((boundary_pts(ibeta,2)-beta_pts(nbeta_pts-1,2)).^2));
                    %    dist_from_prevb=dist_2_ibeta+dist_from_ibeta;
                    %end
                    dist_from_prevb=dist_from_ibeta;
                else
                    % x coord of alpha point
                    beta_pts(nbeta_pts,1)=beta_pts(nbeta_pts-1,1)+x_inc;
                    % y coord of alpha point
                    beta_pts(nbeta_pts,2)=beta_pts(nbeta_pts-1,2)+y_inc;
                    dist_from_ibeta=sqrt(((beta_pts(nbeta_pts,1)-boundary_pts(ibeta,1)).^2)+((beta_pts(nbeta_pts,2)-boundary_pts(ibeta,2)).^2));
                    dis_from_prevb=sqrt(((beta_pts(nbeta_pts,1)-beta_pts(nbeta_pts-1,1)).^2)+((beta_pts(nbeta_pts,2)-beta_pts(nbeta_pts-1,2)).^2));
                end
                avail_dist=beta_dists(jbeta-j)-dist_from_ibeta;
                dist_2go=dist_per_bpt;
                %dist_traveled=dist_traveled+dist_from_prevb;
                dist_traveled=dist_traveled+dist_per_bpt;
            else
                %keyboard
                %dist_traveled=dist_traveled+avail_dist;
                dist_traveled=dist_traveled+dist_per_bpt;
                nbeta_pts=nbeta_pts+1;
                beta_pts(nbeta_pts,1)=boundary_pts(jbeta,1);
                beta_pts(nbeta_pts,2)=boundary_pts(jbeta,2);
                dist_2go=dist_per_bpt;
                avail_dist=beta_dists(jbeta-j);            
                ibeta=ibeta+1;
                jbeta=ibeta+1;
            end
        end
        %keyboard
        % alpha points array (alpha_pts) and beta points array (beta_pts) may contain critical 
        % point at its end so delete
        if (nalfa_pts > ntielines)
            alpha_pts(end,:)=[];
            %error('number alpha points not equal to number beta points')
        end
        if (nbeta_pts > ntielines)
            beta_pts(end,:)=[];
            %error('number beta points not equal to number tie lines')
        end
        na_points=size(alpha_pts,1);
        nb_points=size(beta_pts,1);
        if (na_points ~= ntielines)
            error('number alpha points not equal to number tie lines')
        end
        if (nb_points ~= ntielines)
            error('number beta points not equal to number tie lines')
        end
        % connect each alpha pt with appropriate beta point with a line
        % this line is a tie line, and the collection of these lines is the tie line set
        % for these two critical points
        ab_pts=zeros(ntielines,bxy,bxy);
        %ab_pts=cat(3,alpha_pts,beta_pts);% in 3rd dim, 1 = alpha, 2 = beta
        ab_pts=cat(3,alpha_pts,flipud(beta_pts));% in 3rd dim, 1 = alpha, 2 = beta
        alpha_beta_pts{i,(j-1)}=ab_pts; 
        %clear alpha_pts;
        %clear beta_pts;
        tylines=zeros(ntielines,npts,2);
        for k=1:ntielines
            xline=linspace(ab_pts(k,1,1),ab_pts(k,1,2),npts);
            yline=linspace(ab_pts(k,2,1),ab_pts(k,2,2),npts);
            %xline=linspace(ab_pts(k,1,1),ab_pts(ntielines-(k-1),1,2),npts);
            %yline=linspace(ab_pts(k,2,1),ab_pts(ntielines-(k-1),2,2),npts);
            tylines(k,:,1)=xline;
            tylines(k,:,2)=yline;
        end
        tielines{i,(j-1)}=tylines;
    end %for j=1:nbpts
    
    % shift numbering of boundary points to right (clockwise)
    rep_bounds(2:nbpts,:)=boundary_pts(1:(nbpts-1),:);
    rep_bounds(1,:)=boundary_pts(nbpts,:);
    boundary_pts=rep_bounds;
end



% below is the tieline creation code that creates a tieline between 1st and last boundary pt, 2nd and 2nd last, etc.
%tielines=zeros(ntielines,npts,bxy,nbpts);
%for i=1:nbpts
%    itie=nbpts;
%    for j=1:ntielines
%        x=linspace(boundary_pts(j,1),boundary_pts(itie-(j-1),1),npts);
%        y=linspace(boundary_pts(j,2),boundary_pts(itie-(j-1),2),npts);
%        tielines(j,:,1,i)=x;
%        tielines(j,:,2,i)=y;
%       itie=itie-1;
%    end
%   shift numbering of boundary points to right (clockwise)
%    rep_bounds(2:nbpts,:)=boundary_pts(1:(nbpts-1),:);
%    rep_bounds(1,:)=boundary_pts(nbpts,:);
%    boundary_pts=rep_bounds;
%end